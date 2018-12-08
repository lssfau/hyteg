
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSmoother.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesPressureBlockPreconditioner.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

void setRightBFSBoundaryNeumannPoiseuille( SetupPrimitiveStorage& setupStorage )
{
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   const real_t eps = 0.001;

   for( const auto& it : setupStorage.getVertices() )
   {
      if( std::fabs( it.second->getCoordinates()[0] - 1.0 ) < eps && it.second->getCoordinates()[1] > -1.0 + eps &&
          it.second->getCoordinates()[1] < 1.0 - eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }

   for( const auto& it : setupStorage.getEdges() )
   {
      const auto edgeCoordinates = it.second->getCoordinates();
      if( std::fabs( edgeCoordinates[0][0] - 1.0 ) < eps && std::fabs( edgeCoordinates[1][0] - 1.0 ) < eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t minLevel  = 2;
   const uint_t maxLevel  = 3;
   const real_t tolerance = 1e-17;
   const uint_t maxIter   = 10000;
   const bool   writeVTK  = false;

   //create a Rectangle as mesh with 4 triangles
   auto meshInfo =
       hhg::MeshInfo::meshRectangle( hhg::Point2D( {-1, -1} ), hhg::Point2D( {1, 1} ), hhg::MeshInfo::CRISSCROSS, 1, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );


   setRightBFSBoundaryNeumannPoiseuille( setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   storage->setTimingTree( timingTree );

   hhg::writeDomainPartitioningVTK( storage, "../../output", "UzawaConvergenceTestDomain" );

   hhg::P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > Au( "Au", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > err( "err", storage, minLevel, maxLevel );

   const auto setUVelocityBC = []( const hhg::Point3D& x ) -> real_t {
      if( x[0] < -1.0 + 1e-8 )
      {
         return real_c( 1 - x[1] * x[1] );
      } else
      {
         return real_c( 0 );
      }
   };

   const auto solutionU = []( const hhg::Point3D& x ) -> real_t { return real_c( 1 - x[1] * x[1] ); };

   const auto solutionP = []( const hhg::Point3D& x ) -> real_t { return real_c( -2.0 * x[0] ); };

   hhg::P2P1TaylorHoodStokesOperator L( storage, minLevel, maxLevel );

   u.u.interpolate( setUVelocityBC, maxLevel, hhg::DirichletBoundary );
   u_exact.u.interpolate( solutionU, maxLevel );
   u_exact.p.interpolate( solutionP, maxLevel );

   hhg::communication::syncP2FunctionBetweenPrimitives( u_exact.u, maxLevel );
   hhg::communication::syncFunctionBetweenPrimitives( u_exact.p, maxLevel );

   auto pressurePreconditioner = std::make_shared< hhg::StokesPressureBlockPreconditioner< hhg::P2P1TaylorHoodStokesOperator, hhg::P1LumpedInvMassOperator > >(storage, minLevel, maxLevel);
   auto smoother = std::make_shared< hhg::UzawaSmoother<hhg::P2P1TaylorHoodStokesOperator>  >(storage, minLevel, maxLevel, storage->hasGlobalCells(), 0.37);
   auto coarseGridSolver = std::make_shared< hhg::MinResSolver< hhg::P2P1TaylorHoodStokesOperator > >( storage, minLevel, minLevel, maxIter, tolerance, pressurePreconditioner );
   auto restrictionOperator = std::make_shared< hhg::P2P1StokesToP2P1StokesRestriction>();
   auto prolongationOperator = std::make_shared< hhg::P2P1StokesToP2P1StokesProlongation >();

   auto gmgSolver = hhg::GeometricMultigridSolver< hhg::P2P1TaylorHoodStokesOperator >(
      storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );



   const uint_t npoints = hhg::numberOfGlobalDoFs< hhg::P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   real_t       discr_l2_err, currRes, oldRes = 0;

   L.apply( u, Au, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   oldRes = std::sqrt( r.dotGlobal( r, maxLevel, hhg::All ) ) / real_c( npoints );

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel, hhg::All );
   discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel, hhg::All ) / real_c( npoints ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Totalpoints      = " << npoints );
   WALBERLA_LOG_INFO_ON_ROOT( "initial Residual = " << oldRes );
   WALBERLA_LOG_INFO_ON_ROOT( "initial L2 error = " << discr_l2_err );

   hhg::VTKOutput vtkOutput( "../../output", "P2P1UzawaConvergence", storage );
   vtkOutput.add( u.u );
   vtkOutput.add( u.v );
   vtkOutput.add( u.p );
   vtkOutput.add( u_exact.u );
   vtkOutput.add( u_exact.v );
   vtkOutput.add( u_exact.p );
   vtkOutput.add( err.u );
   vtkOutput.add( err.v );
   vtkOutput.add( err.p );

   if( writeVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   for( int j = 0; j < 8; ++j )
   {
      gmgSolver.solve(L, u, f, maxLevel);

      L.apply( u, Au, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      currRes = std::sqrt( r.dotGlobal( r, maxLevel, hhg::All ) ) / real_c( npoints );

      WALBERLA_LOG_INFO_ON_ROOT( "current Residual = " << currRes );
      WALBERLA_CHECK_LESS( currRes / oldRes, 0.6 );
      oldRes = currRes;
   }

   if( writeVTK )
   {
      vtkOutput.write( maxLevel, 1 );
   }

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
   discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel, hhg::Inner ) ) / real_c( npoints );

   WALBERLA_CHECK_LESS( discr_l2_err, 2e-2 );

   walberla::WcTimingTree tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt.getCopyWithRemainder() );

   return 0;
}
