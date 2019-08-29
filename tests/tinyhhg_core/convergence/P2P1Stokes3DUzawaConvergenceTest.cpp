#include <cmath>

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSmoother.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesPressureBlockPreconditioner.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager manager;

   const std::string meshFileName  = "../../data/meshes/3D/cube_24el.msh";
   const uint_t      minLevel      = 2;
   const uint_t      maxLevel      = 3;
   const uint_t      maxIterations = 3;
   const bool        writeVTK      = false;

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::writeDomainPartitioningVTK( storage, "../../output", "P2P1_Stokes_3D_Uzawa_convergence_partitioning" );

   hhg::P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > uExact( "uExact", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > err( "err", storage, minLevel, maxLevel );
   hhg::P2P1TaylorHoodFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );

   hhg::VTKOutput vtkOutput( "../../output", "P2P1_Stokes_3D_Uzawa_convergence", storage );
   vtkOutput.add( u.u );
   vtkOutput.add( u.v );
   vtkOutput.add( u.w );
   vtkOutput.add( u.p );
   vtkOutput.add( uExact.u );
   vtkOutput.add( uExact.v );
   vtkOutput.add( uExact.w );
   vtkOutput.add( uExact.p );
   vtkOutput.add( err.u );
   vtkOutput.add( err.v );
   vtkOutput.add( err.w );
   vtkOutput.add( err.p );

   hhg::P2P1TaylorHoodStokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > inflowPoiseuille = []( const hhg::Point3D& x ) {
      if ( x[2] < 1e-8 )
      {
         return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
      }
      else
      {
         return 0.0;
      }
   };

   std::function< real_t( const hhg::Point3D& ) > solutionPoiseuille = []( const hhg::Point3D& x ) {
      return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
   };

   std::function< real_t( const hhg::Point3D& ) > collidingFlow_x = []( const hhg::Point3D& x ) {
      return real_c( 20 ) * x[0] * x[1] * x[1] * x[1];
   };

   std::function< real_t( const hhg::Point3D& ) > collidingFlow_y = []( const hhg::Point3D& x ) {
      return real_c( 5 ) * x[0] * x[0] * x[0] * x[0] - real_c( 5 ) * x[1] * x[1] * x[1] * x[1];
   };

   std::function< real_t( const hhg::Point3D& ) > collidingFlow_p = []( const hhg::Point3D& xx ) {
      return real_c( 60 ) * std::pow( xx[0], 2.0 ) * xx[1] - real_c( 20 ) * std::pow( xx[1], 3.0 );
   };

   std::function< real_t( const hhg::Point3D& ) > rhs  = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   u.u.interpolate( collidingFlow_x, maxLevel, hhg::DirichletBoundary );
   u.v.interpolate( collidingFlow_y, maxLevel, hhg::DirichletBoundary );

   uExact.u.interpolate( collidingFlow_x, maxLevel );
   uExact.v.interpolate( collidingFlow_y, maxLevel );
   uExact.p.interpolate( collidingFlow_p, maxLevel );

   if ( writeVTK )
      vtkOutput.write( maxLevel, 0 );

   typedef hhg::StokesPressureBlockPreconditioner< hhg::P2P1TaylorHoodStokesOperator, hhg::P1LumpedInvMassOperator >
        PressurePreconditioner_T;
   auto pressurePrec = std::make_shared< PressurePreconditioner_T >( storage, minLevel, maxLevel );

   auto smoother =
       std::make_shared< hhg::UzawaSmoother< hhg::P2P1TaylorHoodStokesOperator > >( storage, minLevel, maxLevel, 0.4 );
   auto restriction      = std::make_shared< hhg::P2P1StokesToP2P1StokesRestriction >( true );
   auto prolongation     = std::make_shared< hhg::P2P1StokesToP2P1StokesProlongation >();
   auto coarseGridSolver = std::make_shared< hhg::PETScLUSolver< hhg::P2P1TaylorHoodStokesOperator > >( storage, minLevel );
   hhg::GeometricMultigridSolver< hhg::P2P1TaylorHoodStokesOperator > solver(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, 3, 3, 0 );

   const uint_t globalDoFsVelocity = hhg::numberOfGlobalDoFs< hhg::P2FunctionTag >( *storage, maxLevel );
   const uint_t globalDoFsPressure = hhg::numberOfGlobalDoFs< hhg::P1FunctionTag >( *storage, maxLevel );

   L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   err.assign( {1.0, -1.0}, {u, uExact}, maxLevel );
   real_t lastResidual =
       std::sqrt( r.dotGlobal( r, maxLevel ) / ( 3 * (real_t) globalDoFsVelocity + real_c( globalDoFsPressure ) ) );

   real_t discr_l2_err_1_u;
   real_t discr_l2_err_1_v;
   real_t discr_l2_err_1_w;
   real_t discr_l2_err_1_p;
   real_t residuum_l2_1;

   for ( uint_t i = 1; i <= maxIterations; i++ )
   {
      solver.solve( L, u, f, maxLevel );

      hhg::vertexdof::projectMean( u.p, maxLevel );
      hhg::vertexdof::projectMean( uExact.p, maxLevel );

      L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );

      err.assign( {1.0, -1.0}, {u, uExact}, maxLevel );

      if ( writeVTK )
      {
         vtkOutput.write( maxLevel, i );
      }

      discr_l2_err_1_u = std::sqrt( err.u.dotGlobal( err.u, maxLevel ) / (real_t) globalDoFsVelocity );
      discr_l2_err_1_v = std::sqrt( err.v.dotGlobal( err.v, maxLevel ) / (real_t) globalDoFsVelocity );
      discr_l2_err_1_w = std::sqrt( err.w.dotGlobal( err.w, maxLevel ) / (real_t) globalDoFsVelocity );
      discr_l2_err_1_p = std::sqrt( err.p.dotGlobal( err.p, maxLevel ) / (real_t) globalDoFsPressure );
      residuum_l2_1 =
          std::sqrt( r.dotGlobal( r, maxLevel ) / ( 3 * (real_t) globalDoFsVelocity + real_c( globalDoFsPressure ) ) );

      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error w = " << discr_l2_err_1_w );
      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
      WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

      const real_t discrResConvRate = residuum_l2_1 / lastResidual;

      lastResidual = residuum_l2_1;

      WALBERLA_LOG_INFO_ON_ROOT( "After " << std::setw( 3 ) << i << " VCycles: Residual: " << std::scientific << residuum_l2_1
                                          << " | convRate: " << discrResConvRate << " | Error L2 u: " << discr_l2_err_1_u
                                          << " | Error L2 p: " << discr_l2_err_1_p );
   }

//   auto tt        = storage->getTimingTree();
//   auto ttreduced = tt->getReduced().getCopyWithRemainder();
//   WALBERLA_LOG_INFO_ON_ROOT( ttreduced );
//
//   nlohmann::json ttjson = nlohmann::json( ttreduced );
//   std::ofstream  o( "/tmp/uzawa.json" );
//   o << ttjson;
//   o.close();

   WALBERLA_LOG_INFO_ON_ROOT( discr_l2_err_1_u + discr_l2_err_1_v + discr_l2_err_1_w )
   WALBERLA_LOG_INFO_ON_ROOT( discr_l2_err_1_p )
   WALBERLA_LOG_INFO_ON_ROOT( residuum_l2_1 )
   WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v + discr_l2_err_1_w, 3e-03 );
   WALBERLA_CHECK_LESS( discr_l2_err_1_p, 0.8 );
   WALBERLA_CHECK_LESS( residuum_l2_1, 5.0e-05 );

   return EXIT_SUCCESS;
}
