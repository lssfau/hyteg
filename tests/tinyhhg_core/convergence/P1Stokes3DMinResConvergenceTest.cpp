#include <cmath>

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/GaussSeidelSmoother.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesPressureBlockPreconditioner.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;


int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   const std::string meshFileName = "../../data/meshes/3D/cube_24el.msh";
   const uint_t minLevel         =  2;
   const uint_t maxLevel         =  3;
   const uint_t maxIterations    =  5;
   const real_t tolerance = 1e-16;

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

#if 1
   // Get all primitive IDs of primitives at the outflow boundary (z == 1)
   const real_t eps = 1e-8;
   const real_t zBoundary = 1.0;
   for ( const auto & it : setupStorage.getVertices() )
   {
     if ( std::fabs( it.second->getCoordinates()[2] - zBoundary ) < eps &&
          std::fabs( it.second->getCoordinates()[0] ) > eps &&
          std::fabs( it.second->getCoordinates()[1] ) > eps &&
          std::fabs( it.second->getCoordinates()[0] - 1.0 ) > eps &&
          std::fabs( it.second->getCoordinates()[1] - 1.0 ) > eps )
       setupStorage.setMeshBoundaryFlag( it.first, 2 );
   }
  for ( const auto & it : setupStorage.getEdges() )
  {
    if ( std::fabs( it.second->getCoordinates()[0][2] - zBoundary ) < eps &&
         std::fabs( it.second->getCoordinates()[1][2] - zBoundary ) < eps )

      setupStorage.setMeshBoundaryFlag( it.first, 2 );

    if ( std::fabs( it.second->getCoordinates()[0][0] ) < eps &&
         std::fabs( it.second->getCoordinates()[1][0] ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );

    if ( std::fabs( it.second->getCoordinates()[0][0] - 1.0 ) < eps &&
         std::fabs( it.second->getCoordinates()[1][0] - 1.0 ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );

    if ( std::fabs( it.second->getCoordinates()[0][1] ) < eps &&
         std::fabs( it.second->getCoordinates()[1][1] ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );

    if ( std::fabs( it.second->getCoordinates()[0][1] -1.0 ) < eps &&
         std::fabs( it.second->getCoordinates()[1][1] -1.0 ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );
  }
  for ( const auto & it : setupStorage.getFaces() )
  {
    if ( std::fabs( it.second->getCoordinates()[0][2] - zBoundary ) < eps &&
         std::fabs( it.second->getCoordinates()[1][2] - zBoundary ) < eps &&
         std::fabs( it.second->getCoordinates()[2][2] - zBoundary ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 2 );
  }
#endif

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::writeDomainPartitioningVTK( storage, "../../output", "P1_Stokes_3D_MinRes_convergence_partitioning" );

   hhg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > uExact( "uExact", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );

   hhg::VTKOutput vtkOutput( "../../output", "P1_Stokes_3D_MinRes_convergence", storage );
   vtkOutput.add( u.u );
   vtkOutput.add( u.v );
   vtkOutput.add( u.w );
   vtkOutput.add( u.p );
   vtkOutput.add( uExact.u );
   vtkOutput.add( uExact.v );
   vtkOutput.add( uExact.w );
   vtkOutput.add( uExact.p );

   hhg::P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > inflowPoiseuille = []( const hhg::Point3D& x )
   {
      if( x[2] < 1e-8 )
      {
        return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
      }
      else
      {
         return 0.0;
      }
   };


  std::function< real_t( const hhg::Point3D& ) > solutionPoiseuille = []( const hhg::Point3D& x )
  {
      return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
  };

  std::function< real_t( const hhg::Point3D& ) > collidingFlow_x = []( const hhg::Point3D& x )
  {
    return real_c(20) * x[0] * x[1] * x[1] * x[1];
  };

  std::function< real_t( const hhg::Point3D& ) > collidingFlow_y = []( const hhg::Point3D& x )
  {
      return real_c(5) * x[0] * x[0] * x[0] * x[0] - real_c(5) * x[1] * x[1] * x[1] * x[1];
  };

   std::function< real_t( const hhg::Point3D& ) > rhs  = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

#if 1
   u.w.interpolate( inflowPoiseuille, maxLevel, hhg::DirichletBoundary );
#else
   u.u.interpolate( collidingFlow_x, maxLevel, hhg::DirichletBoundary );
   u.v.interpolate( collidingFlow_y, maxLevel, hhg::DirichletBoundary );

   uExact.u.interpolate( collidingFlow_x, maxLevel );
   uExact.v.interpolate( collidingFlow_y, maxLevel );
#endif

   vtkOutput.write( maxLevel, 0 );
#if 1

   typedef hhg::CGSolver< hhg::P1ConstantLaplaceOperator > CoarseGridSolver_T;
   typedef hhg::GeometricMultigridSolver< hhg::P1ConstantLaplaceOperator  > GMGSolver_T;
   typedef hhg::StokesBlockDiagonalPreconditioner< hhg::P1StokesOperator, hhg::P1LumpedInvMassOperator > Preconditioner_T;

   auto coarseGridSolver = std::make_shared< CoarseGridSolver_T  >( storage, minLevel, maxLevel );
   auto smoother = std::make_shared< hhg::GaussSeidelSmoother<hhg::P1ConstantLaplaceOperator>  >();
   auto prolongationOperator = std::make_shared< hhg::P1toP1LinearProlongation >();
   auto restrictionOperator = std::make_shared< hhg::P1toP1LinearRestriction >();
   auto gmgSolver = std::make_shared< GMGSolver_T >( storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2 );
   //hhg::P1LumpedInvMassOperator massOperator( storage, minLevel, maxLevel );
   //Preconditioner_T prec( storage, minLevel, maxLevel, 2, gmgSolver );
   auto prec = std::make_shared< Preconditioner_T >( storage, minLevel, maxLevel, 2, gmgSolver );

   auto solver = hhg::MinResSolver< hhg::P1StokesOperator >( storage, minLevel, maxLevel, maxIterations, tolerance, prec );
   // auto solver = hhg::MinResSolver< hhg::P1StokesFunction< real_t >, hhg::P1StokesOperator, PressurePreconditioner_T >( storage, minLevel, maxLevel, pressurePrec );
   // auto solver = hhg::MinResSolver< hhg::P1StokesFunction< real_t >, hhg::P1StokesOperator >( storage, minLevel, maxLevel );

   solver.solve( L, u, f, maxLevel );
#else
   auto numerator = std::make_shared< hhg::P1StokesFunction< PetscInt > >( "numerator", storage, level, level );
   uint_t globalSize = 0;
   const uint_t localSize = numerator->enumerate(level, globalSize);
   PETScManager petscManager;
   PETScLUSolver< real_t, hhg::P1StokesFunction, hhg::P1StokesOperator > petScLUSolver( numerator, localSize, globalSize );
   f.u.assign( {1.0}, {u.u}, level, DirichletBoundary );
   f.v.assign( {1.0}, {u.v}, level, DirichletBoundary );
   f.w.assign( {1.0}, {u.w}, level, DirichletBoundary );
   petScLUSolver.solve( L, u, f, r, level, tolerance, maxIterations, Inner | NeumannBoundary );
#endif
   vtkOutput.write( maxLevel, 1 );

   L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   real_t final_residual = r.dotGlobal( r, maxLevel, hhg::Inner ) / real_c( hhg::numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, maxLevel ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residual: " << final_residual )
   WALBERLA_CHECK_LESS( final_residual, 3.1e-12 );

   return EXIT_SUCCESS;
}
