#include <cmath>

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultiGrid.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSolver.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesPressureBlockPreconditioner.hpp"

using walberla::real_c;
using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t                ntan   = 3;
   const std::vector< double > layers = {1.0, 2.0, 3.0};

   const double rmin = layers.front();
   const double rmax = layers.back();

   const Point3D sourcePoint  = Point3D( {rmin, 0, 0} ) + 0.5 * Point3D( {rmax - rmin, 0, 0} );
   const real_t  sourceRadius = 0.5;

   const uint_t minLevel            = 2;
   const uint_t maxLevel            = 6;
   const uint_t numVCycles          = 1;
   const uint_t maxMinResIterations = 100;
   const real_t tolerance           = 1e-16;

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::meshSphericalShell( ntan, layers );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< hhg::PrimitiveStorage >  storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   storage->setTimingTree( timingTree );

   hhg::writeDomainPartitioningVTK( storage, "../../output", "StokesSphere_domain" );

   hhg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );

   uint_t totalGlobalDofsStokes = 0;
   for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      uint_t tmpDofStokes = numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
      totalGlobalDofsStokes += tmpDofStokes;
   }
   WALBERLA_LOG_INFO_ON_ROOT( "Total Stokes DoFs on all level :" << totalGlobalDofsStokes );

   hhg::VTKOutput vtkOutput( "../output", "StokesSphere" );
   vtkOutput.set3D();
   vtkOutput.add( &u.u );
   vtkOutput.add( &u.v );
   vtkOutput.add( &u.w );
   vtkOutput.add( &u.p );
   vtkOutput.add( &f.u );
   vtkOutput.add( &f.v );
   vtkOutput.add( &f.w );
   vtkOutput.add( &f.p );

   hhg::P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > rhsPlumeX = [rmin, rmax, sourcePoint, sourceRadius]( const hhg::Point3D& x ) {
      const real_t distToSourcePoint = ( x - sourcePoint ).norm();
      if( distToSourcePoint < sourceRadius )
         return x[0] * ( sourceRadius - distToSourcePoint );
      else
         return 0.0;
   };

   std::function< real_t( const hhg::Point3D& ) > rhsPlumeY = [rmin, rmax, sourcePoint, sourceRadius]( const hhg::Point3D& x ) {
      const real_t distToSourcePoint = ( x - sourcePoint ).norm();
      if( distToSourcePoint < sourceRadius )
         return x[1] * ( sourceRadius - distToSourcePoint );
      else
         return 0.0;
   };

   std::function< real_t( const hhg::Point3D& ) > rhsPlumeZ = [rmin, rmax, sourcePoint, sourceRadius]( const hhg::Point3D& x ) {
      const real_t distToSourcePoint = ( x - sourcePoint ).norm();
      if( distToSourcePoint < sourceRadius )
         return x[2] * ( sourceRadius - distToSourcePoint );
      else
         return 0.0;
   };

   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   f.u.interpolate( rhsPlumeX, maxLevel );
   f.v.interpolate( rhsPlumeY, maxLevel );
   f.w.interpolate( rhsPlumeZ, maxLevel );

   vtkOutput.write( maxLevel, 0 );
#if 1
   typedef CGSolver< hhg::P1Function< real_t >, hhg::P1ConstantLaplaceOperator > CoarseGridSolver_T;
   typedef GMultigridSolver< hhg::P1Function< real_t >,
                             hhg::P1ConstantLaplaceOperator,
                             CoarseGridSolver_T,
                             hhg::P1toP1LinearRestriction,
                             hhg::P1toP1LinearProlongation >
       GMGSolver_T;
   typedef StokesBlockDiagonalPreconditioner< hhg::P1StokesFunction< real_t >,
                                              hhg::P1ConstantLaplaceOperator,
                                              GMGSolver_T,
                                              hhg::P1LumpedInvMassOperator >
       Preconditioner_T;
   typedef StokesPressureBlockPreconditioner< hhg::P1StokesFunction< real_t >, hhg::P1LumpedInvMassOperator >
       PressurePreconditioner_T;

   auto                          coarseGridSolver = std::make_shared< CoarseGridSolver_T >( storage, minLevel, maxLevel );
   hhg::P1toP1LinearProlongation prolongationOperator;
   hhg::P1toP1LinearRestriction  restrictionOperator;
   GMGSolver_T gmgSolver( storage, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2 );
   P1LumpedInvMassOperator  massOperator( storage, minLevel, maxLevel );
   PressurePreconditioner_T pressurePrec( massOperator, storage, minLevel, maxLevel );
   Preconditioner_T         prec( L.A, gmgSolver, massOperator, storage, minLevel, maxLevel, 2 );

   typedef hhg::MinResSolver< hhg::P1StokesFunction< real_t >, hhg::P1StokesOperator, Preconditioner_T > PreconditionedMinRes_T;
   typedef hhg::MinResSolver< hhg::P1StokesFunction< real_t >, hhg::P1StokesOperator, PressurePreconditioner_T >
       PressurePreconditionedMinRes_T;

   auto preconditionedMinResSolver         = PreconditionedMinRes_T( storage, minLevel, maxLevel, prec );
   auto pressurePreconditionedMinResSolver = PressurePreconditionedMinRes_T( storage, minLevel, maxLevel, pressurePrec );
   // auto solver = hhg::MinResSolver< hhg::P1StokesFunction< real_t >, hhg::P1StokesOperator >( storage, minLevel, maxLevel );

#if 0
  preconditionedMinResSolver.solve( L, u, f, r, maxLevel, tolerance, maxMinResIterations, hhg::Inner | hhg::NeumannBoundary, true );
#else

   typedef UzawaSolver< hhg::P1StokesFunction< real_t >,
                        hhg::P1StokesOperator,
                        PressurePreconditionedMinRes_T,
                        P1P1StokesToP1P1StokesRestriction,
                        P1P1StokesToP1P1StokesProlongation,
                        false >
                                      UzawaSolver_T;
   P1P1StokesToP1P1StokesRestriction  stokesRestriction;
   P1P1StokesToP1P1StokesProlongation stokesProlongation;
   UzawaSolver_T                      uzawaSolver(
       storage, pressurePreconditionedMinResSolver, stokesRestriction, stokesProlongation, minLevel, maxLevel, 2, 2, 2 );

   for( uint_t i = 0; i < numVCycles; i++ )
   {
      uzawaSolver.solve(
          L, u, f, r, maxLevel, tolerance, 10000, hhg::Inner | hhg::NeumannBoundary, UzawaSolver_T::CycleType::VCYCLE, false );
      L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      real_t residualMG =
          r.dot( r, maxLevel, hhg::Inner ) / real_c( hhg::numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, maxLevel ) );
      WALBERLA_LOG_INFO_ON_ROOT( "after it " << i << ": " << std::scientific << residualMG );
   }

#endif
#endif

#if 0
  auto numerator = std::make_shared< hhg::P1StokesFunction< PetscInt > >( "numerator", storage, level, level );
   uint_t globalSize = 0;
   const uint_t localSize = numerator->enumerate(level, globalSize);
   PETScManager petscManager;
   PETScLUSolver< real_t, hhg::P1StokesFunction, hhg::P1StokesOperator > petScLUSolver( numerator, localSize, globalSize );
   f.u.assign( {1.0}, {&u.u}, level, DirichletBoundary );
   f.v.assign( {1.0}, {&u.v}, level, DirichletBoundary );
   f.w.assign( {1.0}, {&u.w}, level, DirichletBoundary );
   petScLUSolver.solve( L, u, f, r, level, tolerance, maxIterations, Inner | NeumannBoundary );
#endif
   vtkOutput.write( maxLevel, 1 );

   L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   real_t final_residual =
       r.dot( r, maxLevel, hhg::Inner ) / real_c( hhg::numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, maxLevel ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residual: " << final_residual );
   walberla::WcTimingTree tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt );

   return EXIT_SUCCESS;
}
