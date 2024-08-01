/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/controlflow/AgglomerationWrapper.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void AgglomerationConvergenceTest( const std::string& meshFile,
                                   const uint_t&      minLevel,
                                   const uint_t&      maxLevel,
                                   const real_t&      targetError,
                                   const bool&        superman )
{
   WALBERLA_CHECK_LESS( minLevel, maxLevel );

   const uint_t numIterations = 5;
   const bool vtk = false;
   const int rank = walberla::mpi::MPIManager::instance()->rank();

   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   const auto numberOfProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hyteg::P2ConstantLaplaceOperator L( storage, minLevel, maxLevel );

   hyteg::P2Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, maxLevel, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );

   hyteg::VTKOutput vtkOutput( "../../output", "AgglomerationConvergenceTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( f );
   vtkOutput.add( r );
   vtkOutput.add( err );
   vtkOutput.add( npoints_helper );
   if ( vtk )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );

   WALBERLA_LOG_INFO_ON_ROOT( "DoFs:" )
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      const auto globalDoFs = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
      WALBERLA_LOG_INFO_ON_ROOT( "level " << level << ": " << globalDoFs );
   }


   std::shared_ptr< Solver< P2ConstantLaplaceOperator > > coarseGridSolver;
   std::shared_ptr< PETScLUSolver< P2ConstantLaplaceOperator > > coarseGridSolverInner;
   std::shared_ptr< PrimitiveStorage > agglomerationStorage;
   std::shared_ptr< AgglomerationWrapper< P2ConstantLaplaceOperator > > agglomerationWrapper;

   const uint_t numberOfSupermanProcesses = 2;

   if ( superman )
   {
      // Setup of the agglomeration based solver.
      // Apart from the coarse grid, everything is performed in parallel on all but the superman process.
      // If the superman strategy is active, we use one process for the factorization.
      WALBERLA_LOG_INFO_ON_ROOT( "Total number of processes (incl. superman): "
                                 << numberOfProcesses << ". Number of superman processes: " << numberOfSupermanProcesses )

      loadbalancing::distributed::roundRobin( *storage, numberOfProcesses - numberOfSupermanProcesses );

      bool solveOnEmptyProcesses = false;

      agglomerationWrapper =
          std::make_shared< AgglomerationWrapper< P2ConstantLaplaceOperator > >( storage, minLevel, solveOnEmptyProcesses );
      agglomerationWrapper->setStrategyContinuousProcesses( numberOfProcesses - numberOfSupermanProcesses,
                                                            numberOfProcesses - 1 );
      agglomerationStorage = agglomerationWrapper->getAgglomerationStorage();

      auto numprimitivesagglomerationstorage = agglomerationStorage->getNumberOfLocalPrimitives();
      auto numprimitivesnormalstorage = storage->getNumberOfLocalPrimitives();
      WALBERLA_LOG_DEVEL( "normal primitives: " << numprimitivesnormalstorage )
      WALBERLA_LOG_DEVEL( "agg primitives: " << numprimitivesagglomerationstorage )

      coarseGridSolverInner =
          std::make_shared< PETScLUSolver< P2ConstantLaplaceOperator > >( agglomerationStorage, minLevel );
      coarseGridSolverInner->setManualAssemblyAndFactorization( true );
      coarseGridSolverInner->setReassembleMatrix( true );

      // now wrap the solver
      agglomerationWrapper->setSolver( coarseGridSolverInner );

      coarseGridSolver = agglomerationWrapper;
   }
   else
   {
      coarseGridSolverInner = std::make_shared< PETScLUSolver< P2ConstantLaplaceOperator > >( storage, minLevel );
      coarseGridSolverInner->setReassembleMatrix( true );
      coarseGridSolver = coarseGridSolverInner;
   }

   auto smoother     = std::make_shared< GaussSeidelSmoother< P2ConstantLaplaceOperator > >();
   auto prolongation = std::make_shared< P2toP2QuadraticProlongation >();
   auto restriction  = std::make_shared< P2toP2QuadraticRestriction >();

   auto solver = std::make_shared< GeometricMultigridSolver< P2ConstantLaplaceOperator > >(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel );

   real_t discr_l2_err, discr_l2_residual;
   for ( uint_t iteration = 0; iteration < numIterations; iteration++ )
   {
      if ( superman && uint_c( rank ) >= ( numberOfProcesses - numberOfSupermanProcesses ) )
      {
         coarseGridSolverInner->assembleAndFactorize( *agglomerationWrapper->getAgglomerationOperator() );
      }

      solver->solve( L, u, f, maxLevel );

      err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
      npoints_helper.interpolate( ones, maxLevel );

      const real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );
      discr_l2_err         = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

      L.apply( u, err, maxLevel, Inner | NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, err}, maxLevel, Inner | NeumannBoundary );
      discr_l2_residual = std::sqrt( r.dotGlobal( r, maxLevel, Inner | NeumannBoundary ) / npoints );

      WALBERLA_LOG_INFO_ON_ROOT( "residual " << discr_l2_residual << ", error " << discr_l2_err );

      if ( vtk )
      {
         vtkOutput.write( maxLevel, iteration + 1 );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err << ", (mesh: " << meshFile << ")" );
   WALBERLA_CHECK_LESS( discr_l2_err, targetError );

//   printTimingTree( *storage->getTimingTree() );
//   if ( superman )
//   {
//      printTimingTree( *agglomerationStorage->getTimingTree() );
//   }

}

} // namespace hyteg

void runTests( bool superman )
{
   hyteg::AgglomerationConvergenceTest( hyteg::prependHyTeGMeshDir( "annulus_coarse.msh" ), 4, 8, 3e-7, superman );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );
   WALBERLA_LOG_INFO_ON_ROOT(
       "MUMPS direct solver, no superman - factorization in parallel on all processes right before coarse grid solution phase." )
   runTests( false );

   WALBERLA_LOG_INFO_ON_ROOT( "MUMPS direct solver, with superman - factorization on dedicated process during v-cycle." )
   runTests( true );

}
