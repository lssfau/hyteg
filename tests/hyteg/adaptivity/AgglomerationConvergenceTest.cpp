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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/controlflow/AgglomerationWrapper.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void AgglomerationConvergenceTest( const std::string& meshFile,
                                   const uint_t &     minLevel,
                                   const uint_t &     maxLevel,
                                   const real_t &     targetError,
                                   const bool   &     localMPI,
                                   const bool   &     usePetscOnCoarseGrid)
{
   WALBERLA_CHECK_LESS( minLevel, maxLevel );

   const uint_t numIterations = 3;

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

   if ( localMPI )
   {
      u.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
   }

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, maxLevel, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );

//   hyteg::VTKOutput vtkOutput( "../../output", "AgglomerationConvergenceTest", storage );
//   vtkOutput.add( u );
//   vtkOutput.add( u_exact );
//   vtkOutput.add( f );
//   vtkOutput.add( r );
//   vtkOutput.add( err );
//   vtkOutput.add( npoints_helper );
//   vtkOutput.write( maxLevel, 0 );

#if 1
   // Setup of the agglomeration based solver.
   // Apart from the coarse grid, everything is performed in parallel.
   // The coarse grid problem is, however, solved on a subset of processes.
   const uint_t numberOfSubsetProcesses = ( numberOfProcesses + 1 ) / 2;
   WALBERLA_LOG_INFO_ON_ROOT( "Agglomeration from " << numberOfProcesses << " to " << numberOfSubsetProcesses << " processes." )

   bool solveOnEmptyProcesses = true;
#ifdef HYTEG_BUILD_WITH_PETSC
   // currently the HyTeG solvers rely on global communication, only the PETSc solve() calls can be skipped on empty processes
   if ( usePetscOnCoarseGrid )
      solveOnEmptyProcesses = false;
#endif

   auto agglomerationWrapper =
       std::make_shared< AgglomerationWrapper< P2ConstantLaplaceOperator > >( storage, minLevel, solveOnEmptyProcesses );
   agglomerationWrapper->setStrategyContinuousProcesses( 0, numberOfSubsetProcesses - 1 );
   auto agglomerationStorage = agglomerationWrapper->getAgglomerationStorage();

   // pass agglomeration storage to coarse grid solver
   std::shared_ptr< Solver< P2ConstantLaplaceOperator > > coarseGridSolver;
   auto cgSolver = std::make_shared< CGSolver< P2ConstantLaplaceOperator > >( agglomerationStorage, minLevel, minLevel );
   cgSolver->setPrintInfo( false );
   coarseGridSolver = cgSolver;
#ifdef HYTEG_BUILD_WITH_PETSC
   if ( usePetscOnCoarseGrid )
   {
      coarseGridSolver = std::make_shared< PETScLUSolver< P2ConstantLaplaceOperator > >( agglomerationStorage, minLevel );
   }
#endif

   WALBERLA_UNUSED( usePetscOnCoarseGrid );

   // now wrap the solver
   agglomerationWrapper->setSolver( coarseGridSolver );

   auto smoother     = std::make_shared< GaussSeidelSmoother< P2ConstantLaplaceOperator > >();
   auto prolongation = std::make_shared< P2toP2QuadraticProlongation >();
   auto restriction  = std::make_shared< P2toP2QuadraticRestriction >();
#else
  auto agglomerationWrapper = std::make_shared< CGSolver< P2ConstantLaplaceOperator > >( storage, minLevel, minLevel );
#endif

   auto solver = std::make_shared< GeometricMultigridSolver< P2ConstantLaplaceOperator > >(
       storage, smoother, agglomerationWrapper, restriction, prolongation, minLevel, maxLevel );

   real_t discr_l2_err, discr_l2_residual;
   for ( uint_t iteration = 0; iteration < numIterations; iteration++ )
   {
      solver->solve( L, u, f, maxLevel );

      err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
      npoints_helper.interpolate( ones, maxLevel );

      const real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );
      discr_l2_err         = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

      L.apply( u, err, maxLevel, Inner | NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, err}, maxLevel, Inner | NeumannBoundary );
      discr_l2_residual = std::sqrt( r.dotGlobal( r, maxLevel, Inner | NeumannBoundary ) / npoints );

      WALBERLA_LOG_INFO_ON_ROOT( "residual " << discr_l2_residual << ", error " << discr_l2_err );

//      vtkOutput.write( maxLevel, iteration + 1 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err << ", (mesh: " << meshFile << ")" );
   WALBERLA_CHECK_LESS( discr_l2_err, targetError );
}

} // namespace hyteg

void runTests( bool usePetsc )
{
   auto dp = std::is_same< real_t, double >();
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/tri_1el.msh", 0, 3, real_c( dp ? 3e-7 : 3e-5 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/quad_4el.msh", 0, 3, real_c( dp ? 2e-6 : 5e-5 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest(
       "../../data/meshes/annulus_coarse.msh", 0, 3, real_c( dp ? 8e-7 : 5e-5 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/3D/tet_1el.msh", 0, 3, real_c( dp ? 2e-6 : 6e-5 ), true, usePetsc );
   hyteg::AgglomerationConvergenceTest(
       "../../data/meshes/3D/pyramid_2el.msh", 0, 3, real_c( dp ? 2e-5 : 5e-4 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest(
       "../../data/meshes/3D/regular_octahedron_8el.msh", 0, 3, real_c( dp ? 1e-5 : 2e-4 ), true, usePetsc );

   hyteg::AgglomerationConvergenceTest( "../../data/meshes/tri_1el.msh", 1, 3, real_c( dp ? 3e-7 : 3e-5 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/quad_4el.msh", 1, 3, real_c( dp ? 2e-6 : 5e-5 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/3D/tet_1el.msh", 1, 3, real_c( dp ? 2e-6 : 6e-5 ), true, usePetsc );
   hyteg::AgglomerationConvergenceTest(
       "../../data/meshes/3D/regular_octahedron_8el.msh", 1, 3, real_c( dp ? 1e-5 : 2e-4 ), true, usePetsc );

   hyteg::AgglomerationConvergenceTest( "../../data/meshes/tri_1el.msh", 2, 3, real_c( dp ? 3e-7 : 3e-5 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/quad_4el.msh", 2, 3, real_c( dp ? 2e-6 : 5e-5 ), false, usePetsc );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/3D/tet_1el.msh", 2, 3, real_c( dp ? 2e-6 : 6e-5 ), true, usePetsc );
   hyteg::AgglomerationConvergenceTest(
       "../../data/meshes/3D/regular_octahedron_8el.msh", 2, 3, real_c( dp ? 1e-5 : 2e-4 ), true, usePetsc );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   runTests( false );

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
   runTests( true );
#endif
}