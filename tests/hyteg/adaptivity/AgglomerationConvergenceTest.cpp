/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/controlflow/AgglomerationWrapper.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void AgglomerationConvergenceTest( const std::string& meshFile,
                                   const uint_t       level,
                                   const real_t       targetError,
                                   const bool         localMPI )
{
   WALBERLA_CHECK_GREATER( level, 2 );

   const uint_t numIterations = 5;

   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   const auto numberOfProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hyteg::P2ConstantLaplaceOperator L( storage, 2, level );

   hyteg::P2Function< real_t > r( "r", storage, 2, level );
   hyteg::P2Function< real_t > f( "f", storage, 2, level );
   hyteg::P2Function< real_t > u( "u", storage, 2, level );
   hyteg::P2Function< real_t > u_exact( "u_exact", storage, 2, level );
   hyteg::P2Function< real_t > err( "err", storage, 2, level );
   hyteg::P2Function< real_t > npoints_helper( "npoints_helper", storage, 2, level );

   if ( localMPI )
   {
      u.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
   }

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   hyteg::VTKOutput vtkOutput( "../../output", "AgglomerationConvergenceTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( f );
   vtkOutput.add( r );
   vtkOutput.add( err );
   vtkOutput.add( npoints_helper );
   vtkOutput.write( level, 0 );

   // Setup of the agglomeration based solver.
   // Apart from the coarse grid, everything is performed in parallel.
   // The coarse grid problem is, however, solved on a subset of processes.
   const uint_t numberOfSubsetProcesses = ( numberOfProcesses + 1 ) / 2;
   WALBERLA_LOG_INFO_ON_ROOT( "Agglomeration from " << numberOfProcesses << " to " << numberOfSubsetProcesses << " processes." )
   auto agglomerationStorage = storage->createCopy();

   auto smoother     = std::make_shared< GaussSeidelSmoother< P2ConstantLaplaceOperator > >();
   auto prolongation = std::make_shared< P2toP2QuadraticProlongation >();
   auto restriction  = std::make_shared< P2toP2QuadraticRestriction >();
#if 1
   // pass agglomeration storage to coarse grid solver
   auto coarseGridSolver = std::make_shared< CGSolver< P2ConstantLaplaceOperator > >( agglomerationStorage, 2, 2 );
   // now wrap the solver
   auto coarseGridSolverAgglomeration =
       std::make_shared< AgglomerationWrapper< P2ConstantLaplaceOperator > >( coarseGridSolver, agglomerationStorage, 2, numberOfSubsetProcesses );
#else
  auto coarseGridSolverAgglomeration = std::make_shared< CGSolver< P2ConstantLaplaceOperator > >( storage, 2, 2 );
#endif

   auto solver = std::make_shared< GeometricMultigridSolver< P2ConstantLaplaceOperator > >(
       storage, smoother, coarseGridSolverAgglomeration, restriction, prolongation, 2, level );

   real_t discr_l2_err, discr_l2_residual;
   for ( uint_t iteration = 0; iteration < numIterations; iteration++ )
   {
      solver->solve( L, u, f, level );

      err.assign( {1.0, -1.0}, {u, u_exact}, level );
      npoints_helper.interpolate( ones, level );

      const real_t npoints = npoints_helper.dotGlobal( npoints_helper, level );
      discr_l2_err         = std::sqrt( err.dotGlobal( err, level ) / npoints );

      L.apply( u, err, level, Inner | NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, err}, level, Inner | NeumannBoundary );
      discr_l2_residual = std::sqrt( r.dotGlobal( r, level, Inner | NeumannBoundary ) / npoints );

      WALBERLA_LOG_INFO_ON_ROOT( "residual " << discr_l2_residual << ", error " << discr_l2_err );

      vtkOutput.write( level, iteration + 1 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err << ", (mesh: " << meshFile << ")" );
   WALBERLA_CHECK_LESS( discr_l2_err, targetError );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::AgglomerationConvergenceTest( "../../data/meshes/tri_1el.msh", 3, 1e-7, false );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/quad_4el.msh", 3, 1e-7, false );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/annulus_coarse.msh", 3, 3e-7, false );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/3D/tet_1el.msh", 3, 4e-7, true );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/3D/pyramid_2el.msh", 3, 3e-6, false );
   hyteg::AgglomerationConvergenceTest( "../../data/meshes/3D/regular_octahedron_8el.msh", 3, 1.8e-6, true );
}