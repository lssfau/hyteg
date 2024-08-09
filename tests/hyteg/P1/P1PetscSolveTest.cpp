/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( int solverType, const uint_t& level, const std::string& meshFileName, const real_t& errEps )
{
   WALBERLA_LOG_INFO_ON_ROOT( "##### Mesh file: " << meshFileName << " / level: " << level << " #####" )

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P1PetscSolve_Domain" );

   hyteg::P1Function< real_t > x( "x", storage, level, level );
   hyteg::P1Function< real_t > x_exact( "x_exact", storage, level, level );
   hyteg::P1Function< real_t > b( "x", storage, level, level );
   hyteg::P1Function< real_t > err( "err", storage, level, level );
   hyteg::P1Function< real_t > residuum( "err", storage, level, level );

   hyteg::P1ConstantLaplaceOperator A( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return sin( xx[0] ) * sinh( xx[1] );
   };
   walberla::math::seedRandomGenerator( 0 );
   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) {
      return real_c( walberla::math::realRandom( 0.0, 1.0 ) );
   };

   x.interpolate( exact, level, hyteg::DirichletBoundary );
   x.interpolate( rand, level, hyteg::Inner );
   b.interpolate( exact, level, hyteg::DirichletBoundary );
   x_exact.interpolate( exact, level );

   const uint_t globalDoFs = hyteg::numberOfGlobalDoFs< P1FunctionTag >( *storage, level );

   std::shared_ptr< Solver< hyteg::P1ConstantLaplaceOperator > > solver;

   if ( solverType == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solving with LU (MUMPS)" );
      auto actualSolver = std::make_shared< PETScLUSolver< hyteg::P1ConstantLaplaceOperator > >( storage, level );
      solver            = actualSolver;
   }
   else if ( solverType == 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solving with CG" );
      auto actualSolver = std::make_shared< PETScCGSolver< hyteg::P1ConstantLaplaceOperator > >( storage, level );
      solver            = actualSolver;
   }
   else if ( solverType == 2 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solving with MINRES" );
      auto actualSolver = std::make_shared< PETScMinResSolver< hyteg::P1ConstantLaplaceOperator > >( storage, level );
      solver            = actualSolver;
   }

   walberla::WcTimer timer;
   solver->solve( A, x, b, level );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, residuum, level, hyteg::Inner );

   err.assign( { 1.0, -1.0 }, { x, x_exact }, level );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs );
   real_t residuum_l2  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum = " << residuum_l2 );

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( residuum_l2, 0.0, 1e-12 );

   WALBERLA_CHECK_LESS( discr_l2_err, errEps );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   for ( int solver = 1; solver < 3; solver++ )
   {
      petscSolveTest( solver, 3, prependHyTeGMeshDir( "2D/quad_2el.msh" ), 2.5e-05 );
      petscSolveTest( solver, 3, prependHyTeGMeshDir( "2D/quad_4el.msh" ), 1.4e-05 );
      petscSolveTest( solver, 3, prependHyTeGMeshDir( "3D/tet_1el.msh" ), 4.9e-07 );
      petscSolveTest( solver, 3, prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 4.9e-05 );
      petscSolveTest( solver, 3, prependHyTeGMeshDir( "3D/pyramid_4el.msh" ), 1.5e-05 );
      petscSolveTest( solver, 3, prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 2.3e-04 );
   }

   return EXIT_SUCCESS;
}
