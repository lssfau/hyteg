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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( const uint_t & level, const std::string & meshFileName, const real_t & errEps )
{
   WALBERLA_LOG_INFO_ON_ROOT( "##### Mesh file: " << meshFileName << " / level: " << level << " #####" )

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P1PetscSolve_Domain" );

   hyteg::P1Function< real_t >                      x( "x", storage, level, level + 1 );
   hyteg::P1Function< real_t >                      x_exact( "x_exact", storage, level, level + 1 );
   hyteg::P1Function< real_t >                      b( "x", storage, level, level + 1 );
   hyteg::P1Function< real_t >                      err( "err", storage, level, level + 1 );
   hyteg::P1Function< real_t >                      residuum( "err", storage, level, level + 1 );

   hyteg::P1ConstantLaplaceOperator A( storage, level, level + 1 );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) { return sin( xx[0] ) * sinh( xx[1] ); };
   walberla::math::seedRandomGenerator( 0 );
   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) { return walberla::math::realRandom( 0.0, 1.0 ); };

   x.interpolate( exact, level, hyteg::DirichletBoundary );
   x.interpolate( rand, level, hyteg::Inner );
   b.interpolate( exact, level, hyteg::DirichletBoundary );
   x_exact.interpolate( exact, level );

//   x.interpolate( exact, level + 1, hyteg::DirichletBoundary );
//   x.interpolate( rand, level + 1, hyteg::Inner );
//   b.interpolate( exact, level + 1, hyteg::DirichletBoundary );
//   x_exact.interpolate( exact, level + 1 );

   uint_t localDoFs1 = hyteg::numberOfLocalDoFs< P1FunctionTag >( *storage, level );
//   uint_t localDoFs2 = hyteg::numberOfLocalDoFs< P1FunctionTag >( *storage, level + 1 );
   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P1FunctionTag >( *storage, level );
//   uint_t globalDoFs2 = hyteg::numberOfGlobalDoFs< P1FunctionTag >( *storage, level + 1 );

   WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );
//   WALBERLA_LOG_INFO( "localDoFs2: " << localDoFs2 << " globalDoFs2: " << globalDoFs2 );

   PETScLUSolver< hyteg::P1ConstantLaplaceOperator > solver_1( storage, level );
//   PETScLUSolver< real_t, hyteg::P1Function, hyteg::P1ConstantLaplaceOperator > solver_2( numerator, localDoFs2, globalDoFs2 );

   walberla::WcTimer timer;
   solver_1.solve( A, x, b, level );
//   solver_2.solve( A, x, b, x, level + 1, 0, 0 );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, residuum, level, hyteg::Inner );
//   A.apply( x, residuum, level + 1, hyteg::Inner );

   err.assign( {1.0, -1.0}, {x, x_exact}, level );
//   err.assign( {1.0, -1.0}, {x, x_exact}, level + 1 );

   real_t discr_l2_err_1 = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs1 );
//   real_t discr_l2_err_2 = std::sqrt( err.dotGlobal( err, level + 1 ) / (real_t) globalDoFs2 );
   real_t residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );
//   real_t residuum_l2_2  = std::sqrt( residuum.dotGlobal( residuum, level + 1 ) / (real_t) globalDoFs2 );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error 1 = " << discr_l2_err_1 );
//   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error 2 = " << discr_l2_err_2 );
//   WALBERLA_LOG_INFO_ON_ROOT( "error ratio = " << ( discr_l2_err_1 / discr_l2_err_2 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );
//   WALBERLA_LOG_INFO_ON_ROOT( "residuum 2 = " << residuum_l2_2 );

//   VTKOutput vtkOutput("../../output", "P1PetscSolve", storage);
//   vtkOutput.add( x );
//   vtkOutput.add( x_exact );
//   vtkOutput.add( err );
//   vtkOutput.add( residuum );
//   vtkOutput.write( level );

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( residuum_l2_1, 0.0, 1e-15 );
   //WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( residuum_l2_2, 0.0, 1e-15 );

   WALBERLA_CHECK_LESS( discr_l2_err_1, errEps );

   //WALBERLA_CHECK_LESS( 8.0, ( discr_l2_err_1 / discr_l2_err_2 ) );

}

}

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   petscSolveTest( 3, "../../data/meshes/quad_2el.msh",                  2.5e-05 );
   petscSolveTest( 3, "../../data/meshes/quad_4el.msh",                  1.4e-05 );
   petscSolveTest( 3, "../../data/meshes/3D/tet_1el.msh",                4.9e-07 );
   petscSolveTest( 3, "../../data/meshes/3D/pyramid_2el.msh",            4.9e-05 );
   petscSolveTest( 3, "../../data/meshes/3D/pyramid_4el.msh",            1.5e-05 );
   petscSolveTest( 3, "../../data/meshes/3D/regular_octahedron_8el.msh", 2.3e-04 );

   return EXIT_SUCCESS;
}
