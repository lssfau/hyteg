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
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( const uint_t& level, const std::string& meshFileName, const real_t& errEps )
{
   WALBERLA_LOG_INFO_ON_ROOT( "##### Mesh file: " << meshFileName << " / level: " << level << " #####" )

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hyteg::P2Function< real_t > x( "x", storage, level, level );
   hyteg::P2Function< real_t > x_exact( "x_exact", storage, level, level );
   hyteg::P2Function< real_t > b( "b", storage, level, level );
   hyteg::P2Function< real_t > err( "err", storage, level, level );
   hyteg::P2Function< real_t > residuum( "residuum", storage, level, level );

   hyteg::P2ConstantLaplaceOperator A( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return sin( xx[0] ) * sinh( xx[1] );
   };
   walberla::math::seedRandomGenerator( 0 );
   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) { return real_c( walberla::math::realRandom( 0.0, 1.0 ) ); };

   x.interpolate( exact, level, hyteg::DirichletBoundary );
   x.interpolate( rand, level, hyteg::Inner );
   b.interpolate( exact, level, hyteg::DirichletBoundary );
   x_exact.interpolate( exact, level );

   uint_t localDoFs1  = hyteg::numberOfLocalDoFs< P2FunctionTag >( *storage, level );
   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );

   PETScLUSolver< hyteg::P2ConstantLaplaceOperator > solver_1( storage, level );

   walberla::WcTimer timer;
   solver_1.solve( A, x, b, level );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, residuum, level, hyteg::Inner );

   err.assign( {1.0, -1.0}, {x, x_exact}, level );

   real_t discr_l2_err_1 = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs1 );
   real_t residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error 1 = " << discr_l2_err_1 );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

   VTKOutput vtkOutput( "../../output", "P2PetscSolve", storage );
   vtkOutput.add( x );
   vtkOutput.add( x_exact );
   vtkOutput.add( err );
   vtkOutput.add( residuum );
   vtkOutput.write( level );

   WALBERLA_CHECK_LESS( residuum_l2_1, 4e-15 );
   WALBERLA_CHECK_LESS( discr_l2_err_1, errEps );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   petscSolveTest( 0, "../../data/meshes/quad_4el.msh", 3.0e-04 );
   petscSolveTest( 0, "../../data/meshes/3D/tet_1el.msh", 1.0e-15 );
   petscSolveTest( 0, "../../data/meshes/3D/pyramid_2el.msh", 1.0e-15 );
   petscSolveTest( 0, "../../data/meshes/3D/pyramid_4el.msh", 2.0e-04 );
   petscSolveTest( 0, "../../data/meshes/3D/pyramid_tilted_4el.msh", 2.0e-04 );
   petscSolveTest( 0, "../../data/meshes/3D/regular_octahedron_8el.msh", 3.0e-04 );
   petscSolveTest( 0, "../../data/meshes/3D/cube_24el.msh", 5.0e-04 );

   petscSolveTest( 1, "../../data/meshes/quad_4el.msh", 2.0e-05 );
   petscSolveTest( 1, "../../data/meshes/3D/tet_1el.msh", 3.0e-06 );
   petscSolveTest( 1, "../../data/meshes/3D/pyramid_2el.msh", 2.0e-04 );
   petscSolveTest( 1, "../../data/meshes/3D/pyramid_4el.msh", 2.0e-05 );
   petscSolveTest( 1, "../../data/meshes/3D/pyramid_tilted_4el.msh", 5.0e-05 );
   petscSolveTest( 1, "../../data/meshes/3D/regular_octahedron_8el.msh", 3.0e-04 );
   petscSolveTest( 1, "../../data/meshes/3D/cube_24el.msh", 5.0e-04 );

   petscSolveTest( 3, "../../data/meshes/quad_4el.msh", 3.0e-07 );
   petscSolveTest( 3, "../../data/meshes/3D/tet_1el.msh", 3.0e-07 );
   petscSolveTest( 3, "../../data/meshes/3D/pyramid_2el.msh", 2.7e-06 );
   petscSolveTest( 3, "../../data/meshes/3D/pyramid_4el.msh", 3.2e-07 );
   petscSolveTest( 2, "../../data/meshes/3D/pyramid_tilted_4el.msh", 7.3e-06 );
   petscSolveTest( 3, "../../data/meshes/3D/regular_octahedron_8el.msh", 1.7e-06 );
   petscSolveTest( 2, "../../data/meshes/3D/cube_24el.msh", 1.4e-05 );

   return EXIT_SUCCESS;
}
