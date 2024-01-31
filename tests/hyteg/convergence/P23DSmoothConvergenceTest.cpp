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
#include "core/OpenMP.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/all.h"

#include "hyteg/OpenMPManager.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"

namespace hyteg {

static void testP2SmoothConvergence( const uint_t & level, const std::string & meshFile, const uint_t & numIterations, const real_t & expectedL2Error )
{
  MeshInfo mesh  = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  P2Function< real_t > x  ( "x",   storage, level, level );
  P2Function< real_t > rhs( "rhs", storage, level, level );
  P2ConstantLaplaceOperator A( storage, level, level );

  VTKOutput vtkOutput("../../output", "P23DSmoothConvergenceTest", storage);
  vtkOutput.add( x );

  std::function< real_t( const Point3D & )> zeros = []( const Point3D & ) { return 0; };
  walberla::math::seedRandomGenerator( 0 );
  std::function< real_t( const Point3D & )> rand  = []( const Point3D & ) { return real_c( walberla::math::realRandom( 0.0, 1.0 ) ); };

  rhs.interpolate( zeros, level, All );
  x.interpolate( zeros, level, DirichletBoundary );

  hyteg::OpenMPManager::instance()->forceSerial();
  x.interpolate( rand,  level, Inner );
  hyteg::OpenMPManager::instance()->resetToParallel();

  real_t discreteL2Norm;

  for ( uint_t step = 0; step < numIterations; step++ )
  {
#if 0
   // if'd out to speed up test
   discreteL2Norm = sqrt( x.dotGlobal( x, level, All ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Iteration " << std::setw(10) << step << " - Discrete L2 Norm: " << std::scientific << discreteL2Norm );
   vtkOutput.write( level, step );
#endif
   A.smooth_gs( x, rhs, level, Inner );
  }

  discreteL2Norm = sqrt( x.dotGlobal( x, level, All ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Discrete L2 norm after " << numIterations << " Gauss-Seidel iterations: " << std::scientific << discreteL2Norm );
  WALBERLA_CHECK_LESS( discreteL2Norm, expectedL2Error );
}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2SmoothConvergence( 3, "../../data/meshes/3D/tet_1el.msh", 50, 1.2e-02 );
   hyteg::testP2SmoothConvergence( 2, "../../data/meshes/3D/pyramid_2el.msh", 50, 9.3e-07 );
   hyteg::testP2SmoothConvergence( 2, "../../data/meshes/3D/pyramid_4el.msh", 50, 1.6e-03 );
   hyteg::testP2SmoothConvergence( 2, "../../data/meshes/3D/regular_octahedron_8el.msh", 50, 7.7e-02 );

   return EXIT_SUCCESS;
}
