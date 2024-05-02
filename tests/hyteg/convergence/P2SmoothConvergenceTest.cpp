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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/all.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

namespace hyteg {

static void testP2SmoothConvergence()
{
  const uint_t level = 2;

  MeshInfo mesh  = MeshInfo::fromGmshFile( "../../meshes/quad_16el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto p2Function      = std::make_shared< P2Function< real_t > >( "p2Function", storage, level, level );
  auto rhs             = std::make_shared< P2Function< real_t > >( "rhs",        storage, level, level );
  auto laplaceOperator = std::make_shared< P2ConstantLaplaceOperator >( storage, level, level );

  VTKOutput vtkOutput("../../output", "P2SmoothConvergenceTest", storage);
  vtkOutput.add( *p2Function );

  std::function< real_t( const Point3D & )> zeros = []( const Point3D & ) { return 0; };

  walberla::math::seedRandomGenerator( 0 );
  std::function< real_t( const Point3D & )> rand  = []( const Point3D & ) { return real_c( walberla::math::realRandom( 0.0, 1.0 ) ); };

  rhs->interpolate( zeros, level, All );
  p2Function->interpolate( zeros, level, DirichletBoundary );
  p2Function->interpolate( rand,  level, Inner );

  const uint_t smootherSteps = 2500;
        real_t discreteL2Norm;

  for ( uint_t step = 0; step < smootherSteps; step++ )
  {
#if 0
   // if'd out to speed up test
   discreteL2Norm = sqrt( p2Function->dotGlobal( *p2Function, level, All ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Iteration " << std::setw(10) << step << " - Discrete L2 Norm: " << std::scientific << discreteL2Norm );
   vtkOutput.write( level, step );
#endif
   laplaceOperator->smooth_gs( *p2Function, *rhs, level, Inner );
  }

  discreteL2Norm = sqrt( p2Function->dotGlobal( *p2Function, level, All ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Discrete L2 norm after " << smootherSteps << " Gauss-Seidel iterations: " << std::scientific << discreteL2Norm );

  WALBERLA_CHECK_LESS( discreteL2Norm, 7e-13 );
}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   // walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2SmoothConvergence();

   return EXIT_SUCCESS;
}
