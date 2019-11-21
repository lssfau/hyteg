/*
 * Copyright (c) 2017-2019 Marcus Mohr.
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

#include "hyteg/HytegDefinitions.hpp"

#ifndef HYTEG_BUILD_WITH_EIGEN
#error "This app only works with Eigen enabled. Please enable it via -DHYTEG_BUILD_WITH_EIGEN=ON"
#endif

#include "core/math/Random.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/numerictools/SpectrumEstimation.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

int main( int argc, char* argv[] )
{

  // -------
  //  SETUP
  // -------

  // MPI
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  // mesh, primitives, ...
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0.0, 0.0} ),
                                               Point2D( {1.0, 1.0} ),
                                               MeshInfo::CROSS, 1, 1 );
  SetupPrimitiveStorage setupStorage = SetupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::greedy( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // operator, functions, ...
  uint_t level = 6;
  P1ConstantLaplaceOperator op( storage, level, level );

  std::function<real_t(const Point3D&)> randFunc = []( const Point3D & ) { return walberla::math::realRandom<real_t>(); };
  std::function<real_t(const Point3D&)> zeros = [](const Point3D&) { return 0.0; };

  // commenting in changes results between different runs
  // std::random_device rdev{};
  // walberla::math::seedRandomGenerator( rdev() );

  P1Function< real_t > uRand( "uRand", storage, level, level );
  P1Function< real_t > rhs( "rhs", storage, level, level );

  uRand.interpolate( randFunc, level, Inner );
  uRand.interpolate( zeros, level, DirichletBoundary );
  rhs.interpolate( zeros, level, All );

  // -------
  //  TESTS
  // -------
  uint_t numIts = 50;

  // CG/Lanczos approach
  real_t lowerBound, upperBound;
  estimateSpectralBoundsWithCG( op, uRand2, rhs, numIts, storage, level, lowerBound, upperBound );
  WALBERLA_LOG_INFO_ON_ROOT( "Estimation for spectrum with CG = [" << std::scientific << lowerBound << ", " << upperBound << "]" );
  WALBERLA_CHECK_LESS( std::abs( lambdaMax - upperBound ), 1e-2 );
  WALBERLA_CHECK_FLOAT_EQUAL( upperBound, 7.98704774334941625e+00 );
  WALBERLA_CHECK_LESS( std::abs( lambdaMin - lowerBound ), 1e-2 );
  WALBERLA_CHECK_FLOAT_EQUAL( lowerBound, 5.756730e-03 );

  // If we got up to here, everything's fine
  return EXIT_SUCCESS;

}
