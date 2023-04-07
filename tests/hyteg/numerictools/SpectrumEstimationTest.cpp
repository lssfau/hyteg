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


/// \file
/// This file contains tests for the spectrum estimation functions
/// - estimateSpectralRadiusWithPowerIteration()
/// - estimateSpectralBoundsWithCG()
/// The tests try to estimate eigenvalues for the Finite Element discretisation
/// of the Laplacian with P1 elements on a cross mesh. Besides scaling with the
/// mesh width we get the standard 5-point-stencil from Finite Differences. The
/// resulting eigenvalues are known analytically.

#include "hyteg/HytegDefinitions.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "hyteg/OpenMPManager.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/SpectrumEstimation.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

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
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ),
                                               Point2D( 1.0, 1.0 ),
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

  P1Function< real_t > uRand1( "uRand", storage, level, level );
  P1Function< real_t > uRand2( "uRand", storage, level, level );
  P1Function< real_t > rhs( "rhs", storage, level, level );
  P1Function< real_t > tmp( "tmp", storage, level, level );

  OpenMPManager::instance()->forceSerial();
  uRand1.interpolate( randFunc, level, Inner );
  uRand2.interpolate( randFunc, level, Inner );
  OpenMPManager::instance()->resetToParallel();

  uRand1.interpolate( zeros, level, DirichletBoundary );
  uRand2.interpolate( zeros, level, DirichletBoundary );
  rhs.interpolate( zeros, level, All );

  // analytic spectrum
  uint_t N = 1u << level;
  WALBERLA_LOG_INFO_ON_ROOT( "Analytic spectral bounds:" );

  real_t argOne = real_c( 0.5 ) * pi / static_cast< real_t >( N );
  real_t argTwo = real_c( 0.5 ) * pi * ( static_cast< real_t >( N - 1 ) / static_cast< real_t >( N ) );

  real_t lambdaMin = real_c( 8.0 ) * sin( argOne ) * sin( argOne );
  real_t lambdaMax = real_c( 8.0 ) * sin( argTwo ) * sin( argTwo );
  WALBERLA_LOG_INFO_ON_ROOT( "lambdaMin = " << std::scientific << lambdaMin );
  WALBERLA_LOG_INFO_ON_ROOT( "lambdaMax = " << std::scientific << lambdaMax );

  // -------
  //  TESTS
  // -------
  uint_t numIts = 50;

  // Power Iteration
  real_t radius = estimateSpectralRadiusWithPowerIteration( op, uRand1, tmp, numIts, storage, level );
  WALBERLA_LOG_INFO_ON_ROOT( "Estimation for spectral radius with Power Iteration = " << std::scientific << radius );
  WALBERLA_CHECK_LESS( std::abs( 7.9 - radius ), 5e-2 );
  auto dp = std::is_same< real_t, double >();
  WALBERLA_CHECK_FLOAT_EQUAL( radius, dp ? 7.92126990394241570e+00 : 7.91265965e+00 );

  // CG/Lanczos approach
  real_t lowerBound, upperBound;
  estimateSpectralBoundsWithCG( op, uRand2, rhs, numIts, storage, level, lowerBound, upperBound );
  WALBERLA_LOG_INFO_ON_ROOT( "Estimation for spectrum with CG = [" << std::scientific << lowerBound << ", " << upperBound
                                                                   << "]" );
  WALBERLA_CHECK_LESS( std::abs( lambdaMax - upperBound ), 1e-2 );
  WALBERLA_CHECK_FLOAT_EQUAL( upperBound, dp ? 7.98704774334941625e+00 : 7.99281073e+00 );
  WALBERLA_CHECK_LESS( std::abs( lambdaMin - lowerBound ), 1e-2 );
  WALBERLA_CHECK_FLOAT_EQUAL( lowerBound, 5.756730e-03 );

  // If we got up to here, everything's fine
  return EXIT_SUCCESS;
}
