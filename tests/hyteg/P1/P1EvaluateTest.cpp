/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr.
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
#include <hyteg/communication/Syncing.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroEdge.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::math::pi;

using namespace hyteg;

void test2D()
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t minLevel = 2;
   const size_t maxLevel = 5;

   walberla::math::seedRandomGenerator( 12345678 );

   const uint_t numRandomEvaluations = 1000;

   auto testFunc            = []( const Point3D& x ) { return real_c( 10.0 ) * x[0] + real_c( 3.0 ) * x[1] + real_c( 1.0 ); };
   auto testFuncDerivativeX = []( const Point3D& ) { return real_c( 10.0 ); };
   auto testFuncDerivativeY = []( const Point3D& ) { return real_c( 3.0 ); };

   P1Function< real_t > x( "x", storage, minLevel, maxLevel );
   x.interpolate( testFunc, maxLevel );
   // IMPORTANT
   communication::syncFunctionBetweenPrimitives( x, maxLevel );

   Point3D coordinates(  0.0, 0.5, 0.0  );
   real_t  eval;
   auto    success = x.evaluate( coordinates, maxLevel, eval );
   WALBERLA_CHECK( success );
   Point3D gradient( Point3D::Zero() );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 2.0;
   coordinates[1] = 0.0;
   success        = x.evaluate( coordinates, maxLevel, eval );
   WALBERLA_CHECK( success );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 2.0;
   coordinates[1] = 1.0;
   success        = x.evaluate( coordinates, maxLevel, eval );
   WALBERLA_CHECK( success );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 0.0;
   coordinates[1] = 1.0;
   success        = x.evaluate( coordinates, maxLevel, eval );
   WALBERLA_CHECK( success );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 0.5;
   coordinates[1] = 0.5;
   success        = x.evaluate( coordinates, maxLevel, eval );
   WALBERLA_CHECK( success );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      coordinates[0] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
      coordinates[1] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );

      if ( coordinates[0] < 0.5 && coordinates[1] < 0.5 )
      {
         continue;
      }

      success = x.evaluate( coordinates, maxLevel, eval );
      WALBERLA_CHECK( success );
      x.evaluateGradient( coordinates, maxLevel, gradient );
      WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
      WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
      WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );
   }
}

void test3D()
{
   MeshInfo              meshInfo = MeshInfo::meshSymmetricCuboid( Point3D(  0, 0, 0  ), Point3D(  1, 1, 1  ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t minLevel = 2;
   const size_t maxLevel = 5;

   walberla::math::seedRandomGenerator( 12345678 );

   const uint_t numRandomEvaluations = 1000;

   auto testFunc = []( const Point3D& x ) { return 10.0 * x[0] + 3.0 * x[1] + x[2] + 1.0; };

   P1Function< real_t > x( "x", storage, minLevel, maxLevel );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      x.interpolate( testFunc, level );
      communication::syncFunctionBetweenPrimitives( x, level );

      for ( uint_t i = 0; i < numRandomEvaluations; ++i )
      {
         Point3D coordinates( Point3D::Zero() );
         coordinates[0] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         coordinates[1] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         coordinates[2] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );

         real_t eval;
         auto   success = x.evaluate( coordinates, level, eval );
         WALBERLA_CHECK( success );
         WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ), "Test3D: wrong value at " << coordinates << "." );

         WALBERLA_LOG_DETAIL_ON_ROOT( "Passed: " << coordinates );
      }
   }
}

void testEvaluateWithBlending( uint_t numSamples, uint_t mapType )
{
   std::vector< Point3D > samples;
   samples.reserve( numSamples );

   std::shared_ptr< PrimitiveStorage > storage;
   walberla::math::seedRandomGenerator( 12345678 );

   real_t tolerance = real_c( -1 );

   // annulus map
   if ( mapType == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Testing with AnnulusMap" );
      real_t rMin = real_c( 1 );
      real_t rMax = real_c( 2 );

      MeshInfo              meshInfo = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CRISS, 12, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      AnnulusMap::setMap( setupStorage );

      for ( uint_t k = 0; k < numSamples; ++k )
      {
         real_t radius = walberla::math::realRandom( rMin, rMax );
         real_t angle  = walberla::math::realRandom( real_c( 0 ), real_c( 2 ) * pi );
         samples.push_back( Point3D(  radius * std::cos( angle ), radius * std::sin( angle ), real_c( 0 )  ) );
      }

      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      if constexpr ( std::is_same_v< real_t, double > )
      {
         tolerance = real_c( 1e-13 );
      }
      else
      {
         tolerance = real_c( 5e-7 );
      }
   }

   // 2D affine map
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Testing with AffineMap2D" );
      MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      Point2D  shift(  real_c( 4 ), real_c( -1.0 / 3.0 )  );
      real_t   angle = pi / real_c( 180.0 * 25.0 );
      Matrix2r mat;
      mat( 0, 0 ) = +std::cos( angle );
      mat( 0, 1 ) = -std::sin( angle );
      mat( 1, 0 ) = +std::sin( angle );
      mat( 1, 1 ) = +std::cos( angle );

      AffineMap2D::setMap( setupStorage, mat, shift );
      for ( uint_t k = 0; k < numSamples; ++k )
      {
         real_t xPos = walberla::math::realRandom( real_c( 0 ), real_c( 1 ) );
         real_t yPos = walberla::math::realRandom( real_c( 0 ), real_c( 1 ) );
         real_t x    = mat( 0, 0 ) * xPos + mat( 0, 1 ) * yPos + shift[0];
         real_t y    = mat( 1, 0 ) * xPos + mat( 1, 1 ) * yPos + shift[1];
         samples.emplace_back( x, y, real_c( 0 ) );
      }

      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      if constexpr ( std::is_same_v< real_t, double > )
      {
         tolerance = real_c( 1e-13 );
      }
      else
      {
         tolerance = real_c( 5e-6 );
      }
   }

   const size_t minLevel = 2;
   const size_t maxLevel = 4;

   auto testFunc = [mapType]( const Point3D& x ) {
      // use 1st order polynomial for affine map
      if ( mapType != 0 )
      {
         return real_c( 1 ) + real_c( 5 ) * x[0] + real_c( 3 ) * x[1];
      }

      // for annulus map
      else
      {
         return std::sqrt( x[0] * x[0] + x[1] * x[1] );
      }
   };

   P1Function< real_t > func( "func", storage, minLevel, maxLevel );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      func.interpolate( testFunc, level );
      real_t testValue = real_c( 0 );
      real_t ctrlValue = real_c( 0 );
      communication::syncFunctionBetweenPrimitives( func, level );
      for ( uint_t idx = 0; idx < numSamples; ++idx )
      {
         ctrlValue          = testFunc( samples[idx] );
         bool   coordsFound = func.evaluate( samples[idx], level, testValue, real_c( tolerance ) );
         real_t error       = std::abs( testValue - ctrlValue );
         WALBERLA_LOG_INFO_ON_ROOT( "sampling Point: " << std::scientific << std::showpos << samples[idx] << ", coordsFound = "
                                                       << ( coordsFound ? "TRUE" : "FALSE" ) << ", ctrl = " << ctrlValue
                                                       << ", test = " << testValue << ", error = " << error );
         WALBERLA_CHECK( coordsFound );
         WALBERLA_CHECK_LESS( error, tolerance );
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test2D();
   test3D();
   testEvaluateWithBlending( 10, 1 );
   testEvaluateWithBlending( 10, 0 );

   return EXIT_SUCCESS;
}
