/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

class Result
{
 public:
   static Result zero() { return Result( true ); };
   static Result notZero() { return Result( false ); };

   void check( const real_t dof ) const
   {
      if ( isZero_ )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( dof, 0.0 )
      }
      else
      {
         WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON( dof, 0.0, 10e-8 )
      }
   };

 private:
   Result( const bool isZero )
   : isZero_( isZero ){};
   bool isZero_;
};

void test( const std::function< Point3D( const Point3D& ) >& f, const Result& result )
{
   MeshInfo              meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t minLevel = 2;
   const size_t maxLevel = 4;

   walberla::math::seedRandomGenerator( 42 );

   const uint_t numRandomEvaluations = 1000;

   n1e1::N1E1VectorFunction< real_t >    x( "x", storage, minLevel, maxLevel );
   n1e1::N1E1VectorFunction< real_t >    b( "b", storage, minLevel, maxLevel );
   n1e1::N1E1ElementwiseCurlCurlOperator curlCurl( storage, minLevel, maxLevel );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      x.interpolate( f, level );
      curlCurl.apply( x, b, level, DoFType::All );

      for ( uint_t i = 0; i < numRandomEvaluations; ++i )
      {
         Point3D coordinates;
         coordinates[0] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         coordinates[1] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         coordinates[2] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );

         Point3D eval;
         auto    success = b.evaluate( coordinates, level, eval );

         WALBERLA_CHECK( success );
         result.check( eval[0] );
         result.check( eval[1] );
         result.check( eval[2] );
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test( []( const Point3D& ) { return Point3D{ 0.0, 0.0, 0.0 }; }, Result::zero() );
   WALBERLA_LOG_INFO_ON_ROOT( "Passed (1)" )
   test( []( const Point3D& p ) { return Point3D{ p[0], p[1], p[2] }; }, Result::zero() );
   WALBERLA_LOG_INFO_ON_ROOT( "Passed (2)" )
   test( []( const Point3D& p ) { return Point3D{ p[0] * p[0], p[0] * p[0], p[0] * p[0] }; }, Result::notZero() );
   WALBERLA_LOG_INFO_ON_ROOT( "Passed (3)" )

   return EXIT_SUCCESS;
}
