/*
 * Copyright (c) 2025 Andreas Burkhart.
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
#include "core/mpi/MPIManager.h"

#include "hyteg/geometry/ClosestPoint.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using namespace hyteg;

using tetType = std::tuple< hyteg::Point3D, hyteg::Point3D, hyteg::Point3D, hyteg::Point3D, real_t, hyteg::Point3D >;

Point3D
    closestPointBruteForce( const Point3D& P, const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D, uint_t n )
{
   Point3D closestPoint = A;
   real_t  minDist      = ( P - A ).norm();

   real_t  c1;
   real_t  c2;
   real_t  c3;
   real_t  c4;
   real_t  dist;
   Point3D DP;

   for ( uint_t i = 0; i < n; i++ )
   {
      c1 = real_c( i ) / real_c( n - 1 );
      for ( uint_t k = 0; k < n - i; k++ )
      {
         c2 = real_c( k ) / real_c( n - 1 );
         for ( uint_t p = 0; p < n - i - k; p++ )
         {
            c3 = real_c( p ) / real_c( n - 1 );
            c4 = real_c( 1 ) - c3 - c2 - c1;

            DP   = c1 * A + c2 * B + c3 * C + c4 * D;
            dist = ( DP - P ).norm();

            if ( dist < minDist )
            {
               closestPoint = DP;
               minDist      = dist;
            }
         }
      }
   }

   return closestPoint;
}

void refineTetrahedron( tetType old, hyteg::Point3D P, uint_t n, std::vector< tetType >& out )
{
   std::array< hyteg::Point3D, 6 > newPoints = { real_c( 0.5 ) * ( std::get< 0 >( old ) + std::get< 1 >( old ) ),
                                                 real_c( 0.5 ) * ( std::get< 0 >( old ) + std::get< 2 >( old ) ),
                                                 real_c( 0.5 ) * ( std::get< 1 >( old ) + std::get< 2 >( old ) ),
                                                 real_c( 0.5 ) * ( std::get< 0 >( old ) + std::get< 3 >( old ) ),
                                                 real_c( 0.5 ) * ( std::get< 1 >( old ) + std::get< 3 >( old ) ),
                                                 real_c( 0.5 ) * ( std::get< 2 >( old ) + std::get< 3 >( old ) ) };

   std::array< std::array< hyteg::Point3D, 4 >, 8 > newTets = {
       std::array< hyteg::Point3D, 4 >( { std::get< 0 >( old ), newPoints[0], newPoints[1], newPoints[3] } ),
       std::array< hyteg::Point3D, 4 >( { newPoints[0], std::get< 1 >( old ), newPoints[2], newPoints[4] } ),
       std::array< hyteg::Point3D, 4 >( { newPoints[1], newPoints[2], std::get< 2 >( old ), newPoints[5] } ),
       std::array< hyteg::Point3D, 4 >( { newPoints[3], newPoints[4], newPoints[5], std::get< 3 >( old ) } ),
       std::array< hyteg::Point3D, 4 >( { newPoints[0], newPoints[1], newPoints[2], newPoints[4] } ),
       std::array< hyteg::Point3D, 4 >( { newPoints[0], newPoints[1], newPoints[3], newPoints[4] } ),
       std::array< hyteg::Point3D, 4 >( { newPoints[1], newPoints[2], newPoints[4], newPoints[5] } ),
       std::array< hyteg::Point3D, 4 >( { newPoints[1], newPoints[3], newPoints[4], newPoints[5] } ) };

   for ( auto& a : newTets )
   {
      hyteg::Point3D DP = closestPointBruteForce( P, a[0], a[1], a[2], a[3], n );
      out.push_back( tetType( a[0], a[1], a[2], a[3], ( DP - P ).norm(), DP ) );
   }
}

bool sortByDist( const tetType& a, const tetType& b )
{
   return ( std::get< 4 >( a ) < std::get< 4 >( b ) );
}

Point3D closestPointSubdivision( const Point3D& P,
                                 const Point3D& A,
                                 const Point3D& B,
                                 const Point3D& C,
                                 const Point3D& D,
                                 uint_t         n,
                                 uint_t         it,
                                 uint_t         keepN = 32 )
{
   std::vector< tetType > tetList;

   hyteg::Point3D DP = closestPointBruteForce( P, A, B, C, D, n );
   tetList.push_back( tetType( A, B, C, D, ( DP - P ).norm(), DP ) );

   for ( uint_t k = 0; k < it; k++ )
   {
      std::vector< tetType > tetListNew;

      for ( auto& t : tetList )
      {
         refineTetrahedron( t, P, n, tetListNew );
      }

      std::sort( tetListNew.begin(), tetListNew.end(), sortByDist );

      tetList.clear();

      // we need to keep more tets if n is small
      for ( uint_t i = 0; i < std::min( uint_c( keepN ), uint_c( tetListNew.size() ) ); i++ )
      {
         tetList.push_back( tetListNew[i] );
      }
   }

   return std::get< 5 >( tetList[0] );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   if ( sizeof( real_t ) < 8 )
   {
      // The tolerances in this test are designed for at least double precision.
      WALBERLA_LOG_INFO_ON_ROOT( "Single precision or lower detected. Aborting test." );

      return EXIT_SUCCESS;
   }

   // Define a tetrahedron
   const Point3D A( 1, 1, 1 );
   const Point3D B( 3, 1, 1 );
   const Point3D C( 1, 3, 1 );
   const Point3D D( 1, 1, 3 );

   // Seed random generator for reproducable results
   walberla::math::seedRandomGenerator( 123456 );

   // Now check random points and compare the result to a brute force subdivision approach.
   Point3D P;
   Point3D closestPoint;
   Point3D closestPointCompare;
   real_t  compDist;
   real_t  maxDist = std::numeric_limits< real_t >::min();

   for ( uint_t k = 0; k < 5000; k++ )
   {
      P = { walberla::math::realRandom( real_c( 0 ), real_c( 4 ) ),
            walberla::math::realRandom( real_c( 0 ), real_c( 4 ) ),
            walberla::math::realRandom( real_c( 0 ), real_c( 4 ) ) };

      // inherently also tests closestPointTriangle3D, since it is used in closestPointTetrahedron3D
      closestPoint        = closestPointTetrahedron3D( P, A, B, C, D );
      closestPointCompare = closestPointSubdivision( P, A, B, C, D, 5, 24 );

      compDist = ( closestPoint - closestPointCompare ).norm();

      if ( compDist > maxDist )
      {
         maxDist = compDist;
      }
   }

   WALBERLA_LOG_INFO( "Maximum Distance: " << maxDist );
   WALBERLA_CHECK( maxDist < real_c( 1e-7 ) )

   return EXIT_SUCCESS;
}