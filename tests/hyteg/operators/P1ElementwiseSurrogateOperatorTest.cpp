/*
 * Copyright (c) 2025 Benjamin Mann.
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

#include <core/DataTypes.h>
#include <core/mpi/MPIManager.h>
#include <hyteg/elementwiseoperators/P1ElementwiseOperator.hpp>
#include <hyteg/elementwiseoperators/P1ElementwiseSurrogateOperator.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp>

/* This test checks whether the P1 elementwise surrogate
   operator works correctly. In particular
   we verify the following theorem:
   Let A be the discrete operator associated with
      -div(k grad(u)) = f,
   where k is a polynomial of degree p.
   Then, for the surrogate operator A_q, defined
   by polynomials of degree q, it holds
      q>=p => A_q = A.
*/
using walberla::real_t;
using namespace hyteg;

void P1SurrogateOperatorTest( const std::shared_ptr< PrimitiveStorage >&        storage,
                              std::function< real_t( const hyteg::Point3D& ) >& k,
                              const uint8_t                                     q,
                              const uint_t                                      level )
{
   const real_t epsilon = real_c( std::is_same< real_t, double >() ? 2e-12 : 5e-4 );

   // operators
   forms::p1_div_k_grad_affine_q3               form( k, k );
   P1ElementwiseAffineDivKGradOperator          A( storage, level, level, form );
   P1ElementwiseSurrogateAffineDivKGradOperator A_q( storage, level, level, form );

   A_q.init( q, 1, "", false );

   // functions
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > Au( "Au", storage, level, level );
   hyteg::P1Function< real_t > Aqu( "(A_q)u", storage, level, level );
   hyteg::P1Function< real_t > err( "(A-A_q)u", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > initialU = []( const hyteg::Point3D& x ) {
      return cos( 2 * M_PI * x[0] ) * cos( 2 * M_PI * x[1] ) * cos( 2 * M_PI * x[2] );
   };
   u.interpolate( initialU, level );

   // apply operators
   A.apply( u, Au, level, All, Replace );
   A_q.apply( u, Aqu, level, All, Replace );

   // compute error
   err.assign( { 1.0, -1.0 }, { Au, Aqu }, level, All );
   auto errorMax = err.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "error ||(A - A_q)u||_inf = %e", errorMax ) )
   WALBERLA_CHECK_LESS( errorMax, epsilon, "||(A - A_q)u||_inf" );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // coefficients
   std::map< uint_t, std::function< real_t( const hyteg::Point3D& x ) > > k;

   k[1] = []( const hyteg::Point3D& x ) { return 2 * x[0] + 3 * x[1] + 4 * x[2] + 1; };
   k[2] = []( const hyteg::Point3D& x ) {
      return 5 * pow( x[0], 2 ) + 6 * x[0] * x[1] + 7 * x[0] * x[2] + 2 * x[0] + 8 * pow( x[1], 2 ) + 9 * x[1] * x[2] + 3 * x[1] +
             10 * pow( x[2], 2 ) + 4 * x[2] + 1;
   };
   k[3] = []( const hyteg::Point3D& x ) {
      return 11 * pow( x[0], 3 ) + 12 * pow( x[0], 2 ) * x[1] + 13 * pow( x[0], 2 ) * x[2] + 5 * pow( x[0], 2 ) +
             14 * x[0] * pow( x[1], 2 ) + 15 * x[0] * x[1] * x[2] + 6 * x[0] * x[1] + 16 * x[0] * pow( x[2], 2 ) +
             7 * x[0] * x[2] + 2 * x[0] + 17 * pow( x[1], 3 ) + 18 * pow( x[1], 2 ) * x[2] + 8 * pow( x[1], 2 ) +
             19 * x[1] * pow( x[2], 2 ) + 9 * x[1] * x[2] + 3 * x[1] + 20 * pow( x[2], 3 ) + 10 * pow( x[2], 2 ) + 4 * x[2] + 1;
   };
   k[4] = []( const hyteg::Point3D& x ) {
      return 21 * pow( x[0], 4 ) + 22 * pow( x[0], 3 ) * x[1] + 23 * pow( x[0], 3 ) * x[2] + 11 * pow( x[0], 3 ) +
             24 * pow( x[0], 2 ) * pow( x[1], 2 ) + 25 * pow( x[0], 2 ) * x[1] * x[2] + 12 * pow( x[0], 2 ) * x[1] +
             26 * pow( x[0], 2 ) * pow( x[2], 2 ) + 13 * pow( x[0], 2 ) * x[2] + 5 * pow( x[0], 2 ) + 27 * x[0] * pow( x[1], 3 ) +
             28 * x[0] * pow( x[1], 2 ) * x[2] + 14 * x[0] * pow( x[1], 2 ) + 29 * x[0] * x[1] * pow( x[2], 2 ) +
             15 * x[0] * x[1] * x[2] + 6 * x[0] * x[1] + 30 * x[0] * pow( x[2], 3 ) + 16 * x[0] * pow( x[2], 2 ) +
             7 * x[0] * x[2] + 2 * x[0] + 31 * pow( x[1], 4 ) + 32 * pow( x[1], 3 ) * x[2] + 17 * pow( x[1], 3 ) +
             33 * pow( x[1], 2 ) * pow( x[2], 2 ) + 18 * pow( x[1], 2 ) * x[2] + 8 * pow( x[1], 2 ) + 34 * x[1] * pow( x[2], 3 ) +
             19 * x[1] * pow( x[2], 2 ) + 9 * x[1] * x[2] + 3 * x[1] + 35 * pow( x[2], 4 ) + 20 * pow( x[2], 3 ) +
             10 * pow( x[2], 2 ) + 4 * x[2] + 1;
   };

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------

   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // ----------------------------
   //  Prepare setup for 3D tests
   // ----------------------------
   MeshInfo              meshInfo3d = MeshInfo::meshCuboid( Point3D( 0.0, 0.0, 0.0 ), Point3D( 1.0, 1.0, 1.0 ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage3d( meshInfo3d, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage3d );
   std::shared_ptr< PrimitiveStorage > storage3d = std::make_shared< PrimitiveStorage >( setupStorage3d );

   // -------------------
   //  Run tests
   // -------------------
   for ( uint_t lvl = 3; lvl <= 5; ++lvl )
   {
      for ( uint8_t q = 1; q <= 4; ++q )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "2d, q=p=%d, level=%d", q, lvl ) );
         P1SurrogateOperatorTest( storage, k[q], q, lvl );
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "3d, q=p=%d, level=%d", q, lvl ) );
         P1SurrogateOperatorTest( storage3d, k[q], q, lvl );
      }
   }
   return 0;
}
