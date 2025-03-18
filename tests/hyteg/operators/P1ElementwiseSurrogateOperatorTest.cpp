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
#include <core/math/Random.h>
#include <core/mpi/MPIManager.h>
#include <hyteg/elementwiseoperators/P1ElementwiseOperator.hpp>
#include <hyteg/elementwiseoperators/P1ElementwiseSurrogateOperator.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp>

using walberla::real_t;
using namespace hyteg;

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
template < uint8_t DIM, uint8_t DEGREE >
void P1SurrogateOperatorTest( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t level )
{
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "level=%d, degree=%d, %dD", level, DEGREE, DIM ) );

   double epsilon, errorMax;

   // setup pde coefficient k ∈ P_q
   hyteg::surrogate::polynomial::Polynomial< real_t, DIM, DEGREE > k_poly;
   for ( auto& c : k_poly )
   {
      c = walberla::math::realRandom();
   }
   auto k = [&]( const hyteg::Point3D& x ) { return k_poly.eval_naive( x ); };

   // operators
   forms::p1_div_k_grad_affine_q3                         form( k, k );
   P1ElementwiseAffineDivKGradOperator                    A( storage, level, level, form );
   P1ElementwiseSurrogateAffineDivKGradOperator< DEGREE > A_q( storage, level, level, form );

   A_q.init( 1, "", false );

   // functions
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > Au( "Au", storage, level, level );
   hyteg::P1Function< real_t > Aqu( "(A_q)u", storage, level, level );
   hyteg::P1Function< real_t > err( "(A-A_q)u", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > initialU = []( const hyteg::Point3D& x ) {
      return cos( 2 * M_PI * x[0] ) * cos( 2 * M_PI * x[1] ) * cos( 2 * M_PI * x[2] );
   };
   u.interpolate( initialU, level );

   epsilon = std::is_same< real_t, double >() ? 2e-12 : 1e-5;

   // apply operators
   A.apply( u, Au, level, All, Replace );
   A_q.apply( u, Aqu, level, All, Replace );
   err.assign( { 1.0, -1.0 }, { Au, Aqu }, level, All );
   errorMax = err.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%37s = %e", "||(A - A_q)u||_inf", errorMax ) )
   WALBERLA_CHECK_LESS( errorMax, epsilon, "||(A - A_q)u||_inf" );

   /* for some reason the diagonals of the 2d operators are quite different for Surrogate and
      conventional operator. The difference seems to vanish for level->∞, though.
      This should probably be investigated!
    */
   if ( storage->hasGlobalCells() )
   {
      epsilon = std::is_same< real_t, double >() ? 1e-8 : 3e-4;
   }
   else
   {
      epsilon = std::is_same< real_t, double >() ? 1e-5 : 3e-4;
   }

   // diagonal values
   A.computeDiagonalOperatorValues();
   A_q.computeDiagonalOperatorValues();
   err.assign( { 1.0, -1.0 }, { *A.getDiagonalValues(), *A_q.getDiagonalValues() }, level, All );
   errorMax = err.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%37s = %e", "||diag(A) - diag(A_q)||_inf", errorMax ) )
   WALBERLA_CHECK_LESS( errorMax, epsilon, "||diag(A) - diag(A_q)||_inf" );

   // inverse diagonal values
   A.computeInverseDiagonalOperatorValues();
   A_q.computeInverseDiagonalOperatorValues();
   err.assign( { 1.0, -1.0 }, { *A.getInverseDiagonalValues(), *A_q.getInverseDiagonalValues() }, level, All );
   errorMax = err.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%37s = %e", "||inv(diag(A)) - inv(diag(A_q))||_inf", errorMax ) )
   WALBERLA_CHECK_LESS( errorMax, epsilon, "||inv(diag(A)) - inv(diag(A_q))||_inf" );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage2d = std::make_shared< PrimitiveStorage >( setupStorage );

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
      P1SurrogateOperatorTest< 2, 1 >( storage2d, lvl );
      P1SurrogateOperatorTest< 3, 1 >( storage3d, lvl );
      P1SurrogateOperatorTest< 2, 2 >( storage2d, lvl );
      P1SurrogateOperatorTest< 3, 2 >( storage3d, lvl );
      P1SurrogateOperatorTest< 2, 3 >( storage2d, lvl );
      P1SurrogateOperatorTest< 3, 3 >( storage3d, lvl );
      P1SurrogateOperatorTest< 2, 4 >( storage2d, lvl );
      P1SurrogateOperatorTest< 3, 4 >( storage3d, lvl );
      P1SurrogateOperatorTest< 2, 5 >( storage2d, lvl );
      P1SurrogateOperatorTest< 3, 5 >( storage3d, lvl );
   }
   return 0;
}
