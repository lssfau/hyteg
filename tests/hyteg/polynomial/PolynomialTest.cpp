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

#include <core/Environment.h>
#include <core/Format.hpp>
#include <core/config/Create.h>
#include <core/logging/Logging.h>
#include <core/math/Constants.h>
#include <core/math/Random.h>
#include <core/timing/Timer.h>
#include <filesystem>
#include <hyteg/polynomial/elementwise/polynomial.hpp>

using hyteg::idx_t;
using hyteg::real_t;
using hyteg::uint_t;
using walberla::math::pi;
using walberla::math::realRandom;

// test correctness of basis functions
void MonomialBasisTest( int i, int j, int k )
{
   // ---------------------------------------------------------
   /// initialize polynomial p(x,y,z) = x^i * y^j * z^k
   // ---------------------------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Monomial(%d, %d, %d)", i, j, k ) );
   auto p = hyteg::surrogate::polynomial::Monomial( uint_t( i ), uint_t( j ), uint_t( k ) );

   // ---------------------------------------------------------
   /// test compression
   // ---------------------------------------------------------
   WALBERLA_ASSERT_EQUAL( i, p.i(), "p.i" );
   WALBERLA_ASSERT_EQUAL( j, p.j(), "p.j" );
   WALBERLA_ASSERT_EQUAL( k, p.k(), "p.k" );
   WALBERLA_ASSERT_EQUAL( i, p.expand()[0], "p.expand" );
   WALBERLA_ASSERT_EQUAL( j, p.expand()[1], "p.expand" );
   WALBERLA_ASSERT_EQUAL( k, p.expand()[2], "p.expand" );
   WALBERLA_ASSERT_EQUAL( i + j + k, p.degree(), "p.degree" );

   // ---------------------------------------------------------
   /// test p.eval
   // ---------------------------------------------------------
   // choose random point x in R^3
   hyteg::Point3D x{ realRandom(), realRandom(), realRandom() };
   // evaluate p manually
   auto px_manual = std::pow( x[0], i ) * std::pow( x[1], j ) * std::pow( x[2], k );
   // evaluate p using built-in function
   auto p_eval_x = p.eval( x );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "px_manual - p.eval(x) = %e", px_manual - p_eval_x ) );
   WALBERLA_ASSERT_FLOAT_EQUAL( px_manual, p_eval_x, "p.eval" );
}

// test different algorithms to evaluate polynomials
void PolynomialTest( uint8_t d, uint8_t q )
{
   // ---------------------------------------------------------
   /// initialize
   // ---------------------------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Polynomial(d=%d, q=%d)", d, q ) );
   // choose random element p from P_q(R^d)
   hyteg::surrogate::polynomial::Polynomial< real_t > p( d, q );
   for ( auto& c : p )
   {
      c = realRandom();
   }
   // choose random point x from R^3
   hyteg::Point3D x{ realRandom(), realRandom(), realRandom() };

   // ---------------------------------------------------------
   /// evaluate p(x) manually
   // ---------------------------------------------------------
   real_t px_manual = 0.0;
   // initialize powers of x,y,z
   Eigen::Vector< real_t, -1 > x_pow( q + 1 ); // x^0, x^1, ...
   Eigen::Vector< real_t, -1 > y_pow( q + 1 ); // y^0, y^1, ...
   Eigen::Vector< real_t, -1 > z_pow( q + 1 ); // z^0, z^1, ...
   x_pow[0] = y_pow[0] = z_pow[0] = 1.0;
   for ( int i = 1; i <= q; ++i )
   {
      x_pow[i] = x_pow[i - 1] * x[0];
      y_pow[i] = y_pow[i - 1] * x[1];
      z_pow[i] = z_pow[i - 1] * x[2];
   }
   // sum up contributions of each basis function φ_n, i.e., p(x) = ∑_n c_n φ_n(x)
   for ( idx_t n = 0; n < p.n_coefficients(); ++n )
   {
      auto [i, j, k] = p.phi( n ).expand(); // φ_n = x^i y^j z^k
      px_manual += p.c( n ) * x_pow[i] * y_pow[j] * z_pow[k];
   }

   // ---------------------------------------------------------
   /// use naive evaluation (should be equivalent to the above)
   // ---------------------------------------------------------
   real_t px_naive = p.eval_naive( x );

   // ---------------------------------------------------------
   /// evaluate p(x) using Horner's method
   // ---------------------------------------------------------
   if ( d == 3 )
   {
      p.fix_z( x[2] ); // restrict p to P_q(R^2)
   }
   if ( d >= 2 )
   {
      p.fix_y( x[1] ); // restrict to P_q(R)
   }
   // evaluate 1d polynomial by Horner's method
   real_t px_horner = p.eval( x[0] );

   // ---------------------------------------------------------
   /// compare solutions
   // ---------------------------------------------------------
   double epsilon = std::is_same< real_t, double >() ? 1e-15 : 1e-6;
   auto   check   = [&]( const std::string& method, double px ) {
      double diff   = std::abs( px_manual - px );
      double px_abs = std::abs( px_manual );
      WALBERLA_LOG_INFO_ON_ROOT( "   " << method << ":" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "      |px_manual - p(x)|/|px_manual| = %e", diff / px_abs ) );
      WALBERLA_CHECK_LESS( diff, epsilon * px_abs, "accuracy " << method );
   };
   check( "Naive evaluation", px_naive );
   check( "Horner's method", px_horner );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // -------------------
   //  Run tests
   // -------------------
   for ( int i = 0; i < 5; ++i )
   {
      for ( int j = 0; j < 5; ++j )
      {
         for ( int k = 0; k < 5; ++k )
         {
            MonomialBasisTest( i, j, k );
         }
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "=======================================================" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   for ( uint8_t d = 1; d <= 3; ++d )
   {
      for ( uint8_t deg = 0; deg <= 12; ++deg )
      {
         PolynomialTest( d, deg );
      }
   }

   return 0;
}
