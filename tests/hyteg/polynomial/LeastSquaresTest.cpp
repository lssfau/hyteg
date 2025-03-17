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
#include <core/math/Constants.h>
#include <core/math/Random.h>
#include <core/timing/Timer.h>
#include <filesystem>
#include <hyteg/polynomial/elementwise/leastSquares.hpp>

using hyteg::idx_t;
using hyteg::real_t;
using hyteg::uint_t;
using walberla::math::pi;
using walberla::math::realRandom;

// test accuracy of least squares fit
template < uint_t D >
double LeastSquaresTest( uint_t q, uint_t lvl, uint_t downsampling, bool use_precomputed )
{
   // ---------------------------------------------------------
   /// initialize least squares system
   // ---------------------------------------------------------
   std::unique_ptr< hyteg::surrogate::LeastSquares< real_t > > lsq;
   if ( use_precomputed )
   {
      lsq = std::make_unique< hyteg::surrogate::LeastSquares< real_t > >( "svd", D, q, lvl, downsampling );
   }
   else
   {
      lsq = std::make_unique< hyteg::surrogate::LeastSquares< real_t > >( D, q, lvl, downsampling );
      // write svd to file to be used when calling the function again with `use_precomputed=true`
      lsq->write_to_file( "svd" );
   }
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "         n_samples = %d; dim(P) = %d", lsq->rows, lsq->cols ) );

   // ---------------------------------------------------------
   /// initialize data
   // ---------------------------------------------------------
   // define smooth function
   auto f = []( const hyteg::Point3D& x ) {
      return real_t( std::sin( pi / 4.0 * x[0] ) * std::cos( pi / 4.0 * x[1] ) * std::sin( pi / 4.0 * x[2] ) );
   };
   // coordinates on tetrahedron
   hyteg::surrogate::polynomial::Domain< real_t > X( lvl );
   // fill right-hand side
   auto it = lsq->samplingIterator();
   while ( it != it.end() )
   {
      auto x = X( it.ijk() );
      lsq->setRHS( it(), f( x ) );
      ++it;
   }

   // ---------------------------------------------------------
   /// polynomial approximation
   // ---------------------------------------------------------
   auto& coeffs = lsq->solve();
   // initialize polynomial
   hyteg::surrogate::polynomial::Polynomial< real_t, D > p( uint8_t( q ), coeffs );

   // ---------------------------------------------------------
   /// evaluate polynomial and compute error
   // ---------------------------------------------------------
   idx_t  len_edge = ( idx_t( 1 ) << lvl );
   double error    = 0.0;
   idx_t  k_max    = ( D == 3 ) ? len_edge : 1;
   idx_t  j_max    = ( D >= 2 ) ? len_edge : 1;
   idx_t  i_max    = len_edge;
   for ( idx_t k = 0; k < k_max; ++k )
   {
      auto z = X[k];
      if constexpr ( D == 3 )
      {
         p.fix_z( z );
      }
      for ( idx_t j = 0; j < j_max - k; ++j )
      {
         auto y = X[j];
         if constexpr ( D >= 2 )
         {
            p.fix_y( y );
         }
         for ( idx_t i = 0; i < i_max - k - j; ++i )
         {
            auto x = X[i];
            auto e = f( { x, y, z } ) - p.eval( x );
            error += e * e;
         }
      }
   }
   auto h  = 1.0 / double( len_edge );
   auto dV = std::pow( h, D );
   error   = std::sqrt( error * dV );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "         discrete L2 error: ||f-p|| = %e", error ) );

   return error;
}

template < uint8_t D >
void all_tests()
{
   const uint_t          lvl = 4;
   std::vector< double > epsilon;
   if constexpr ( D == 2 )
   {
      epsilon = { 3e-1, 2e-1, 1e-1, 2e-2, 2e-2, 2e-3, 8e-4, 8e-5, 2e-5, 4e-7 };
   }
   else
   {
      epsilon = { 2e-1, 1e-1, 5e-2, 3e-2, 8e-3, 3e-3, 4e-4, 2e-4, 9e-6, 2e-6 };
   }

   for ( uint_t q = 0; q <= 9; ++q )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%dD, level %d, polynomial degree %d", D, lvl, q ) );
      for ( uint_t ds = 1; ds <= 3; ++ds )
      {
         // with single precision, accuracy won't get better than 5e-4
         if ( std::is_same< real_t, float >() && epsilon[q] < 5e-4 )
         {
            epsilon[q] = 5e-4;
         }
         if ( q <= hyteg::surrogate::LeastSquares< real_t >::max_degree( lvl, ds ) )
         {
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "   downsampling factor %d", ds ) );

            WALBERLA_LOG_INFO_ON_ROOT( "      compute svd" );
            auto err1 = LeastSquaresTest< D >( q, lvl, ds, 0 );
            WALBERLA_CHECK_LESS( err1, epsilon[q], "accuracy LSQ" );

            WALBERLA_LOG_INFO_ON_ROOT( "      use precomputed svd" );
            auto   err2  = LeastSquaresTest< D >( q, lvl, ds, 1 );
            double eps_m = std::is_same< real_t, double >() ? 1e-14 : 1e-6;
            WALBERLA_CHECK_LESS( std::abs( err1 - err2 ), err1 * eps_m, "precomputed svd" );
         }
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   // cleanup
   std::filesystem::remove_all( "svd" );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::filesystem::create_directory( "svd" );

   // -------------------
   //  Run tests
   // -------------------
   WALBERLA_LOG_INFO_ON_ROOT( "LeastSquaresTest" );

   all_tests< 2 >();
   all_tests< 3 >();

   return 0;
}
