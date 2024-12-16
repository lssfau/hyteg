/*
 * Copyright (c) 2024 Benjamin Mann.
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

#pragma once

#include <core/logging/Logging.h>
#include <hyteg/eigen/EigenWrapper.hpp>

#include "polynomial.hpp"

namespace hyteg {
namespace surrogate {

namespace interpolation {

/**
 * @brief Computes the number of vertices in a triangle with edge length n
 *
 * @param n number of vertices along an edge
 * @return The triangular number of n.
 */
static constexpr inline uint_t tri( uint_t n )
{
   return n * ( n + 1 ) / 2;
}

/**
 * @brief Computes the number of vertices in a tetrahedron with edge length n
 *
 * @param n number of vertices along an edge
 * @return The tetrahedral number of n.
 */
static constexpr inline uint_t tet( uint_t n )
{
   return n * ( n + 1 ) * ( n + 2 ) / 6;
}

/**
 * @brief Computes the downsampled number of vertices along an edge on a given level
 *
 * @param lvl The level of refinement.
 * @param downsampling The downsampling factor.
 * @return The number of downsampled vertices along the edge.
 */
static constexpr inline uint_t n_edge( uint_t lvl, uint_t downsampling )
{
   return ( ( ( 1 << lvl ) - 1 ) / downsampling ) + 1; // ceil(n/d) with n = 2^lvl, d = downsampling
}

/**
 * @brief Computes the downsampled number of vertices in a volume element
 *
 * @tparam D The dimension (2 for 2D, 3 for 3D).
 * @param lvl The level of refinement.
 * @param downsampling The downsampling factor.
 * @return The number of vertices in the volume.
 */
template < uint8_t D >
static constexpr inline uint_t n_volume( uint_t n_edge )
{
   return ( D == 2 ) ? tri( n_edge ) : tet( n_edge );
}

/**
 * @brief Computes the number of vertices in a volume element with given edge length
 *
 * @tparam D The dimension (2 for 2D, 3 for 3D).
 * @param n_edge The number of vertices along an edge.
 * @return The number of vertices in the volume.
 */
template < uint8_t D >
static constexpr inline uint_t n_volume( uint_t lvl, uint_t downsampling )
{
   return n_volume< D >( n_edge( lvl, downsampling ) );
}

} // namespace interpolation

using walberla::uint_t;

// todo change lvl and l_sample to runtime args
template < uint8_t D, uint8_t Q, uint8_t lvl, uint8_t downsampling = 1, bool precomputed = true >
class LeastSquares
{
 private:
   class Iterator
   {
    public:
      inline Iterator( uint_t mylvl, uint_t d )
      : _n( 0 )
      , _i( 0 )
      , _j( 0 )
      , _k( 0 )
      , _ijk_max( interpolation::n_edge( mylvl, 1 ) )
      , _n_max( interpolation::n_volume< D >( mylvl, d ) )
      , _stride( d )
      {}

      inline Iterator& operator++()
      {
         _i += _stride;
         if ( _i >= _ijk_max - _j - _k )
         {
            _i = 0;
            _j += _stride;
            if ( _j >= _ijk_max - _k )
            {
               _j = 0;
               _k += _stride;
            }
         }
         ++_n;
         return *this;
      }

      inline constexpr uint_t i() const { return _i; }
      inline constexpr uint_t j() const { return _j; }
      inline constexpr uint_t k() const { return _k; }
      inline constexpr uint_t operator()() const { return _n; }
      inline constexpr uint_t end() const { return _n_max; }
      inline constexpr bool   operator!=( const uint_t other ) const { return _n != other; }
      inline constexpr bool   operator==( const uint_t other ) const { return _n == other; }

    private:
      uint_t _n; // index of sample point (row in A)
      uint_t _i; // index of x-coordinate
      uint_t _j; // index of y-coordinate
      uint_t _k; // index of z-coordinate

      uint_t _ijk_max; // number of microedges along an edge (only coincides with sample points along edge if downsampling=1)
      uint_t _n_max;   // number of sample points
      uint_t _stride;  // downsampling factor
   };

 public:
   // number of interpolation points
   static constexpr int rows = interpolation::n_volume< D >( lvl, downsampling );
   // dimension of polynomial space
   static constexpr int cols = polynomial::dimP< D, Q >;

   LeastSquares()
   : A( rows, cols )
   , Uh( cols, rows )
   , Si( cols )
   , V( cols, cols )
   , b( rows )
   , c( cols )
   {
      // todo change to runtime check
      if constexpr ( precomputed )
      {
         WALBERLA_LOG_WARNING_ON_ROOT(
             "For the specified template arguments, surrogate::LeasSquares has not been precomputed.\nComputing SVD now." );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Setup Vandermonde matrix" );
      constexpr polynomial::Basis< D, Q > phi;
      constexpr polynomial::Coordinates   coords( lvl );

      auto it = samplingIterator();
      while ( it != it.end() )
      {
         auto x = coords.x( it.i(), it.j(), it.k() );

         for ( uint_t col = 0; col < cols; ++col )
         {
            A( it(), col ) = phi.eval( col, x );
         }
         ++it;
      }

      // SVD
      WALBERLA_LOG_INFO_ON_ROOT( "Compute SVD" );
      Eigen::JacobiSVD< Matrix > svd( A, Eigen::ComputeThinU | Eigen::ComputeThinV );
      WALBERLA_LOG_INFO_ON_ROOT( "Store SVD" );
      Uh = svd.matrixU().adjoint();
      Si = svd.singularValues().cwiseInverse();
      V  = svd.matrixV();
   }

   Iterator samplingIterator() const { return Iterator( lvl, downsampling ); }
   void     setRHS( const uint_t n, const double val ) { b( n ) = val; }

   const polynomial::Polynomial< D, Q > solve()
   {
      c = V * ( Si.cwiseProduct( Uh * b ) );
      return polynomial::Polynomial< D, Q >( c );
   }

   using Matrix = Eigen::Matrix< double, -1, -1, Eigen::RowMajor >;
   using Vector = Eigen::Matrix< double, -1, 1, Eigen::ColMajor >;

 private:
   // Vandermonde matrix
   Matrix A;
   // SVD of Vandermonde matrix
   Matrix Uh; // left singular vectors
   Vector Si; // inverse of singular values
   Matrix V;  // right singular vectors
   // rhs of least squares problem
   Vector b;
   // coefficients of the polynomial
   Vector c;
};

} // namespace surrogate

} // namespace hyteg