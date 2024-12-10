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
template < uint8_t lvl >
static constexpr uint_t n_edge = 1 << lvl;
template < uint8_t lvl >
static constexpr uint_t n_edge_int = ( lvl > 1 ) ? n_edge< lvl > - 2 : 0;
template < uint8_t lvl >
static constexpr uint_t n_face = n_edge< lvl > * ( n_edge< lvl > + 1 ) / 2;
template < uint8_t lvl >
static constexpr uint_t n_face_int = ( lvl > 1 ) ? n_face< lvl > - 3 * n_edge_int< lvl > - 3 : 0;
template < uint8_t lvl >
static constexpr uint_t n_cell = polynomial::dimP< 3, n_edge< lvl > - 1 >;
template < uint8_t lvl >
static constexpr uint_t n_cell_int = ( lvl > 2 ) ? n_cell< lvl > - 4 * n_face_int< lvl > - 6 * n_edge_int< lvl > - 4 : 0;
// use different level for each primitive type
template < uint8_t l_face, uint8_t l_edge >
static constexpr uint_t n_face_var = n_face_int< l_face > + 3 * n_edge_int< l_edge > + 3;
template < uint8_t l_cell, uint8_t l_face, uint8_t l_edge >
static constexpr uint_t n_cell_var = n_cell_int< l_cell > + 4 * n_face_int< l_face > + 6 * n_edge_int< l_edge > + 4;
} // namespace interpolation

// todo remove
// static constexpr uint_t n = interpolation::n_cell<5>;
// static constexpr uint_t m = interpolation::n_cell_var<4, 5, 5>;

using walberla::uint_t;

template < uint8_t dim, uint8_t degree, uint8_t lvl, bool reduced_sample_size = false, bool precomputed = true >
class LeastSquares
{
 private:
   static constexpr uint_t n_edge = interpolation::n_edge< lvl >;
   // reduced number of interpolation points for higher levels
   static constexpr uint8_t face_lvl_red = ( lvl > 5 ) ? 5 : lvl;
   static constexpr uint8_t cell_lvl_red = ( lvl > 4 ) ? 4 : lvl;
   static constexpr int     rows_reduced = ( dim == 2 ) ? interpolation::n_face_var< face_lvl_red, lvl > :
                                                          interpolation::n_cell_var< cell_lvl_red, face_lvl_red, lvl >;

 public:
   // number of interpolation points
   static constexpr int rows = reduced_sample_size ? rows_reduced : interpolation::n_cell< lvl >;
   // dimension of polynomial space
   static constexpr int cols = polynomial::dimP< dim, degree >;

   struct Iterator
   {
      uint_t n; // index of interpolation point
      uint_t i; // index of x-coordinate
      uint_t j; // index of y-coordinate
      uint_t k; // index of z-coordinate

      inline Iterator()
      : n( 0 )
      , i( 0 )
      , j( 0 )
      , k( 0 )
      , di( 1 )
      , dj( 1 )
      {}

      inline Iterator& operator++()
      {
         if constexpr ( reduced_sample_size )
         {
            increment_ijk_reduced();
         }
         else
         {
            increment_ijk();
         }
         ++n;
         return *this;
      }

    private:
      inline constexpr void increment_ijk()
      {
         ++i;
         if ( i >= n_edge - j - k )
         {
            i = 0;
            ++j;
            if ( j >= n_edge - k )
            {
               j = 0;
               ++k;
            }
         }
      }

      inline constexpr void increment_ijk_reduced()
      {
         i += di;
         if ( i >= n_edge - j - k ) // leaving triangle in x direction
         {
            if ( i - di < n_edge - j - k - 1 ) // jumped over 1-2-3 face
            {
               // add point on the 1-2-3 face (or 1-3 edge)
               i = n_edge - j - k - 1;
            }
            else
            {
               // increment j
               i = 0;
               j += dj;
               // adjust di
               if ( k == 0 && j % f_stride == 0 )
               {
                  // 0-1-2 face interior
                  di = f_stride; // reduced number of points on face
               }
               else if ( j % c_stride == 0 && k % c_stride == 0 )
               {
                  // cell interior
                  di = c_stride; // reduced number of points in cell
               }
               else
               {
                  // skip these lines
                  di = n_edge;
               }
            }
            if ( j >= n_edge - k ) // leaving triangle in y direction
            {
               if ( j - dj < n_edge - k - 1 ) // jumped over 2-3 edge
               {
                  // add point on the 2-3 edge
                  j = n_edge - k - 1;
               }
               else
               {
                  // increment k
                  j = 0;
                  ++k;
                  // adjust di and dj
                  if ( k % f_stride == 0 )
                  {
                     // 0-1-2 face interior
                     di = f_stride; // reduced number of points on face
                     // 0-1-3 face interior
                     dj = f_stride; // reduced number of points on face
                  }
                  else
                  {
                     // skip these planes
                     di = n_edge;
                     dj = n_edge;
                  }
               }
            }
         }
      }

      int di;
      int dj;

      static constexpr uint_t f_stride = 1 << ( lvl - face_lvl_red );
      static constexpr uint_t c_stride = 1 << ( lvl - cell_lvl_red );
   };

   LeastSquares()
   : A( rows, cols )
   , Uh( cols, rows )
   , Si( cols )
   , V( cols, cols )
   , b( rows )
   , c( cols )
   {
      if constexpr ( precomputed )
      {
         WALBERLA_LOG_WARNING_ON_ROOT(
             "For the specified template arguments, surrogate::LeasSquares has not been precomputed.\nComputing SVD now." );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Setup Vandermonde matrix" );
      constexpr polynomial::Basis< dim, degree > phi;
      constexpr polynomial::Coordinates          coords( lvl );

      Iterator it;
      while ( it.n < rows )
      {
         auto x = coords.x( it.i, it.j, it.k );

         for ( uint_t col = 0; col < cols; ++col )
         {
            A( it.n, col ) = phi.eval( col, x );
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

   void     setRHS( const Eigen::Matrix< double, rows, 1 >& rhs ) { b = rhs; }
   Iterator samplingIterator() const { return Iterator(); }
   void     setRHS( const uint_t n, const double val ) { b( n ) = val; }

   const auto& solve()
   {
      c = V * ( Si * ( Uh * b ) );
      return c;
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