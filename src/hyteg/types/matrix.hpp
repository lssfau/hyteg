/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include <array>
#include <cmath>
#include <iomanip>

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/types/pointnd.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief  NxM Matrix
/// \author Daniel Drzisga (drzisga@ma.tum.de)
/// \date   September, 2017
///
/// The Matrix class represents an MxN-dimensional matrix with basic support for algebraic operations
/// \tparam T Matrix value data type
/// \tparam M Number of rows
/// \tparam N Number of columns
template < typename T, uint_t M, uint_t N >
class Matrix
{
 public:
   static const uint_t Size = M * N;

   /// Default constructor setting all components to zero
   Matrix()
   {
      for ( uint_t i = 0; i < Size; ++i )
      {
         x[i] = (T) 0;
      }
   }

   /// Sets all values to the given constant
   Matrix( const T& constant )
   {
      for ( uint_t i = 0; i < Size; ++i )
      {
         x[i] = constant;
      }
   }

   /// Constructs the matrix using values from M*N-dimensional array \p _x in row-major order
   /// \param _x Pointer to M*N-dimensional array
   Matrix( T _x[Size] )
   {
      for ( uint_t i = 0; i < Size; ++i )
      {
         x[i] = _x[i];
      }
   }

   /// Copy constructor
   /// \param b Reference to another instance of Matrix
   Matrix( const Matrix& b )
   {
      for ( uint_t i = 0; i < Size; ++i )
      {
         x[i] = b.x[i];
      }
   }

   /// Get reference to a single matrix component
   /// \param row Row index
   /// \param col Column index
   /// \returns Reference to component at position [row,col] in matrix
   T& operator()( uint_t row, uint_t col )
   {
      WALBERLA_ASSERT( row < M, "Matrix row index out of bounds: row = " << row << " but M = " << M );
      WALBERLA_ASSERT( col < N, "Matrix column index out of bounds: col = " << col << " but N = " << N );
      return x[N * row + col];
   }

   /// Get const reference to a single matrix component
   /// \param row Row index
   /// \param col Column index
   /// \returns Const reference to component at position [row,col] in matrix
   const T& operator()( uint_t row, uint_t col ) const
   {
      WALBERLA_ASSERT( row < M, "Matrix row index out of bounds: row = " << row << " but M = " << M );
      WALBERLA_ASSERT( col < N, "Matrix column index out of bounds: col = " << col << " but N = " << N );
      return x[N * row + col];
   }

   /// Get raw pointer to underlying matrix data
   /// \returns Pointer to first element of underlying matrix data
   T* data() { return &x[0]; }

   /// Get const raw pointer to underlying matrix data
   /// \returns Constant pointer to first element of underlying matrix data
   const T* data() const { return &x[0]; }

   Matrix< T, M, N >& operator+=( const Matrix< T, M, N >& rhs )
   {
      for ( uint_t i = 0; i < Size; ++i )
      {
         x[i] += rhs.x[i];
      }
      return *this;
   }

   Matrix< T, M, N >& operator*=( T scalar )
   {
      for ( uint_t i = 0; i < Size; ++i )
      {
         x[i] *= scalar;
      }
      return *this;
   }

   Matrix< T, N, M > transpose()
   {
      Matrix< T, N, M > out;
      for ( uint_t i = 0; i < M; ++i )
      {
         for ( uint_t j = 0; j < N; ++j )
         {
            out( j, i ) = ( *this )( i, j );
         }
      }
      return out;
   }

   template < uint_t N_rhs >
   Matrix< T, M, N_rhs > mul( const Matrix< T, N, N_rhs >& rhs )
   {
      Matrix< T, M, N_rhs > out;
      for ( uint_t i = 0; i < M; ++i )
      {
         for ( uint_t j = 0; j < N_rhs; ++j )
         {
            for ( uint_t k = 0; k < N; ++k )
            {
               out( i, j ) += ( *this )( i, k ) * rhs( k, j );
            }
         }
      }
      return out;
   }

   PointND< T, M > mul( const PointND< T, N >& rhs ) const
   {
      PointND< T, M > out;
      for ( uint_t i = 0; i < M; ++i )
      {
         for ( uint_t j = 0; j < N; ++j )
         {
            out[i] += ( *this )( i, j ) * rhs[j];
         }
      }
      return out;
   }

   Matrix< T, M, N > inverse() const
   {
      if constexpr ( N == M && N == 2 )
      {
         Matrix< T, M, N > out;

         T det       = ( *this )( 0, 0 ) * ( *this )( 1, 1 ) - ( *this )( 0, 1 ) * ( *this )( 1, 0 );
         det         = static_cast< T >( 1 ) / det;
         out( 0, 0 ) = +( *this )( 1, 1 ) * det;
         out( 0, 1 ) = -( *this )( 0, 1 ) * det;
         out( 1, 0 ) = -( *this )( 1, 0 ) * det;
         out( 1, 1 ) = +( *this )( 0, 0 ) * det;

         return out;
      }
      WALBERLA_ABORT( "Adjugate computation not implemented for matrix dimensions " << M << " x " << N )
   }

   T det() const
   {
      if constexpr ( N == M && N == 2 )
      {
         return ( *this )( 0, 0 ) * ( *this )( 1, 1 ) - ( *this )( 0, 1 ) * ( *this )( 1, 0 );
      }
      else if constexpr ( N == M && N == 3 )
      {
         return ( *this )( 0, 0 ) * ( *this )( 1, 1 ) * ( *this )( 2, 2 ) +
                ( *this )( 0, 1 ) * ( *this )( 1, 2 ) * ( *this )( 2, 0 ) +
                ( *this )( 0, 2 ) * ( *this )( 1, 0 ) * ( *this )( 2, 1 ) -
                ( *this )( 2, 0 ) * ( *this )( 1, 1 ) * ( *this )( 0, 2 ) -
                ( *this )( 2, 1 ) * ( *this )( 1, 2 ) * ( *this )( 0, 0 ) -
                ( *this )( 2, 2 ) * ( *this )( 1, 0 ) * ( *this )( 0, 1 );
      }

      WALBERLA_ABORT( "Determinant computation not implemented for matrix dimensions " << M << " x " << N )
   }

   void setAll( const T& constant )
   {
      for ( uint_t i = 0; i < Size; ++i )
      {
         x[i] = constant;
      }
   }

 private:
   T x[Size];
};

template < typename T, uint_t M, uint_t N >
inline Matrix< T, M, N > operator*( T scalar, Matrix< T, M, N > rhs )
{
   rhs *= scalar;
   return rhs;
}

template < typename T, uint_t M, uint_t N >
inline Matrix< T, M, N > operator*( Matrix< T, M, N > lhs, T scalar )
{
   lhs *= scalar;
   return lhs;
}

template < typename T, uint_t M, uint_t N >
inline std::ostream& operator<<( std::ostream& os, const Matrix< T, M, N >& matrix )
{
   os << "[\n";

   for ( uint_t i = 0; i < M; ++i )
   {
      os << "[";
      for ( uint_t j = 0; j < N; ++j )
      {
         os << std::scientific << std::setw( 13 ) << matrix( i, j );
         if ( j != N - 1 )
         {
            os << ", ";
         }
      }
      os << "],\n";
   }

   os << "]";

   return os;
}

template < uint_t M, uint_t N >
using Matrixr = Matrix< real_t, M, N >;
typedef Matrix< real_t, 2, 2 >   Matrix2r;
typedef Matrix< real_t, 3, 3 >   Matrix3r;
typedef Matrix< real_t, 4, 4 >   Matrix4r;
typedef Matrix< real_t, 6, 6 >   Matrix6r;
typedef Matrix< real_t, 10, 10 > Matrix10r;

} // namespace hyteg
