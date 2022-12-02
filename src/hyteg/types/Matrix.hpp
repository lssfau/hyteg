/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_t;
using walberla::uint_t;

/// \brief  Dense NxM Matrix
///
/// The Matrix class represents a dense MxN-dimensional matrix with basic support for algebraic operations.
/// Entries are stored in row-major ordering. Internally the class uses a dense Eigen::Matrix to store
/// and perform algebraic ops.
/// \tparam T Matrix value data type
/// \tparam M Number of rows
/// \tparam N Number of columns
template < typename T, uint_t M, uint_t N >
class Matrix
{
   /// Eigen does not allow to create column-vectors in row-major ordering
   static constexpr Eigen::StorageOptions storageType()
   {
      if constexpr ( N == 1 )
      {
         return Eigen::ColMajor;
      }
      else
      {
         return Eigen::RowMajor;
      }
   }

 public:
   // static const uint_t Size = M * N;

   /// Default constructor setting all components to zero
   Matrix() { matrix_.array() = static_cast< T >( 0 ); }

   /// Sets all values to the given constant
   explicit Matrix( const T& constant ) { matrix_.array() = constant; }

   /// Construct Matrix from an EigenMatrix
   explicit Matrix( const Eigen::Matrix< T, M, N, Eigen::RowMajor >& eigenMat ) { matrix_.array() = eigenMat.matrix_.array(); }

   /// Get reference to a single matrix component
   /// \param rIdx row index
   /// \param cIdx column index
   /// \returns Reference to component at position [rIdx,cIdx] in matrix
   T& operator()( uint_t rIdx, uint_t cIdx )
   {
      WALBERLA_ASSERT( rIdx < M, "Matrix row index out of bounds: rIdx = " << rIdx << " but M = " << M );
      WALBERLA_ASSERT( cIdx < N, "Matrix column index out of bounds: cIdx = " << cIdx << " but N = " << N );
      return matrix_( int_c( rIdx ), int_c( cIdx ) );
   }

   /// Get const reference to a single matrix component
   /// \param rIdx row index
   /// \param cIdx column index
   /// \returns Const reference to component at position [rIdx,cIdx] in matrix
   const T& operator()( uint_t rIdx, uint_t cIdx ) const
   {
      WALBERLA_ASSERT( rIdx < M, "Matrix row index out of bounds: rIdx = " << rIdx << " but M = " << M );
      WALBERLA_ASSERT( cIdx < N, "Matrix column index out of bounds: cIdx = " << cIdx << " but N = " << N );
      return matrix_( int_c( rIdx ), int_c( cIdx ) );
   }

   /// Get raw pointer to underlying matrix data
   /// \returns Pointer to first element of underlying matrix data
   T* data() { return matrix_.data(); }

   /// Get const raw pointer to underlying matrix data
   /// \returns Constant pointer to first element of underlying matrix data
   const T* data() const { return matrix_.data(); }

   Matrix< T, M, N >& operator+=( const Matrix< T, M, N >& rhs )
   {
      matrix_ += rhs.matrix_;
      return *this;
   }

   Matrix< T, M, N >& operator*=( T scalar )
   {
      matrix_.array() *= scalar;
      return *this;
   }

   Matrix< T, N, M > transpose()
   {
      Matrix< T, N, M > out;
      out.matrix_ = matrix_.transpose();
      return out;
   }

   template < uint_t N_rhs >
   Matrix< T, M, N_rhs > mul( const Matrix< T, N, N_rhs >& rhs )
   {
      Matrix< T, M, N_rhs > out;
      out.matrix_ = matrix_ * rhs.matrix_;
      return out;
   }

   PointND< T, M > mul( const PointND< T, N >& rhs ) const
   {
      PointND< T, M > out;
      out.vector_ = matrix_ * rhs.vector_;
      return out;
   }

   Matrix< T, M, N > inverse() const
   {
      if constexpr ( N != M )
      {
         WALBERLA_ABORT( "Inverse not available for matrix dimensions " << M << " x " << N );
      }
      else
      {
// #define MATRIX_CLASS_USES_INVERSE_FROM_EIGEN
#ifndef MATRIX_CLASS_USES_INVERSE_FROM_EIGEN
         if constexpr ( N == 2 )
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
         WALBERLA_ABORT( "Inverse computation not implemented for matrix dimensions " << M << " x " << N );
#else
         Matrix< T, M, N > out;
         out.matrix_ = matrix_.inverse();
#endif
      }
   }

   T det() const { return matrix_.determinant(); }

   void setAll( const T& constant ) { matrix_.array() = constant; }

   /// Internally the matrix is stored as a dense Eigen::Matrix
   Eigen::Matrix< T, M, N, Matrix< T, M, N >::storageType() > matrix_;
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
