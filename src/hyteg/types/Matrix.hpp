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

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/eigen/EigenWrapper.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType, int M, int N >
class Matrix : public Eigen::Matrix< ValueType, M, N, N == 1 ? Eigen::ColMajor : Eigen::RowMajor >
{
   using Eigen::Matrix< ValueType, M, N, N == 1 ? Eigen::ColMajor : Eigen::RowMajor >::Matrix;
};

// template < typename ValueType, int M, int N >
// using Matrix = Eigen::Matrix< ValueType, M, N, N == 1 ? Eigen::ColMajor : Eigen::RowMajor >;

template < int M, int N >
using Matrixr = Matrix< real_t, M, N >;

using Matrix2r  = Matrixr< 2, 2 >;
using Matrix3r  = Matrixr< 3, 3 >;
using Matrix4r  = Matrixr< 4, 4 >;
using Matrix6r  = Matrixr< 6, 6 >;
using Matrix10r = Matrixr< 10, 10 >;
using MatrixXr  = Matrixr< Eigen::Dynamic, Eigen::Dynamic >;
using VectorXr  = Matrixr< Eigen::Dynamic, 1 >;

extern template class Matrix< real_t, 2, 2 >;
extern template class Matrix< real_t, 3, 3 >;
extern template class Matrix< real_t, 4, 4 >;
extern template class Matrix< real_t, 6, 6 >;
extern template class Matrix< real_t, 10, 10 >;
extern template class Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >;
extern template class Matrix< real_t, Eigen::Dynamic, 1 >;

} // namespace hyteg

namespace walberla::mpi {

template < typename T, // Element type of SendBuffer
           typename G, // Growth policy of SendBuffer
           typename MatrixDataType,
           int M,
           int N >
GenericSendBuffer< T, G >& operator<<( GenericSendBuffer< T, G >& buf, const hyteg::Matrix< MatrixDataType, M, N >& matrix )
{
   for ( int i = 0; i < M; ++i )
   {
      for ( int j = 0; j < N; j++ )
      {
         buf << matrix( i, j );
      }
   }

   return buf;
}

template < typename T, // Element type  of RecvBuffer
           typename MatrixDataType,
           int M,
           int N >
GenericRecvBuffer< T >& operator>>( GenericRecvBuffer< T >& buf, hyteg::Matrix< MatrixDataType, M, N >& matrix )
{
   for ( int i = 0; i < M; ++i )
   {
      for ( int j = 0; j < N; j++ )
      {
         buf >> matrix( i, j );
      }
   }
   return buf;
}


} // namespace walberla::mpi
