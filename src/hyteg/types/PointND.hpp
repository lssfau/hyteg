/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl,
 * Marcus Mohr.
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

namespace walberla::math {
template < typename Type >
class Vector3;
}

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief  N-dimensional vector
///
/// The PointND class represents an N-dimensional vector with support for algebraic operations.
/// Internally the class uses a dense Eigen::Matrix for this and at some point we might want to
/// switch to that representation altogether
/// \tparam ValueType Vector data type
/// \tparam N Dimension of vector
template < typename ValueType, int N >
class PointND : public Eigen::Matrix< ValueType, N, 1 >
{
   using Eigen::Matrix< ValueType, N, 1 >::Matrix;

 public:
   PointND()
   : Eigen::Matrix< ValueType, N, 1 >( Eigen::Matrix< ValueType, N, 1 >::Zero() )
   {}
};

using Point2D  = PointND< walberla::real_t, 2 >;
using Point3D  = PointND< walberla::real_t, 3 >;
using Point4D  = PointND< walberla::real_t, 4 >;
using Point6D  = PointND< walberla::real_t, 6 >;
using Point10D = PointND< walberla::real_t, 10 >;

extern template class PointND< walberla::real_t, 2 >;
extern template class PointND< walberla::real_t, 3 >;
extern template class PointND< walberla::real_t, 4 >;
extern template class PointND< walberla::real_t, 6 >;
extern template class PointND< walberla::real_t, 10 >;

inline Point3D crossProduct( const Point3D& a, const Point3D& b )
{
   return a.cross( b );
}

walberla::math::Vector3< real_t > toVec3( const Point3D& p );

Point3D toPoint3D( const walberla::math::Vector3< real_t >& v );

} // namespace hyteg

namespace walberla::mpi {

template < typename T, // Element type of SendBuffer
           typename G, // Growth policy of SendBuffer
           typename PointNDDataType,
           int PointNDDimension >
GenericSendBuffer< T, G >& operator<<( GenericSendBuffer< T, G >&                                 buf,
                                       const hyteg::PointND< PointNDDataType, PointNDDimension >& pointND )
{
   for ( int i = 0; i < PointNDDimension; ++i )
   {
      buf << pointND[i];
   }

   return buf;
}

template < typename T, // Element type  of RecvBuffer
           typename PointNDDataType,
           int PointNDDimension >
GenericRecvBuffer< T >& operator>>( GenericRecvBuffer< T >& buf, hyteg::PointND< PointNDDataType, PointNDDimension >& pointND )
{
   for ( int i = 0; i < PointNDDimension; ++i )
   {
      buf >> pointND[i];
   }
   return buf;
}

} // namespace walberla
