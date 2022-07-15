/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <core/math/Vector3.h>
#include <iostream>

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/eigen/EigenWrapper.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// \brief  N-dimensional vector
/// \author Daniel Drzisga (drzisga@ma.tum.de)
/// \date   March, 2017
///
/// The PointND class represents an N-dimensional vector with support for algebraic operations
/// \tparam T Vector data type
/// \tparam N Dimension of vector
template < typename T, size_t N >
class PointND
{
 public:
   /// Constructor setting all components to zero
   PointND()
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] = (T) 0;
      }
   }

   /// Constructs the vector using values from n-dimensional array \p _x
   /// \param _x Pointer to N-dimensional array
   PointND( T _x[N] )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] = _x[i];
      }
   }

   /// Constructs the vector using values from n-dimensional array \p list, required for list initializer
   /// \param list N-dimensional array
   PointND( std::array< T, N > list )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] = list[i];
      }
   }

   /// Copy constructor
   /// \param b Reference to another instance of PointND
   PointND( const PointND& b )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] = b.x_[i];
      }
   }

   /// Computes the dot product between two PointND vectors
   /// \param b Right hand side of dot operator
   /// \returns Dot product between \p this and \p b
   T dot( const PointND& b ) const
   {
      T tmp = 0.0;
      for ( size_t i = 0; i < N; ++i )
      {
         tmp += x_[i] * b.x_[i];
      }
      return tmp;
   }

   /// Computes the 2D normal of this PointND
   /// \returns 2D Point normal to this PointND
   PointND< T, 2 > normal2D() { return PointND< T, 2 >( { this->x_[1], -this->x_[0] } ); }

   /// Computes the squared Euclidean norm of \p this
   /// \returns Squared Euclidean norm of \p this
   T normSq() const { return dot( *this ); }

   /// Computes the Euclidean norm of \p this
   /// \returns Euclidean norm of \p this
   T norm() const { return std::sqrt( normSq() ); }

   /// Add another PointND component wise to \p this
   /// \param rhs Reference to PointND that will be added to \p this
   /// \returns Reference to \p this after addition
   PointND& operator+=( const PointND& rhs )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] += rhs.x_[i];
      }
      return *this;
   }

   /// Subtract another PointND component wise from \p this
   /// \param rhs Reference to PointND that will be subtracted from \p this
   /// \returns Reference to \p this after subtraction
   PointND& operator-=( const PointND& rhs )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] -= rhs.x_[i];
      }
      return *this;
   }

   /// Multiply \p this with scalar value
   /// \param scalar Scalar value that \p this gets multiplied with
   /// \returns Reference to \p this after multiplication with \p scalar
   PointND& operator*=( T scalar )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] *= scalar;
      }
      return *this;
   }

   /// Divide \p this with scalar value
   /// \param scalar Scalar value that \p gets divided by
   /// \returns Reference to \p this after division by \p scalar
   PointND& operator/=( T scalar )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] /= scalar;
      }
      return *this;
   }

   /// Return negated copy of \p this
   /// \returns Copy of \p this with negated components
   PointND operator-() const { return static_cast< T >( -1 ) * ( *this ); }

   /// Reference to component of vector at position \p index
   /// \param index The index of the component to access
   /// \returns Reference to component at position \p index
   T& operator[]( const uint_t index )
   {
      WALBERLA_ASSERT( index < N, "PointND index out of bounds: index = " << index << " but N = " << N );
      return x_[index];
   }

   /// Value of component of vector at position \p index
   /// \param index The index of the component to access
   /// \returns Value of component at position \p index
   T operator[]( const uint_t index ) const
   {
      WALBERLA_ASSERT( index < N, "PointND index out of bounds: index = " << index << " but N = " << N );
      return x_[index];
   }

   void setAll( const T& scalar )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         x_[i] = scalar;
      }
   }

   void serialize( walberla::mpi::SendBuffer& sendBuffer ) const
   {
      for ( size_t index = 0; index < N; index++ )
      {
         sendBuffer << x_[index];
      }
   }

   void deserialize( walberla::mpi::RecvBuffer& recvBuffer )
   {
      for ( size_t index = 0; index < N; index++ )
      {
         recvBuffer >> x_[index];
      }
   }

 protected:
   T x_[N];
};

template < typename T, size_t N >
inline PointND< T, N > operator+( PointND< T, N > lhs, const PointND< T, N >& rhs )
{
   return lhs += rhs;
}

template < typename T, size_t N >
inline PointND< T, N > operator-( PointND< T, N > lhs, const PointND< T, N >& rhs )
{
   return lhs -= rhs;
}

template < typename T, size_t N >
inline PointND< T, N > operator*( real_t scalar, PointND< T, N > rhs )
{
   rhs *= scalar;
   return rhs;
}

template < typename T, size_t N >
inline PointND< T, N > operator*( PointND< T, N > lhs, real_t scalar )
{
   lhs *= scalar;
   return lhs;
}

template < typename T, size_t N >
inline PointND< T, N > operator/( PointND< T, N > lhs, real_t scalar )
{
   lhs /= scalar;
   return lhs;
}

template < typename T, size_t N >
inline std::ostream& operator<<( std::ostream& os, const PointND< T, N >& pointnd )
{
   os << "[";

   for ( size_t i = 0; i < N; ++i )
   {
      os << pointnd[i];
      if ( i != N - 1 )
      {
         os << ", ";
      }
   }

   os << "]";

   return os;
}

typedef PointND< real_t, 2 >  Point2D;
typedef PointND< real_t, 3 >  Point3D;
typedef PointND< real_t, 4 >  Point4D;
typedef PointND< real_t, 6 >  Point6D;
typedef PointND< real_t, 10 > Point10D;

template < size_t N >
using PointNDr = PointND< real_t, N >;

inline Point3D crossProduct( const Point3D& a, const Point3D& b )
{
   return Point3D( { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] } );
}

inline walberla::math::Vector3< real_t > toVec3( const Point3D& p )
{
   return walberla::math::Vector3< real_t >( p[0], p[1], p[2] );
}

inline Point3D toPoint3D( const walberla::math::Vector3< real_t >& v )
{
   return Point3D( { v[0], v[1], v[2] } );
}

template < typename T, size_t N >
inline Eigen::Matrix< T, N, 1 > toEigen( const PointND< T, N >& p )
{
   Eigen::Matrix< T, N, 1 > eigen;
   for ( size_t i = 0; i < N; ++i )
   {
      eigen[i] = p[i];
   }
   return eigen;
}

} // namespace hyteg

namespace walberla {
namespace mpi {

template < typename T, // Element type of SendBuffer
           typename G, // Growth policy of SendBuffer
           typename PointNDDataType,
           size_t PointNDDimension >
GenericSendBuffer< T, G >& operator<<( GenericSendBuffer< T, G >&                                 buf,
                                       const hyteg::PointND< PointNDDataType, PointNDDimension >& pointND )
{
   pointND.serialize( buf );
   return buf;
}

template < typename T, // Element type  of RecvBuffer
           typename PointNDDataType,
           size_t PointNDDimension >
GenericRecvBuffer< T >& operator>>( GenericRecvBuffer< T >& buf, hyteg::PointND< PointNDDataType, PointNDDimension >& pointND )
{
   pointND.deserialize( buf );
   return buf;
}

} // namespace mpi
} // namespace walberla

