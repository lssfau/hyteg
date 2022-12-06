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

#include <array>
#include <cmath>
#include <iostream>

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
/// \tparam T Vector data type
/// \tparam N Dimension of vector
template < typename T, size_t N >
class PointND
{
 public:
   /// Constructor setting all components to zero
   PointND() { vector_.setZero(); }

   /// Constructs the vector using values from n-dimensional array \p _x
   /// \param _x Pointer to N-dimensional array
   PointND( T _x[N] )
   {
      for ( size_t i = 0; i < N; ++i )
      {
         vector_[i] = _x[i];
      }
   }

   /// Constructs the vector using values from n-dimensional array \p list, required for list initializer
   /// \param list N-dimensional array
   PointND( std::array< T, N > list )
   {
      for ( uint_t i = 0; i < N; ++i )
      {
         vector_[walberla::int_c( i )] = list[i];
      }
   }

   /// Computes the dot product between two PointND vectors
   /// \param b Right hand side of dot operator
   /// \returns Dot product between \p this and \p b
   T dot( const PointND& other ) const { return vector_.dot( other.vector_ ); }

   /// Computes the 2D normal of this PointND
   /// \returns 2D Point normal to this PointND
   PointND< T, 2 > normal2D() { return PointND< T, 2 >( { this->vector_[1], -this->vector_[0] } ); }

   /// Computes the squared Euclidean norm of \p this
   /// \returns Squared Euclidean norm of \p this
   T normSq() const { return vector_.dot( vector_ ); }

   /// Computes the Euclidean norm of \p this
   /// \returns Euclidean norm of \p this
   T norm() const { return vector_.norm(); }

   /// Add another PointND component wise to \p this
   /// \param rhs Reference to PointND that will be added to \p this
   /// \returns Reference to \p this after addition
   PointND& operator+=( const PointND& rhs )
   {
      vector_ += rhs.vector_;
      return *this;
   }

   /// Subtract another PointND component wise from \p this
   /// \param rhs Reference to PointND that will be subtracted from \p this
   /// \returns Reference to \p this after subtraction
   PointND& operator-=( const PointND& rhs )
   {
      vector_ -= rhs.vector_;
      return *this;
   }

   /// Multiply \p this with scalar value
   /// \param scalar Scalar value that \p this gets multiplied with
   /// \returns Reference to \p this after multiplication with \p scalar
   PointND& operator*=( T scalar )
   {
      vector_ *= scalar;
      return *this;
   }

   /// Divide \p this with scalar value
   /// \param scalar Scalar value that \p gets divided by
   /// \returns Reference to \p this after division by \p scalar
   PointND& operator/=( T scalar )
   {
      vector_ /= scalar;
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
      return vector_( walberla::int_c( index ) );
   }

   /// Value of component of vector at position \p index
   /// \param index The index of the component to access
   /// \returns Value of component at position \p index
   T operator[]( const uint_t index ) const
   {
      WALBERLA_ASSERT( index < N, "PointND index out of bounds: index = " << index << " but N = " << N );
      return vector_( walberla::int_c( index ) );
   }

   void setAll( const T& scalar ) { vector_.array() = scalar; }

#ifndef _MSC_VER
   void serialize( walberla::mpi::SendBuffer& sendBuffer ) const { sendBuffer << vector_; }

   void deserialize( walberla::mpi::RecvBuffer& recvBuffer ) { recvBuffer >> vector_; }
#else
   void serialize( walberla::mpi::SendBuffer& sendBuffer ) const
   {
      for ( int k = 0; k < N; ++k )
      {
         sendBuffer << vector_( k );
      }
   }

   void deserialize( walberla::mpi::RecvBuffer& recvBuffer )
   {
      for ( int k = 0; k < N; ++k )
      {
         recvBuffer >> vector_( k );
      }
   }
#endif

   Eigen::Matrix< T, static_cast< int >(N), 1 > vector_;
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
   Point3D cross;
   cross.vector_ = a.vector_.cross( b.vector_ );
   return cross;
}

walberla::math::Vector3< real_t > toVec3( const Point3D& p );

Point3D toPoint3D( const walberla::math::Vector3< real_t >& v );

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
