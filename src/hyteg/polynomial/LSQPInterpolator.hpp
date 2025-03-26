/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Benjamin Mann.
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

#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/types/Matrix.hpp"

#include "Polynomial.hpp"

namespace hyteg {

using walberla::uint_c;

enum class LSQPType
{
   EDGE,    // one interpolation point on every edge of one type (HO,VE or DI) connected to a vertex.
   VERTEX,  // one interpolation point on every vertex.
   EDGE_ALL // one interpolation point on every edge of one type (HO,VE or DI).
};

template < LSQPType Type >
constexpr uint_t GetNumInterpolationPoints( uint_t level )
{
   switch ( Type )
   {
   case LSQPType::EDGE:
      return ( levelinfo::num_microedges_per_face( level ) - 3 * levelinfo::num_microedges_per_edge( level ) - 3 ) / 3;

   case LSQPType::VERTEX:
      return levelinfo::num_microvertices_per_face( level ) - 3 * ( levelinfo::num_microvertices_per_edge( level ) - 2 ) - 3;

   case LSQPType::EDGE_ALL:
      return ( levelinfo::num_microedges_per_face( level ) - 3 * levelinfo::num_microedges_per_edge( level ) ) / 3;

   default:
      return 0;
   }
}

template < LSQPType Type >
constexpr uint_t GetNumInterpolationPoints3D( uint_t level )
{
   const uint_t ip_on_edge = levelinfo::num_microvertices_per_edge( level ) - 2;
   const uint_t ip_on_face = levelinfo::num_microvertices_per_face( level ) - 3 * ip_on_edge - 3;

   switch ( Type )
   {
   case LSQPType::VERTEX:
      return levelinfo::num_microvertices_per_cell( level ) - 4 * ip_on_face - 6 * ip_on_edge - 4;

   // todo other node types
   default:
      return 0;
      // WALBERLA_ABORT("3D Polynomial Interpolation only implemented for P1 elements");
   }
}

template < typename Basis, LSQPType Type >
class LSQPInterpolator
{
 public:
   LSQPInterpolator( uint_t degree )
   : degree_( degree )
   , numCoefficients_( Polynomial2D< Basis >::getNumCoefficients( degree ) )
   , offset_( 0 )
   , A( 0, Polynomial2D< Basis >::getNumCoefficients( degree ) )
   , rhs( 0, 1 )
   {}

   void addInterpolationPoint( const Point2D& x, real_t value )
   {
      A.conservativeResize( A.rows() + 1, Eigen::NoChange );
      rhs.conservativeResize( rhs.rows() + 1, Eigen::NoChange );

      for ( uint_t k = 0; k < numCoefficients_; ++k )
      {
         A( offset_, k ) = Basis::eval( k, x );
      }

      rhs( offset_ ) = value;

      ++offset_;
   }

   void interpolate( Polynomial2D< Basis >& poly )
   {
      WALBERLA_ASSERT( degree_ == poly.getDegree(), "Polynomial degrees don't match!" );

      MatrixXr coeffs;
      coeffs = A.colPivHouseholderQr().solve( rhs );

      for ( uint_t i = 0; i < numCoefficients_; ++i )
      {
         poly.setCoefficient( i, coeffs( i ) );
      }
   }

   uint_t numInterpolationPoints() const { return uint_c( A.rows() ); }

 private:
   uint_t   degree_;
   uint_t   numCoefficients_;
   uint_t   offset_;
   MatrixXr A;
   MatrixXr rhs;
};

template < typename Basis, LSQPType Type >
class LSQPInterpolator3D
{
 public:
   LSQPInterpolator3D( uint_t degree )
   : degree_( degree )
   , numCoefficients_( Polynomial3D< Basis >::getNumCoefficients( degree ) )
   , offset_( 0 )
   , A( 0, numCoefficients_ )
   , rhs( 0, 1 )
   {}

   void addInterpolationPoint( const Point3D& x, const real_t& value )
   {
      A.conservativeResize( A.rows() + 1, Eigen::NoChange );
      rhs.conservativeResize( rhs.rows() + 1, Eigen::NoChange );

      for ( uint_t k = 0; k < numCoefficients_; ++k )
      {
         A( offset_, k ) = Basis::eval( k, x );
      }

      rhs( offset_ ) = value;

      ++offset_;
   }

   void interpolate( Polynomial3D< Basis >& poly )
   {
      WALBERLA_ASSERT( degree_ == poly.getDegree(), "Polynomial degrees don't match!" );

      MatrixXr coeffs;
      coeffs = A.colPivHouseholderQr().solve( rhs );

      for ( uint_t i = 0; i < numCoefficients_; ++i )
      {
         poly.setCoefficient( i, coeffs( i ) );
      }
   }

   uint_t numInterpolationPoints() const { return uint_c( A.rows() ); }

 private:
   uint_t   degree_;
   uint_t   numCoefficients_;
   uint_t   offset_;
   MatrixXr A;
   //todo factor out QR-decomposition
   //  Eigen::ColPivHouseholderQR<typename MatrixXr::PlainObject> QR;
   MatrixXr rhs;
};
} // namespace hyteg
