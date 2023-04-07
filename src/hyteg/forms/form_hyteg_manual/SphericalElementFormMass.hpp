/*
 * Copyright (c) 2022 Marcus Mohr.
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

#include "hyteg/forms/P1Form.hpp"

using walberla::real_c;
using walberla::real_t;

namespace hyteg {

/// Form for computing the element mass matrix for elements on spherical triangles
///
/// See hyteg::SphericalElementFunction< ValueType > for info on these elements.
class SphericalElementFormMass : public P1Form
{
 public:
   SphericalElementFormMass( real_t radius = real_c( 1 ) )
   : radius_( radius )
   {
      WALBERLA_ASSERT( radius > real_c( 0 ) );
   };

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix3r& elMat ) const override
   {
      Point3D v0, v1, v2;
      geometryMap_->evalF( coords[0], v0 );
      geometryMap_->evalF( coords[1], v1 );
      geometryMap_->evalF( coords[2], v2 );

      real_t factor = real_c( 1 ) / ( radius_ * radius_ );

      // compute edge lengths of spherical triangle
      real_t a0 = std::acos( v0.dot( v1 ) * factor );
      real_t a1 = std::acos( v0.dot( v2 ) * factor );
      real_t a2 = std::acos( v1.dot( v2 ) * factor );

      // compute spherical excess using L'Huilier's formula
      a0 *= real_c( 0.5 );
      a1 *= real_c( 0.5 );
      a2 *= real_c( 0.5 );
      real_t aux    = real_c( 0.5 ) * ( a0 + a1 + a2 );
      aux           = std::tan( aux ) * std::tan( aux - a0 ) * std::tan( aux - a1 ) * std::tan( aux - a2 );
      real_t excess = real_c( 4 ) * std::atan( std::sqrt( aux ) );

      // compute factor with is the scaled area of spherical triangle (using Girard's theorem)
      factor = excess * ( radius_ * radius_ ) / real_c( 12 );

      // assemble element mass matrix
      elMat( 0, 0 ) = factor * real_c( 2 );
      elMat( 1, 0 ) = factor * real_c( 1 );
      elMat( 2, 0 ) = factor * real_c( 1 );

      elMat( 0, 1 ) = factor * real_c( 1 );
      elMat( 1, 1 ) = factor * real_c( 2 );
      elMat( 2, 1 ) = factor * real_c( 1 );

      elMat( 0, 2 ) = factor * real_c( 1 );
      elMat( 1, 2 ) = factor * real_c( 1 );
      elMat( 2, 2 ) = factor * real_c( 2 );
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const { WALBERLA_ABORT( "Not implemented." ); }

   void setRadius( real_t radius )
   {
      WALBERLA_ASSERT( radius > real_c( 0 ) );
      radius_ = radius;
   }

 private:
   real_t radius_;
};

} // namespace hyteg
