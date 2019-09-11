/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
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

#include "core/math/Utility.h"

#include "GeometryMap.hpp"

namespace hyteg {

class AffineMap : public GeometryMap
{
 public:
   AffineMap( const std::array< Point3D, 3 >& fromCoords, const std::array< Point3D, 3 >& toCoords )
   {
      x1_ = fromCoords[0];
      x2_ = fromCoords[1];
      x3_ = fromCoords[2];

      xb1_ = toCoords[0];
      xb2_ = toCoords[1];
      xb3_ = toCoords[2];
   }

   AffineMap( walberla::mpi::RecvBuffer& recvBuffer )
   {
      recvBuffer >> x1_;
      recvBuffer >> x2_;
      recvBuffer >> x3_;
      recvBuffer >> xb1_;
      recvBuffer >> xb2_;
      recvBuffer >> xb3_;
   }

   void evalF( const Point3D& x, Point3D& Fx ) const
   {
      real_t tmp0 = x1_[0] - x2_[0];
      real_t tmp1 = x1_[1] - x3_[1];
      real_t tmp2 = x1_[0] - x3_[0];
      real_t tmp3 = x1_[1] - x2_[1];
      real_t tmp4 = tmp0 * tmp1 - tmp2 * tmp3;
      real_t tmp5 = 1.0 / tmp4;
      real_t tmp6 = x[0] - x1_[0];
      real_t tmp7 = x[1] - x1_[1];
      real_t tmp8 = tmp1 * tmp6 - tmp2 * tmp7;
      real_t tmp9 = -tmp0 * tmp7 + tmp3 * tmp6;
      Fx[0]       = tmp5 * ( tmp4 * xb1_[0] + tmp8 * ( xb1_[0] - xb2_[0] ) - tmp9 * ( xb1_[0] - xb3_[0] ) );
      Fx[1]       = tmp5 * ( tmp4 * xb1_[1] + tmp8 * ( xb1_[1] - xb2_[1] ) - tmp9 * ( xb1_[1] - xb3_[1] ) );
   }

   void evalDF( const Point3D&, Matrix2r& DFx ) const
   {
      real_t tmp0 = xb1_[0] - xb2_[0];
      real_t tmp1 = x1_[1] - x3_[1];
      real_t tmp2 = xb1_[0] - xb3_[0];
      real_t tmp3 = x1_[1] - x2_[1];
      real_t tmp4 = x1_[0] - x2_[0];
      real_t tmp5 = x1_[0] - x3_[0];
      real_t tmp6 = 1.0 / ( tmp1 * tmp4 - tmp3 * tmp5 );
      real_t tmp7 = xb1_[1] - xb2_[1];
      real_t tmp8 = xb1_[1] - xb3_[1];
      DFx( 0, 0 ) = tmp6 * ( tmp0 * tmp1 - tmp2 * tmp3 );
      DFx( 0, 1 ) = tmp6 * ( -tmp0 * tmp5 + tmp2 * tmp4 );
      DFx( 1, 0 ) = tmp6 * ( tmp1 * tmp7 - tmp3 * tmp8 );
      DFx( 1, 1 ) = tmp6 * ( tmp4 * tmp8 - tmp5 * tmp7 );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFxInv ) const
   {
      Matrix2r tmp;
      evalDF( x, tmp );
      real_t invDet = 1.0 / tmp.det();
      DFxInv        = tmp.adj();
      DFxInv *= invDet;
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const
   {
      sendBuffer << Type::AFFINE;
      sendBuffer << x1_;
      sendBuffer << x2_;
      sendBuffer << x3_;
      sendBuffer << xb1_;
      sendBuffer << xb2_;
      sendBuffer << xb3_;
   }

 private:
   Point3D x1_;
   Point3D x2_;
   Point3D x3_;
   Point3D xb1_;
   Point3D xb2_;
   Point3D xb3_;
};

} // namespace hyteg
