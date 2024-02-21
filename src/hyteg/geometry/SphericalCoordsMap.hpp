/*
 * Copyright (c) 2023 Andreas Burkhart.
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

#include <cmath>

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

/// Class providing geometry mapping based on spherical coordinates, x[0] = r, x[1] = theta, x[2] = phi
class SphericalCoordsMap : public GeometryMap
{
 public:
   SphericalCoordsMap() {}

   static void setMap( SetupPrimitiveStorage& setupStorage )
   {
      auto blend = std::make_shared< SphericalCoordsMap >();

      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(), blend );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), blend );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), blend );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap( vertex.getID(), blend );
      }
   }

   void evalDF( const Point3D& x, Matrix2r& DFx ) const final override
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFx );
      WALBERLA_ABORT( "SphericalCoordsMap::evalDF unimplemented for 2D!" );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const final override
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvx );
      WALBERLA_ABORT( "SphericalCoordsMap::evalDFinv unimplemented for 2D!" );
   }

   //x[0] = r, x[1] = theta, x[2] = phi
   void evalF( const Point3D& x, Point3D& Fx ) const override final
   {
      Fx[0] = x[0] * std::sin( x[1] ) * std::cos( x[2] );
      Fx[1] = x[0] * std::sin( x[1] ) * std::sin( x[2] );
      Fx[2] = x[0] * std::cos( x[1] );
   }

   real_t evalDF( const Point3D& x, Matrix3r& DFx ) const override final
   {
      DFx( 0, 0 ) = std::sin( x[1] ) * std::cos( x[2] );
      DFx( 0, 1 ) = x[0] * std::cos( x[1] ) * std::cos( x[2] );
      DFx( 0, 2 ) = -x[0] * std::sin( x[1] ) * std::sin( x[2] );

      DFx( 1, 0 ) = std::sin( x[1] ) * std::sin( x[2] );
      DFx( 1, 1 ) = x[0] * std::cos( x[1] ) * std::sin( x[2] );
      DFx( 1, 2 ) = x[0] * std::sin( x[1] ) * std::cos( x[2] );

      DFx( 2, 0 ) = std::cos( x[1] );
      DFx( 2, 1 ) = -x[0] * std::sin( x[1] );
      DFx( 2, 2 ) = 0;

      return DFx( 0, 0 ) * DFx( 1, 1 ) * DFx( 2, 2 ) - DFx( 0, 0 ) * DFx( 2, 1 ) * DFx( 1, 2 ) -
             DFx( 1, 0 ) * DFx( 0, 1 ) * DFx( 2, 2 ) + DFx( 1, 0 ) * DFx( 2, 1 ) * DFx( 0, 2 ) +
             DFx( 2, 0 ) * DFx( 0, 1 ) * DFx( 1, 2 ) - DFx( 2, 0 ) * DFx( 1, 1 ) * DFx( 0, 2 );
   }

   void evalDFinvDF( const Point3D& x, Matrixr<3,9>& DFinvDFx ) const override final
   {
      const real_t tmp0 = sin(x[2]);
      const real_t tmp1 = (tmp0*tmp0);
      const real_t tmp2 = (x[0]*x[0]);
      const real_t tmp3 = sin(x[1]);
      const real_t tmp4 = (tmp3*tmp3*tmp3);
      const real_t tmp5 = tmp2*tmp4;
      const real_t tmp6 = cos(x[2]);
      const real_t tmp7 = (tmp6*tmp6);
      const real_t tmp8 = cos(x[1]);
      const real_t tmp9 = (tmp8*tmp8);
      const real_t tmp10 = tmp2*tmp9;
      const real_t tmp11 = tmp1*tmp10*tmp3 + tmp1*tmp5 + tmp10*tmp3*tmp7 + tmp5*tmp7;
      const real_t tmp12 = 1.0 / (tmp11);
      const real_t tmp13 = (tmp3*tmp3);
      const real_t tmp14 = tmp13*tmp6;
      const real_t tmp15 = tmp14*x[0];
      const real_t tmp16 = tmp12*tmp15;
      const real_t tmp17 = 1.0 / (tmp11*tmp11);
      const real_t tmp18 = 2*x[0];
      const real_t tmp19 = tmp18*tmp4;
      const real_t tmp20 = tmp3*tmp9;
      const real_t tmp21 = tmp1*tmp18;
      const real_t tmp22 = tmp18*tmp7;
      const real_t tmp23 = tmp17*(-tmp1*tmp19 - tmp19*tmp7 - tmp20*tmp21 - tmp20*tmp22);
      const real_t tmp24 = tmp2*tmp23;
      const real_t tmp25 = tmp0*tmp13;
      const real_t tmp26 = tmp25*x[0];
      const real_t tmp27 = tmp12*tmp26;
      const real_t tmp28 = tmp3*tmp8;
      const real_t tmp29 = tmp21*tmp28 + tmp22*tmp28;
      const real_t tmp30 = tmp2*tmp28;
      const real_t tmp31 = tmp1*tmp30 + tmp30*tmp7;
      const real_t tmp32 = tmp12*tmp28;
      const real_t tmp33 = tmp32*tmp6;
      const real_t tmp34 = tmp28*x[0];
      const real_t tmp35 = tmp23*tmp34;
      const real_t tmp36 = tmp0*tmp32;
      const real_t tmp37 = tmp1*tmp13;
      const real_t tmp38 = tmp13*tmp7;
      const real_t tmp39 = -tmp37*x[0] - tmp38*x[0];
      const real_t tmp40 = tmp0*tmp9;
      const real_t tmp41 = tmp40*x[0];
      const real_t tmp42 = -tmp26 - tmp41;
      const real_t tmp43 = tmp6*tmp9;
      const real_t tmp44 = tmp43*x[0];
      const real_t tmp45 = tmp15 + tmp44;
      const real_t tmp46 = 2*tmp12*tmp30;
      const real_t tmp47 = tmp2*(tmp8*tmp8*tmp8);
      const real_t tmp48 = tmp2*tmp37;
      const real_t tmp49 = tmp2*tmp38;
      const real_t tmp50 = tmp17*(-tmp1*tmp47 - tmp47*tmp7 - tmp48*tmp8 - tmp49*tmp8);
      const real_t tmp51 = tmp2*tmp50;
      const real_t tmp52 = tmp34*tmp50;
      const real_t tmp53 = tmp12*tmp2;
      DFinvDFx(0,0) = tmp14*tmp24 + 2*tmp16;
      DFinvDFx(0,1) = tmp24*tmp25 + 2*tmp27;
      DFinvDFx(0,2) = tmp12*tmp29 + tmp23*tmp31;
      DFinvDFx(1,0) = tmp33 + tmp35*tmp6;
      DFinvDFx(1,1) = tmp0*tmp35 + tmp36;
      DFinvDFx(1,2) = tmp12*(-tmp37 - tmp38) + tmp23*tmp39;
      DFinvDFx(2,0) = tmp12*(-tmp25 - tmp40) + tmp23*tmp42;
      DFinvDFx(2,1) = tmp12*(tmp14 + tmp43) + tmp23*tmp45;
      DFinvDFx(2,2) = 0;
      DFinvDFx(0,3) = tmp14*tmp51 + tmp46*tmp6;
      DFinvDFx(0,4) = tmp0*tmp46 + tmp25*tmp51;
      DFinvDFx(0,5) = tmp12*(tmp1*tmp2*tmp9 + tmp2*tmp7*tmp9 - tmp48 - tmp49) + tmp31*tmp50;
      DFinvDFx(1,3) = tmp12*tmp44 - tmp16 + tmp52*tmp6;
      DFinvDFx(1,4) = tmp0*tmp52 + tmp12*tmp41 - tmp27;
      DFinvDFx(1,5) = -tmp12*tmp29 + tmp39*tmp50;
      DFinvDFx(2,3) = tmp42*tmp50;
      DFinvDFx(2,4) = tmp45*tmp50;
      DFinvDFx(2,5) = 0;
      DFinvDFx(0,6) = -tmp25*tmp53;
      DFinvDFx(0,7) = tmp14*tmp53;
      DFinvDFx(0,8) = 0;
      DFinvDFx(1,6) = -tmp36*x[0];
      DFinvDFx(1,7) = tmp33*x[0];
      DFinvDFx(1,8) = 0;
      DFinvDFx(2,6) = -tmp12*tmp45;
      DFinvDFx(2,7) = tmp12*tmp42;
      DFinvDFx(2,8) = 0;
   }; 

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override final { sendBuffer << Type::SPHERICAL_COORDS; }
};

} // namespace hyteg
