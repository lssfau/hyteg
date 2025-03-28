/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Andreas Burkhart.
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
#include "hyteg/types/PointND.hpp"
#include "hyteg/types/Matrix.hpp"

namespace hyteg {

/// Class providing geometry mapping based on polar coordinates; convention is x[0] = \f$\rho\f$, x[1] = \f$\varphi\f$
class PolarCoordsMap : public GeometryMap
{
 public:
   PolarCoordsMap() {}

   static void setMap( SetupPrimitiveStorage& setupStorage )
   {
      auto blend = std::make_shared< PolarCoordsMap >();

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

   void evalF( const Point3D& x, Point3D& Fx ) const override final
   {
      Fx[0] = x[0] * std::cos( x[1] );
      Fx[1] = x[0] * std::sin( x[1] );
   }

   void evalDF( const Point3D& x, Matrix2r& DFx ) const override final
   {
      DFx( 0, 0 ) = std::cos( x[1] );
      DFx( 0, 1 ) = -x[0] * std::sin( x[1] );
      DFx( 1, 0 ) = std::sin( x[1] );
      DFx( 1, 1 ) = x[0] * std::cos( x[1] );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const override final
   {
      DFinvx( 0, 0 ) = std::cos( x[1] );
      DFinvx( 0, 1 ) = std::sin( x[1] );
      DFinvx( 1, 0 ) = -std::sin( x[1] ) / x[0];
      DFinvx( 1, 1 ) = std::cos( x[1] ) / x[0];
   }

   void evalDFinvDF( const Point3D& x, Matrixr< 2, 4 >& DFinvDFx ) const override final
   {
      const real_t tmp0 = std::sin( x[1] );
      const real_t tmp1 = ( tmp0 * tmp0 );
      const real_t tmp2 = std::cos( x[1] );
      const real_t tmp3 = ( tmp2 * tmp2 );
      const real_t tmp4 = tmp1 * x[0] + tmp3 * x[0];
      const real_t tmp5 = tmp1 + tmp3;
      const real_t tmp6 = -tmp5 / ( tmp4 * tmp4 );
      const real_t tmp7 = real_c( 1 ) / ( tmp5 );
      const real_t tmp8 = real_c( 1 ) / ( tmp4 );
      DFinvDFx( 0, 0 )  = real_c( 0 );
      DFinvDFx( 0, 1 )  = real_c( 0 );
      DFinvDFx( 1, 0 )  = -tmp0 * tmp6;
      DFinvDFx( 1, 1 )  = tmp2 * tmp6;
      DFinvDFx( 0, 2 )  = -tmp0 * tmp7;
      DFinvDFx( 0, 3 )  = tmp2 * tmp7;
      DFinvDFx( 1, 2 )  = -tmp2 * tmp8;
      DFinvDFx( 1, 3 )  = -tmp0 * tmp8;
   };

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override final { sendBuffer << Type::POLAR_COORDS; }
};

} // namespace hyteg
