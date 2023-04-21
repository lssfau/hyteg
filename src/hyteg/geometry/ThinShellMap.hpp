/*
 * Copyright (c) 2022-2023 Marcus Mohr.
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

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "GeometryMap.hpp"

namespace hyteg {

using walberla::real_c;

/// Class providing geometry mapping for a thin spherical shell
///
/// This geometry map provides a blending operation for any kind of base mesh
/// representing a thin spherical shell. Geometric nodes on refined hierarchy
/// levels are projected onto the thin spherical shell.
class ThinShellMap : public GeometryMap
{
 public:
   ThinShellMap( const Face& face, real_t radius )
   : radius_( radius )
   {
      primitiveVertices.reserve( 3 );
      const auto coords = face.getCoordinates();
      for ( const auto& point : coords )
      {
         primitiveVertices.push_back( point );
      }
   }

   ThinShellMap( const Edge& edge, real_t radius )
   : radius_( radius )
   {
      primitiveVertices.reserve( 2 );
      const auto coords = edge.getCoordinates();
      for ( const auto& point : coords )
      {
         primitiveVertices.push_back( point );
      }
   }

   ThinShellMap( walberla::mpi::RecvBuffer& recvBuffer ) { recvBuffer >> radius_ >> primitiveVertices; }

   void evalF( const Point3D& xComp, Point3D& xPhys ) const override
   {
      real_t oldRad = std::sqrt( xComp[0] * xComp[0] + xComp[1] * xComp[1] + xComp[2] * xComp[2] );

      xPhys[0] = xComp[0] * radius_ / oldRad;
      xPhys[1] = xComp[1] * radius_ / oldRad;
      xPhys[2] = xComp[2] * radius_ / oldRad;
   }

   void evalFinv( const Point3D& xPhys, Point3D& xComp ) const override
   {
      WALBERLA_UNUSED( xPhys );
      WALBERLA_UNUSED( xComp );
      WALBERLA_ABORT( "ThinShellMap::evalFinv() not implemented, yet!" );
   }

   void evalDF( const Point3D&, Matrix2r& ) const final override
   {
      WALBERLA_ABORT( "ThinShellMap::evalDF unimplemented for 2D!" );
   }

   void evalDFinv( const Point3D&, Matrix2r& ) const final override
   {
      WALBERLA_ABORT( "ThinShellMap::evalDFinv unimplemented for 2D!" );
   }

   real_t evalDF( const Point3D& x, Matrix3r& DFx ) const final override
   {
      real_t rad    = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      real_t factor = radius_ / ( rad * rad * rad );

      DFx( 0, 0 ) = x[1] * x[1] + x[2] * x[2];
      DFx( 0, 1 ) = -x[0] * x[1];
      DFx( 0, 2 ) = -x[0] * x[2];

      DFx( 1, 0 ) = -x[0] * x[1];
      DFx( 1, 1 ) = x[0] * x[0] + x[2] * x[2];
      DFx( 1, 2 ) = -x[1] * x[2];

      DFx( 2, 0 ) = -x[0] * x[2];
      DFx( 2, 1 ) = -x[1] * x[2];
      DFx( 2, 2 ) = x[0] * x[0] + x[1] * x[1];

      DFx *= factor;

      return DFx( 0, 0 ) * DFx( 1, 1 ) * DFx( 2, 2 ) - DFx( 0, 0 ) * DFx( 2, 1 ) * DFx( 1, 2 ) -
             DFx( 1, 0 ) * DFx( 0, 1 ) * DFx( 2, 2 ) + DFx( 1, 0 ) * DFx( 2, 1 ) * DFx( 0, 2 ) +
             DFx( 2, 0 ) * DFx( 0, 1 ) * DFx( 1, 2 ) - DFx( 2, 0 ) * DFx( 1, 1 ) * DFx( 0, 2 );
   };

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override
   {
      sendBuffer << Type::THIN_SHELL << radius_ << primitiveVertices;
   }

   static void setMap( SetupPrimitiveStorage& setupStorage, real_t radius )
   {
      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< ThinShellMap >( face, radius ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< ThinShellMap >( edge, radius ) );
      }
   }

 private:
   real_t                 radius_;
   std::vector< Point3D > primitiveVertices; // macro-vertex coordinates of the map's primitive for evalFinv()
};

} // namespace hyteg
