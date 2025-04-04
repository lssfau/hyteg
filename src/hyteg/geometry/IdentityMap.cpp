/*
 * Copyright (c) 2024 Marcus Mohr.
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

#include "hyteg/geometry/IdentityMap.hpp"

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/types/Matrix.hpp"

namespace hyteg {
void IdentityMap::evalF( const Point3D& x, Point3D& Fx ) const
{
   Fx = x;
}

void IdentityMap::evalFinv( const Point3D& xPhys, Point3D& xComp ) const
{
   xComp = xPhys;
}

void IdentityMap::evalDF( const Point3D&, Matrix2r& DFx ) const
{
   DFx( 0, 0 ) = real_c( 1.0 );
   DFx( 0, 1 ) = real_c( 0.0 );
   DFx( 1, 0 ) = real_c( 0.0 );
   DFx( 1, 1 ) = real_c( 1.0 );
}

real_t IdentityMap::evalDF( const Point3D& point_3d, Matrix3r& DFx ) const
{
   DFx( 0, 0 ) = real_c( 1.0 );
   DFx( 0, 1 ) = real_c( 0.0 );
   DFx( 0, 2 ) = real_c( 0.0 );

   DFx( 1, 0 ) = real_c( 0.0 );
   DFx( 1, 1 ) = real_c( 1.0 );
   DFx( 1, 2 ) = real_c( 0.0 );

   DFx( 2, 0 ) = real_c( 0.0 );
   DFx( 2, 1 ) = real_c( 0.0 );
   DFx( 2, 2 ) = real_c( 1.0 );

   return real_c( 1.0 );
}

void IdentityMap::evalDFinv( const Point3D&, Matrix2r& DFinvx ) const
{
   DFinvx( 0, 0 ) = real_c( 1.0 );
   DFinvx( 0, 1 ) = real_c( 0.0 );
   DFinvx( 1, 0 ) = real_c( 0.0 );
   DFinvx( 1, 1 ) = real_c( 1.0 );
}

void IdentityMap::evalDFinvDF( const Point3D& x, Matrixr< 2, 4 >& DFinvDFx ) const
{
   WALBERLA_UNUSED( x );
   DFinvDFx.setZero();
}

void IdentityMap::evalDFinvDF( const Point3D& x, Matrixr< 3, 9 >& DFinvDFx ) const
{
   WALBERLA_UNUSED( x );
   DFinvDFx.setZero();
}

void IdentityMap::serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const
{
   sendBuffer << Type::IDENTITY;
}

void IdentityMap::setMap( SetupPrimitiveStorage& setupStorage )
{
   auto blendingMap = std::make_shared< IdentityMap >();

   for ( auto it : setupStorage.getCells() )
   {
      Cell& cell = *it.second;
      setupStorage.setGeometryMap( cell.getID(), blendingMap );
   }

   for ( auto it : setupStorage.getFaces() )
   {
      Face& face = *it.second;
      setupStorage.setGeometryMap( face.getID(), blendingMap );
   }

   for ( auto it : setupStorage.getEdges() )
   {
      Edge& edge = *it.second;
      setupStorage.setGeometryMap( edge.getID(), blendingMap );
   }

   for ( auto it : setupStorage.getVertices() )
   {
      Vertex& vertex = *it.second;
      setupStorage.setGeometryMap( vertex.getID(), blendingMap );
   }
}
} // namespace hyteg
