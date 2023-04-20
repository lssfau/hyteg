/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
#include "AffineMap2D.hpp"
#include "AffineMap3D.hpp"
#include "AnnulusMap.hpp"
#include "CircularMap.hpp"
#include "IcosahedralShellMap.hpp"
#include "IdentityMap.hpp"
#include "PolarCoordsMap.hpp"
#include "ThinShellMap.hpp"

namespace hyteg {

real_t GeometryMap::evalDetDF( const Point3D& x )
{
   Matrix2r DF;
   evalDF( x, DF );
   return DF.determinant();
}

void GeometryMap::serialize( const std::shared_ptr< GeometryMap >& map, walberla::mpi::SendBuffer& sendBuffer )
{
   map->serializeSubClass( sendBuffer );
}

std::shared_ptr< GeometryMap > GeometryMap::deserialize( walberla::mpi::RecvBuffer& recvBuffer )
{
   Type type;
   recvBuffer >> type;

   switch ( type )
   {
   case Type::AFFINE:
   case Type::AFFINE_2D:
      return std::make_shared< AffineMap2D >( recvBuffer );
   case Type::AFFINE_3D:
      return std::make_shared< AffineMap3D >( recvBuffer );
   case Type::ANNULUS:
      return std::make_shared< AnnulusMap >( recvBuffer );
   case Type::CIRCULAR:
      return std::make_shared< CircularMap >( recvBuffer );
   case Type::IDENTITY:
      return std::make_shared< IdentityMap >();
   case Type::ICOSAHEDRAL_SHELL:
      return std::make_shared< IcosahedralShellMap >( recvBuffer );
   case Type::THIN_SHELL:
      return std::make_shared< ThinShellMap >( recvBuffer );
   case Type::POLAR_COORDS:
      return std::make_shared< PolarCoordsMap >();
   default:
      WALBERLA_ABORT( "Error in deserializing GeometryMap: Unsupported Map Type" )
   }
}

bool GeometryMap::verifyPointPairing( const Point3D& computationalCoordinates,
                                      const Point3D& physicalCoordinates,
                                      real_t         threshold ) const
{
   Point3D mapped;
   this->evalF( computationalCoordinates, mapped );
   return ( mapped - physicalCoordinates ).norm() < threshold;
}

} // namespace hyteg
