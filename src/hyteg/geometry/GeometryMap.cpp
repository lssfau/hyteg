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

#include "hyteg/geometry/GeometryMap.hpp"

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "AffineMap2D.hpp"
#include "AffineMap3D.hpp"
#include "AnnulusAlignedMap.hpp"
#include "AnnulusMap.hpp"
#include "CircularMap.hpp"
#include "IcosahedralShellAlignedMap.hpp"
#include "IcosahedralShellMap.hpp"
#include "IdentityMap.hpp"
#include "PolarCoordsMap.hpp"
#include "SphericalCoordsMap.hpp"
#include "ThinShellMap.hpp"
#include "TokamakMap.hpp"
#include "TorusMap.hpp"

namespace hyteg {
void GeometryMap::evalFinv( const Point3D& xPhys, Point3D& xComp ) const
{
   WALBERLA_UNUSED( xPhys );
   WALBERLA_UNUSED( xComp );
   WALBERLA_ABORT( "GeometryMap::evalFinv() not implemented for current child map" );
}

real_t GeometryMap::evalDF( const Point3D& x, Matrix3r& DFx ) const
{
   WALBERLA_UNUSED( x );
   WALBERLA_UNUSED( DFx );
   WALBERLA_ABORT( "GeometryMap::evalDF() not implemented for 3D in child class!" );
}

void GeometryMap::evalDFinvDF( const Point3D& x, Matrixr< 2, 4 >& DFinvDFx ) const
{
   WALBERLA_UNUSED( x );
   WALBERLA_UNUSED( DFinvDFx );
   WALBERLA_ABORT( "GeometryMap::evalDFinvDF() not implemented for 2D in child class!" );
}

void GeometryMap::evalDFinvDF( const Point3D& x, Matrixr< 3, 9 >& DFinvDFx ) const
{
   WALBERLA_UNUSED( x );
   WALBERLA_UNUSED( DFinvDFx );
   WALBERLA_ABORT( "GeometryMap::evalDFinvDF() not implemented for 3D in child class!" );
}

void GeometryMap::evalDFinv( const Point3D& x, Matrix3r& DFinvx ) const
{
   WALBERLA_UNUSED( x );
   WALBERLA_UNUSED( DFinvx );
   WALBERLA_ABORT( "GeometryMap::evalDFinv() not implemented for 3D in child class!" );
}

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
   case Type::ANNULUS_ALIGNED:
      return std::make_shared< AnnulusAlignedMap >( recvBuffer );
   case Type::CIRCULAR:
      return std::make_shared< CircularMap >( recvBuffer );
   case Type::IDENTITY:
      return std::make_shared< IdentityMap >();
   case Type::ICOSAHEDRAL_SHELL:
      return std::make_shared< IcosahedralShellMap >( recvBuffer );
   case Type::ICOSAHEDRAL_SHELL_ALIGNED:
      return std::make_shared< IcosahedralShellAlignedMap >( recvBuffer );
   case Type::THIN_SHELL:
      return std::make_shared< ThinShellMap >( recvBuffer );
   case Type::POLAR_COORDS:
      return std::make_shared< PolarCoordsMap >();
   case Type::SPHERICAL_COORDS:
      return std::make_shared< SphericalCoordsMap >();
   case Type::TOKAMAK:
      return std::make_shared< TokamakMap >( recvBuffer );
   case Type::TORUS:
      return std::make_shared< TorusMap >( recvBuffer );
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
