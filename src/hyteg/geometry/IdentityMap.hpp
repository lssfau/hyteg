/*
 * Copyright (c) 2017-2024 Daniel Drzisga, Dominik Thoennes, Andreas Burkhart, Marcus Mohr.
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

#include "hyteg/geometry/GeometryMap.hpp"

using walberla::real_c;

namespace hyteg {

// cannot include the header here, but neet a declaration for the
// interface of IdentityMap::setMap()
class SetupPrimitiveStorage;

class IdentityMap : public GeometryMap
{
 public:
   IdentityMap() {}

   void evalF( const Point3D& x, Point3D& Fx ) const override final;

   void evalFinv( const Point3D& xPhys, Point3D& xComp ) const override final;

   void evalDF( const Point3D&, Matrix2r& DFx ) const override final;

   real_t evalDF( const Point3D&, Matrix3r& DFx ) const override final;

   void evalDFinv( const Point3D&, Matrix2r& DFinvx ) const override final;

   void evalDFinvDF( const Point3D& x, Matrixr< 2, 4 >& DFinvDFx ) const override final;

   void evalDFinvDF( const Point3D& x, Matrixr< 3, 9 >& DFinvDFx ) const override final;

   bool isIdentity() const final { return true; }
   bool isAffine() const final { return true; }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override final;

   static void setMap( SetupPrimitiveStorage& setupStorage );
};

} // namespace hyteg
