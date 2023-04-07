/*
* Copyright (c) 2017-2022 Dominik Thoennes.
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

#include <core/math/Vector3.h>

#include "PointND.hpp"

walberla::math::Vector3< walberla::real_t > hyteg::toVec3( const Point3D& p )
{
   return walberla::math::Vector3< real_t >( p[0], p[1], p[2] );
}

hyteg::Point3D hyteg::toPoint3D( const walberla::math::Vector3< real_t >& v )
{
   return Point3D( v[0], v[1], v[2] );
}

template class hyteg::PointND< walberla::real_t, 2 >;
template class hyteg::PointND< walberla::real_t, 3 >;
template class hyteg::PointND< walberla::real_t, 4 >;
template class hyteg::PointND< walberla::real_t, 6 >;
template class hyteg::PointND< walberla::real_t, 10 >;
