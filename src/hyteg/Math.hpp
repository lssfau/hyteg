/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#ifndef MATH_HPP
#define MATH_HPP
#ifdef _MSC_VER
  #define _USE_MATH_DEFINES //for M_PI
#endif
#include <cmath>

#include "core/DataTypes.h"

#include "types/pointnd.hpp"

namespace hyteg {
namespace math {

using walberla::real_t;

template < size_t M, size_t N >
inline real_t det2( const std::array< PointND< real_t, M >, N >& m )
{
   return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

inline real_t faceOrientation2D( const Point3D& a, const Point3D& b, const Point3D& c )
{
   std::array< Point3D, 2 > jac;
   jac[0] = b - a;
   jac[1] = c - a;

   return std::copysign( 1.0, det2( jac ) );
}

constexpr uint_t binomialCoefficient(uint_t n, uint_t k)
{
  uint_t out = 1;

  if (k > n - k) {
    k = n - k;
  }

  for (uint_t i = 0; i < k; ++i)
  {
    out *= (n - i);
    out /= (i + 1);
  }

  return out;
}

/// Returns the given point from cartesian coordinates to spherical coordinates.
/// (x, y, z) -> (r, theta, phi) / (radius, inclination, azimuth)
inline Point3D toSpherical( const Point3D & pointInCartesianCoordinates )
{
   const real_t r = pointInCartesianCoordinates.norm();
   const real_t theta = std::acos( pointInCartesianCoordinates[2] / r );
   const real_t phi = std::atan2( pointInCartesianCoordinates[1], pointInCartesianCoordinates[0] );
   return Point3D( { r, theta, phi } );
}

/// Returns the given point from spherical coordinates to cartesian coordinates.
/// (r, theta, phi) / (radius, inclination, azimuth) -> (x, y, z)
inline Point3D toCartesian( const Point3D & pointInSphericalCoordinates )
{
    const real_t r     = pointInSphericalCoordinates[0];
    const real_t theta = pointInSphericalCoordinates[1];
    const real_t phi   = pointInSphericalCoordinates[2];
    const real_t x = r * std::sin( theta ) * std::cos( phi );
    const real_t y = r * std::sin( theta ) * std::sin( phi );
    const real_t z = r * std::cos( theta );
    return Point3D( { x, y, z } );
}

} // namespace math
} // namespace hyteg

#endif /* MATH_HPP */
