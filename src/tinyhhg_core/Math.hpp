#ifndef MATH_HPP
#define MATH_HPP

#include <cmath>

#include "core/DataTypes.h"

#include "types/pointnd.hpp"

namespace hhg {
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
} // namespace hhg

#endif /* MATH_HPP */
