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

} // namespace math
} // namespace hhg

#endif /* MATH_HPP */
