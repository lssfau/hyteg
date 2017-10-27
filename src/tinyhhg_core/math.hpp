#ifndef MATH_HPP
#define MATH_HPP

#include <cmath>
#include "types/pointnd.hpp"
#include "core/DataTypes.h"

namespace hhg
{
namespace math
{

using walberla::real_t;

template<size_t M, size_t N>
inline real_t det2(const std::array<PointND<walberla::real_t, M>, N>& m)
{
  return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

inline real_t faceOrientation2D(const Point3D& a, const Point3D& b, const Point3D& c)
{
  std::array<Point3D, 2> jac;
  jac[0] = b-a;
  jac[1] = c-a;

  return std::copysign(1.0, det2(jac));
}

}
}

#endif /* MATH_HPP */