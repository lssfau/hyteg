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
  walberla::real_t det2(const std::array<PointND<walberla::real_t, M>, N>& m)
  {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
  }

}
}

#endif /* MATH_HPP */