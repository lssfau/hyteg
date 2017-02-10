#ifndef MATH_HPP
#define MATH_HPP

#include <cmath>
#include "types/pointnd.hpp"

namespace hhg
{
namespace math
{

  template<size_t M, size_t N>
  double det2(const std::array<PointND<double, M>, N>& m)
  {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
  }

}
}

#endif /* MATH_HPP */