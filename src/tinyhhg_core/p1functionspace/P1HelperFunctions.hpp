#pragma once

#include <array>

#include "core/DataTypes.h"

namespace hyteg {
namespace P1 {

template<uint_t N>
inline real_t getSupremumNorm(std::array<P1Function<real_t>*, N> functions, uint_t level) {
  WALBERLA_ASSERT(N >= 1, "getSupremumNorm requires at least one function as argument");

  real_t maxValue = functions[0]->getMaxValue(level);

  for (uint_t k = 1; k < functions.size(); ++k) {
    maxValue = std::max(maxValue, functions[k]->getMaxValue(level));
  }

  return maxValue;
}

template<uint_t N>
inline real_t getApproximateEuclideanNorm(std::array<P1Function<real_t>*, N> functions, uint_t level) {
  return std::sqrt(walberla::real_c(N)) * getSupremumNorm<N>(functions, level);
}

}
}