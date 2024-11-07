/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
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

#include <array>

#include "core/DataTypes.h"

namespace hyteg {
namespace P1 {

template<uint_t N>
inline real_t getSupremumNorm(std::array<P1Function<real_t>*, N> functions, uint_t level) {
  WALBERLA_ASSERT(N >= 1, "getSupremumNorm requires at least one function as argument");

  real_t maxValue = functions[0]->getMaxDoFValue(level);

  for (uint_t k = 1; k < functions.size(); ++k) {
    maxValue = std::max(maxValue, functions[k]->getMaxDoFValue(level));
  }

  return maxValue;
}

template<uint_t N>
inline real_t getApproximateEuclideanNorm(std::array<P1Function<real_t>*, N> functions, uint_t level) {
  return std::sqrt(walberla::real_c(N)) * getSupremumNorm<N>(functions, level);
}

}
}
