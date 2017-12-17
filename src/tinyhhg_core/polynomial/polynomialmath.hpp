#pragma once

#include "hierarchicalbasis.hpp"

namespace hhg {

namespace PolyMath {

template<class Polynomial>
real_t scalarProduct2D(const Polynomial& lhs, const Polynomial& rhs, uint_t level) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  real_t h = 1.0 / (rowsize-1);

  Point2D x;

  real_t sp = real_c(0);

  for (uint_t i = 0; i < rowsize; ++i) {

    x[0] = i * h;

    for (uint_t j = 0; j < inner_rowsize; ++j) {

      x[1] = j * h;
      sp += lhs.eval(x) * rhs.eval(x);
    }

    --inner_rowsize;
  }

  return sp;
}

real_t scalarProduct2DHierarchical(uint_t lhs, uint_t rhs, uint_t level) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  real_t h = 1.0 / (rowsize-1);

  Point2D x;

  real_t sp = real_c(0);

  for (uint_t i = 0; i < rowsize; ++i) {

    x[0] = i * h;

    for (uint_t j = 0; j < inner_rowsize; ++j) {

      x[1] = j * h;
      sp += HierarchicalBasis::eval(level, lhs, x) * HierarchicalBasis::eval(level, rhs, x);
    }

    --inner_rowsize;
  }

  return sp;
}

}

}