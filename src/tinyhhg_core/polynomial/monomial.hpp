#pragma once

#include <tinyhhg_core/types/pointnd.hpp>

namespace hhg {

class Monomial2D {
 public:

  static real_t eval(uint_t xExp, uint_t yExp, const Point2D &x) {
    return std::pow(x[0], xExp)*std::pow(x[1], yExp);
  }

};

}