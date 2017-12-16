#pragma once

#include "polynomial.hpp"
#include "polynomialmath.hpp"

namespace hhg {

template<uint_t Degree>
class Polynomial2DBasis
{
public:

  Polynomial2DBasis() {
    coeffs_.resize(Polynomial2D<Degree>::getNumCoefficients());
    polys_.resize(Polynomial2D<Degree>::getNumCoefficients());

    for (uint_t d = 0; d < Polynomial2D<Degree>::getNumCoefficients(); ++d) {
      polys_[d].setCoefficient(d, 1.0);
    }
  }

  void orthogonalize(uint_t level) {

    real_t sp = PolyMath::scalarProduct2D(polys_[0], polys_[0], level);
    polys_[0].scale(1.0 / std::sqrt(sp));

    for (uint_t i = 1; i < polys_.size(); ++i) {

      for (uint_t j = 0; j < i; ++j) {

        real_t c = PolyMath::scalarProduct2D(polys_[j], polys_[i], level);

        polys_[i].scaleAdd(-c, polys_[j]);

      }

      sp = PolyMath::scalarProduct2D(polys_[i], polys_[i], level);
      polys_[i].scale(1.0 / std::sqrt(sp));
    }
  }

  real_t eval(const Point2D &x) const {

    real_t eval = coeffs_[0] * polys_[0](x);

    for (uint_t d = 1; d < Polynomial2D<Degree>::getNumCoefficients(); ++d) {
      eval += coeffs_[d] * polys_[d](x);
    }

  }

  std::vector<real_t> coeffs_;
  std::vector<Polynomial2D<Degree>> polys_;

};

}