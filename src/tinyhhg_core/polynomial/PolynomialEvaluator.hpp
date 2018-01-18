#pragma once

#include "Polynomial1D.hpp"
#include "Polynomial2D.hpp"

namespace hhg {

template<uint_t Degree, uint_t InterpolationLevel>
class Polynomial2DEvaluator {
public:

  typedef Polynomial1D<Degree, MonomialBasis1D> Polynomial1;
  typedef Polynomial2D<Degree, InterpolationLevel, MonomialBasis2D> Polynomial2;

  Polynomial2DEvaluator(const Polynomial2& poly)
    : poly2_(poly) {
    static_assert(Degree <= 2, "Polynomial2DEvaluator not implemented for degree larger than 2");
  }

  real_t eval(const Point2D &x) const {
    return poly2_.eval(x);
  }

  void setY(real_t y) {
    // set 1D Polynomial coefficients

    if (Degree == 0) {
      poly1_.setCoefficient(0, poly2_.getCoefficient(0));
    }

    if (Degree == 1) {
      poly1_.setCoefficient(0, poly2_.getCoefficient(0) + poly2_.getCoefficient(2) * y);
      poly1_.setCoefficient(1, poly2_.getCoefficient(1));
    }

    if (Degree == 2) {
      poly1_.setCoefficient(0, poly2_.getCoefficient(0) + poly2_.getCoefficient(2) * y + poly2_.getCoefficient(5) * y * y);
      poly1_.setCoefficient(1, poly2_.getCoefficient(1) + poly2_.getCoefficient(4) * y);
      poly1_.setCoefficient(2, poly2_.getCoefficient(3));
    }
  }

  real_t evalX(real_t x) const {
    return poly1_.eval(x);
  }

  real_t setStartX(real_t x, real_t h) {
    WALBERLA_ASSERT(Degree == 2, "Incremental evaluation only impelemented for Degree 2");
    deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*pow(x, 2);
    deltas[1] = -poly1_.getCoefficient(1)*x + poly1_.getCoefficient(1)*(h + x) - poly1_.getCoefficient(2)*pow(x, 2) + poly1_.getCoefficient(2)*pow(h + x, 2);
    deltas[2] = poly1_.getCoefficient(1)*x - 2*poly1_.getCoefficient(1)*(h + x) + poly1_.getCoefficient(1)*(2*h + x) + poly1_.getCoefficient(2)*pow(x, 2) - 2*poly1_.getCoefficient(2)*pow(h + x, 2) + poly1_.getCoefficient(2)*pow(2*h + x, 2);
    return deltas[0];
  }

  real_t incrementEval() {
    WALBERLA_ASSERT(Degree == 2, "Incremental evaluation only impelemented for Degree 2");
    deltas[0] += deltas[1];
    deltas[1] += deltas[2];
    return deltas[0];
  }

private:
  const Polynomial2& poly2_;
  Polynomial1 poly1_;

  std::array<real_t, Degree+1> deltas;

};

}