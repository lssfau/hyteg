#pragma once

#include "Polynomial1D.hpp"
#include "Polynomial2D.hpp"

namespace hhg {

class Polynomial2DEvaluator {
public:

  typedef Polynomial1D<MonomialBasis1D> Polynomial1;
  typedef Polynomial2D<MonomialBasis2D> Polynomial2;

  Polynomial2DEvaluator(const Polynomial2& poly)
    : degree_(poly.getDegree()),
      poly2_(poly),
      poly1_(poly.getDegree()),
      deltas(poly.getDegree() + 1)
  {
  }

  real_t eval(const Point2D &x) const {
    return poly2_.eval(x);
  }

  void setY(real_t y) {

    for (uint_t degree = 0; degree <= degree_; ++degree) {
      poly1_.setCoefficient(degree, 0.0);
    }

    int start = 0;
    real_t y_;

    for (uint_t coeff = 0; coeff <= degree_; ++coeff) {

      int idx = start;
      y_ = walberla::real_c(1.0);

      for(uint_t degree = 0; degree <= degree_-coeff; ++degree) {

        poly1_.addToCoefficient(coeff, poly2_.getCoefficient(idx) * y_);

        idx += coeff + degree + 2;
        y_ *= y;
      }

      start += coeff + 1;
    }
  }

  real_t evalX(real_t x) const {
    return poly1_.eval(x);
  }

  template<uint_t Degree>
  real_t setStartX(real_t x, real_t h) {
    static_assert(Degree <= 7, "Polynomial2DEvaluator not implemented for degree larger than 7");
    if (Degree == 0) {
      deltas[0] = poly1_.getCoefficient(0);
    }
    if (Degree == 1) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x;
      deltas[1] = h*poly1_.getCoefficient(1);
    }
    if (Degree == 2) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x;
      deltas[1] = 2*h*poly1_.getCoefficient(2)*x + h*(h*poly1_.getCoefficient(2) + poly1_.getCoefficient(1));
      deltas[2] = 2*h*h*poly1_.getCoefficient(2);
    }
    if (Degree == 3) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x;
      deltas[1] = h*(h*(h*poly1_.getCoefficient(3) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(3*h*poly1_.getCoefficient(3)*x + h*(3*h*poly1_.getCoefficient(3) + 2*poly1_.getCoefficient(2)));
      deltas[2] = 6*h*h*poly1_.getCoefficient(3)*x + h*h*(6*h*poly1_.getCoefficient(3) + 2*poly1_.getCoefficient(2));
      deltas[3] = 6*h*h*h*poly1_.getCoefficient(3);
    }
    if (Degree == 4) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x;
      deltas[1] = h*(h*(h*(h*poly1_.getCoefficient(4) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(4*h*poly1_.getCoefficient(4) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(4*h*poly1_.getCoefficient(4)*x + h*(6*h*poly1_.getCoefficient(4) + 3*poly1_.getCoefficient(3))));
      deltas[2] = h*h*(h*(14*h*poly1_.getCoefficient(4) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(12*h*h*poly1_.getCoefficient(4)*x + h*h*(24*h*poly1_.getCoefficient(4) + 6*poly1_.getCoefficient(3)));
      deltas[3] = 24*h*h*h*poly1_.getCoefficient(4)*x + h*h*h*(36*h*poly1_.getCoefficient(4) + 6*poly1_.getCoefficient(3));
      deltas[4] = 24*h*h*h*h*poly1_.getCoefficient(4);
    }
    if (Degree == 5) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*poly1_.getCoefficient(5) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(5*h*poly1_.getCoefficient(5) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(10*h*poly1_.getCoefficient(5) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(5*h*poly1_.getCoefficient(5)*x + h*(10*h*poly1_.getCoefficient(5) + 4*poly1_.getCoefficient(4)))));
      deltas[2] = h*h*(h*(h*(30*h*poly1_.getCoefficient(5) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(70*h*poly1_.getCoefficient(5) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(20*h*h*poly1_.getCoefficient(5)*x + h*h*(60*h*poly1_.getCoefficient(5) + 12*poly1_.getCoefficient(4))));
      deltas[3] = h*h*h*(h*(150*h*poly1_.getCoefficient(5) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(60*h*h*h*poly1_.getCoefficient(5)*x + h*h*h*(180*h*poly1_.getCoefficient(5) + 24*poly1_.getCoefficient(4)));
      deltas[4] = 120*h*h*h*h*poly1_.getCoefficient(5)*x + h*h*h*h*(240*h*poly1_.getCoefficient(5) + 24*poly1_.getCoefficient(4));
      deltas[5] = 120*h*h*h*h*h*poly1_.getCoefficient(5);
    }
    if (Degree == 6) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*poly1_.getCoefficient(6) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(6*h*poly1_.getCoefficient(6) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(15*h*poly1_.getCoefficient(6) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(20*h*poly1_.getCoefficient(6) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(6*h*poly1_.getCoefficient(6)*x + h*(15*h*poly1_.getCoefficient(6) + 5*poly1_.getCoefficient(5))))));
      deltas[2] = h*h*(h*(h*(h*(62*h*poly1_.getCoefficient(6) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(180*h*poly1_.getCoefficient(6) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(210*h*poly1_.getCoefficient(6) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(30*h*h*poly1_.getCoefficient(6)*x + h*h*(120*h*poly1_.getCoefficient(6) + 20*poly1_.getCoefficient(5)))));
      deltas[3] = h*h*h*(h*(h*(540*h*poly1_.getCoefficient(6) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(900*h*poly1_.getCoefficient(6) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(120*h*h*h*poly1_.getCoefficient(6)*x + h*h*h*(540*h*poly1_.getCoefficient(6) + 60*poly1_.getCoefficient(5))));
      deltas[4] = h*h*h*h*(h*(1560*h*poly1_.getCoefficient(6) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(360*h*h*h*h*poly1_.getCoefficient(6)*x + h*h*h*h*(1440*h*poly1_.getCoefficient(6) + 120*poly1_.getCoefficient(5)));
      deltas[5] = 720*h*h*h*h*h*poly1_.getCoefficient(6)*x + h*h*h*h*h*(1800*h*poly1_.getCoefficient(6) + 120*poly1_.getCoefficient(5));
      deltas[6] = 720*h*h*h*h*h*h*poly1_.getCoefficient(6);
    }
    if (Degree == 7) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x + poly1_.getCoefficient(7)*x*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*(h*poly1_.getCoefficient(7) + poly1_.getCoefficient(6)) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(h*(7*h*poly1_.getCoefficient(7) + 6*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(h*(21*h*poly1_.getCoefficient(7) + 15*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(h*(35*h*poly1_.getCoefficient(7) + 20*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(h*(h*(35*h*poly1_.getCoefficient(7) + 15*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + x*(7*h*poly1_.getCoefficient(7)*x + h*(21*h*poly1_.getCoefficient(7) + 6*poly1_.getCoefficient(6)))))));
      deltas[2] = h*h*(h*(h*(h*(h*(126*h*poly1_.getCoefficient(7) + 62*poly1_.getCoefficient(6)) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(h*(434*h*poly1_.getCoefficient(7) + 180*poly1_.getCoefficient(6)) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(h*(630*h*poly1_.getCoefficient(7) + 210*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(h*h*(h*(490*h*poly1_.getCoefficient(7) + 120*poly1_.getCoefficient(6)) + 20*poly1_.getCoefficient(5)) + x*(42*h*h*poly1_.getCoefficient(7)*x + h*h*(210*h*poly1_.getCoefficient(7) + 30*poly1_.getCoefficient(6))))));
      deltas[3] = h*h*h*(h*(h*(h*(1806*h*poly1_.getCoefficient(7) + 540*poly1_.getCoefficient(6)) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(h*(3780*h*poly1_.getCoefficient(7) + 900*poly1_.getCoefficient(6)) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*(h*(3150*h*poly1_.getCoefficient(7) + 540*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + x*(210*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*(1260*h*poly1_.getCoefficient(7) + 120*poly1_.getCoefficient(6)))));
      deltas[4] = h*h*h*h*(h*(h*(8400*h*poly1_.getCoefficient(7) + 1560*poly1_.getCoefficient(6)) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*h*(h*(10920*h*poly1_.getCoefficient(7) + 1440*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(840*h*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*h*(5040*h*poly1_.getCoefficient(7) + 360*poly1_.getCoefficient(6))));
      deltas[5] = h*h*h*h*h*(h*(16800*h*poly1_.getCoefficient(7) + 1800*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(2520*h*h*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*h*h*(12600*h*poly1_.getCoefficient(7) + 720*poly1_.getCoefficient(6)));
      deltas[6] = 5040*h*h*h*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*h*h*h*(15120*h*poly1_.getCoefficient(7) + 720*poly1_.getCoefficient(6));
      deltas[7] = 5040*h*h*h*h*h*h*h*poly1_.getCoefficient(7);
    }
    return deltas[0];
  }

  template<uint_t Degree>
  real_t incrementEval() {
    if (Degree >= 1) {
      deltas[0] += deltas[1];
    }
    if (Degree >= 2) {
      deltas[1] += deltas[2];
    }
    if (Degree >= 3) {
      deltas[2] += deltas[3];
    }
    if (Degree >= 4) {
      deltas[3] += deltas[4];
    }
    if (Degree >= 5) {
      deltas[4] += deltas[5];
    }
    if (Degree >= 6) {
      deltas[5] += deltas[6];
    }
    if (Degree >= 7) {
      deltas[6] += deltas[7];
    }
    return deltas[0];
  }

private:
  uint_t degree_;
  const Polynomial2& poly2_;
  Polynomial1 poly1_;

  std::vector<real_t> deltas;

};

}