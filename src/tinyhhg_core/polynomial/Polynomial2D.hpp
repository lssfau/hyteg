#pragma once

#include "MonomialBasis2D.hpp"

namespace hhg {

template<uint_t Degree, typename Basis>
class Polynomial2D {
 public:

  static constexpr uint_t getNumCoefficientsForDegree(uint_t degree) {
    return math::binomialCoefficient(2 + degree - 1, degree);
  }

  static constexpr uint_t getNumCoefficients(uint_t degree) {
    return math::binomialCoefficient(2 + degree, degree);
  }

  static const uint_t NumCoefficients_ = getNumCoefficients(Degree);

  Polynomial2D() {
  }

  real_t eval(const Point2D &x) const {

    real_t eval = coeffs_[0] * Basis::eval(0, x);

    for (uint_t c = 1; c < NumCoefficients_; ++c) {
      eval += coeffs_[c] * Basis::eval(c, x);
    }

    return eval;
  }

  void setCoefficient(uint_t idx, real_t value) {
    WALBERLA_ASSERT(idx < NumCoefficients_);
    coeffs_[idx] = value;
  }

  real_t getCoefficient(uint_t idx) const {
    WALBERLA_ASSERT(idx < NumCoefficients_);
    return coeffs_[idx];
  }

  void scale(real_t scalar) {
    for (uint_t i = 0; i < NumCoefficients_; ++i) {
      coeffs_[i] *= scalar;
    }
  }

  void scaleAdd(real_t scalar, const Polynomial2D<Degree, Basis>& rhs) {
    for (uint_t i = 0; i < NumCoefficients_; ++i) {
      coeffs_[i] += scalar * rhs.coeffs_[i];
    }
  }

private:
  std::array<real_t, NumCoefficients_> coeffs_;

};

template<uint_t Degree, uint_t InterpolationLevel, typename Basis>
inline std::ostream& operator<<(std::ostream &os, const Polynomial2D<Degree, Basis> &poly)
{
  os << "[";

  for (size_t i = 0; i < poly.NumCoefficients_; ++i)
  {
    os << poly.getCoefficient(i);
    if (i != poly.NumCoefficients_-1)
    {
      os << ", ";
    }
  }

  os << "]";

  return os;
}

template<uint_t Degree>
using GeneralPolynomial2D = Polynomial2D<Degree, MonomialBasis2D>;

}