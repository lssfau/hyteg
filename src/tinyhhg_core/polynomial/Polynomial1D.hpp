#pragma once

#include "MonomialBasis1D.hpp"

namespace hhg {

template<uint_t Degree, typename Basis>
class Polynomial1D {
 public:

  static constexpr uint_t getNumCoefficientsForDegree(uint_t degree) {
    return 1;
  }

  static constexpr uint_t getNumCoefficients(uint_t degree) {
    return degree + 1;
  }

  static const uint_t NumCoefficients_ = getNumCoefficients(Degree);

  Polynomial1D() {
  }

  real_t eval(real_t x) const {

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

  void scaleAdd(real_t scalar, const Polynomial1D<Degree, Basis>& rhs) {
    for (uint_t i = 0; i < NumCoefficients_; ++i) {
      coeffs_[i] += scalar * rhs.coeffs_[i];
    }
  }

private:
  std::array<real_t, NumCoefficients_> coeffs_;

};

template<uint_t Degree, uint_t InterpolationLevel, typename Basis>
inline std::ostream& operator<<(std::ostream &os, const Polynomial1D<Degree, Basis> &poly)
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

template<uint_t Degree, uint_t InterpolationLevel>
using GeneralPolynomial1D = Polynomial1D<Degree, MonomialBasis1D>;

}