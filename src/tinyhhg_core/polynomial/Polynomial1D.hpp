#pragma once

#include "MonomialBasis1D.hpp"

namespace hhg {

template<typename Basis>
class Polynomial1D {
 public:

  static constexpr uint_t getNumCoefficientsForDegree(uint_t degree) {
    return 1;
  }

  static constexpr uint_t getNumCoefficients(uint_t degree) {
    return degree + 1;
  }

  Polynomial1D(uint_t degree)
    : degree_(degree),
      numCoefficients_(getNumCoefficients(degree)),
      coeffs_(getNumCoefficients(degree))
  {
  }

  real_t eval(real_t x) const {

    real_t eval = coeffs_[0] * Basis::eval(0, x);

    for (uint_t c = 1; c < numCoefficients_; ++c) {
      eval = std::fma(coeffs_[c], Basis::eval(c, x), eval);
    }

    return eval;
  }

  void setCoefficient(uint_t idx, real_t value) {
    WALBERLA_ASSERT(idx < numCoefficients_);
    coeffs_[idx] = value;
  }

  void addToCoefficient(uint_t idx, real_t value) {
    WALBERLA_ASSERT(idx < numCoefficients_);
    coeffs_[idx] += value;
  }

  real_t getCoefficient(uint_t idx) const {
    WALBERLA_ASSERT(idx < numCoefficients_);
    return coeffs_[idx];
  }

  void scale(real_t scalar) {
    for (uint_t i = 0; i < numCoefficients_; ++i) {
      coeffs_[i] *= scalar;
    }
  }

  void scaleAdd(real_t scalar, const Polynomial1D<Basis>& rhs) {
    for (uint_t i = 0; i < numCoefficients_; ++i) {
      coeffs_[i] += scalar * rhs.coeffs_[i];
    }
  }

private:
  uint_t degree_;
  uint_t numCoefficients_;
  std::vector<real_t> coeffs_;

};

template<uint_t Degree, uint_t InterpolationLevel, typename Basis>
inline std::ostream& operator<<(std::ostream &os, const Polynomial1D<Basis> &poly)
{
  os << "[";

  for (size_t i = 0; i < poly.numCoefficients_; ++i)
  {
    os << poly.getCoefficient(i);
    if (i != poly.numCoefficients_-1)
    {
      os << ", ";
    }
  }

  os << "]";

  return os;
}

template<uint_t InterpolationLevel>
using GeneralPolynomial1D = Polynomial1D<MonomialBasis1D>;

}