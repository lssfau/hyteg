#pragma once

#include "MonomialBasis2D.hpp"

namespace hyteg {

template<typename Basis>
class Polynomial2D {
 public:

  static constexpr uint_t getNumCoefficientsForDegree(uint_t degree) {
    return math::binomialCoefficient(2 + degree - 1, degree);
  }

  static constexpr uint_t getNumCoefficients(uint_t degree) {
    return math::binomialCoefficient(2 + degree, degree);
  }

  Polynomial2D(uint_t degree)
    : degree_(degree),
      numCoefficients_(getNumCoefficients(degree)),
      coeffs_(getNumCoefficients(degree))
  {
  }

  uint_t getDegree() const {
    return degree_;
  }

  real_t eval(const Point2D &x) const {

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

  real_t getCoefficient(uint_t idx) const {
    WALBERLA_ASSERT(idx < numCoefficients_);
    return coeffs_[idx];
  }

  void scale(real_t scalar) {
    for (uint_t i = 0; i < numCoefficients_; ++i) {
      coeffs_[i] *= scalar;
    }
  }

  void scaleAdd(real_t scalar, const Polynomial2D<Basis>& rhs) {
    for (uint_t i = 0; i < numCoefficients_; ++i) {
      coeffs_[i] += scalar * rhs.coeffs_[i];
    }
  }

  real_t lInfinityError(const Polynomial2D& rhs) {
    real_t error = std::numeric_limits<real_t>::min();

    uint_t i = 0;

    for (; i < std::min(numCoefficients_, rhs.numCoefficients_); ++i) {
      error = std::max(error, std::abs(coeffs_[i] - rhs.coeffs_[i]));
    }

    for (; i < numCoefficients_; ++i) {
      error = std::max(error, std::abs(coeffs_[i]));
    }

    for (; i < rhs.numCoefficients_; ++i) {
      error = std::max(error, std::abs(rhs.coeffs_[i]));
    }

    return error;
  }

private:
  uint_t degree_;
  uint_t numCoefficients_;
  std::vector<real_t> coeffs_;

};

template<typename Basis>
inline std::ostream& operator<<(std::ostream &os, const Polynomial2D<Basis> &poly)
{
  os << "[";

  uint_t numCoefficients = poly.getNumCoefficients(poly.getDegree());

  for (size_t i = 0; i < numCoefficients; ++i)
  {
    os << poly.getCoefficient(i);
    if (i != numCoefficients-1)
    {
      os << ", ";
    }
  }

  os << "]";

  return os;
}

using GeneralPolynomial2D = Polynomial2D<MonomialBasis2D>;

}