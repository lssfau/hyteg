#pragma once

#include <boost/math/special_functions/binomial.hpp>

#include "monomial.hpp"

namespace hhg {

template<uint_t Degree>
class Polynomial2D {
 public:

  Polynomial2D() {
//    WALBERLA_LOG_DEVEL("numCoefficients = " << getNumCoefficients());
    coeffs_.resize(getNumCoefficients());
  }

  real_t eval(const Point2D &x) const {

    real_t eval = real_c(0);

    uint_t offset = 0;
    for (uint_t d = 0; d <= Degree; ++d) {

      uint_t i = d;
      uint_t j = 0;

      for (uint_t k = 0; k < getNumCoefficientsForDegree(d); ++k) {

//        WALBERLA_LOG_DEVEL("(d,i,j) = (" << d << ", "  << i << ", " << j << ")");

        eval += coeffs_[offset] * Monomial2D::eval(i, j, x);

        ++offset;
        --i;
        ++j;
      }
    }

    return eval;
  }

  void setCoefficient(uint_t idx, real_t value) {
    WALBERLA_ASSERT(idx < getNumCoefficients());
    coeffs_[idx] = value;
  }

  real_t getCoefficient(uint_t idx) const {
    WALBERLA_ASSERT(idx < getNumCoefficients());
    return coeffs_[idx];
  }

  void scale(real_t scalar) {
    for (uint_t i = 0; i < coeffs_.size(); ++i) {
      coeffs_[i] *= scalar;
    }
  }

  void scaleAdd(real_t scalar, const Polynomial2D<Degree>& rhs) {
    for (uint_t i = 0; i < coeffs_.size(); ++i) {
      coeffs_[i] += scalar * rhs.coeffs_[i];
    }
  }

  static uint_t getNumCoefficientsForDegree(uint_t degree) {
    return static_cast<uint_t>(boost::math::binomial_coefficient<real_t>(2 + degree - 1, degree));
  }

  static uint_t getNumCoefficients() {
    return static_cast<uint_t>(boost::math::binomial_coefficient<real_t>(2 + Degree, Degree));
  }

private:
  std::vector<real_t> coeffs_;

};

template<uint_t Degree>
inline std::ostream& operator<<(std::ostream &os, const Polynomial2D<Degree> &poly)
{
  os << "[";

  for (size_t i = 0; i < poly.getNumCoefficients(); ++i)
  {
    os << poly.getCoefficient(i);
    if (i != poly.getNumCoefficients()-1)
    {
      os << ", ";
    }
  }

  os << "]";

  return os;
}

}