#pragma once

#include "hierarchicalbasis.hpp"

namespace hhg {

constexpr uint_t getNumCoefficientsForDegree(uint_t degree) {
  return math::binomialCoefficient(2 + degree - 1, degree);
}

constexpr uint_t getNumCoefficients(uint_t degree) {
  return math::binomialCoefficient(2 + degree, degree);
}

template<uint_t Degree, uint_t InterpolationLevel, typename HierarchicalBasis>
class Polynomial2D {
 public:

  static const uint_t NumCoefficients_ = getNumCoefficients(Degree);

  Polynomial2D() {
    coeffs_.resize(NumCoefficients_);
  }

  real_t eval(const Point2D &x) const {

    real_t eval = coeffs_[0] * HierarchicalBasis::eval(InterpolationLevel, 0, x);

    for (uint_t c = 1; c < NumCoefficients_; ++c) {
      eval += coeffs_[c] * HierarchicalBasis::eval(InterpolationLevel, c, x);
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

  void scaleAdd(real_t scalar, const Polynomial2D<Degree, InterpolationLevel, HierarchicalBasis>& rhs) {
    for (uint_t i = 0; i < NumCoefficients_; ++i) {
      coeffs_[i] += scalar * rhs.coeffs_[i];
    }
  }

private:
  std::vector<real_t> coeffs_;

};

template<uint_t Degree, uint_t InterpolationLevel, typename HierarchicalBasis>
inline std::ostream& operator<<(std::ostream &os, const Polynomial2D<Degree, InterpolationLevel, HierarchicalBasis> &poly)
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
using HorizontalPolynomial2D = Polynomial2D<Degree, InterpolationLevel, HorizontalEdgeBasis>;

template<uint_t Degree, uint_t InterpolationLevel>
using VerticalPolynomial2D = Polynomial2D<Degree, InterpolationLevel, VerticalEdgeBasis>;

}