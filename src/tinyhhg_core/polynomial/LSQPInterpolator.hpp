#pragma once

#include "tinyhhg_core/eigen/EigenWrapper.hpp"

#ifdef HHG_BUILD_WITH_EIGEN

#include "Polynomial2D.hpp"

namespace hhg {

template<uint_t Degree, typename Basis>
class LSQPInterpolator {
public:

  static const uint_t NumCoefficients;

  LSQPInterpolator(uint_t interpolationLevel)
      : interpolationLevel_(interpolationLevel),
        numInterpolationPoints_((levelinfo::num_microedges_per_face(interpolationLevel) - 3 * levelinfo::num_microedges_per_edge(interpolationLevel) - 3) / 3),
        offset_(0),
        A(numInterpolationPoints_, NumCoefficients),
        rhs(numInterpolationPoints_, 1) {
  }

  void addInterpolationPoint(const Point2D& x, real_t value) {
    WALBERLA_ASSERT(offset_ < numInterpolationPoints_, "Added too many interpolation points");

    for (uint_t k = 0; k < NumCoefficients; ++k) {
      A(offset_, k) = Basis::eval(k, x);
    }

    rhs(offset_) = value;

    ++offset_;
  }

  void interpolate(Polynomial2D<Degree, Basis>& poly) {
    WALBERLA_ASSERT(offset_ == numInterpolationPoints_, "Not enough interpolation points were added");

    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> coeffs;
    coeffs = A.colPivHouseholderQr().solve(rhs);

    for (uint_t i = 0; i < NumCoefficients; ++i) {
      poly.setCoefficient(i, coeffs(i));
    }
  }

private:
  uint_t interpolationLevel_;
  uint_t numInterpolationPoints_;
  uint_t offset_;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> rhs;
};

template<uint_t Degree, typename Basis>
const uint_t LSQPInterpolator<Degree, Basis>::NumCoefficients = Polynomial2D<Degree, Basis>::NumCoefficients_;

}

#endif