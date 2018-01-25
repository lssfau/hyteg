#pragma once

#include "tinyhhg_core/eigen/EigenWrapper.hpp"

#ifdef HHG_BUILD_WITH_EIGEN

#include "Polynomial2D.hpp"

namespace hhg {

template<typename Basis>
class LSQPInterpolator {
public:

  LSQPInterpolator(uint_t degree, uint_t interpolationLevel)
      : degree_(degree),
        numCoefficients_(Polynomial2D<Basis>::getNumCoefficients(degree)),
        interpolationLevel_(interpolationLevel),
        numInterpolationPoints_((levelinfo::num_microedges_per_face(interpolationLevel) - 3 * levelinfo::num_microedges_per_edge(interpolationLevel) - 3) / 3),
        offset_(0),
        A(numInterpolationPoints_, Polynomial2D<Basis>::getNumCoefficients(degree)),
        rhs(numInterpolationPoints_, 1) {
  }

  void addInterpolationPoint(const Point2D& x, real_t value) {
    WALBERLA_ASSERT(offset_ < numInterpolationPoints_, "Added too many interpolation points");

    for (uint_t k = 0; k < numCoefficients_; ++k) {
      A(offset_, k) = Basis::eval(k, x);
    }

    rhs(offset_) = value;

    ++offset_;
  }

  void interpolate(Polynomial2D<Basis>& poly) {
    WALBERLA_ASSERT(offset_ == numInterpolationPoints_, "Not enough interpolation points were added");

    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> coeffs;
    coeffs = A.colPivHouseholderQr().solve(rhs);

    for (uint_t i = 0; i < numCoefficients_; ++i) {
      poly.setCoefficient(i, coeffs(i));
    }
  }

private:
  uint_t degree_;
  uint_t numCoefficients_;
  uint_t interpolationLevel_;
  uint_t numInterpolationPoints_;
  uint_t offset_;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> rhs;
};

}

#endif