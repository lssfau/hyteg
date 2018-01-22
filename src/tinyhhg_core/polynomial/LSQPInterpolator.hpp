#pragma once

#include "tinyhhg_core/eigen/EigenWrapper.hpp"

#ifdef HHG_BUILD_WITH_EIGEN

#include "Polynomial2D.hpp"

namespace hhg {

template<uint_t Degree, uint_t InterpolationLevel, typename Basis>
class LSQPInterpolator {
public:

  static const uint_t NumVertices = (levelinfo::num_microedges_per_face(InterpolationLevel) - 3 * levelinfo::num_microedges_per_edge(InterpolationLevel) - 3) / 3;
  static const uint_t NumCoefficients = Polynomial2D<Degree, InterpolationLevel, Basis>::NumCoefficients_;

  LSQPInterpolator()
      : offset_(0), A(NumVertices, NumCoefficients), rhs(NumVertices, 1) {
  }

  void addInterpolationPoint(const Point2D& x, real_t value) {
    WALBERLA_ASSERT(offset_ < NumVertices, "Added too many interpolation points");

    for (uint_t k = 0; k < NumCoefficients; ++k) {
      A(offset_, k) = Basis::eval(InterpolationLevel, k, x);
    }

    rhs(offset_) = value;

    ++offset_;
  }

  void interpolate(Polynomial2D<Degree, InterpolationLevel, Basis>& poly) {
    WALBERLA_ASSERT(values.size() == NumVertices, "values vector must have the same size as the number of vertices on the interpolation level");
    WALBERLA_ASSERT(offset_ == NumVertices, "Not enough interpolation points were added");

    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> coeffs;
    coeffs = A.colPivHouseholderQr().solve(rhs);

    for (uint_t i = 0; i < NumCoefficients; ++i) {
      poly.setCoefficient(i, coeffs(i));
    }
  }

private:
  uint_t offset_;
  Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<real_t, Eigen::Dynamic, 1> rhs;
};

}

#endif