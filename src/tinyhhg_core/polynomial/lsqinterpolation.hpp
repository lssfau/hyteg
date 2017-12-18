#pragma once

#include "tinyhhg_core/eigen/EigenWrapper.hpp"

#ifdef HHG_BUILD_WITH_EIGEN

#include "hierarchicalbasis.hpp"
#include "polynomial.hpp"

namespace hhg {

template<uint_t Degree, uint_t InterpolationLevel>
class LSQInterpolator {
public:

  static const uint_t NumVertices = levelinfo::num_microvertices_per_face(InterpolationLevel);
  static const uint_t NumCoefficients = Polynomial2D<Degree, InterpolationLevel>::NumCoefficients_;

  LSQInterpolator() {
    uint_t rowsize = levelinfo::num_microvertices_per_edge(InterpolationLevel);

    uint_t inner_rowsize = rowsize;
    real_t h = 1.0 / (rowsize-1);

    Point2D x;

    uint_t offset = 0;

    for (uint_t i = 0; i < rowsize; ++i) {
      x[0] = i * h;

      for (uint_t j = 0; j < inner_rowsize; ++j) {
        x[1] = j * h;

        for (uint_t k = 0; k < NumCoefficients; ++k) {
          A(offset, k) = HierarchicalBasis::eval(InterpolationLevel, k, x);
        }

        ++offset;
      }

      --inner_rowsize;
    }
  }

  void interpolate(std::vector<real_t> values, Polynomial2D<Degree, InterpolationLevel>& poly) {
    WALBERLA_ASSERT(values.size() == NumVertices, "values vector must have the same size as the number of vertices on the interpolation level")

    Eigen::Map<Eigen::Matrix<real_t, NumVertices, 1>> rhsVector(values.data());

    // We do not need to solve the system since the basis is orthonormal and a simple multiplication is sufficient
    // but solving the system is numerically more stable
    Eigen::Matrix<real_t, NumCoefficients, 1> coeffs = A.colPivHouseholderQr().solve(rhsVector);
//    Eigen::Matrix<real_t, NumCoefficients, 1> coeffs = A.transpose() * rhsVector;

    for (uint_t i = 0; i < NumCoefficients; ++i) {
      poly.setCoefficient(i, coeffs(i));
    }
  }

private:
  uint_t interpolationLevel_;
  Eigen::Matrix<real_t, NumVertices, NumCoefficients> A;

};

}

#endif