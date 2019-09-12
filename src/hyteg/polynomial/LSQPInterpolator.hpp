/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "hyteg/eigen/EigenWrapper.hpp"

#ifdef HYTEG_BUILD_WITH_EIGEN

#include "Polynomial2D.hpp"

namespace hyteg {

enum class LSQPType
{
  EDGE,
  VERTEX
};

template<LSQPType Type>
constexpr uint_t GetNumInterpolationPoints(uint_t level)
{
  switch(Type)
  {
    case LSQPType::EDGE:
      return (levelinfo::num_microedges_per_face(level) - 3 * levelinfo::num_microedges_per_edge(level) - 3) / 3;
    case LSQPType::VERTEX:
      return levelinfo::num_microvertices_per_face(level) - 3 * (levelinfo::num_microvertices_per_edge(level) - 2) - 3;
  }
}

template<typename Basis, LSQPType Type>
class LSQPInterpolator {
public:

  LSQPInterpolator(uint_t degree, uint_t interpolationLevel)
      : degree_(degree),
        numCoefficients_(Polynomial2D<Basis>::getNumCoefficients(degree)),
        interpolationLevel_(interpolationLevel),
        numInterpolationPoints_(GetNumInterpolationPoints<Type>(interpolationLevel)),
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

    if (numInterpolationPoints_ < poly.getNumCoefficients(poly.getDegree())) {
      WALBERLA_LOG_WARNING("Polynomial interpolation may have poor quality since there are less interpolation points "
                           "than coefficients. Please try to increase the interpolation level to fix this.");
    }

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