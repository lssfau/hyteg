#pragma once

#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/matrix.hpp"

namespace hhg {

class FaceMap {
public:
  virtual ~FaceMap() {}

  virtual void evalF(const Point3D&, Point3D&) = 0;
  virtual void evalDF(const Point3D&, Matrix2r&) = 0;
  virtual void evalDFinv(const Point3D&, Matrix2r&) = 0;

  void evalTensorCoeff(const Point3D& x, Matrix2r& coeff) {
    evalDFinv(x, coeff);
    real_t invDet = coeff.det();
    Matrix2r coeffT = coeff.transpose();
    coeff = coeff.mul(coeffT);
    coeff *= 1.0 / std::abs(invDet);
  }
};

}
