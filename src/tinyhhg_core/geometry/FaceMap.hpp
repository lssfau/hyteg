#pragma once

#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/matrix.hpp"

namespace hhg {

class IdentityMap;
class AffineMap;
class CircularMap;

class FaceMap {
public:

  enum Type {
    IDENTITY = 0,
    AFFINE = 1,
    CIRCULAR = 2
  };

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

  virtual void serialize(walberla::mpi::SendBuffer& sendBuffer) = 0;
  static std::shared_ptr<FaceMap> deserialize(walberla::mpi::RecvBuffer& recvBuffer);
};

}
