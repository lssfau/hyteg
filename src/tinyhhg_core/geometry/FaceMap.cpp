#include "IdentityMap.hpp"
#include "CircularMap.hpp"
#include "AffineMap.hpp"

namespace hhg {

std::shared_ptr<FaceMap> FaceMap::deserialize(walberla::mpi::RecvBuffer& recvBuffer) {
  Type type;
  recvBuffer >> type;

  switch (type) {
    case IDENTITY:
      return std::shared_ptr<FaceMap>(new IdentityMap());
    case AFFINE:
      return std::shared_ptr<FaceMap>(new AffineMap(recvBuffer));
    case CIRCULAR:
      return std::shared_ptr<FaceMap>(new CircularMap(recvBuffer));
  }
}

void FaceMap::evalTensorCoeff(const Point3D& x, Matrix2r& coeff) {
  Point3D xPhy;
  evalF(x, xPhy);
  Matrix2r Khat;
  // TODO: remove hardcoded TENSOR
  Khat(0,0) = 3.0 * xPhy[0]*xPhy[0] + 2.0 * xPhy[1]*xPhy[1] + 1.0;
  Khat(0,1) = -xPhy[0]*xPhy[0] - xPhy[1]*xPhy[1];
  Khat(1,0) = Khat(0,1);
  Khat(1,1) = 4.0 * xPhy[0]*xPhy[0] + 5.0 * xPhy[1]*xPhy[1] + 1.0;
  evalDFinv(x, coeff);
  real_t invDet = coeff.det();
  Matrix2r coeffT = coeff.transpose();
  coeff = coeff.mul(Khat.mul(coeffT));
  coeff *= 1.0 / std::abs(invDet);
}

real_t FaceMap::evalDetDF(const Point3D& x) {
  Matrix2r DF;
  evalDF(x, DF);
  return DF.det();
}

}