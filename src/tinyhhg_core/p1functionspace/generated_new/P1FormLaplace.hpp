#pragma once

#include "tinyhhg_core/geometry/FaceMap.hpp"

namespace hhg {

class P1Form_laplace {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    faceMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = coords[0][1] - coords[1][1];
    real_t tmp1 = coords[0][0] - coords[1][0];
    real_t tmp2 = DFinv(0,0)*tmp0 - DFinv(1,0)*tmp1;
    real_t tmp3 = coords[0][0] - coords[2][0];
    real_t tmp4 = DFinv(1,0)*tmp3;
    real_t tmp5 = coords[0][1] - coords[2][1];
    real_t tmp6 = DFinv(0,0)*tmp5;
    real_t tmp7 = tmp2 + tmp4 - tmp6;
    real_t tmp8 = DFinv(0,1)*tmp0 - DFinv(1,1)*tmp1;
    real_t tmp9 = DFinv(1,1)*tmp3;
    real_t tmp10 = DFinv(0,1)*tmp5;
    real_t tmp11 = -tmp10 + tmp8 + tmp9;
    real_t tmp12 = 1.0/fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0));
    real_t tmp13 = -tmp0*tmp3 + tmp1*tmp5;
    real_t tmp14 = pow(tmp13, -2);
    real_t tmp15 = 1.0/fabs(1.0/tmp13);
    real_t tmp16 = 0.5*tmp12*tmp14*tmp15;
    out[0] = tmp16*(pow(tmp11, 2) + pow(tmp7, 2));
    out[1] = tmp16*(tmp11*(tmp10 - tmp9) + tmp7*(-tmp4 + tmp6));
    out[2] = -tmp12*tmp14*tmp15*(0.5*tmp11*tmp8 + 0.5*tmp2*tmp7);
  }

  std::shared_ptr<FaceMap> faceMap;
};

}