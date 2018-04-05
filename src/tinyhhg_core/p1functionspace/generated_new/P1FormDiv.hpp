#pragma once

#include "tinyhhg_core/geometry/FaceMap.hpp"

namespace hhg {

class P1Form_div_1 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    faceMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][0];
    real_t tmp1 = -coords[2][1];
    real_t tmp2 = coords[0][0] - coords[1][0];
    real_t tmp3 = coords[0][1] + tmp1;
    real_t tmp4 = coords[0][0] + tmp0;
    real_t tmp5 = coords[0][1] - coords[1][1];
    real_t tmp6 = 1.0/(tmp2*tmp3 - tmp4*tmp5);
    real_t tmp7 = 0.166666666666667*tmp6/(fabs(tmp6)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    out[0] = tmp7*(-DFinv(0,0)*(coords[1][1] + tmp1) + DFinv(1,0)*(coords[1][0] + tmp0));
    out[1] = tmp7*(DFinv(0,0)*tmp3 - DFinv(1,0)*tmp4);
    out[2] = tmp7*(-DFinv(0,0)*tmp5 + DFinv(1,0)*tmp2);
  }

  std::shared_ptr<FaceMap> faceMap;
};

class P1Form_div_2 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    faceMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][0];
    real_t tmp1 = -coords[2][1];
    real_t tmp2 = coords[0][0] - coords[1][0];
    real_t tmp3 = coords[0][1] + tmp1;
    real_t tmp4 = coords[0][0] + tmp0;
    real_t tmp5 = coords[0][1] - coords[1][1];
    real_t tmp6 = 1.0/(tmp2*tmp3 - tmp4*tmp5);
    real_t tmp7 = 0.166666666666667*tmp6/(fabs(tmp6)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    out[0] = tmp7*(-DFinv(0,1)*(coords[1][1] + tmp1) + DFinv(1,1)*(coords[1][0] + tmp0));
    out[1] = tmp7*(DFinv(0,1)*tmp3 - DFinv(1,1)*tmp4);
    out[2] = tmp7*(-DFinv(0,1)*tmp5 + DFinv(1,1)*tmp2);
  }

  std::shared_ptr<FaceMap> faceMap;
};

}