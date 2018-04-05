#pragma once

#include "tinyhhg_core/geometry/FaceMap.hpp"

namespace hhg {

class P1Form_divT_1 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    faceMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][0];
    real_t tmp1 = -coords[2][1];
    real_t tmp2 = 1.0/((coords[0][0] - coords[1][0])*(coords[0][1] + tmp1) - (coords[0][0] + tmp0)*(coords[0][1] - coords[1][1]));
    real_t tmp3 = tmp2*(-DFinv(0,0)*(coords[1][1] + tmp1) + DFinv(1,0)*(coords[1][0] + tmp0))/(fabs(tmp2)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    real_t tmp4 = 0.166666666666667*tmp3;
    out[0] = 0.166666666666667*tmp3;
    out[1] = tmp4;
    out[2] = tmp4;
  }

  std::shared_ptr<FaceMap> faceMap;
};

class P1Form_divT_2 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    faceMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][0];
    real_t tmp1 = -coords[2][1];
    real_t tmp2 = 1.0/((coords[0][0] - coords[1][0])*(coords[0][1] + tmp1) - (coords[0][0] + tmp0)*(coords[0][1] - coords[1][1]));
    real_t tmp3 = tmp2*(-DFinv(0,1)*(coords[1][1] + tmp1) + DFinv(1,1)*(coords[1][0] + tmp0))/(fabs(tmp2)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    real_t tmp4 = 0.166666666666667*tmp3;
    out[0] = 0.166666666666667*tmp3;
    out[1] = tmp4;
    out[2] = tmp4;
  }

  std::shared_ptr<FaceMap> faceMap;
};

}
