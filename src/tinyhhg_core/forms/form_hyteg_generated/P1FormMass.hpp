#pragma once

#include "tinyhhg_core/geometry/GeometryMap.hpp"

namespace hhg {

class P1Form_mass {
public:
  void evalQuadraturePoint(const Point3D& x_hat, const std::array<Point3D,3>& coords, const Matrix2r DFinv, real_t w, Point3D& out) const
  {
     real_t tmp0 = x_hat[0] + x_hat[1] - 1;
     real_t tmp1 = 1.0/fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0));
     real_t tmp2 = 1.0/fabs(1.0/((coords[0][0] - coords[1][0])*(coords[0][1] - coords[2][1]) - (coords[0][0] - coords[2][0])*(coords[0][1] - coords[1][1])));
     real_t tmp3 = tmp0*tmp1*tmp2;
     out[0] += w * pow(tmp0, 2)*tmp1*tmp2;
     out[1] += w * -tmp3*x_hat[0];
     out[2] += w * -tmp3*x_hat[1];
  }

  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
     Point3D x_hat;
     Point3D x_tilde;
     Matrix2r DFinv;
     out[0] = 0;
     out[1] = 0;
     out[2] = 0;
     x_hat[0] = 0.16666666666666666;
     x_hat[1] = 0.16666666666666666;
     x_tilde[0] = 0.666666666666667*coords[0][0] + 0.166666666666667*coords[1][0] + 0.166666666666667*coords[2][0];
     x_tilde[1] = 0.666666666666667*coords[0][1] + 0.166666666666667*coords[1][1] + 0.166666666666667*coords[2][1];
     geometryMap->evalDFinv(x_tilde, DFinv);
     evalQuadraturePoint(x_hat, coords, DFinv, 0.16666666666666666, out);
     x_hat[0] = 0.6666666666666666;
     x_hat[1] = 0.16666666666666666;
     x_tilde[0] = 0.166666666666667*coords[0][0] + 0.666666666666667*coords[1][0] + 0.166666666666667*coords[2][0];
     x_tilde[1] = 0.166666666666667*coords[0][1] + 0.666666666666667*coords[1][1] + 0.166666666666667*coords[2][1];
     geometryMap->evalDFinv(x_tilde, DFinv);
     evalQuadraturePoint(x_hat, coords, DFinv, 0.16666666666666666, out);
     x_hat[0] = 0.16666666666666666;
     x_hat[1] = 0.6666666666666666;
     x_tilde[0] = 0.166666666666667*coords[0][0] + 0.166666666666667*coords[1][0] + 0.666666666666667*coords[2][0];
     x_tilde[1] = 0.166666666666667*coords[0][1] + 0.166666666666667*coords[1][1] + 0.666666666666667*coords[2][1];
     geometryMap->evalDFinv(x_tilde, DFinv);
     evalQuadraturePoint(x_hat, coords, DFinv, 0.16666666666666666, out);
  }

  std::shared_ptr<GeometryMap> geometryMap;
};

}