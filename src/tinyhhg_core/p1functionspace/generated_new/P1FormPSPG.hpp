#pragma once

#include "tinyhhg_core/geometry/GeometryMap.hpp"

namespace hhg {

class P1Form_pspg {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
     Point3D x_hat({0.333333333333333, 0.333333333333333});
     Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
     Matrix2r DFinv;
     geometryMap->evalDFinv(x_tilde, DFinv);
     real_t tmp0 = coords[0][0] - coords[1][0];
     real_t tmp1 = coords[0][1] - coords[2][1];
     real_t tmp2 = coords[0][0] - coords[2][0];
     real_t tmp3 = coords[0][1] - coords[1][1];
     real_t tmp4 = tmp0*tmp1 - tmp2*tmp3;
     real_t tmp5 = 1/(pow(tmp4, 2)*pow(fabs(1.0/tmp4), 2)*pow(fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)), 2));
     real_t tmp6 = DFinv(0,0)*tmp3 - DFinv(1,0)*tmp0;
     real_t tmp7 = DFinv(1,0)*tmp2;
     real_t tmp8 = DFinv(0,0)*tmp1;
     real_t tmp9 = tmp6 + tmp7 - tmp8;
     real_t tmp10 = DFinv(0,1)*tmp3 - DFinv(1,1)*tmp0;
     real_t tmp11 = DFinv(1,1)*tmp2;
     real_t tmp12 = DFinv(0,1)*tmp1;
     real_t tmp13 = tmp10 + tmp11 - tmp12;
     out[0] = -tmp5*(0.1*pow(tmp13, 2) + 0.1*pow(tmp9, 2));
     out[1] = -tmp5*(0.1*tmp13*(-tmp11 + tmp12) + 0.1*tmp9*(-tmp7 + tmp8));
     out[2] = tmp5*(0.1*tmp10*tmp13 + 0.1*tmp6*tmp9);
  }

  std::shared_ptr<GeometryMap> geometryMap;
};

}