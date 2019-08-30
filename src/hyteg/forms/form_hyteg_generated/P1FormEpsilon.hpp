#pragma once

#include "hyteg/geometry/GeometryMap.hpp"

namespace hyteg {

class P1Form_epsilon_11 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    geometryMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][1];
    real_t tmp1 = coords[1][1] + tmp0;
    real_t tmp2 = -coords[2][0];
    real_t tmp3 = coords[1][0] + tmp2;
    real_t tmp4 = coords[0][0] - coords[1][0];
    real_t tmp5 = coords[0][1] + tmp0;
    real_t tmp6 = coords[0][0] + tmp2;
    real_t tmp7 = coords[0][1] - coords[1][1];
    real_t tmp8 = tmp4*tmp5 - tmp6*tmp7;
    real_t tmp9 = 1/(pow(tmp8, 2)*fabs(1.0/tmp8)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    real_t tmp10 = 1.0*DFinv(0,0)*tmp1 - 1.0*DFinv(1,0)*tmp3;
    real_t tmp11 = 0.5*DFinv(0,1)*tmp1 - 0.5*DFinv(1,1)*tmp3;
    out[0] = tmp9*(1.0*pow(DFinv(0,0)*tmp1 - DFinv(1,0)*tmp3, 2) + 0.5*pow(DFinv(0,1)*tmp1 - DFinv(1,1)*tmp3, 2));
    out[1] = -tmp9*(tmp10*(DFinv(0,0)*tmp5 - DFinv(1,0)*tmp6) + tmp11*(DFinv(0,1)*tmp5 - DFinv(1,1)*tmp6));
    out[2] = tmp9*(tmp10*(DFinv(0,0)*tmp7 - DFinv(1,0)*tmp4) + tmp11*(DFinv(0,1)*tmp7 - DFinv(1,1)*tmp4));
  }

  std::shared_ptr<GeometryMap> geometryMap;
};

class P1Form_epsilon_12 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    geometryMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][1];
    real_t tmp1 = coords[1][1] + tmp0;
    real_t tmp2 = -coords[2][0];
    real_t tmp3 = coords[1][0] + tmp2;
    real_t tmp4 = coords[0][0] - coords[1][0];
    real_t tmp5 = coords[0][1] + tmp0;
    real_t tmp6 = coords[0][0] + tmp2;
    real_t tmp7 = coords[0][1] - coords[1][1];
    real_t tmp8 = tmp4*tmp5 - tmp6*tmp7;
    real_t tmp9 = 0.5*(DFinv(0,1)*tmp1 - DFinv(1,1)*tmp3)/(pow(tmp8, 2)*fabs(1.0/tmp8)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    out[0] = tmp9*(DFinv(0,0)*tmp1 - DFinv(1,0)*tmp3);
    out[1] = -tmp9*(DFinv(0,0)*tmp5 - DFinv(1,0)*tmp6);
    out[2] = tmp9*(DFinv(0,0)*tmp7 - DFinv(1,0)*tmp4);
  }

  std::shared_ptr<GeometryMap> geometryMap;
};

class P1Form_epsilon_21 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    geometryMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][1];
    real_t tmp1 = coords[1][1] + tmp0;
    real_t tmp2 = -coords[2][0];
    real_t tmp3 = coords[1][0] + tmp2;
    real_t tmp4 = coords[0][0] - coords[1][0];
    real_t tmp5 = coords[0][1] + tmp0;
    real_t tmp6 = coords[0][0] + tmp2;
    real_t tmp7 = coords[0][1] - coords[1][1];
    real_t tmp8 = tmp4*tmp5 - tmp6*tmp7;
    real_t tmp9 = 0.5*(DFinv(0,0)*tmp1 - DFinv(1,0)*tmp3)/(pow(tmp8, 2)*fabs(1.0/tmp8)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    out[0] = tmp9*(DFinv(0,1)*tmp1 - DFinv(1,1)*tmp3);
    out[1] = -tmp9*(DFinv(0,1)*tmp5 - DFinv(1,1)*tmp6);
    out[2] = tmp9*(DFinv(0,1)*tmp7 - DFinv(1,1)*tmp4);
  }

  std::shared_ptr<GeometryMap> geometryMap;
};

class P1Form_epsilon_22 {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    geometryMap->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][1];
    real_t tmp1 = coords[1][1] + tmp0;
    real_t tmp2 = -coords[2][0];
    real_t tmp3 = coords[1][0] + tmp2;
    real_t tmp4 = coords[0][0] - coords[1][0];
    real_t tmp5 = coords[0][1] + tmp0;
    real_t tmp6 = coords[0][0] + tmp2;
    real_t tmp7 = coords[0][1] - coords[1][1];
    real_t tmp8 = tmp4*tmp5 - tmp6*tmp7;
    real_t tmp9 = 1/(pow(tmp8, 2)*fabs(1.0/tmp8)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    real_t tmp10 = 0.5*DFinv(0,0)*tmp1 - 0.5*DFinv(1,0)*tmp3;
    real_t tmp11 = 1.0*DFinv(0,1)*tmp1 - 1.0*DFinv(1,1)*tmp3;
    out[0] = tmp9*(0.5*pow(DFinv(0,0)*tmp1 - DFinv(1,0)*tmp3, 2) + 1.0*pow(DFinv(0,1)*tmp1 - DFinv(1,1)*tmp3, 2));
    out[1] = -tmp9*(tmp10*(DFinv(0,0)*tmp5 - DFinv(1,0)*tmp6) + tmp11*(DFinv(0,1)*tmp5 - DFinv(1,1)*tmp6));
    out[2] = tmp9*(tmp10*(DFinv(0,0)*tmp7 - DFinv(1,0)*tmp4) + tmp11*(DFinv(0,1)*tmp7 - DFinv(1,1)*tmp4));
  }

  std::shared_ptr<GeometryMap> geometryMap;
};

}