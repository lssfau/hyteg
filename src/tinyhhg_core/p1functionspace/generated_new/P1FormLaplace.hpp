#pragma once

#include "tinyhhg_core/geometry/FaceMap.hpp"

namespace hhg {

class P1FormLaplace {

public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D xRef = 1.0/3.0 * (coords[0] + coords[1] + coords[2]);
    Point3D x;
    Matrix2r DFinvT;
    faceMap->evalF(xRef, x);
    faceMap->evalDFinv(xRef, DFinvT);
    real_t tmp0 = 1.0/fabs(DFinvT(0,0)*DFinvT(1,1) - DFinvT(0,1)*DFinvT(1,0));
    real_t tmp1 = 0.5*tmp0;
    real_t tmp2 = 0.5*DFinvT(0,0) + 0.5*DFinvT(0,1);
    real_t tmp3 = 0.5*DFinvT(1,0) + 0.5*DFinvT(1,1);
    real_t tmp4 = -tmp0*(DFinvT(0,0)*tmp2 + DFinvT(1,0)*tmp3);
    real_t tmp5 = -tmp0*(DFinvT(0,1)*tmp2 + DFinvT(1,1)*tmp3);
    real_t tmp6 = tmp1*(DFinvT(0,0)*DFinvT(0,1) + DFinvT(1,0)*DFinvT(1,1));
    out[0] = tmp1*(pow(DFinvT(0,0) + DFinvT(0,1), 2) + pow(DFinvT(1,0) + DFinvT(1,1), 2));
    out[1] = tmp4;
    out[2] = tmp5;
  }

  std::shared_ptr<FaceMap> faceMap;

};

}