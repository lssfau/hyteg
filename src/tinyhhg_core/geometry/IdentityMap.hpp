#pragma once

#include "FaceMap.hpp"

namespace hhg {

class IdentityMap : public FaceMap {
public:

  IdentityMap() {}

  void evalF(const Point3D& x, Point3D& Fx) {
    Fx = x;
  }

  void evalDF(const Point3D& x, Matrix2r& DFx) {
    DFx(0,0) = 1.0;
    DFx(0,1) = 0.0;
    DFx(1,0) = 0.0;
    DFx(1,1) = 1.0;
  }

  real_t detDF(const Point2D& x) {
    return 1.0;
  }
};

}
