#pragma once

#include "FaceMap.hpp"

namespace hhg {

class IdentityMap : public FaceMap2D {
public:

  IdentityMap(const Face& face) : FaceMap2D(face) {}

  void evalF(const Point2D& x, Point2D& Fx) {
    Fx = x;
  }

  void evalDF(const Point2D& x, Matrix2r& DFx) {
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
