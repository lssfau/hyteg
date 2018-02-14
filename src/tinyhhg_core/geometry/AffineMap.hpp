#pragma once

#include "FaceMap.hpp"

namespace hhg {

class AffineMap : public FaceMap {
public:

  AffineMap(const std::array<Point3D, 3>& coords) {
    c = coords[0];
    dx = coords[1] - coords[0];
    dy = coords[2] - coords[0];
  }

  void evalF(const Point3D& xRef, Point3D& xPhy) {
    xPhy[0] = c[0] + dx[0] * xRef[0] + dy[0] * xRef[1];
    xPhy[1] = c[1] + dx[1] * xRef[0] + dy[1] * xRef[1];
  }

  void evalDF(const Point2D&, Matrix2r& DFx) {
    DFx(0,0) = dx[0];
    DFx(0,1) = dy[0];
    DFx(1,0) = dx[1];
    DFx(1,1) = dy[1];
  }

  real_t detDF(const Point2D& x) {
    return dx[0] * dy[1] - dy[0] * dx[1];
  }

private:
  Point3D c;
  Point3D dx;
  Point3D dy;
};

}
