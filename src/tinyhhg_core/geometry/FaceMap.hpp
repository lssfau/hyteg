#pragma once

#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/matrix.hpp"

namespace hhg {

class FaceMap {
public:
  virtual ~FaceMap() {}

  virtual void evalF(const Point3D&, Point3D&) = 0;
  virtual void evalDF(const Point3D&, Matrix2r&) = 0;
  virtual real_t detDF(const Point2D&) = 0;
};

}
