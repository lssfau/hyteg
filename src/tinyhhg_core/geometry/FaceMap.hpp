#pragma once

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/matrix.hpp"

namespace hhg {

class FaceMap2D {
public:
  FaceMap2D(const Face& face) : face_(face) {}

  virtual ~FaceMap2D() {}

  virtual void evalF(const Point2D&, Point2D&) = 0;
  virtual void evalDF(const Point2D&, Matrix2r&) = 0;
  virtual real_t detDF(const Point2D&) = 0;

protected:
  const Face& face_;
};

}
