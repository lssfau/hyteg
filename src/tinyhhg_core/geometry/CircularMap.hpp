#pragma once

#include "FaceMap.hpp"

namespace hhg {

class CircularMap : public FaceMap {
public:

  CircularMap(const Face& face, const std::shared_ptr<PrimitiveStorage>& storage, const Point2D& center, real_t radius)
      : center_(center), radius_(radius)
  {
    // Get edge on boundary
    // TODO: ATTENTION - ONLY ONE EDGE MAY LIE ON THE BOUNDARY
    WALBERLA_ASSERT_EQUAL(face.hasBoundaryEdge(), true);
    WALBERLA_ASSERT_EQUAL(face.edgesOnBoundary.size(), 1);

    const Edge& edge = *storage->getEdge(face.edgesOnBoundary[0]);
    const Vertex& vertex = *storage->getVertex(face.get_vertex_opposite_to_edge(edge.getID()));

    x1_[0] = edge.getCoordinates()[0][0];
    x1_[1] = edge.getCoordinates()[0][1];
    x2_[0] = vertex.getCoordinates()[0];
    x2_[1] = vertex.getCoordinates()[1];
    x3_[0] = edge.getCoordinates()[1][0];
    x3_[1] = edge.getCoordinates()[1][1];

    s1_ = std::atan2((x1_-center)[1], (x1_-center)[0]);
    real_t s3 = std::atan2((x3_-center)[1], (x3_-center)[0]);
    s3bar_ = s3 - s1_;

    if (s3bar_ < -0.5 * walberla::math::PI) {
      s3bar_ += 2.0 * walberla::math::PI;
    } else if (s3bar_ > 0.5 * walberla::math::PI) {
      s3bar_ -= 2.0 * walberla::math::PI;
    }
  }

  void evalF(const Point3D& xRef, Point3D& Fx) {
    Fx[0] =  -xRef[1]*(x1_[0] - x3_[0]) + x1_[0] - xRef[0]*(x1_[0] - x2_[0]);
    Fx[1] =  (-xRef[1]*(x1_[1] - x3_[1]) - xRef[0]*(x1_[1] - x2_[1]) + x1_[1]);

    if (std::abs(xRef[1]-1.0) > 1e-12) {
      Fx[0] += ((xRef[1] + xRef[0] - 1)*(xRef[1]*(x1_[0] - x3_[0]) + radius_*cos(xRef[1]*s3bar_ + s1_) - x1_[0] + center_[0]))/(xRef[1] - 1);
      Fx[1] += ((xRef[1] + xRef[0] - 1)*(xRef[1]*(x1_[1] - x3_[1]) + radius_*sin(xRef[1]*s3bar_ + s1_) - x1_[1] + center_[1]))/(xRef[1] - 1);
    }
  }

  void evalDF(const Point2D& xRef, Matrix2r& DFx) {
    DFx(0,0) =  (xRef[1]*x2_[0] - xRef[1]*x3_[0] + radius_*cos(xRef[1]*s3bar_ + s1_) - x2_[0] + center_[0])/(xRef[1] - 1);
    DFx(0,1) =  ((xRef[1] - 1)*(xRef[1]*(x1_[0] - x3_[0]) - x1_[0] + xRef[0]*(x1_[0] - x2_[0])) + (xRef[1] - 1)*(radius_*cos(xRef[1]*s3bar_ + s1_) + center_[0] - xRef[0]*(x1_[0] - x2_[0]) - (xRef[1] - 1)*(x1_[0] - x3_[0]) + (-xRef[1] - xRef[0] + 1)*(radius_*s3bar_*sin(xRef[1]*s3bar_ + s1_) - x1_[0] + x3_[0])) + (-xRef[1] - xRef[0] + 1)*(xRef[1]*(x1_[0] - x3_[0]) + radius_*cos(xRef[1]*s3bar_ + s1_) - x1_[0] + center_[0]))/pow(xRef[1] - 1, 2);
    DFx(1,0) =  (xRef[1]*x2_[1] - xRef[1]*x3_[1] + radius_*sin(xRef[1]*s3bar_ + s1_) - x2_[1] + center_[1])/(xRef[1] - 1);
    DFx(1,1) =  ((xRef[1] - 1)*(xRef[1]*(x1_[1] - x3_[1]) + xRef[0]*(x1_[1] - x2_[1]) - x1_[1]) + (xRef[1] - 1)*(radius_*sin(xRef[1]*s3bar_ + s1_) - xRef[0]*(x1_[1] - x2_[1]) + center_[1] - (xRef[1] - 1)*(x1_[1] - x3_[1]) + (xRef[1] + xRef[0] - 1)*(radius_*s3bar_*cos(xRef[1]*s3bar_ + s1_) + x1_[1] - x3_[1])) + (-xRef[1] - xRef[0] + 1)*(xRef[1]*(x1_[1] - x3_[1]) + radius_*sin(xRef[1]*s3bar_ + s1_) - x1_[1] + center_[1]))/pow(xRef[1] - 1, 2);
  }

  real_t detDF(const Point2D& x) {
    WALBERLA_ABORT("not implemented")
    return 1.0;
  }

private:

  Point2D x1_;
  Point2D x2_;
  Point2D x3_;

  Point2D center_;
  real_t radius_;
  real_t s1_;
  real_t s3bar_;
};

}
