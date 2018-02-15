#pragma once

#include "FaceMap.hpp"

namespace hhg {

class CircularMap : public FaceMap {
 public:

  CircularMap(const Face& face, const std::shared_ptr<PrimitiveStorage>& storage, const Point3D& center, real_t radius)
      : center_(center), radius_(radius)
  {
    // Get edge on boundary
    // TODO: ATTENTION - ONLY ONE EDGE MAY LIE ON THE BOUNDARY
    WALBERLA_ASSERT_EQUAL(face.hasBoundaryEdge(), true);
    WALBERLA_ASSERT_EQUAL(face.edgesOnBoundary.size(), 1);

    const Edge& edge = *storage->getEdge(face.edgesOnBoundary[0]);
    const Vertex& vertex = *storage->getVertex(face.get_vertex_opposite_to_edge(edge.getID()));

    Point3D x2, x3;

    x1_ = edge.getCoordinates()[0];
    x2 = vertex.getCoordinates();
    x3 = edge.getCoordinates()[1];

    x2bar_ = x2 - x1_;
    x3bar_ = x3 - x1_;

    s1_ = std::atan2((x1_-center)[1], (x1_-center)[0]);
    real_t s3 = std::atan2((x3-center)[1], (x3-center)[0]);
    s3bar_ = s3 - s1_;

    if (s3bar_ < -0.5 * walberla::math::PI) {
      s3bar_ += 2.0 * walberla::math::PI;
    } else if (s3bar_ > 0.5 * walberla::math::PI) {
      s3bar_ -= 2.0 * walberla::math::PI;
    }
  }

  void evalF(const Point3D& x, Point3D& Fx) {
    Fx[0] =  x1_[0] + x2bar_[0]*(-x3bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])) + x3bar_[0]*(x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])) + (radius_*cos(s1_ + s3bar_*(x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[0] - x3bar_[0]*(x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])) + center_[0])*(-x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/(-x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1);
    Fx[1] =  x1_[1] + x2bar_[1]*(-x3bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])) + x3bar_[1]*(x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])) + (radius_*sin(s1_ + s3bar_*(x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[1] - x3bar_[1]*(x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])) + center_[1])*(-x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/(-x2bar_[0]*(x[1] - x1_[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(x[0] - x1_[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1);
  }

  void evalDF(const Point3D& x, Matrix2r& DFx) {
    DFx(0,0) =  x2bar_[0]*x3bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[0]*x2bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(center_[0] + radius_*cos(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[0] - x3bar_[0]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/((x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])*pow(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1, 2)) + (x2bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(center_[0] + radius_*cos(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[0] - x3bar_[0]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1) + (radius_*s3bar_*x2bar_[1]*sin(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*x2bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1);
    DFx(0,1) =  x2bar_[0]*(center_[0] + radius_*cos(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[0] - x3bar_[0]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/((x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])*pow(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1, 2)) + (-x2bar_[0]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(center_[0] + radius_*cos(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[0] - x3bar_[0]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1) + (-radius_*s3bar_*x2bar_[0]*sin(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[0]*x3bar_[0]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1);
    DFx(1,0) =  -x2bar_[1]*(center_[1] + radius_*sin(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[1] - x3bar_[1]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/((x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])*pow(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1, 2)) + (x2bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(center_[1] + radius_*sin(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[1] - x3bar_[1]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1) + (-radius_*s3bar_*x2bar_[1]*cos(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*x3bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1);
    DFx(1,1) =  x2bar_[0]*x3bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[0]*(center_[1] + radius_*sin(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[1] - x3bar_[1]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/((x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])*pow(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1, 2)) - x3bar_[0]*x2bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + (-x2bar_[0]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(center_[1] + radius_*sin(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))) - x1_[1] - x3bar_[1]*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1) + (radius_*s3bar_*x2bar_[0]*cos(s1_ + s3bar_*(x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1])))/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x2bar_[0]*x3bar_[1]/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]))*(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x3bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) - x3bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1)/(-x2bar_[0]*(-x1_[1] + x[1])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + x2bar_[1]*(-x1_[0] + x[0])/(x2bar_[0]*x3bar_[1] - x3bar_[0]*x2bar_[1]) + 1);
  }

  real_t detDF(const Point2D& x) {
    return 1.0;
  }

 private:

  Point3D x1_;
  Point3D x2bar_;
  Point3D x3bar_;

  Point3D center_;
  real_t radius_;
  real_t s1_;
  real_t s3bar_;
};

}
