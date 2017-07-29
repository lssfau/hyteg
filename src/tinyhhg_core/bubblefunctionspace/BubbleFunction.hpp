#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {

class VertexBubbleFunctionMemory;
class EdgeBubbleFunctionMemory;
class FaceBubbleFunctionMemory;

class BubbleFunction : public Function {
 public:
  BubbleFunction(const std::string &name,
                 const std::shared_ptr<PrimitiveStorage> &storage,
                 uint_t minLevel,
                 uint_t maxLevel);

  ~BubbleFunction();

  /// Interpolates a given expression to a P1Function
  void interpolate(std::function<real_t(const Point3D &)> &expr, uint_t level, DoFType flag = All);

  void assign(const std::vector<walberla::real_t> scalars,
              const std::vector<BubbleFunction *> functions,
              size_t level,
              DoFType flag = All);

  void add(const std::vector<walberla::real_t> scalars,
           const std::vector<BubbleFunction *> functions,
           size_t level,
           DoFType flag = All);

  real_t dot(BubbleFunction &rhs, size_t level, DoFType flag = All);

//  void prolongate(size_t level, DoFType flag = All);

//  void restrict(size_t level, DoFType flag = All);

  const PrimitiveDataID<VertexBubbleFunctionMemory, Vertex> &getVertexDataID() const { return vertexDataID_; }

  const PrimitiveDataID<EdgeBubbleFunctionMemory, Edge> &getEdgeDataID() const { return edgeDataID_; }

  const PrimitiveDataID<FaceBubbleFunctionMemory, Face> &getFaceDataID() const { return faceDataID_; }

 private:
  PrimitiveDataID<VertexBubbleFunctionMemory, Vertex> vertexDataID_;
  PrimitiveDataID<EdgeBubbleFunctionMemory, Edge> edgeDataID_;
  PrimitiveDataID<FaceBubbleFunctionMemory, Face> faceDataID_;
};
}