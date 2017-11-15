#pragma once

#include "VertexDoFToEdgeDoFMemory.hpp"
#include "VertexDoFToEdgeDoFDataHandling.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics.hpp"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif


namespace hhg
{

class VertexDoFToEdgeDoFOperator : public Operator<P1Function< real_t >, EdgeDoFFunction< real_t > >
{
public:
VertexDoFToEdgeDoFOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
  : Operator(storage, minLevel, maxLevel)
{
  auto faceVertexDoFToEdgeDoFStencilMemoryDataHandling = std::make_shared< FaceVertexDoFToEdgeDoFStencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);
  auto edgeVertexDoFToEdgeDoFStencilMemoryDataHandling = std::make_shared< EdgeVertexDoFToEdgeDoFStencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);
  auto vertexeVertexDoFToEdgeDoFStencilMemoryDataHandling = std::make_shared< VertexVertexDoFToEdgeDoFStencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);

  storage->addFaceData(faceStencilID_, faceVertexDoFToEdgeDoFStencilMemoryDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addEdgeData(edgeStencilID_, edgeVertexDoFToEdgeDoFStencilMemoryDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addVertexData(vertexStencilID_, vertexeVertexDoFToEdgeDoFStencilMemoryDataHandling, "VertexDoFToEdgeDoFOperatorVertexStencil");

  for (uint_t level = minLevel_; level <= maxLevel_; ++level)
  {
    WALBERLA_ABORT("implement me");
  }
}

~VertexDoFToEdgeDoFOperator()
{
}

void apply_impl(P1Function< real_t > & src, BubbleFunction< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace)
{
  WALBERLA_ABORT("implement me");
}

  const PrimitiveDataID<VertexVertexDoFToEdgeDoFStencilMemory< real_t >, Vertex> &getVertexStencilID() const { return vertexStencilID_; }

  const PrimitiveDataID<EdgeVertexDoFToEdgeDoFStencilMemory< real_t >, Edge> &getEdgeStencilID() const { return edgeStencilID_; }

  const PrimitiveDataID<FaceVertexDoFToEdgeDoFStencilMemory< real_t >, Face> &getFaceStencilID() const { return faceStencilID_; }


private:

  PrimitiveDataID<VertexVertexDoFToEdgeDoFStencilMemory< real_t >, Vertex> vertexStencilID_;
  PrimitiveDataID<EdgeVertexDoFToEdgeDoFStencilMemory< real_t >, Edge> edgeStencilID_;
  PrimitiveDataID<FaceVertexDoFToEdgeDoFStencilMemory< real_t >, Face> faceStencilID_;

void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[1][3], fenics::ElementType element_type) {

}
};

}
