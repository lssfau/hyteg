#pragma once

#include "VertexDoFToEdgeDoFMemory.hpp"
#include "VertexDoFToEdgeDoFDataHandling.hpp"
#include "VertexDoFToEdgeDoFFace.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

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
  /// since the Vertex does not own any EdgeDoFs only edge and face are needed
    auto faceVertexDoFToEdgeDoFDataHandling = std::make_shared< MacroFaceVertexDoFToEdgeDoFDataHandling >(minLevel_, maxLevel_);
    auto edgeVertexDoFToEdgeDoFDataHandling = std::make_shared< MacroEdgeVertexDoFToEdgeDoFDataHandling >(minLevel_, maxLevel_);

    storage->addEdgeData(edgeStencilID_, edgeVertexDoFToEdgeDoFDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
    storage->addFaceData(faceStencilID_, faceVertexDoFToEdgeDoFDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil");
  }

  ~VertexDoFToEdgeDoFOperator() = default;

  void apply_impl(P1Function< real_t > & src, EdgeDoFFunction< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    WALBERLA_ABORT("implement me");
  }

  /// since the Vertex does not own any EdgeDoFs only edge and face are needed
  const PrimitiveDataID<EdgeVertexDoFToEdgeDoFStencilMemory< real_t >, Edge> &getEdgeStencilID() const { return edgeStencilID_; }
  const PrimitiveDataID<FaceVertexDoFToEdgeDoFStencilMemory< real_t >, Face> &getFaceStencilID() const { return faceStencilID_; }

private:
  PrimitiveDataID<EdgeVertexDoFToEdgeDoFStencilMemory< real_t >, Edge> edgeStencilID_;
  PrimitiveDataID<FaceVertexDoFToEdgeDoFStencilMemory< real_t >, Face> faceStencilID_;

};

}
