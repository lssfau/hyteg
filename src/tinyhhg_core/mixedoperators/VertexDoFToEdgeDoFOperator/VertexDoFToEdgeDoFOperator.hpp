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
    auto faceVertexDoFToEdgeDoFStencilMemoryDataHandling = std::make_shared< FaceVertexDoFToEdgeDoFStencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);
    auto edgeVertexDoFToEdgeDoFStencilMemoryDataHandling = std::make_shared< EdgeVertexDoFToEdgeDoFStencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);

    storage->addEdgeData(edgeStencilID_, edgeVertexDoFToEdgeDoFStencilMemoryDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
    storage->addFaceData(faceStencilID_, faceVertexDoFToEdgeDoFStencilMemoryDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil");

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {
      for(auto& it : storage_->getFaces()) {
        Face &face = *it.second;
        auto face_stencil = face.getData(faceStencilID_)->getPointer(level);
        for(uint_t i = 1; i <= face.getData(faceStencilID_)->getSize(level); ++i){
          WALBERLA_LOG_DEVEL("this is wrong!");
          face_stencil[i] = i;
        }
      }
      for(auto& it : storage_->getEdges()) {
        Edge &edge = *it.second;
        auto edge_stencil = edge.getData(edgeStencilID_)->getPointer(level);
        for(uint_t i = 1; i <= edge.getData(edgeStencilID_)->getSize(level); ++i){
          WALBERLA_LOG_DEVEL("this is wrong!");
          edge_stencil[i] = i;
        }
      }


      WALBERLA_LOG_DEVEL("implement me");
    }
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

void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[1][3], fenics::ElementType element_type) {

}
};

}
