#include "EdgeDoFOperator.hpp"

namespace hhg {

EdgeDoFOperator::EdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                 size_t minLevel,
                                 size_t maxLevel)
  : Operator(storage, minLevel, maxLevel)
{
  auto edgeVertexDoFToEdgeDoFDataHandling   = std::make_shared< MacroEdgeEdgeDoFToEdgeDoFDataHandling   >(minLevel_, maxLevel_);
  auto faceVertexDoFToEdgeDoFDataHandling   = std::make_shared< MacroFaceEdgeDoFToEdgeDoFDataHandling   >(minLevel_, maxLevel_);

  storage->addEdgeData(edgeStencilID_, edgeVertexDoFToEdgeDoFDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencilID_, faceVertexDoFToEdgeDoFDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
}

void
EdgeDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src, EdgeDoFFunction<real_t> &dst, uint_t level, DoFType flag, UpdateType updateType) {

  src.getCommunicator(level)->startCommunication<Face, Edge>();
  src.getCommunicator(level)->startCommunication<Edge, Face>();
  src.getCommunicator(level)->endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      edgedof::macroedge::apply(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
    }
  }

  src.getCommunicator(level)->endCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      edgedof::macroface::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
    }
  }


  WALBERLA_ABORT("implement me");
}


const PrimitiveDataID<StencilMemory<real_t>, Edge> &EdgeDoFOperator::getEdgeStencilID_() const {
  return edgeStencilID_;
}

const PrimitiveDataID<StencilMemory<real_t>, Face> &EdgeDoFOperator::getFaceStencilID_() const {
  return faceStencilID_;
}

}

