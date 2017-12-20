#include "EdgeDoFToVertexDoFOperator.hpp"
#include "EdgeDoFToVertexDoFDataHandling.hpp"
#include "EdgeDoFToVertexDoFApply.hpp"

namespace hhg{

EdgeDoFToVertexDoFOperator::EdgeDoFToVertexDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                                       size_t minLevel,
                                                       size_t maxLevel)
  :Operator(storage,minLevel,maxLevel)
{

  auto vertexVertexDoFToEdgeDoFDataHandling = std::make_shared< MacroVertexEdgeDoFToVertexDoFDataHandling >(minLevel_, maxLevel_);
  auto edgeVertexDoFToEdgeDoFDataHandling   = std::make_shared< MacroEdgeEdgeDoFToVertexDoFDataHandling   >(minLevel_, maxLevel_);
  auto faceVertexDoFToEdgeDoFDataHandling   = std::make_shared< MacroFaceEdgeDoFToVertexDoFDataHandling   >(minLevel_, maxLevel_);

  storage->addVertexData(vertexStencilID_, vertexVertexDoFToEdgeDoFDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addEdgeData(edgeStencilID_, edgeVertexDoFToEdgeDoFDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencilID_, faceVertexDoFToEdgeDoFDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  /// the stencil assembly will be done in the P2-Operator
}

void EdgeDoFToVertexDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src,
                                            P1Function<real_t> &dst,
                                            uint_t level,
                                            DoFType flag,
                                            UpdateType updateType)
{

  src.getCommunicator(level)->startCommunication<Edge, Vertex>();
  src.getCommunicator(level)->startCommunication<Face, Edge>();
  src.getCommunicator(level)->endCommunication<Edge, Vertex>();

  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      EdgeDoFToVertexDoFVertex::apply(level, vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), updateType);
    }
  }

  dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

  // end pulling edge halos
  src.getCommunicator(level)->endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      EdgeDoFToVertexDoFEdge::apply(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
    }
  }

  dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

  dst.getCommunicator(level)->startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      EdgeDoFToVertexDoFFace::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
    }
  }

  dst.getCommunicator(level)->endCommunication<Edge, Face>();
}

const PrimitiveDataID<StencilMemory< real_t >, Vertex > &EdgeDoFToVertexDoFOperator::getVertexStencilID() const {
  return vertexStencilID_;
}

const PrimitiveDataID<StencilMemory< real_t >, Edge > &EdgeDoFToVertexDoFOperator::getEdgeStencilID() const {
  return edgeStencilID_;
}

const PrimitiveDataID<StencilMemory< real_t >, Face > &EdgeDoFToVertexDoFOperator::getFaceStencilID() const {
  return faceStencilID_;
}


}
