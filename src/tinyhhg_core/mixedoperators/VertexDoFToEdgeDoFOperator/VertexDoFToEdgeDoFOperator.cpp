#include "VertexDoFToEdgeDoFOperator.hpp"

namespace hhg {

VertexDoFToEdgeDoFOperator::VertexDoFToEdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel)
  : Operator(storage, minLevel, maxLevel) {
  /// since the Vertex does not own any EdgeDoFs only edge and face are needed
  auto faceVertexDoFToEdgeDoFDataHandling = std::make_shared<MacroFaceVertexDoFToEdgeDoFDataHandling>(minLevel_, maxLevel_);
  auto edgeVertexDoFToEdgeDoFDataHandling = std::make_shared<MacroEdgeVertexDoFToEdgeDoFDataHandling>(minLevel_, maxLevel_);

  storage->addEdgeData(edgeStencilID_, edgeVertexDoFToEdgeDoFDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addFaceData(faceStencilID_, faceVertexDoFToEdgeDoFDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil");
}

void VertexDoFToEdgeDoFOperator::apply_impl(P1Function<real_t> &src, EdgeDoFFunction<real_t> &dst, size_t level, DoFType flag,
                                            UpdateType updateType) {

  src.getCommunicator(level)->startCommunication<Face, Edge>();
  src.getCommunicator(level)->startCommunication<Edge, Face>();
  src.getCommunicator(level)->endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      VertexDoFToEdgeDoF::applyEdge(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
    }
  }

  src.getCommunicator(level)->endCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      VertexDoFToEdgeDoF::applyFace(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
    }
  }

}

}/// namespace hhg