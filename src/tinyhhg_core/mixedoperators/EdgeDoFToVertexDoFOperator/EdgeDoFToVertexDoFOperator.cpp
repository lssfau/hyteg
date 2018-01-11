#include "EdgeDoFToVertexDoFOperator.hpp"
#include "EdgeDoFToVertexDoFApply.hpp"

namespace hhg{

EdgeDoFToVertexDoFOperator::EdgeDoFToVertexDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                                       const size_t & minLevel,
                                                       const size_t & maxLevel)
  :Operator(storage,minLevel,maxLevel)
{

  using namespace EdgeDoFToVertexDoF;

  auto vertexDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Vertex >>(minLevel_, maxLevel_, macroVertexEdgeDoFToVertexDoFStencilSize);

  auto edgeDataHandling   =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge   >>(minLevel_, maxLevel_, macroEdgeEdgeDoFToVertexDoFStencilSize);

  auto faceDataHandling   =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face   >>(minLevel_, maxLevel_, macroFaceEdgeDoFToVertexDoFStencilSize);

  storage->addVertexData(vertexStencilID_, vertexDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addEdgeData(edgeStencilID_, edgeDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencilID_, faceDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  /// the stencil assembly will be done in the P2-Operator
}

real_t* EdgeDoFToVertexDoFOperator::getFaceStencil(const PrimitiveID& faceId, uint_t level) {
  return storage_->getFace(faceId)->getData(faceStencilID_)->getPointer( level );
}

void EdgeDoFToVertexDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src,
                                            P1Function<real_t> &dst,
                                            uint_t level,
                                            DoFType flag,
                                            UpdateType updateType)
{
  using namespace EdgeDoFToVertexDoF;
  src.getCommunicator(level)->startCommunication<Edge, Vertex>();
  src.getCommunicator(level)->startCommunication<Face, Edge>();
  src.getCommunicator(level)->endCommunication<Edge, Vertex>();

  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      applyVertex(level, vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), updateType);
    }
  }

  dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

  // end pulling edge halos
  src.getCommunicator(level)->endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      applyEdge(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
    }
  }

  dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

  dst.getCommunicator(level)->startCommunication<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      applyFace(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
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


namespace EdgeDoFToVertexDoF {
////////// Stencil sizes //////////
uint_t macroVertexEdgeDoFToVertexDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 * numDependencies;
}

uint_t macroEdgeEdgeDoFToVertexDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 + 5 * numDependencies;
}

uint_t macroFaceEdgeDoFToVertexDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(numDependencies);
  return 12;
}

}/// namespace EdgeDoFToVertexDoF
}/// namespace hhg
