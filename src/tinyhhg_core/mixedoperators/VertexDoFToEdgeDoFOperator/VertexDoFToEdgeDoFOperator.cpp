#include "VertexDoFToEdgeDoFOperator.hpp"

namespace hhg {

VertexDoFToEdgeDoFOperator::VertexDoFToEdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel)
  : Operator(storage, minLevel, maxLevel) {
  /// since the Vertex does not own any EdgeDoFs only edge and face are needed

  auto edgeDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge >>(minLevel_,
                                                                        maxLevel_,
                                                                        VertexDoFToEdgeDoF::macroEdgeVertexDoFToEdgeDoFStencilSize);

  auto faceDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face >>(minLevel_,
                                                                        maxLevel_,
                                                                        VertexDoFToEdgeDoF::macroFaceVertexDoFToEdgeDoFStencilSize);

  storage->addEdgeData(edgeStencilID_, edgeDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addFaceData(faceStencilID_, faceDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil");
}

real_t* VertexDoFToEdgeDoFOperator::getEdgeStencil(const PrimitiveID& edgeId, uint_t level) {
  return storage_->getEdge(edgeId)->getData(edgeStencilID_)->getPointer( level );
}

real_t* VertexDoFToEdgeDoFOperator::getFaceStencil(const PrimitiveID& faceId, uint_t level) {
  return storage_->getFace(faceId)->getData(faceStencilID_)->getPointer( level );
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

namespace VertexDoFToEdgeDoF {

///
/// \param level stencil size is independent of level
/// \param numDependencies number of adjacent faces of the edge
/// \return number of the stencil entries
uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  return 2 + numDependencies;
}

uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( numDependencies );
  return 4 + 4 + 4;
}
}


}/// namespace hhg
