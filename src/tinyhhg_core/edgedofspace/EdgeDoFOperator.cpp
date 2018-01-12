#include "EdgeDoFOperator.hpp"

namespace hhg {

EdgeDoFOperator::EdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                 size_t minLevel,
                                 size_t maxLevel)
  : Operator(storage, minLevel, maxLevel)
{
  auto edgeDataHandling   =
      std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge   >>(minLevel_, maxLevel_, macroEdgeEdgeDoFToEdgeDoFStencilSize);

  auto faceDataHandling   =
      std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face   >>(minLevel_, maxLevel_, macroFaceEdgeDoFToEdgeDoFStencilSize);


  storage->addEdgeData(edgeStencilID_, edgeDataHandling  , "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addFaceData(faceStencilID_, faceDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
}

real_t* EdgeDoFOperator::getFaceStencil(const PrimitiveID& faceId, uint_t level) {
  return storage_->getFace(faceId)->getData(faceStencilID_)->getPointer( level );
}

void
EdgeDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src, EdgeDoFFunction<real_t> &dst, uint_t level, DoFType flag, UpdateType updateType) {

  src.getCommunicator(level)->startCommunication<Face, Edge>();
  src.getCommunicator(level)->startCommunication<Edge, Face>();
  src.getCommunicator(level)->endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

//    if (testFlag(edge.getDoFType(), flag))
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
}


const PrimitiveDataID<StencilMemory<real_t>, Edge> &EdgeDoFOperator::getEdgeStencilID_() const {
  return edgeStencilID_;
}

const PrimitiveDataID<StencilMemory<real_t>, Face> &EdgeDoFOperator::getFaceStencilID_() const {
  return faceStencilID_;
}

/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 1 + 2 * numDependencies;
}
/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(numDependencies);
  return 5 + 5 + 5;
}


}

