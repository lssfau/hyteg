#include "EdgeDoFOperator.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

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

void
EdgeDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src, EdgeDoFFunction<real_t> &dst, uint_t level, DoFType flag, UpdateType updateType) {

  this->startTiming( "EdgeDoFOperator - Apply" );
  src.startCommunication<Face, Edge>( level );
  src.startCommunication<Edge, Face>( level );
  src.endCommunication<Face, Edge>( level );

  for (auto& it : storage_->getEdges())
  {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if ( testFlag( edgeBC, flag ) )
    {
      edgedof::macroedge::apply(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
    }
  }

  src.endCommunication<Edge, Face>( level );

  for (auto& it : storage_->getFaces())
  {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if ( testFlag( faceBC, flag ) )
    {
      edgedof::macroface::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
    }
  }
  this->stopTiming( "EdgeDoFOperator - Apply" );
}


const PrimitiveDataID<StencilMemory<real_t>, Edge> &EdgeDoFOperator::getEdgeStencilID() const {
  return edgeStencilID_;
}

const PrimitiveDataID<StencilMemory<real_t>, Face> &EdgeDoFOperator::getFaceStencilID() const {
  return faceStencilID_;
}

/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return 1 + 2 * primitive.getNumNeighborFaces();
}

/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 5 + 5 + 5;
}


}

