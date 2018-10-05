#include "tinyhhg_core/HHGDefinitions.hpp"
#include "EdgeDoFOperator.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "generatedKernels/generatedKernels.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroCell.hpp"

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

  auto face3DDataHandling   =
      std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >>( minLevel_, maxLevel_ );

  auto cellDataHandling   =
      std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >>( minLevel_, maxLevel_ );

  storage->addEdgeData(edgeStencilID_, edgeDataHandling  , "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addFaceData(faceStencilID_, faceDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencil3DID_, face3DDataHandling  , "VertexDoFToEdgeDoFOperatorFace3DStencil");
  storage->addCellData(cellStencilID_, cellDataHandling  , "VertexDoFToEdgeDoFOperatorCellStencil");
}

void
EdgeDoFOperator::apply_impl(EdgeDoFFunction<real_t> &src, EdgeDoFFunction<real_t> &dst, uint_t level, DoFType flag, UpdateType updateType) {

  this->startTiming( "EdgeDoFOperator - Apply" );

  src.startCommunication<Edge, Face>( level );
  src.endCommunication<Edge, Face>( level );
  src.startCommunication<Face, Cell>( level );
  src.endCommunication<Face, Cell>( level );
  src.communicate< Cell, Face >( level );
  src.startCommunication<Face, Edge>( level );

  for (auto& it : storage_->getCells())
  {
    Cell & cell = *it.second;

    const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if ( testFlag( cellBC, flag ) )
    {
      edgedof::macrocell::apply(level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType);
    }
  }



  for (auto& it : storage_->getFaces())
  {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if ( testFlag( faceBC, flag ) )
    {
      if ( storage_->hasGlobalCells() )
      {
        edgedof::macroface::apply3D( level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
      }
      else if( hhg::globalDefines::useGeneratedKernels )
      {
        real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
        real_t* src_data = face.getData( src.getFaceDataID() )->getPointer( level );
        real_t*       dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
        if( updateType == hhg::Replace )
        {
          edgedof::macroface::generated::applyReplace( dst_data, src_data, opr_data, level );
        } else if( updateType == hhg::Add )
        {
          edgedof::macroface::generated::applyAdd( dst_data, src_data, opr_data, level );
        }
      }
      else
      {
        edgedof::macroface::apply( level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
      }
    }
  }

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





  this->stopTiming( "EdgeDoFOperator - Apply" );
}


const PrimitiveDataID<StencilMemory<real_t>, Edge> &EdgeDoFOperator::getEdgeStencilID() const {
  return edgeStencilID_;
}

const PrimitiveDataID<StencilMemory<real_t>, Face> &EdgeDoFOperator::getFaceStencilID() const {
  return faceStencilID_;
}

const PrimitiveDataID<LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face> &EdgeDoFOperator::getFaceStencil3DID() const {
  return faceStencil3DID_;
}

const PrimitiveDataID<LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell> &EdgeDoFOperator::getCellStencilID() const {
  return cellStencilID_;
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

uint_t macroCellEdgeDoFToEdgeDoFStencilSize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 7 * 7 * 27;
}


}

