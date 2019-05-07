#include "EdgeDoFOperator.hpp"

#include "tinyhhg_core/HHGDefinitions.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/edgedofspace/generatedKernels/all.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroCell.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroFace.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroEdge.hpp"

namespace hhg {

EdgeDoFOperator::EdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                 size_t minLevel,
                                 size_t maxLevel)
  : Operator(storage, minLevel, maxLevel)
{
  auto edgeDataHandling   =
      std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge   >>(minLevel_, maxLevel_, macroEdgeEdgeDoFToEdgeDoFStencilSize);

  auto edge3DDataHandling   =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >>( minLevel_, maxLevel_ );
  
  auto faceDataHandling   =
      std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face   >>(minLevel_, maxLevel_, macroFaceEdgeDoFToEdgeDoFStencilSize);

  auto face3DDataHandling   =
      std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >>( minLevel_, maxLevel_ );

  auto cellDataHandling   =
      std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >>( minLevel_, maxLevel_ );

  storage->addEdgeData(edgeStencilID_, edgeDataHandling  , "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addEdgeData(edgeStencil3DID_, edge3DDataHandling  , "VertexDoFToEdgeDoFOperatorEdge3DStencil");
  storage->addFaceData(faceStencilID_, faceDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencil3DID_, face3DDataHandling  , "VertexDoFToEdgeDoFOperatorFace3DStencil");
  storage->addCellData(cellStencilID_, cellDataHandling  , "VertexDoFToEdgeDoFOperatorCellStencil");
}

void EdgeDoFOperator::apply(const EdgeDoFFunction<real_t> &src,const EdgeDoFFunction<real_t> &dst, uint_t level, DoFType flag, UpdateType updateType) const {

  this->startTiming( "Apply" );

  src.startCommunication<Edge, Face>( level );
  src.endCommunication<Edge, Face>( level );
  src.startCommunication<Face, Cell>( level );
  src.endCommunication<Face, Cell>( level );
  src.communicate< Cell, Face >( level );
  src.startCommunication<Face, Edge>( level );

  this->timingTree_->start( "Macro-Cell" );

  for (auto& it : storage_->getCells())
  {
    Cell & cell = *it.second;

    const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if ( testFlag( cellBC, flag ) )
    {
      if( hhg::globalDefines::useGeneratedKernels && updateType == Replace )
      {
        typedef edgedof::EdgeDoFOrientation eo;
        auto dstData = cell.getData( dst.getCellDataID() )->getPointer( level );
        auto srcData = cell.getData( src.getCellDataID() )->getPointer( level );
        auto stencilData = cell.getData( cellStencilID_ )->getData( level );
        std::map< eo, uint_t > firstIdx;
        for ( auto e : edgedof::allEdgeDoFOrientations )
            firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );
        edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_replace ( &dstData[firstIdx[eo::X]],
                                                                                   &dstData[firstIdx[eo::XY]],
                                                                                   &dstData[firstIdx[eo::XYZ]],
                                                                                   &dstData[firstIdx[eo::XZ]],
                                                                                   &dstData[firstIdx[eo::Y]],
                                                                                   &dstData[firstIdx[eo::YZ]],
                                                                                   &dstData[firstIdx[eo::Z]],
                                                                                   &srcData[firstIdx[eo::X]],
                                                                                   &srcData[firstIdx[eo::XY]],
                                                                                   &srcData[firstIdx[eo::XYZ]],
                                                                                   &srcData[firstIdx[eo::XZ]],
                                                                                   &srcData[firstIdx[eo::Y]],
                                                                                   &srcData[firstIdx[eo::YZ]],
                                                                                   &srcData[firstIdx[eo::Z]],
                                                                                   stencilData,
                                                                                   static_cast< int64_t >( level ) );
          
      }
      else
      {
        edgedof::macrocell::apply(level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType);
      }
    }
  }

  this->timingTree_->stop( "Macro-Cell" );

  this->timingTree_->start( "Macro-Face" );

  for ( auto& it : storage_->getFaces() )
  {
     Face& face = *it.second;

     const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
     if ( testFlag( faceBC, flag ) )
     {
        if ( storage_->hasGlobalCells() )
        {
           if ( hhg::globalDefines::useGeneratedKernels && face.getNumNeighborCells() == 2 )
           {
             auto opr_data = face.getData( faceStencil3DID_ )->getData( level );
             auto src_data = face.getData( src.getFaceDataID() )->getPointer( level );
             auto dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );

             const uint_t offset_x = edgedof::macroface::index( level, 0, 0, edgedof::EdgeDoFOrientation::X );
             const uint_t offset_xy = edgedof::macroface::index( level, 0, 0, edgedof::EdgeDoFOrientation::XY );
             const uint_t offset_y = edgedof::macroface::index( level, 0, 0, edgedof::EdgeDoFOrientation::Y );

             std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, uint_t > > offset_gl_orientation;
             for ( uint_t gl = 0; gl < 2; gl++ )
             {
               for ( const auto & eo : edgedof::allEdgeDoFOrientations )
               {
                 offset_gl_orientation[gl][eo] = edgedof::macroface::index( level, 0, 0, eo, gl );
               }
             }

             auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );
             auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

             auto neighbor_cell_0_local_vertex_id_0 = static_cast< int64_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at(0) );
             auto neighbor_cell_0_local_vertex_id_1 = static_cast< int64_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at(1) );
             auto neighbor_cell_0_local_vertex_id_2 = static_cast< int64_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at(2) );

             auto neighbor_cell_1_local_vertex_id_0 = static_cast< int64_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at(0) );
             auto neighbor_cell_1_local_vertex_id_1 = static_cast< int64_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at(1) );
             auto neighbor_cell_1_local_vertex_id_2 = static_cast< int64_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at(2) );

             if ( updateType == Replace )
             {
               edgedof::macroface::generated::apply_3D_macroface_one_sided_edgedof_to_edgedof_replace(
               &dst_data[offset_x],
               &dst_data[offset_xy],
               &dst_data[offset_y],
               &src_data[offset_x],
               &src_data[offset_xy],
               &src_data[offset_y],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
               opr_data[0],
               static_cast< int64_t >( level ),
               neighbor_cell_0_local_vertex_id_0,
               neighbor_cell_0_local_vertex_id_1,
               neighbor_cell_0_local_vertex_id_2
               );
             }
             else
             {
               edgedof::macroface::generated::apply_3D_macroface_one_sided_edgedof_to_edgedof_add(
               &dst_data[offset_x],
               &dst_data[offset_xy],
               &dst_data[offset_y],
               &src_data[offset_x],
               &src_data[offset_xy],
               &src_data[offset_y],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
               &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
               opr_data[0],
               static_cast< int64_t >( level ),
               neighbor_cell_0_local_vertex_id_0,
               neighbor_cell_0_local_vertex_id_1,
               neighbor_cell_0_local_vertex_id_2
               );
             }

             edgedof::macroface::generated::apply_3D_macroface_one_sided_edgedof_to_edgedof_add(
             &dst_data[offset_x],
             &dst_data[offset_xy],
             &dst_data[offset_y],
             &src_data[offset_x],
             &src_data[offset_xy],
             &src_data[offset_y],
             &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
             &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
             &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
             &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
             &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
             &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
             &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
             opr_data[1],
             static_cast< int64_t >( level ),
             neighbor_cell_1_local_vertex_id_0,
             neighbor_cell_1_local_vertex_id_1,
             neighbor_cell_1_local_vertex_id_2
             );
           }
           else
           {
              edgedof::macroface::apply3D(
                  level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
           }
        }
        else
        {
           if ( hhg::globalDefines::useGeneratedKernels )
           {
              typedef edgedof::EdgeDoFOrientation eo;
              real_t*                             opr_data = face.getData( faceStencilID_ )->getPointer( level );
              real_t*                             src_data = face.getData( src.getFaceDataID() )->getPointer( level );
              real_t*                             dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
              std::map< eo, uint_t >              firstIdx;
              for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
                 firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

              if ( updateType == hhg::Replace )
              {
                 edgedof::macroface::generated::apply_2D_macroface_edgedof_to_edgedof_replace( &dst_data[firstIdx[eo::X]],
                                                                                               &dst_data[firstIdx[eo::XY]],
                                                                                               &dst_data[firstIdx[eo::Y]],
                                                                                               &src_data[firstIdx[eo::X]],
                                                                                               &src_data[firstIdx[eo::XY]],
                                                                                               &src_data[firstIdx[eo::Y]],
                                                                                               &opr_data[5],
                                                                                               &opr_data[0],
                                                                                               &opr_data[10],
                                                                                               static_cast< int64_t >( level ) );
              }
              else if ( updateType == hhg::Add )
              {
                 edgedof::macroface::generated::apply_2D_macroface_edgedof_to_edgedof_add( &dst_data[firstIdx[eo::X]],
                                                                                           &dst_data[firstIdx[eo::XY]],
                                                                                           &dst_data[firstIdx[eo::Y]],
                                                                                           &src_data[firstIdx[eo::X]],
                                                                                           &src_data[firstIdx[eo::XY]],
                                                                                           &src_data[firstIdx[eo::Y]],
                                                                                           &opr_data[5],
                                                                                           &opr_data[0],
                                                                                           &opr_data[10],
                                                                                           static_cast< int64_t >( level ) );
              }
           }
           else
           {
              edgedof::macroface::apply( level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
           }
        }
     }
  }

  this->timingTree_->stop( "Macro-Face" );

  src.endCommunication<Face, Edge>( level );

  this->timingTree_->start( "Macro-Edge" );

  for (auto& it : storage_->getEdges())
  {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if ( testFlag( edgeBC, flag ) )
    {
      if ( storage_->hasGlobalCells() )
      {
        edgedof::macroedge::apply3D(level, edge, *storage_, edgeStencil3DID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
      else
      {
        edgedof::macroedge::apply(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }

    }
  }

  this->timingTree_->stop( "Macro-Edge" );

  this->stopTiming( "Apply" );
}


const PrimitiveDataID<StencilMemory<real_t>, Edge> &EdgeDoFOperator::getEdgeStencilID() const {
  return edgeStencilID_;
}

const PrimitiveDataID<StencilMemory<real_t>, Face> &EdgeDoFOperator::getFaceStencilID() const {
  return faceStencilID_;
}

const PrimitiveDataID<LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge> &EdgeDoFOperator::getEdgeStencil3DID() const {
  return edgeStencil3DID_;
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

