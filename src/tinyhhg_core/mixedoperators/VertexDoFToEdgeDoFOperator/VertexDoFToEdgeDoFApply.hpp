#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements3D.hpp"
#include "tinyhhg_core/Algorithms.hpp"

namespace hhg{
namespace VertexDoFToEdgeDoF{

/// map[neighborCellID][centerOrientation][indexOffset] = weight
typedef std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > > MacroEdgeStencilMap_T;
/// map[neighborCellID][centerOrientation][indexOffset] = weight
typedef std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > > MacroFaceStencilMap_T;

/// map[centerOrientation][indexOffset] = weight
typedef std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > MacroCellStencilMap_T;

inline void applyEdge(const uint_t & Level, Edge &edge,
                      const PrimitiveDataID<StencilMemory < real_t >, Edge> &operatorId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                      UpdateType update) {

  size_t rowsize = levelinfo::num_microedges_per_edge(Level);

  real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  real_t * src      = edge.getData(srcId)->getPointer( Level );
  real_t * dst      = edge.getData(dstId)->getPointer( Level );

  real_t tmp;

  for(uint_t i = 0; i < rowsize; ++i){
    tmp = 0.0;
    for(uint_t k = 0; k < hhg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF.size(); ++k){
      tmp += opr_data[hhg::vertexdof::stencilIndexFromHorizontalEdge(hhg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF[k])] *
             src[hhg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hhg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF[k] )];
    }
    for(uint_t k = 0; k < hhg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF.size(); ++k){
      tmp += opr_data[hhg::vertexdof::stencilIndexFromHorizontalEdge(hhg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF[k])] *
             src[hhg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hhg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF[k] )];
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < hhg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF.size(); ++k){
        tmp += opr_data[hhg::vertexdof::stencilIndexFromHorizontalEdge(hhg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF[k])] *
               src[hhg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hhg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF[k] )];
      }
    }

    if (update==Replace) {
      dst[hhg::edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] = tmp;
    } else if (update==Add) {
      dst[hhg::edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] += tmp;
    }
  }

}


inline void applyEdge3D( const uint_t & level, const Edge & edge,
                         const PrimitiveStorage & storage,
                         const PrimitiveDataID<LevelWiseMemory< MacroEdgeStencilMap_T >, Edge > &operatorId,
                         const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                         const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                         UpdateType update )
{
  auto opr_data = edge.getData(operatorId)->getData( level );
  real_t * src  = edge.getData(srcId)->getPointer( level );
  real_t * dst  = edge.getData(dstId)->getPointer( level );

  for ( const auto & centerIndexOnEdge : hhg::edgedof::macroedge::Iterator( level, 0 ) )
  {
    const edgedof::EdgeDoFOrientation edgeCenterOrientation = edgedof::EdgeDoFOrientation::X;

    real_t tmp = real_c( 0 );

    for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++  )
    {
      const Cell & neighborCell = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
      auto cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID() );

      const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >( { neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at(cellLocalEdgeID).at(0),
                                                                                  neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at(cellLocalEdgeID).at(1) } );

      const auto centerIndexInCell = indexing::basisConversion( centerIndexOnEdge, basisInCell, {0, 1, 2, 3}, levelinfo::num_microedges_per_edge( level ) );
      const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientation( edgeCenterOrientation, basisInCell.at(0), basisInCell.at(1), basisInCell.at(2) );

      for ( const auto & stencilIt : opr_data[neighborCellID][cellCenterOrientation] )
      {
        const auto stencilOffset = stencilIt.first;
        const auto stencilWeight = stencilIt.second;

        const auto leafIndexInCell = centerIndexInCell + stencilOffset;
        const auto leafIndexOnEdge = indexing::basisConversion( leafIndexInCell, {0, 1, 2, 3}, basisInCell, levelinfo::num_microvertices_per_edge( level ) );

        const auto onCellFacesSet = vertexdof::macrocell::isOnCellFace( leafIndexInCell, level );
        const auto onCellFacesSetOnEdge = vertexdof::macrocell::isOnCellFace( leafIndexOnEdge, level );

        WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() );

        uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

        if ( onCellFacesSet.size() == 0 )
        {
          // leaf in macro-cell
          leafArrayIndexOnEdge = vertexdof::macroedge::indexOnNeighborCell( level, leafIndexOnEdge.x(), neighborCellID, edge.getNumNeighborFaces() );
        }
        else if ( onCellFacesSet.size() == 1 )
        {
          // leaf on macro-face
          const auto faceID = neighborCell.neighborFaces().at( *onCellFacesSet.begin() );
          WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), faceID ) != edge.neighborFaces().end() )
          const auto localFaceIDOnEdge = edge.face_index( faceID );
          leafArrayIndexOnEdge = vertexdof::macroedge::indexOnNeighborFace( level, leafIndexOnEdge.x(), localFaceIDOnEdge );
        }
        else
        {
          // leaf on macro-edge
          WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), 2 );
          leafArrayIndexOnEdge = vertexdof::macroedge::index( level, leafIndexOnEdge.x() );
        }

        tmp += src[ leafArrayIndexOnEdge ] * stencilWeight;
      }

    }

    if ( update == Replace )
    {
      dst[ edgedof::macroedge::index( level, centerIndexOnEdge.x() ) ] = tmp;
    }
    else if ( update == Add )
    {
      dst[ edgedof::macroedge::index( level, centerIndexOnEdge.x() ) ] += tmp;
    }
  }
}



inline void applyFace(const uint_t & Level, Face &face,
                      const PrimitiveDataID<StencilMemory < real_t >, Face> &operatorId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Face> &srcId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstId,
                      UpdateType update){

  real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  real_t * src      = face.getData(srcId)->getPointer( Level );
  real_t * dst      = face.getData(dstId)->getPointer( Level );

  real_t tmp;

  using namespace vertexdof::macroface;

  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromHorizontalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromHorizontalEdge(neighborsFromHorizontalEdge[k])] *
               src[indexFromHorizontalEdge( Level, it.col(), it.row(), neighborsFromHorizontalEdge[k] )];
      }
      if (update==Replace) {
        dst[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] += tmp;
      }
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromDiagonalEdge(neighborsFromDiagonalEdge[k])] *
               src[indexFromDiagonalEdge( Level, it.col(), it.row(), neighborsFromDiagonalEdge[k] )];
      }
      if (update==Replace) {
        dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] += tmp;
      }
    }
    if( it.col() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromVerticalEdge(neighborsFromVerticalEdge[k])] *
               src[indexFromVerticalEdge( Level, it.col(), it.row(), neighborsFromVerticalEdge[k] )];
      }

      if (update==Replace) {
        dst[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] += tmp;
      }
    }
  }
}


inline void applyFace3D( const uint_t & level, Face &face,
                         const PrimitiveStorage & storage,
                         const PrimitiveDataID< LevelWiseMemory< MacroFaceStencilMap_T >, Face > &operatorId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Face > &srcId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Face > &dstId,
                         UpdateType update)
{
  auto opr_data = face.getData(operatorId)->getData( level );
  real_t * src  = face.getData(srcId)->getPointer( level );
  real_t * dst  = face.getData(dstId)->getPointer( level );

  for ( const auto & centerIndexInFace : hhg::edgedof::macroface::Iterator( level, 0 ) )
  {
    std::map< edgedof::EdgeDoFOrientation, real_t > tmpResults = {
    { edgedof::EdgeDoFOrientation::X, real_c(0) },
    { edgedof::EdgeDoFOrientation::Y, real_c(0) },
    { edgedof::EdgeDoFOrientation::XY, real_c(0) },
    };

    for ( const auto & faceCenterOrientation : edgedof::faceLocalEdgeDoFOrientations )
    {
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::X && edgedof::isHorizontalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::Y && edgedof::isVerticalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::XY && edgedof::isDiagonalEdgeOnBoundary( level, centerIndexInFace )  )
        continue;

      for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++  )
      {
        const Cell & neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
        const uint_t localFaceID = neighborCell.getLocalFaceID( face.getID() );

        const std::array< uint_t, 4 > localVertexIDsAtCell = {
        neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(0),
        neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(1),
        neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(2),
        6 - neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(0)
        - neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(1)
        - neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(2)
        };

        const auto centerIndexInCell = indexing::basisConversion( centerIndexInFace, localVertexIDsAtCell, {0, 1, 2, 3}, levelinfo::num_microedges_per_edge( level ) );
        const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientation( faceCenterOrientation, localVertexIDsAtCell.at(0), localVertexIDsAtCell.at(1), localVertexIDsAtCell.at(2) );

        for ( const auto & stencilIt : opr_data[neighborCellID][cellCenterOrientation] )
        {
          const auto stencilOffset = stencilIt.first;
          const auto stencilWeight = stencilIt.second;

          const auto leafIndexInCell = centerIndexInCell + stencilOffset;
          const auto leafIndexInFace = indexing::basisConversion( leafIndexInCell, {0, 1, 2, 3}, localVertexIDsAtCell, levelinfo::num_microvertices_per_edge( level ) );
          WALBERLA_ASSERT_LESS_EQUAL( leafIndexInFace.z(), 1 );
          uint_t leafArrayIndexInFace;
          if ( leafIndexInFace.z() == 0 )
          {
            leafArrayIndexInFace = vertexdof::macroface::index( level, leafIndexInFace.x(), leafIndexInFace.y() );
          }
          else
          {
            leafArrayIndexInFace = vertexdof::macroface::index( level, leafIndexInFace.x(), leafIndexInFace.y(), neighborCellID );
          }

          tmpResults[faceCenterOrientation] += stencilWeight * src[leafArrayIndexInFace];
        }
      }
      const auto dstIdx = edgedof::macroface::index( level, centerIndexInFace.x(), centerIndexInFace.y(), faceCenterOrientation );
      if ( update == Replace )
      {
        dst[dstIdx] = tmpResults[faceCenterOrientation];
      }
      else
      {
        dst[dstIdx] += tmpResults[faceCenterOrientation];
      }
    }
  }
}


inline void applyCell(const uint_t & Level, Cell & cell,
                      const PrimitiveDataID<LevelWiseMemory< MacroCellStencilMap_T >, Cell> & operatorId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Cell> &srcId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Cell> &dstId,
                      UpdateType update){

  auto opr_data = cell.getData(operatorId)->getData( Level );
  real_t * src  = cell.getData(srcId)->getPointer( Level );
  real_t * dst  = cell.getData(dstId)->getPointer( Level );

  for ( const auto & it : hhg::edgedof::macrocell::Iterator( Level, 0 ) )
  {
    std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

    if ( edgedof::macrocell::isInnerXEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
    if ( edgedof::macrocell::isInnerYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
    if ( edgedof::macrocell::isInnerZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
    if ( edgedof::macrocell::isInnerXYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
    if ( edgedof::macrocell::isInnerXZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
    if ( edgedof::macrocell::isInnerYZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

    for ( const auto & centerOrientation : innerOrientations )
    {
      real_t tmp = 0.0;

      const auto vertexDoFNeighbors = P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );
      for ( const auto & neighbor : vertexDoFNeighbors )
      {
        const auto   srcIdx      = it + neighbor;
        const auto   srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
        tmp += opr_data[centerOrientation][neighbor] * src[srcArrayIdx];
      }

      const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );

      if ( update == Replace )
      {
        dst[dstArrayIdx] = tmp;
      }
      else if ( update == Add )
      {
        dst[dstArrayIdx] += tmp;
      }
    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    real_t tmp = 0.0;
    const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;

    const auto vertexDoFNeighbors = P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );
    for ( const auto & neighbor : vertexDoFNeighbors )
    {
      const auto   srcIdx      = it + neighbor;
      const auto   srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
      tmp += opr_data[centerOrientation][neighbor] * src[srcArrayIdx];
    }

    const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );

    if ( update == Replace )
    {
      dst[dstArrayIdx] = tmp;
    }
    else if ( update == Add )
    {
      dst[dstArrayIdx] += tmp;
    }
  }
}


} // VertexDoFToEdgeDoFOperator
} // namespace hhg
