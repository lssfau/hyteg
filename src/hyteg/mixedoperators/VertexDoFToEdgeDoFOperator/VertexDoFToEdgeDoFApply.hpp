/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p2functionspace/P2Elements3D.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
namespace VertexDoFToEdgeDoF{

/// map[neighborCellID][centerOrientation][indexOffset] = weight
typedef std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, std::map< indexing::Index, real_t > > > MacroEdgeStencilMap_T;
/// map[neighborCellID][centerOrientation][indexOffset] = weight
typedef std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, std::map< indexing::Index, real_t > > > MacroFaceStencilMap_T;

/// map[centerOrientation][indexOffset] = weight
typedef std::map< edgedof::EdgeDoFOrientation, std::map< indexing::Index, real_t > > MacroCellStencilMap_T;

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
    for(uint_t k = 0; k < hyteg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF.size(); ++k){
      tmp += opr_data[hyteg::vertexdof::stencilIndexFromHorizontalEdge(
                  hyteg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF[k])] *
             src[hyteg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hyteg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF[k] )];
    }
    for(uint_t k = 0; k < hyteg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF.size(); ++k){
      tmp += opr_data[hyteg::vertexdof::stencilIndexFromHorizontalEdge(
                  hyteg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF[k])] *
             src[hyteg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hyteg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF[k] )];
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < hyteg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF.size(); ++k){
        tmp += opr_data[hyteg::vertexdof::stencilIndexFromHorizontalEdge(
                    hyteg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF[k])] *
               src[hyteg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hyteg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF[k] )];
      }
    }

    if (update==Replace) {
      dst[hyteg::edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] = tmp;
    } else if (update==Add) {
      dst[hyteg::edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] += tmp;
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

  for ( const auto & centerIndexOnEdge : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
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
      const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(edgeCenterOrientation, basisInCell.at(0),
                                                                                      basisInCell.at(1), basisInCell.at(2));

      for ( const auto & stencilIt : opr_data[neighborCellID][cellCenterOrientation] )
      {
        const auto stencilOffset = stencilIt.first;
        const auto stencilWeight = stencilIt.second;

        const auto leafIndexInCell = centerIndexInCell + stencilOffset;
        const auto leafIndexOnEdge = indexing::basisConversion(
            leafIndexInCell, { 0, 1, 2, 3 }, basisInCell, levelinfo::num_microvertices_per_edge( level ) );

        const auto onCellFacesSet = vertexdof::macrocell::isOnCellFace( leafIndexInCell, level );
        const auto onCellFacesSetOnEdge = vertexdof::macrocell::isOnCellFace( leafIndexOnEdge, level );

        WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() );

        const auto cellLocalIDsOfNeighborFaces = indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID );
        std::vector< uint_t > cellLocalIDsOfNeighborFacesWithLeafOnThem;
        std::set_intersection( cellLocalIDsOfNeighborFaces.begin(), cellLocalIDsOfNeighborFaces.end(),
                               onCellFacesSet.begin(), onCellFacesSet.end(), std::back_inserter( cellLocalIDsOfNeighborFacesWithLeafOnThem ) );

        uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

        if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 0 )
        {
          // leaf in macro-cell
          leafArrayIndexOnEdge = vertexdof::macroedge::indexOnNeighborCell( level, leafIndexOnEdge.x(), neighborCellID, edge.getNumNeighborFaces() );
        }
        else if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 1 )
        {
          // leaf on macro-face

          const auto faceID = neighborCell.neighborFaces().at( *cellLocalIDsOfNeighborFacesWithLeafOnThem.begin() );
          WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), faceID ) != edge.neighborFaces().end() );

          const auto localFaceIDOnEdge = edge.face_index( faceID );
          leafArrayIndexOnEdge = vertexdof::macroedge::indexOnNeighborFace( level, leafIndexOnEdge.x(), localFaceIDOnEdge );

        }
        else
        {
          // leaf on macro-edge
          WALBERLA_ASSERT_EQUAL( cellLocalIDsOfNeighborFacesWithLeafOnThem.size(), 2 );
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

  for ( const auto & it : hyteg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.y() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromHorizontalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromHorizontalEdge(neighborsFromHorizontalEdge[k])] *
               src[indexFromHorizontalEdge( Level, it.x(), it.y(), neighborsFromHorizontalEdge[k] )];
      }
      if (update==Replace) {
        dst[edgedof::macroface::indexFromHorizontalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_HO_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromHorizontalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_HO_C )] += tmp;
      }
    }
    if( it.x() + it.y() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromDiagonalEdge(neighborsFromDiagonalEdge[k])] *
               src[indexFromDiagonalEdge( Level, it.x(), it.y(), neighborsFromDiagonalEdge[k] )];
      }
      if (update==Replace) {
        dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_DI_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_DI_C )] += tmp;
      }
    }
    if( it.x() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromVerticalEdge(neighborsFromVerticalEdge[k])] *
               src[indexFromVerticalEdge( Level, it.x(), it.y(), neighborsFromVerticalEdge[k] )];
      }

      if (update==Replace) {
        dst[edgedof::macroface::indexFromVerticalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_VE_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromVerticalEdge( Level, it.x(), it.y(), stencilDirection::EDGE_VE_C )] += tmp;
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

  for ( const auto & centerIndexInFace : hyteg::edgedof::macroface::Iterator( level, 0 ) )
  {
    std::map< edgedof::EdgeDoFOrientation, real_t > tmpResults = {
    { edgedof::EdgeDoFOrientation::X, real_c(0) },
    { edgedof::EdgeDoFOrientation::Y, real_c(0) },
    { edgedof::EdgeDoFOrientation::XY, real_c(0) },
    };

    for ( const auto & faceCenterOrientation : edgedof::faceLocalEdgeDoFOrientations )
    {
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::X && edgedof::macroface::isHorizontalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::Y && edgedof::macroface::isVerticalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::XY && edgedof::macroface::isDiagonalEdgeOnBoundary( level, centerIndexInFace )  )
        continue;

      for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++  )
      {
        const Cell & neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
        const uint_t localFaceID = neighborCell.getLocalFaceID( face.getID() );

        const auto centerIndexInCell = edgedof::macroface::getIndexInNeighboringMacroCell( centerIndexInFace, face, neighborCellID, storage, level );
        const auto cellCenterOrientation = edgedof::macroface::getOrientattionInNeighboringMacroCell( faceCenterOrientation, face, neighborCellID, storage );

        for ( const auto & stencilIt : opr_data[neighborCellID][cellCenterOrientation] )
        {
          const auto stencilOffset = stencilIt.first;
          const auto stencilWeight = stencilIt.second;

          const auto leafIndexInCell = centerIndexInCell + stencilOffset;
          const auto leafIndexInFace =
              vertexdof::macrocell::getIndexInNeighboringMacroFace( leafIndexInCell, neighborCell, localFaceID, storage, level );

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

  for ( const auto & it : hyteg::edgedof::macrocell::Iterator( Level, 0 ) )
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
        const auto srcIdx      = it + neighbor;
        const auto srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
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
      const auto srcIdx      = it + neighbor;
      const auto srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
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
} // namespace hyteg
