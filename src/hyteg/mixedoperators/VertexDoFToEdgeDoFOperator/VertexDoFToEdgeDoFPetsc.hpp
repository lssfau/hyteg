
#pragma once

#include "core/DataTypes.h"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFApply.hpp"

namespace hyteg {
namespace VertexDoFToEdgeDoF {

using walberla::real_t;
using walberla::uint_t;

#ifdef HYTEG_BUILD_WITH_PETSC


inline void saveEdgeOperator( const uint_t & Level, const Edge & edge,
                              const PrimitiveDataID< StencilMemory< real_t >, Edge>    & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & dstId,
                              Mat & mat )
{
  const real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  const PetscInt * src      = edge.getData(srcId)->getPointer( Level );
  const PetscInt * dst      = edge.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetscInt dstInt;

  for( const auto it : edgedof::macroedge::Iterator( Level, 0 ) )
  {
    dstInt = dst[edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C )];

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF )
    {
      srcInt = src[vertexdof::macroedge::indexFromHorizontalEdge( Level, it.col(), neighbor )];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], ADD_VALUES );
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF )
    {
      srcInt = src[vertexdof::macroedge::indexFromHorizontalEdge( Level, it.col(), neighbor )];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], ADD_VALUES );
    }

    if( edge.getNumNeighborFaces() == 2 )
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF )
      {
        srcInt = src[vertexdof::macroedge::indexFromHorizontalEdge( Level, it.col(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], ADD_VALUES );
      }
    }
  }
}


inline void saveEdgeOperator3D( const uint_t & level, const Edge & edge,
                                const PrimitiveStorage                                                 & storage,
                                const PrimitiveDataID<LevelWiseMemory< MacroEdgeStencilMap_T >, Edge > & operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & dstId,
                                Mat & mat )
{
  auto opr_data = edge.getData(operatorId)->getData( level );
  auto src  = edge.getData(srcId)->getPointer( level );
  auto dst  = edge.getData(dstId)->getPointer( level );

  for ( const auto & centerIndexOnEdge : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
  {
    const edgedof::EdgeDoFOrientation edgeCenterOrientation = edgedof::EdgeDoFOrientation::X;

    const auto dstInt = dst[ edgedof::macroedge::index( level, centerIndexOnEdge.x() ) ];

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
        const auto leafIndexOnEdge = indexing::basisConversion( leafIndexInCell, {0, 1, 2, 3}, basisInCell, levelinfo::num_microvertices_per_edge( level ) );

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

        const auto srcInt = src[ leafArrayIndexOnEdge ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &stencilWeight, ADD_VALUES );
      }
    }
  }
}


inline void saveFaceOperator( const uint_t & Level, const Face & face,
                              const PrimitiveDataID< StencilMemory< real_t >, Face>    & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & dstId,
                              Mat & mat )
{
  const real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  const PetscInt * src      = face.getData(srcId)->getPointer( Level );
  const PetscInt * dst      = face.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetscInt dstInt;

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0 )
    {
      dstInt = dst[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromHorizontalEdge )
      {
        srcInt = src[vertexdof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], ADD_VALUES );
      }
    }

    if( it.col() != 0 )
    {
      dstInt = dst[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromVerticalEdge )
      {
        srcInt = src[vertexdof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromVerticalEdge( neighbor ) ], ADD_VALUES );
      }
    }

    if( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1) )
    {
      dstInt = dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromDiagonalEdge )
      {
        srcInt = src[vertexdof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromDiagonalEdge( neighbor ) ], ADD_VALUES );
      }
    }

  }
}


inline void saveFaceOperator3D( const uint_t & level, const Face & face,
                                const PrimitiveStorage                                                  & storage,
                                const PrimitiveDataID< LevelWiseMemory< MacroFaceStencilMap_T >, Face > & operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & dstId,
                                Mat & mat )
{
  auto opr_data = face.getData(operatorId)->getData( level );
  auto src  = face.getData(srcId)->getPointer( level );
  auto dst  = face.getData(dstId)->getPointer( level );

  for ( const auto & centerIndexInFace : hyteg::edgedof::macroface::Iterator( level, 0 ) )
  {
     for ( const auto & faceCenterOrientation : edgedof::faceLocalEdgeDoFOrientations )
    {
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::X && edgedof::macroface::isHorizontalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::Y && edgedof::macroface::isVerticalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::XY && edgedof::macroface::isDiagonalEdgeOnBoundary( level, centerIndexInFace )  )
        continue;

      const auto dstIdx = edgedof::macroface::index( level, centerIndexInFace.x(), centerIndexInFace.y(), faceCenterOrientation );
      const auto dstInt = dst[ dstIdx ];

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
          const auto leafIndexInFace = vertexdof::macrocell::getIndexInNeighboringMacroFace( leafIndexInCell, neighborCell, localFaceID, storage, level );

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

          const auto srcInt = src[ leafArrayIndexInFace ];
          MatSetValues( mat, 1, &dstInt, 1, &srcInt, &stencilWeight, ADD_VALUES );
        }
      }
    }
  }
}


inline void saveCellOperator( const uint_t & Level, const Cell & cell,
                              const PrimitiveDataID<LevelWiseMemory< MacroCellStencilMap_T >, Cell> & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Cell> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Cell> & dstId,
                              Mat & mat )
{
  auto opr_data = cell.getData( operatorId )->getData( Level );
  const PetscInt *src = cell.getData( srcId )->getPointer( Level );
  const PetscInt *dst = cell.getData( dstId )->getPointer( Level );

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
      const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
      const auto dstInt      = dst[ dstArrayIdx ];

      const auto vertexDoFNeighbors = P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );

      for ( const auto & neighbor : vertexDoFNeighbors )
      {
        const auto   srcIdx      = it + neighbor;
        const auto   srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
        const auto   srcInt      = src[ srcArrayIdx ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[centerOrientation][neighbor], ADD_VALUES );
      }
    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;

    const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
    const auto dstInt      = dst[ dstArrayIdx ];

    const auto vertexDoFNeighbors = P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );

    for ( const auto & neighbor : vertexDoFNeighbors )
    {
      const auto   srcIdx      = it + neighbor;
      const auto   srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
      const auto   srcInt      = src[ srcArrayIdx ];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[centerOrientation][neighbor], ADD_VALUES );
    }
  }
}

template < class OperatorType >
inline void createMatrix( const OperatorType&                opr,
                          const P1Function< PetscInt >&      src,
                          const EdgeDoFFunction< PetscInt >& dst,
                          Mat&                               mat,
                          size_t                             level,
                          DoFType                            flag )
{
  const auto storage = src.getStorage();

  for (auto& it : opr.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag))
    {
      if ( storage->hasGlobalCells() )
      {
        saveEdgeOperator3D(level, edge, *storage, opr.getEdgeStencil3DID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
      }
      else
      {
        saveEdgeOperator(level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
      }
    }
  }

  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      if ( storage->hasGlobalCells() )
      {
        saveFaceOperator3D( level, face, *storage, opr.getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
      }
      else
      {
        saveFaceOperator( level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
      }
    }
  }

  for (auto& it : opr.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag))
    {
      saveCellOperator(level, cell, opr.getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat);
    }
  }
}

#endif

}
}
