
#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroFace.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

namespace hhg {
namespace edgedof {

using walberla::real_t;
using walberla::uint_t;

#ifdef HHG_BUILD_WITH_PETSC

inline void createVectorFromFunction(EdgeDoFFunction<PetscScalar> &function,
                                     EdgeDoFFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      macroedge::createVectorFromFunction<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
      macroface::createVectorFromFunction<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag)) {
      macrocell::createVectorFromFunction<PetscScalar>(level, cell, function.getCellDataID(), numerator.getCellDataID(), vec);
    }
  }
}

inline void createFunctionFromVector(EdgeDoFFunction<PetscScalar> &function,
                                     EdgeDoFFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  function.startCommunication<Vertex, Edge>( level );
  function.endCommunication<Vertex, Edge>( level );

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      edgedof::macroedge::createFunctionFromVector<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  function.startCommunication<Edge, Face>( level );
  function.endCommunication<Edge, Face>( level );

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
      edgedof::macroface::createFunctionFromVector<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag)) {
      edgedof::macrocell::createFunctionFromVector<PetscScalar>(level, cell, function.getCellDataID(), numerator.getCellDataID(), vec);
    }
  }
}

inline void applyDirichletBC(EdgeDoFFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {

  for (auto &it : numerator.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = numerator.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, DirichletBoundary)) {
      edgedof::macroedge::applyDirichletBC(level, edge, mat, numerator.getEdgeDataID());
    }
  }

  for (auto &it : numerator.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = numerator.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, DirichletBoundary)) {
      edgedof::macroface::applyDirichletBC(level, face, mat, numerator.getFaceDataID());
    }
  }

}



inline void saveEdgeOperator( const uint_t & Level, const Edge & edge,
                              const PrimitiveDataID< StencilMemory< real_t >, Edge>    & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & dstId,
                              Mat & mat )
{
  using namespace hhg::edgedof::macroedge;
  size_t rowsize = levelinfo::num_microedges_per_edge(Level);

  real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  PetscInt * src      = edge.getData(srcId)->getPointer( Level );
  PetscInt * dst      = edge.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetscInt dstInt;

  for(uint_t i = 0; i < rowsize; ++i){
    dstInt = dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )];

    for(uint_t k = 0; k < neighborsOnEdgeFromHorizontalEdge.size(); ++k){
      srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnEdgeFromHorizontalEdge[k] )];
      MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnEdgeFromHorizontalEdge[k])], INSERT_VALUES);
    }
    for(uint_t k = 0; k < neighborsOnSouthFaceFromHorizontalEdge.size(); ++k){
      srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnSouthFaceFromHorizontalEdge[k] )];
      MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnSouthFaceFromHorizontalEdge[k])], INSERT_VALUES);
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < neighborsOnNorthFaceFromHorizontalEdge.size(); ++k){
        srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnNorthFaceFromHorizontalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnNorthFaceFromHorizontalEdge[k])], INSERT_VALUES);
      }
    }
  }
}


inline void saveEdgeOperator3D( const uint_t & level, const Edge & edge,
                                const PrimitiveStorage                                                   & storage,
                                const PrimitiveDataID<LevelWiseMemory< macroedge::StencilMap_T >, Edge > & operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & dstId,
                                Mat & mat )
{
  auto opr_data = edge.getData(operatorId)->getData( level );
  auto src  = edge.getData(srcId)->getPointer( level );
  auto dst  = edge.getData(dstId)->getPointer( level );

  for ( const auto & centerIndexOnEdge : hhg::edgedof::macroedge::Iterator( level, 0 ) )
  {
    const EdgeDoFOrientation edgeCenterOrientation = EdgeDoFOrientation::X;

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

      for ( const auto & leafOrientationInCell : edgedof::allEdgeDoFOrientations )
      {
        for ( const auto & stencilIt : opr_data[neighborCellID][cellCenterOrientation][leafOrientationInCell] )
        {
          const auto stencilOffset = stencilIt.first;
          const auto stencilWeight = stencilIt.second;

          const auto leafOrientationOnEdge = edgedof::convertEdgeDoFOrientationCellToFace( leafOrientationInCell, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ));
          const auto leafIndexInCell = centerIndexInCell + stencilOffset;

          const auto leafIndexOnEdge = indexing::basisConversion( leafIndexInCell, {0, 1, 2, 3}, basisInCell, levelinfo::num_microedges_per_edge( level ) );

          const auto onCellFacesSet = edgedof::macrocell::isOnCellFaces( level, leafIndexInCell, leafOrientationInCell );
          const auto onCellFacesSetOnEdge = edgedof::macrocell::isOnCellFaces( level, leafIndexOnEdge, leafOrientationOnEdge );

          WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() );

          uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

          const auto cellLocalIDsOfNeighborFaces = indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID );
          std::vector< uint_t > cellLocalIDsOfNeighborFacesWithLeafOnThem;
          std::set_intersection( cellLocalIDsOfNeighborFaces.begin(), cellLocalIDsOfNeighborFaces.end(),
                                 onCellFacesSet.begin(), onCellFacesSet.end(), std::back_inserter( cellLocalIDsOfNeighborFacesWithLeafOnThem ) );

          if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 0 )
          {
            // leaf in macro-cell
            leafArrayIndexOnEdge = edgedof::macroedge::indexOnNeighborCell( level, leafIndexOnEdge.x(), neighborCellID, edge.getNumNeighborFaces(), leafOrientationOnEdge );
          }
          else if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 1 )
          {
            // leaf on macro-face
            WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) );

            const auto cellLocalFaceID = *cellLocalIDsOfNeighborFacesWithLeafOnThem.begin();
            const auto facePrimitiveID = neighborCell.neighborFaces().at( cellLocalFaceID );
            WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), facePrimitiveID ) != edge.neighborFaces().end() );


            // The leaf orientation on the edge must be X, Y or XY since it is located on a neighboring face.
            // Therefore we need to know the three spanning vertex IDs and convert the leaf orientation again.
            const auto spanningCellLocalVertices = indexing::cellLocalFaceIDsToSpanningVertexIDs.at( cellLocalFaceID );
            std::array< uint_t, 4 > faceBasisInCell;
            if ( spanningCellLocalVertices.count( basisInCell.at( 2 ) ) == 1 )
            {
              faceBasisInCell = basisInCell;
            }
            else
            {
              WALBERLA_ASSERT( spanningCellLocalVertices.count( basisInCell.at( 3 ) ) == 1 );
              faceBasisInCell = basisInCell;
              faceBasisInCell[2] = basisInCell.at(3);
              faceBasisInCell[3] = basisInCell.at(2);
            }

            const auto leafIndexOnEdgeGhostLayer = indexing::basisConversion( leafIndexInCell, {0, 1, 2, 3}, faceBasisInCell, levelinfo::num_microedges_per_edge( level ) );
            const auto leafOrientationOnEdgeGhostLayer = edgedof::convertEdgeDoFOrientationCellToFace( leafOrientationInCell, faceBasisInCell.at( 0 ), faceBasisInCell.at( 1 ), faceBasisInCell.at( 2 ));

            const auto localFaceIDOnEdge = edge.face_index( facePrimitiveID );
            leafArrayIndexOnEdge = edgedof::macroedge::indexOnNeighborFace( level, leafIndexOnEdgeGhostLayer.x(), localFaceIDOnEdge, leafOrientationOnEdgeGhostLayer );
          }
          else
          {
            // leaf on macro-edge
            WALBERLA_ASSERT_EQUAL( cellLocalIDsOfNeighborFacesWithLeafOnThem.size(), 2 );
            WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) );
            WALBERLA_ASSERT_EQUAL( leafOrientationOnEdge, EdgeDoFOrientation::X );
            leafArrayIndexOnEdge = edgedof::macroedge::index( level, leafIndexOnEdge.x() );
          }

          const auto srcInt = src[ leafArrayIndexOnEdge ];
          MatSetValues(mat, 1, &dstInt, 1, &srcInt, &stencilWeight, ADD_VALUES);
        }
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
  real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  PetscInt * src      = face.getData(srcId)->getPointer( Level );
  PetscInt * dst      = face.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetscInt dstInt;

  using namespace edgedof::macroface;

  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0) {
      dstInt = dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )];
      for(uint_t k = 0; k < neighborsFromHorizontalEdge.size(); ++k){
        srcInt = src[indexFromHorizontalEdge( Level, it.col(), it.row(), neighborsFromHorizontalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[edgedof::stencilIndexFromHorizontalEdge(neighborsFromHorizontalEdge[k])], INSERT_VALUES);
      }
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      dstInt = dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )];
      for(uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k){
        srcInt = src[indexFromDiagonalEdge( Level, it.col(), it.row(), neighborsFromDiagonalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[edgedof::stencilIndexFromDiagonalEdge(neighborsFromDiagonalEdge[k])], INSERT_VALUES);
      }
    }
    if( it.col() != 0) {
      dstInt = dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )];
      for(uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k){
        srcInt = src[indexFromVerticalEdge( Level, it.col(), it.row(), neighborsFromVerticalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[edgedof::stencilIndexFromVerticalEdge(neighborsFromVerticalEdge[k])], INSERT_VALUES);
      }
    }
  }
}

inline void saveFaceOperator3D( const uint_t & level, const Face & face,
                                const PrimitiveStorage                                                    & storage,
                                const PrimitiveDataID< LevelWiseMemory< macroface::StencilMap_T >, Face > & operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & dstId,
                                Mat & mat )
{
  auto opr_data = face.getData(operatorId)->getData( level );
  auto src      = face.getData(srcId)->getPointer( level );
  auto dst      = face.getData(dstId)->getPointer( level );

  for ( const auto & centerIndexInFace : hhg::edgedof::macroface::Iterator( level, 0 ) )
  {
    for ( const auto & faceCenterOrientation : edgedof::faceLocalEdgeDoFOrientations )
    {
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::X && edgedof::isHorizontalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::Y && edgedof::isVerticalEdgeOnBoundary( level, centerIndexInFace ) )
        continue;
      if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::XY && edgedof::isDiagonalEdgeOnBoundary( level, centerIndexInFace )  )
        continue;

      const auto dstIdx = edgedof::macroface::index( level, centerIndexInFace.x(), centerIndexInFace.y(), faceCenterOrientation );
      const auto dstInt = dst[ dstIdx ];

      for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++  )
      {
        const Cell & neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
        const uint_t localFaceID = neighborCell.getLocalFaceID( face.getID() );

        const auto centerIndexInCell = macroface::getIndexInNeighboringMacroCell( centerIndexInFace, face, neighborCellID, storage, level );
        const auto cellCenterOrientation = macroface::getOrientattionInNeighboringMacroCell( faceCenterOrientation, face, neighborCellID, storage );

        for ( const auto & leafOrientation : edgedof::allEdgeDoFOrientations )
        {
          for ( const auto & stencilIt : opr_data[neighborCellID][cellCenterOrientation][leafOrientation] )
          {
            const auto stencilOffset = stencilIt.first;
            const auto stencilWeight = stencilIt.second;

            const auto leafOrientationInFace = macrocell::getOrientattionInNeighboringMacroFace( leafOrientation, neighborCell, localFaceID, storage );

            const auto leafIndexInCell = centerIndexInCell + stencilOffset;
            const auto leafIndexInFace = macrocell::getIndexInNeighboringMacroFace( leafIndexInCell, neighborCell, localFaceID, storage, level );

            WALBERLA_ASSERT_LESS_EQUAL( leafIndexInFace.z(), 1 );

            uint_t leafArrayIndexInFace;
            if ( algorithms::contains( edgedof::faceLocalEdgeDoFOrientations, leafOrientationInFace ) && leafIndexInFace.z() == 0 )
            {
              leafArrayIndexInFace = edgedof::macroface::index( level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace );
            }
            else
            {
              leafArrayIndexInFace = edgedof::macroface::index( level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace, neighborCellID );
            }

            const auto srcInt = src[ leafArrayIndexInFace ];
            MatSetValues(mat, 1, &dstInt, 1, &srcInt, &stencilWeight, ADD_VALUES);
          }
        }
      }
    }
  }
}


inline void saveCellOperator( const uint_t & Level, const Cell & cell,
                              const PrimitiveDataID< LevelWiseMemory< macrocell::StencilMap_T >, Cell> &operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Cell> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Cell> & dstId,
                              Mat & mat )
{
  auto srcData = cell.getData( srcId )->getPointer( Level );
  auto dstData = cell.getData( dstId )->getPointer( Level );
  auto opr_data = cell.getData( operatorId )->getData( Level );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

    if ( macrocell::isInnerXEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
    if ( macrocell::isInnerYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
    if ( macrocell::isInnerZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
    if ( macrocell::isInnerXYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
    if ( macrocell::isInnerXZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
    if ( macrocell::isInnerYZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

    for ( const auto & centerOrientation : innerOrientations )
    {
      const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
      const auto dstInt      = dstData[ dstArrayIdx ];

      for ( const auto & leafOrientation : edgedof::allEdgeDoFOrientations )
      {
        const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
        for ( const auto & neighbor : edgeDoFNeighbors )
        {
          const auto   srcIdx        = it + neighbor;
          const auto   srcArrayIdx   = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
          const auto   srcInt        = srcData[ srcArrayIdx ];
          const real_t stencilWeight = opr_data[centerOrientation][leafOrientation][neighbor];
          MatSetValues(mat, 1, &dstInt, 1, &srcInt, &stencilWeight, ADD_VALUES);
        }
      }
    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;

    const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
    const auto dstInt      = dstData[ dstArrayIdx ];

    for ( const auto & leafOrientation : edgedof::allEdgeDoFOrientations )
    {
      const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
      for ( const auto & neighbor : edgeDoFNeighbors )
      {
        const auto   srcIdx        = it + neighbor;
        const auto   srcArrayIdx   = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
        const auto   srcInt        = srcData[ srcArrayIdx ];
        const real_t stencilWeight = opr_data[centerOrientation][leafOrientation][neighbor];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &stencilWeight, ADD_VALUES);
      }
    }
  }
}


template<class OperatorType>
inline void createMatrix(OperatorType& opr, EdgeDoFFunction< PetscInt > & src, EdgeDoFFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  auto storage = src.getStorage();

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
        saveFaceOperator3D(level, face, *storage, opr.getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
      }
      else
      {
        saveFaceOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
      }
    }
  }

  for (auto& it : opr.getStorage()->getCells()) {
    Cell& cell = *it.second;

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