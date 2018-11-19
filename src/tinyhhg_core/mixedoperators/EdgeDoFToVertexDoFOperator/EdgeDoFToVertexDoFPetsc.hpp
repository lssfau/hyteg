
#pragma once

#include "core/DataTypes.h"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {

using walberla::real_t;
using walberla::uint_t;

#ifdef HHG_BUILD_WITH_PETSC

inline void saveVertexOperator( const uint_t & level,
                                const Vertex & vertex,
                                const PrimitiveDataID< StencilMemory< real_t >, Vertex> & operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Vertex> & srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Vertex> & dstId,
                                Mat & mat )
{
  auto opr_data = vertex.getData(operatorId)->getPointer( level );
  auto src = vertex.getData(srcId)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );

  WALBERLA_ASSERT_LESS_EQUAL( vertex.getNumNeighborEdges() + vertex.getNumNeighborFaces(), vertex.getData(srcId)->getSize( level ),
                              "Stencil memory size is smaller than it should be." );

  MatSetValues(mat, 1, dst, (PetscInt) ( vertex.getNumNeighborEdges() + vertex.getNumNeighborFaces() ), src, opr_data, INSERT_VALUES );
}

inline void saveEdgeOperator( const uint_t & Level,  const Edge & edge,
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

  for( const auto it : vertexdof::macroedge::Iterator( Level, 1 ) )
  {
    dstInt = dst[vertexdof::macroedge::indexFromVertex( Level, it.col(), stencilDirection::VERTEX_C )];

    for ( const auto & neighbor : edgedof::macroedge::neighborsOnEdgeFromVertex )
    {
      srcInt = src[edgedof::macroedge::indexFromVertex( Level, it.col(), neighbor )];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
    }

    for ( const auto & neighbor : edgedof::macroedge::neighborsOnSouthFaceFromVertex )
    {
      srcInt = src[edgedof::macroedge::indexFromVertex( Level, it.col(), neighbor )];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
    }

    if( edge.getNumNeighborFaces() == 2 )
    {
      for ( const auto & neighbor : edgedof::macroedge::neighborsOnNorthFaceFromVertex )
      {
        srcInt = src[edgedof::macroedge::indexFromVertex( Level, it.col(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
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

  for ( const auto & it : vertexdof::macroface::Iterator( Level, 1 ) )
  {
    dstInt = dst[vertexdof::macroface::indexFromVertex( Level, it.col(), it.row(),
                                                                 stencilDirection::VERTEX_C )];

    for ( const auto & neighbor : edgedof::macroface::neighborsFromVertex )
    {
      srcInt = src[edgedof::macroface::indexFromVertex( Level, it.col(), it.row(), neighbor )];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
    }
  }
}


inline void saveCellOperator( const uint_t & Level, const Cell & cell,
                              const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell> &operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Cell> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Cell> & dstId,
                              Mat & mat )
{
  auto opr_data = cell.getData(operatorId)->getData( Level );
  PetscInt * src  = cell.getData(srcId)->getPointer( Level );
  PetscInt * dst  = cell.getData(dstId)->getPointer( Level );

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
    const auto dstArrayIdx = vertexdof::macrocell::index( Level, it.x(), it.y(), it.z() );
    const auto dstInt      = dst[ dstArrayIdx ];

    for ( const auto & orientation : edgedof::allEdgeDoFOrientations )
    {
      const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromVertexDoFInMacroCell( orientation );
      for ( const auto & neighbor : edgeDoFNeighbors )
      {
        const auto   srcIdx      = it + neighbor;
        const auto   srcArrayIdx = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), orientation );
        const auto   srcInt      = src[ srcArrayIdx ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[orientation][neighbor], ADD_VALUES );
      }
    }
  }
}


template<class OperatorType>
inline void createMatrix(OperatorType& opr, EdgeDoFFunction< PetscInt > & src, P1Function< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  for (auto& it : opr.getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag))
    {
      saveVertexOperator(level, vertex, opr.getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat);
    }
  }

  for (auto& it : opr.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag))
    {
      saveEdgeOperator(level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
    }
  }

  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      saveFaceOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
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
