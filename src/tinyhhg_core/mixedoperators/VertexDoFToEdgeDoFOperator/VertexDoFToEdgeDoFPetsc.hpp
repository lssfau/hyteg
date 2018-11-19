
#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace VertexDoFToEdgeDoF {

using walberla::real_t;
using walberla::uint_t;

#ifdef HHG_BUILD_WITH_PETSC


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
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF )
    {
      srcInt = src[vertexdof::macroedge::indexFromHorizontalEdge( Level, it.col(), neighbor )];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
    }

    if( edge.getNumNeighborFaces() == 2 )
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF )
      {
        srcInt = src[vertexdof::macroedge::indexFromHorizontalEdge( Level, it.col(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
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
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
      }
    }

    if( it.col() != 0 )
    {
      dstInt = dst[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromVerticalEdge )
      {
        srcInt = src[vertexdof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromVerticalEdge( neighbor ) ], INSERT_VALUES );
      }
    }

    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1) )
    {
      dstInt = dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromDiagonalEdge )
      {
        srcInt = src[vertexdof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), neighbor )];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromDiagonalEdge( neighbor ) ], INSERT_VALUES );
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


template<class OperatorType>
inline void createMatrix(OperatorType& opr, P1Function< PetscInt > & src, EdgeDoFFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{

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
