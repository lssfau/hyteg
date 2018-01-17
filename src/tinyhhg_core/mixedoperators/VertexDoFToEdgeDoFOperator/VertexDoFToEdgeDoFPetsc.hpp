
#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace VertexDoFToEdgeDoF {

using walberla::real_t;
using walberla::uint_t;

#ifdef HHG_BUILD_WITH_PETSC

template<size_t Level>
inline void saveEdgeOperatorTmpl( const Edge & edge,
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
    dstInt = dst[ edgedof::macroedge::indexFromHorizontalEdge<Level>( it.col(), stencilDirection::EDGE_HO_C ) ];

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF )
    {
      srcInt = src[ vertexdof::macroedge::indexFromHorizontalEdge< Level >( it.col(), neighbor ) ];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF )
    {
      srcInt = src[ vertexdof::macroedge::indexFromHorizontalEdge< Level >( it.col(), neighbor ) ];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
    }

    if( edge.getNumNeighborFaces() == 2 )
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF )
      {
        srcInt = src[ vertexdof::macroedge::indexFromHorizontalEdge< Level >( it.col(), neighbor ) ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
      }
    }
  }
}

SPECIALIZE(void, saveEdgeOperatorTmpl, saveEdgeOperator);

template<size_t Level>
inline void saveFaceOperatorTmpl( const Face & face,
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
      dstInt = dst[ edgedof::macroface::indexFromHorizontalEdge<Level>( it.col(), it.row(), stencilDirection::EDGE_HO_C ) ];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromHorizontalEdge )
      {
        srcInt = src[ vertexdof::macroface::indexFromHorizontalEdge< Level >( it.col(), it.row(), neighbor ) ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ], INSERT_VALUES );
      }
    }

    if( it.row() != 0 )
    {
      dstInt = dst[ edgedof::macroface::indexFromVerticalEdge<Level>( it.col(), it.row(), stencilDirection::EDGE_VE_C ) ];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromVerticalEdge )
      {
        srcInt = src[ vertexdof::macroface::indexFromVerticalEdge< Level >( it.col(), it.row(), neighbor ) ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromVerticalEdge( neighbor ) ], INSERT_VALUES );
      }
    }

    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1) )
    {
      dstInt = dst[ edgedof::macroface::indexFromDiagonalEdge<Level>( it.col(), it.row(), stencilDirection::EDGE_DI_C ) ];
      for ( const auto & neighbor : vertexdof::macroface::neighborsFromDiagonalEdge )
      {
        srcInt = src[ vertexdof::macroface::indexFromDiagonalEdge< Level >( it.col(), it.row(), neighbor ) ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ vertexdof::stencilIndexFromDiagonalEdge( neighbor ) ], INSERT_VALUES );
      }
    }

  }
}

SPECIALIZE(void, saveFaceOperatorTmpl, saveFaceOperator);

template<class OperatorType>
inline void createMatrix(OperatorType& opr, P1Function< PetscInt > & src, EdgeDoFFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{

  for (auto& it : opr.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      saveEdgeOperator(level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
    }
  }

  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      saveFaceOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
    }
  }
}

#endif

}
}
