
#pragma once

#include "core/DataTypes.h"

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

  for( const auto it : vertexdof::macroedge::Iterator( Level, 1 ) )
  {
    dstInt = dst[vertexdof::macroedge::indexFromVertex( Level, it.col(), stencilDirection::VERTEX_C )];

    for ( const auto & neighbor : edgedof::macroedge::neighborsOnEdgeFromVertex )
    {
      srcInt = src[ edgedof::macroedge::indexFromVertex< Level >( it.col(), neighbor ) ];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
    }

    for ( const auto & neighbor : edgedof::macroedge::neighborsOnSouthFaceFromVertex )
    {
      srcInt = src[ edgedof::macroedge::indexFromVertex< Level >( it.col(), neighbor ) ];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
    }

    if( edge.getNumNeighborFaces() == 2 )
    {
      for ( const auto & neighbor : edgedof::macroedge::neighborsOnNorthFaceFromVertex )
      {
        srcInt = src[ edgedof::macroedge::indexFromVertex< Level >( it.col(), neighbor ) ];
        MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
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

  for ( const auto & it : vertexdof::macroface::Iterator( Level, 1 ) )
  {
    dstInt = dst[vertexdof::macroface::indexFromVertex( Level, it.col(), it.row(),
                                                                 stencilDirection::VERTEX_C )];

    for ( const auto & neighbor : edgedof::macroface::neighborsFromVertex )
    {
      srcInt = src[ edgedof::macroface::indexFromVertex< Level >( it.col(), it.row(), neighbor ) ];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
    }
  }
}

SPECIALIZE(void, saveFaceOperatorTmpl, saveFaceOperator);

template<class OperatorType>
inline void createMatrix(OperatorType& opr, EdgeDoFFunction< PetscInt > & src, P1Function< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  for (auto& it : opr.getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      saveVertexOperator(level, vertex, opr.getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat);
    }
  }

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
