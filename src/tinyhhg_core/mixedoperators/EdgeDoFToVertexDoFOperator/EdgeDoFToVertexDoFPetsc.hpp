
#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace EdgeDoFToVertexDoF {

using walberla::real_t;
using walberla::uint_t;

#ifdef HHG_BUILD_WITH_PETSC

inline void saveOperator( const uint_t & level,
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
inline void saveOperatorTmpl( const Edge & edge,
                              const PrimitiveDataID< StencilMemory< real_t >, Edge>    & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & dstId,
                              Mat & mat )
{
  const real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  const real_t * src      = edge.getData(srcId)->getPointer( Level );
  const real_t * dst      = edge.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetcsInt dstInt;

  for( const auto it : vertexdof::macroedge::Iterator( Level, 1 ) )
  {
    dstInt = dst[ vertexdof::macroedge::indexFromVertex<Level>( it.col(), stencilDirection::VERTEX_C) ];

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

SPECIALIZE(void, saveOperatorTmpl, saveOperator);

template<size_t Level>
inline void saveOperatorTmpl( const Face & face,
                              const PrimitiveDataID< StencilMemory< real_t >, Face>    & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & dstId,
                              Mat & mat )
{
  const real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  const real_t * src      = face.getData(srcId)->getPointer( Level );
  const real_t * dst      = face.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetcsInt dstInt;

  for ( const auto & it : vertexdof::macroface::Iterator( Level, 1 ) )
  {
    dstInt = dst[ vertexdof::macroface::indexFromVertex<Level>( it.col(), it.row(), stencilDirection::VERTEX_C) ];

    for ( const auto & neighbor : edgedof::macroface::neighborsFromVertex )
    {
      srcInt = src[ edgedof::macroface::indexFromVertex< Level >( it.col(), it.row(), neighbor ) ];
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[ edgedof::stencilIndexFromVertex( neighbor ) ], INSERT_VALUES );
    }
  }
}

SPECIALIZE(void, saveOperatorTmpl, saveOperator);

#endif

}
}
