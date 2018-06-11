
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"

namespace hhg {

void P1toP1LinearRestriction::restrictMacroVertex( const real_t *src, real_t *dst, const uint_t & sourceLevel,
                                                   const uint_t & numNeighborEdges ) const
{
  WALBERLA_UNUSED( sourceLevel );
  dst[0] = src[0];

  for (uint_t i = 0; i < numNeighborEdges; ++i) {
    dst[0] += 0.5*src[i+1];
    i += 1;
  }
}

void P1toP1LinearRestriction::restrictMacroEdge( const real_t *src, real_t *dst, const uint_t & sourceLevel,
                                                 const uint_t & numNeighborFaces ) const
{
  size_t rowsize_c = levelinfo::num_microvertices_per_edge( sourceLevel - 1 );

  uint_t i_c;
  for ( i_c = 1; i_c < rowsize_c - 1; ++i_c )
  {
    dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] = 1.0 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                      2 * i_c,
                                                                                                                                                      stencilDirection::VERTEX_C )];

    for ( auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] += 0.5 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                         2 * i_c,
                                                                                                                                                         neighbor )];
    }

    for ( auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] += 0.5 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                         2 * i_c,
                                                                                                                                                         neighbor )];
    }

    if ( numNeighborFaces == 2 )
    {
      for ( auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] += 0.5 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                           2 *
                                                                                                                                                           i_c,
                                                                                                                                                           neighbor )];
      }
    }
  }
}

void P1toP1LinearRestriction::restrictMacroFace( const real_t *src, real_t *dst, const uint_t & sourceLevel,
                                                 const uint_t & numNeighborCells ) const
{
  WALBERLA_UNUSED( numNeighborCells );
  uint_t N_c = levelinfo::num_microvertices_per_edge( sourceLevel - 1 );
  uint_t N_c_i = N_c;

  real_t tmp;

  for ( uint_t j = 1; j < N_c - 2; ++j )
  {
    for ( uint_t i = 1; i < N_c_i - 2; ++i )
    {

      tmp = src[vertexdof::macroface::indexFromVertex( sourceLevel, 2 * i, 2 * j,
                                                       stencilDirection::VERTEX_C )];

      for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter )
      {
        tmp += 0.5 * src[vertexdof::macroface::indexFromVertex( sourceLevel, 2 * i, 2 * j, neighbor )];
      }

      dst[vertexdof::macroface::indexFromVertex( sourceLevel - 1, i, j, stencilDirection::VERTEX_C )] = tmp;
    }

    --N_c_i;
  }
}

}