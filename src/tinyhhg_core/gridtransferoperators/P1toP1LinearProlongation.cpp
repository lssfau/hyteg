
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"

namespace hhg {

void P1toP1LinearProlongation::prolongateMacroVertex( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const
{
  dst[0] = src[0];
}

void P1toP1LinearProlongation::prolongateMacroEdge( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const
{
  uint_t rowsize_c = levelinfo::num_microvertices_per_edge( sourceLevel );
  uint_t i_c;
  
  for ( i_c = 1; i_c < rowsize_c - 1; ++i_c ) {

    dst[vertexdof::macroedge::indexFromVertex( sourceLevel + 1, 2 * i_c, stencilDirection::VERTEX_C )] =
      src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c, stencilDirection::VERTEX_C )];
    
    dst[vertexdof::macroedge::indexFromVertex( sourceLevel + 1, 2 * i_c - 1, stencilDirection::VERTEX_C )]
    = 0.5 * ( src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c - 1, stencilDirection::VERTEX_C )]
              + src[vertexdof::macroedge::indexFromVertex(
    sourceLevel, i_c, stencilDirection::VERTEX_C )] );
  }

  dst[vertexdof::macroedge::indexFromVertex( sourceLevel + 1, 2 * i_c - 1, stencilDirection::VERTEX_C )] =
    0.5 * ( src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c - 1, stencilDirection::VERTEX_C )] + src[vertexdof::macroedge::indexFromVertex(
  sourceLevel, i_c, stencilDirection::VERTEX_C )] );
}

void P1toP1LinearProlongation::prolongateMacroFace( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const
{
  typedef stencilDirection SD;
  using namespace vertexdof::macroface;

  uint_t N_c = levelinfo::num_microvertices_per_edge(sourceLevel);
  uint_t N_c_i = N_c;

  uint_t j;

  for (uint_t i = 1; i < N_c - 1; ++i) {
    for (j = 1; j < N_c_i - 2; ++j) {
      dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j, SD::VERTEX_C )] = src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )];
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] =
      0.5*( src[indexFromVertex( sourceLevel, i - 1, j, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel, i, j - 1, SD::VERTEX_C )]);
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )]
                                                                                  + src[indexFromVertex( sourceLevel, i - 1, j, SD::VERTEX_C )]);
      dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex(
      sourceLevel,
      i, j,
      SD::VERTEX_C )] + src[indexFromVertex(
      sourceLevel, i, j - 1, SD::VERTEX_C )]);
    }

    dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex(
    sourceLevel, i - 1, j, SD::VERTEX_C )] + src[indexFromVertex(
    sourceLevel, i, j - 1, SD::VERTEX_C )]);
    dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex(
    sourceLevel, i,
    j, SD::VERTEX_C )] + src[indexFromVertex(
    sourceLevel, i - 1, j, SD::VERTEX_C )]);
    dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex(
    sourceLevel, i,
    j, SD::VERTEX_C )] + src[indexFromVertex(
    sourceLevel, i, j - 1, SD::VERTEX_C )]);

    --N_c_i;
  }
}

}