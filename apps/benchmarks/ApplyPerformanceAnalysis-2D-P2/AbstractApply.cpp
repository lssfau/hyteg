#include "AbstractApply.hpp"

#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"

/// these are copies of the abstract apply function found in HyTeG
/// adjusted a little for better comparison

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

void apply_2d_vertex_to_vertex( double* WALBERLA_RESTRICT dstPtr,
                                double const* WALBERLA_RESTRICT const srcPtr,
                                double const* WALBERLA_RESTRICT const stencilPtr,
                                const uint_t                          level )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < rowsize -j - 1; ++i )
      {
         auto tmp = real_t( 0 );
         for ( const auto direction : vertexdof::macroface::neighborsWithCenter )
         {
            tmp += stencilPtr[vertexdof::stencilIndexFromVertex( direction )] *
                   srcPtr[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
         }
         dstPtr[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] = tmp;
      }
   }
}