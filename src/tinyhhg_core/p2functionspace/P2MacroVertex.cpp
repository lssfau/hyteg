#include "P2MacroVertex.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/primitives/all.hpp"

using walberla::real_t;

namespace hhg {
namespace P2 {

namespace macrovertex {

void smoothSORVertexDoF( uint_t                                                     level,
                         Vertex&                                                    vertex,
                         const real_t&                                              relax,
                         const PrimitiveDataID< StencilMemory< double >, Vertex >&  vertexDoFStencilID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstVertexDoFID,
                         const PrimitiveDataID< StencilMemory< double >, Vertex >&  edgeDoFStencilID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstEdgeDoFID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& rhsVertexDoFID )
{
   real_t* vertexDoFStencil = vertex.getData( vertexDoFStencilID )->getPointer( level );
   real_t* dstVertexDoF     = vertex.getData( dstVertexDoFID )->getPointer( level );
   real_t* edgeDoFStencil   = vertex.getData( edgeDoFStencilID )->getPointer( level );
   real_t* dstEdgeDoF       = vertex.getData( dstEdgeDoFID )->getPointer( level );
   real_t* rhs              = vertex.getData( rhsVertexDoFID )->getPointer( level );

   real_t tmp = 0;
   tmp        = rhs[0];
   for( uint_t i = 0; i < vertex.getData( edgeDoFStencilID )->getSize( level ); ++i )
   {
      tmp -= dstEdgeDoF[i] * edgeDoFStencil[i];
   }
   for( uint_t i = 1; i < vertex.getData( vertexDoFStencilID )->getSize( level ); ++i )
   {
      tmp -= dstVertexDoF[i] * vertexDoFStencil[i];
   }

   dstVertexDoF[0] = (1.0 - relax) * dstVertexDoF[0] + (relax * tmp) / vertexDoFStencil[0];
}
} // namespace vertex

} // namespace P2
} // namespace hhg
