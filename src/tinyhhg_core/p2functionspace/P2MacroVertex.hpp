#pragma once

#include "core/DataTypes.h"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"

class Vertex;
class Edge;
class Face;

template<typename ValueType>
class StencilMemory;
template< typename ValueType >
class FunctionMemory;
template< typename DataType, typename PrimitiveType >
class PrimitiveDataID;

using walberla::uint_t;
using walberla::real_t;


namespace hyteg {
namespace P2 {

namespace macrovertex {

void smoothSORVertexDoF( uint_t                                                     level,
                         Vertex&                                                    vertex,
                         const real_t&                                              relax,
                         const PrimitiveDataID< StencilMemory< real_t >, Vertex >&  vertexDoFStencilID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstVertexDoFID,
                         const PrimitiveDataID< StencilMemory< real_t >, Vertex >&  edgeDoFStencilID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstEdgeDoFID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& rhsVertexDoFID );

void smoothSOR3D(
    const uint_t&                                                                                    level,
    const PrimitiveStorage&                                                                          storage,
    Vertex&                                                                                          vertex,
    const real_t&                                                                                    relax,
    const PrimitiveDataID< StencilMemory< real_t >, Vertex >&                                        vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex >& edgeToVertexOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                                       vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                                       vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                                       edgeDoFDstId );

} // namespace vertex

} // namespace P2
} // namespace hyteg
