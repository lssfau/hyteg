#pragma once

#include "core/DataTypes.h"

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


namespace hhg {
namespace P2 {

namespace vertex {

void smoothGSvertexDoF( uint_t                                                     level,
                        Vertex&                                                    vertex,
                        const PrimitiveDataID< StencilMemory< real_t >, Vertex >&  vertexDoFStencilID,
                        const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstVertexDoFID,
                        const PrimitiveDataID< StencilMemory< real_t >, Vertex >&  edgeDoFStencilID,
                        const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstEdgeDoFID,
                        const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& rhsVertexDoFID );

} // namespace vertex

//namespace edge {
//
//void smoothGSvertexDoF( const uint_t&                                            Level,
//                        Edge&                                                    edge,
//                        const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexDoFStencilID,
//                        const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
//                        const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeDoFStencilID,
//                        const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
//                        const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsVertexDoFID );
//
//void smoothGSedgeDoF( const uint_t&                                            Level,
//                      Edge&                                                    edge,
//                      const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexDoFStencilID,
//                      const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
//                      const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeDoFStencilID,
//                      const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
//                      const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsEdgeDoFID );
//
//} // namespace edge
//
//namespace face {
//
//void smoothGSvertexDoF( const uint_t&                                            Level,
//                        Face&                                                    face,
//                        const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexDoFStencilID,
//                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstVertexDoFID,
//                        const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeDoFStencilID,
//                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstEdgeDoFID,
//                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsVertexDoFID );
//
//void smoothGSedgeDoF( const uint_t&                                            Level,
//                      Face&                                                    face,
//                      const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexDoFStencilID,
//                      const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstVertexDoFID,
//                      const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeDoFStencilID,
//                      const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstEdgeDoFID,
//                      const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsEdgeDoFID );
//
//} // namespace face

} // namespace P2
} // namespace hhg