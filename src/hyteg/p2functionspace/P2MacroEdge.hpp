#pragma once

#include "hyteg/FunctionMemory.hpp"
#include "hyteg/PrimitiveID.hpp"
#include "hyteg/StencilMemory.hpp"
#include "hyteg/primitives/Edge.hpp"

namespace hyteg {
namespace P2 {
namespace macroedge {

void smoothSOR( const uint_t&                                            level,
                const Edge&                                              edge,
                const real_t&                                            relax,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToVertexStencilID,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToVertexStencilID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToEdgeStencilID,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToEdgeStencilID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsVertexDoFID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsEdgeDoFID );

inline void smoothGaussSeidel( const uint_t&                                            level,
                               const Edge&                                              edge,
                               const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToVertexStencilID,
                               const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToVertexStencilID,
                               const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
                               const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToEdgeStencilID,
                               const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToEdgeStencilID,
                               const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
                               const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsVertexDoFID,
                               const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsEdgeDoFID )
{
     smoothSOR( level,
                edge,
                1.0,
                vertexToVertexStencilID,
                edgeToVertexStencilID,
                dstVertexDoFID,
                vertexToEdgeStencilID,
                edgeToEdgeStencilID,
                dstEdgeDoFID,
                rhsVertexDoFID,
                rhsEdgeDoFID );
}

void smoothSOR3D(
    const uint_t&                                                                                level,
    const PrimitiveStorage&                                                                      storage,
    Edge&                                                                                        edge,
    const real_t&                                                                                relax,
    const PrimitiveDataID< StencilMemory< real_t >, Edge >&                                      vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFRhsId,
    const bool&                                                                                  backwards = false );

} // namespace macroedge
} // namespace P2
} // namespace hyteg
