#pragma once

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"

namespace hhg {
namespace P2 {
namespace macroedge {

void smoothGaussSeidl( const uint_t &                                           level,
                       const Edge &                                             edge,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge > & vertexToVertexStencilID,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge > & edgeToVertexStencilID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge > &dstVertexDoFID,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge > & vertexToEdgeStencilID,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge > & edgeToEdgeStencilID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge > &dstEdgeDoFID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge > &rhsVertexDoFID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge > &rhsEdgeDoFID );

} // namespace macroedge
} // namespace P2
} // namespace hhg