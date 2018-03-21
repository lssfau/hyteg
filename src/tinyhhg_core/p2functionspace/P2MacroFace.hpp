#pragma once

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/primitiveid.hpp"

namespace hhg {
namespace P2 {
namespace macroface {

void smoothJacobiVertexDoF(const uint_t &level, Face &face,
                           const PrimitiveDataID <StencilMemory<real_t>, Face> &vertexDoFStencilID,
                           const PrimitiveDataID <FunctionMemory<real_t>, Face> &srcVertexDoFID,
                           const PrimitiveDataID <FunctionMemory<real_t>, Face> &dstVertexDoFID,
                           const PrimitiveDataID <StencilMemory<real_t>, Face> &edgeDoFStencilID,
                           const PrimitiveDataID <FunctionMemory<real_t>, Face> &srcEdgeDoFID,
                           const PrimitiveDataID <FunctionMemory<real_t>, Face> &rhsVertexDoFID);

void smoothJacobiEdgeDoF(const uint_t &Level, Face &face,
                         const PrimitiveDataID <StencilMemory<real_t>, Face> &vertexDoFStencilID,
                         const PrimitiveDataID <FunctionMemory<real_t>, Face> &srcVertexDoFID,
                         const PrimitiveDataID <StencilMemory<real_t>, Face> &edgeDoFStencilID,
                         const PrimitiveDataID <FunctionMemory<real_t>, Face> &srcEdgeDoFID,
                         const PrimitiveDataID <FunctionMemory<real_t>, Face> &dstEdgeDoFID,
                         const PrimitiveDataID <FunctionMemory<real_t>, Face> &rhsEdgeDoFID);

}
}
}