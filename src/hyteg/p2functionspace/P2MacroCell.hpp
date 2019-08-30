
#pragma once

#include "core/DataTypes.h"

#include "hyteg/LevelWiseMemory.hpp"
#include "hyteg/StencilMemory.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "hyteg/primitives/Cell.hpp"

namespace hyteg {
namespace P2 {
namespace macrocell {

using walberla::real_t;
using walberla::uint_t;

void smoothSOR(
    const uint_t&                                                                                level,
    Cell&                                                                                        cell,
    const real_t&                                                                                relax,
    const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >&        vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     edgeDoFRhsId );

} // namespace macrocell
} // namespace P2
} // namespace hyteg