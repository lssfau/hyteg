/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFApply.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"

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
    const PrimitiveDataID< LevelWiseMemory< vertexdof::macroedge::StencilMap_T >, Edge >&        vertexToVertexOperator3DId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFRhsId,
    const bool&                                                                                  backwards = false );

void smoothJacobi( const uint_t&                                            level,
                   const Edge&                                              edge,
                   const real_t&                                            relax,
                   const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToVertexStencilID,
                   const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToVertexStencilID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcVertexDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
                   const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToEdgeStencilID,
                   const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToEdgeStencilID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcEdgeDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsVertexDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsEdgeDoFID );

} // namespace macroedge
} // namespace P2
} // namespace hyteg
