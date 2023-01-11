/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/PrimitiveID.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "hyteg/primitives/Face.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {

real_t evaluate( const uint_t&                                            level,
                 Face&                                                    face,
                 const Point3D&                                           coordinates,
                 const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcVertexDoFID,
                 const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcEdgeDoFID );

void evaluateGradient( const uint_t&                                            level,
                       Face&                                                    face,
                       const Point3D&                                           coordinates,
                       const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcVertexDoFID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcEdgeDoFID,
                       Point3D&                                                 gradient );

void smoothJacobiVertexDoF( const uint_t&                                            level,
                            const Face&                                              face,
                            const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexDoFStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcVertexDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstVertexDoFID,
                            const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeDoFStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcEdgeDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsVertexDoFID,
                            const real_t                                             dampingFactor = real_c( 2.0 / 3.0 ) );

void smoothJacobiEdgeDoF( const uint_t&                                            Level,
                          const Face&                                              face,
                          const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexDoFStencilID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcVertexDoFID,
                          const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeDoFStencilID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcEdgeDoFID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstEdgeDoFID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsEdgeDoFID,
                          const real_t                                             dampingFactor = real_c( 2.0 / 3.0 ) );

void smoothSOR( const uint_t&                                            level,
                const Face&                                              face,
                const real_t&                                            relax,
                const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexToVertexStencilID,
                const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeToVertexStencilID,
                const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstVertexDoFID,
                const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexToEdgeStencilID,
                const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeToEdgeStencilID,
                const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstEdgeDoFID,
                const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsVertexDoFID,
                const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsEdgeDoFID );

void smoothSOR3D(
    const uint_t&                                                                                level,
    const PrimitiveStorage&                                                                      storage,
    const Face&                                                                                  face,
    const real_t&                                                                                relax,
    const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face >&        vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >&                                     edgeDoFRhsId );

void smoothGaussSeidel( const uint_t&                                            level,
                        const Face&                                              face,
                        const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexToVertexStencilID,
                        const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeToVertexStencilID,
                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstVertexDoFID,
                        const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexToEdgeStencilID,
                        const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeToEdgeStencilID,
                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstEdgeDoFID,
                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsVertexDoFID,
                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsEdgeDoFID );

} // namespace macroface
} // namespace P2
} // namespace hyteg
