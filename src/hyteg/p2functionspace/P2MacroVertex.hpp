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

#include "core/DataTypes.h"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"

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

void smoothJacobiVertexDoF( uint_t                                                     level,
                            Vertex&                                                    vertex,
                            const real_t&                                              relax,
                            const PrimitiveDataID< StencilMemory < real_t >, Vertex >& vertexToVertexStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& srcVertexDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstVertexDoFID,
                            const PrimitiveDataID< StencilMemory < real_t >, Vertex >& edgeToVertexStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& srcEdgeDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& rhsVertexDoFID );

} // namespace vertex

} // namespace P2
} // namespace hyteg
