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