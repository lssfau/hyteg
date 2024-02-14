/*
 * Copyright (c) 2019-2023 Nils Kohl, Dominik Thoennes, Michael Zikeli.
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

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include <map>
#include "core/Macros.h"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

template < typename ValueType >
void apply_3D_macrocell_vertexdof_to_vertexdof_add(ValueType * RESTRICT _data_p1CellDstAdd, ValueType const * RESTRICT const _data_p1CellSrcAdd, int level, std::map< hyteg::indexing::Index, ValueType >& p1CellStencil);

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg