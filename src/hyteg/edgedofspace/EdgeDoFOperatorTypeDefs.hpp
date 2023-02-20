/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include <map>

#include "core/DataTypes.h"

namespace hyteg {


namespace edgedof {

using walberla::real_t;
using walberla::uint_t;

enum class EdgeDoFOrientation : walberla::uint_t;

namespace macroedge {

/// map[neighborCellID][centerOrientation][leafOrientation][indexOffset] = weight
typedef std::map< uint_t,
                  std::map< edgedof::EdgeDoFOrientation,
                            std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > > >
    StencilMap_T;

} // namespace macroedge

namespace macroface {

/// map[neighborCellID][centerOrientation][leafOrientation][indexOffset] = weight
typedef std::map< uint_t,
                  std::map< edgedof::EdgeDoFOrientation,
                            std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > > >
    StencilMap_T;

} // namespace macroface

namespace macrocell {

/// map[centerOrientation][leafOrientation][indexOffset] = weight
typedef std::map< edgedof::EdgeDoFOrientation,
                  std::map< edgedof::EdgeDoFOrientation, std::map< indexing::IndexIncrement, real_t > > >
    StencilMap_T;

} // namespace macrocell

} // namespace edgedof
} // namespace hyteg