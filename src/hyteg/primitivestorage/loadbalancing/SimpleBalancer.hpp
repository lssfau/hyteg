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

namespace hyteg {
namespace loadbalancing {

using walberla::uint_t;

/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on a certain rank
void allPrimitivesOnOneRank( SetupPrimitiveStorage & storage, const uint_t & targetRank );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on root
void allPrimitivesOnRoot( SetupPrimitiveStorage & storage );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributed all primitives in a round robin fashion
void roundRobin( SetupPrimitiveStorage & storage );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributes all primitives in a greedy fashion.
/// It is expected to result in a much better edge-cut ratio than the roundRobin algorithm. But still not optimal.
void greedy( SetupPrimitiveStorage & storage );


} // namespace loadbalancing
} // namespace hyteg
