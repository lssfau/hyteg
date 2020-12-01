/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {
namespace loadbalancing {

using walberla::uint_t;

/// \brief Load balancing with the ParMETIS library.
///
/// The load balancer works in 2 steps:
///
/// First, a graph is built treating only the volume-primitives as nodes and
/// builds edges along the neighboring volume-primitives. Two volume-primitives
/// are in this case neighbors if they share at least one macro-vertex.
/// Each graph-edge is assigned a weight, that is larger for volume-primitive pairs
/// that share many vertices and decreases with fewer shared vertices.
/// The resulting graph is partitioned with ParMETIS
///
/// In a second step, the remaining interface primitives are assigned depending
/// on the neighboring volume-primitives: of all ranks of neighboring volume-primitives,
/// the interface primitive is assigned to the rank that most neighbors are already assigned to.
///
/// ParMETIS is run in parallel. However, for massively parallel runs, the number
/// of processes that execute ParMETIS may be reduced to subCommunicatorSize processes.
/// If specified, a subcommunicator of that size is created.
///
/// \param setupStorage the SetupPrimitiveStorage that shall be re-distributed
/// \param subCommunicatorSize number of processes that call the ParMETIS graph partitioning routine
///
void parmetis( SetupPrimitiveStorage& setupStorage, uint_t subCommunicatorSize );
void parmetis( SetupPrimitiveStorage& setupStorage );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on a certain rank
void allPrimitivesOnOneRank( SetupPrimitiveStorage & storage, const uint_t & targetRank );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that locates all primitives on root
void allPrimitivesOnRoot( SetupPrimitiveStorage & storage );


/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributed all primitives in a round robin fashion
void roundRobin( SetupPrimitiveStorage & storage );
void roundRobin( SetupPrimitiveStorage & storage, uint_t numRanks );

/// \brief First balances the volume primitives by a round robin approach.
///        Then assigns lower-dimensional primitives to neighbor-ranks that carry least amount of similar primitives.
///        This algorithm is optimized towards low communication overhead, lower-dim primitives may not be distributed equally.
void roundRobinVolume( SetupPrimitiveStorage& storage );
void roundRobinVolume( SetupPrimitiveStorage& storage, uint_t numRanks );

/// \brief Load balancing function for \ref SetupPrimitiveStorage that distributes all primitives in a greedy fashion.
/// It is expected to result in a much better edge-cut ratio than the roundRobin algorithm. But still not optimal.
void greedy( SetupPrimitiveStorage & storage );


} // namespace loadbalancing
} // namespace hyteg
