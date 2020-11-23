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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "core/mpi/MPIWrapper.h"

namespace hyteg {
namespace loadbalancing {
namespace distributed {


MigrationMap_T parmetis( PrimitiveStorage & storage );

/// \brief Performs a round robin distribution  in parallel.
///
/// \param storage                 the PrimitiveStorage, the primitives are distributed on
MigrationInfo roundRobin( PrimitiveStorage & storage );

/// \brief Performs a round robin distribution to a subset of processes in parallel.
///
/// \param storage                 the PrimitiveStorage, the primitives are distributed on
/// \param numberOfTargetProcesses if smaller than total number of processes, 
///                                the round robin distributes
///                                among the first numberOfTargetProcesses processes
MigrationInfo roundRobin( PrimitiveStorage & storage, uint_t numberOfTargetProcesses );

/// \brief Performs a round robin distribution to a subset of processes in parallel.
///
/// \param storage the PrimitiveStorage, the primitives are distributed on
/// \param minRank lowest included rank in the distribution
/// \param maxRank highest included rank in the distribution
MigrationInfo roundRobin( PrimitiveStorage & storage, uint_t minRank, uint_t maxRank );

/// \brief Performs a round robin distribution to a subset of processes in parallel.
///
/// \param storage the PrimitiveStorage, the primitives are distributed on
/// \param interval the distribution happens on every interval'th process, e.g. if interval == 24,
///                 the storage is distributed among rank 0, 24, 48, ...
MigrationInfo roundRobinInterval( PrimitiveStorage & storage, uint_t interval );

/// \brief Performs a round robin distribution to a subset of processes in parallel.
///
/// \param storage      the PrimitiveStorage, the primitives are distributed on
/// \param interval     the distribution happens on every interval'th process, e.g. if interval == 24,
///                     the storage is distributed among rank 0, 24, 48, ...
/// \param numProcesses number of processes that will obtain primitives
MigrationInfo roundRobinInterval( PrimitiveStorage & storage, uint_t interval, uint_t numProcesses );

/// \brief Reverses the previous distribution previously performed by some load balancing algorithm.
///
/// \param originalMigrationInfo the MigrationInfo that was calculated by the previous algorithm
/// \param storageToRedistribute this PrimitiveStorage is redistributed
MigrationInfo reverseDistribution( const MigrationInfo & originalMigrationInfo, PrimitiveStorage& storageToRedistribute );

/// \brief Same as reverseDistribution() but does not perform the actual migration.
///
/// \param originalMigrationInfo the MigrationInfo that was calculated by the previous algorithm
MigrationInfo reverseDistributionDry( const MigrationInfo & originalMigrationInfo );

void diffusiveSmooth( PrimitiveStorage& storage, uint_t outerIterations, uint_t smoothingIterations );


}
}
}
