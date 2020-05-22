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

namespace hyteg {
namespace loadbalancing {
namespace distributed {

/// Maps local primitive ID to target rank.
/// Returned by all distributed load balancing implementations.
/// Is always passed to PrimitiveStorage::migratePrimitives().
typedef std::map< PrimitiveID::IDType, uint_t > MigrationMap_T;

MigrationMap_T parmetis( PrimitiveStorage & storage );

/// \brief Performs a round robin distribution  in parallel.
///
/// \param storage                 the PrimitiveStorage, the primitives are distributed on
MigrationMap_T roundRobin( PrimitiveStorage & storage );

/// \brief Performs a round robin distribution to a subset of processes in parallel.
///
/// \param storage                 the PrimitiveStorage, the primitives are distributed on
/// \param numberOfTargetProcesses if smaller than total number of processes, 
///                                the round robin distributes
///                                among the first numberOfTargetProcesses processes
MigrationMap_T roundRobin( PrimitiveStorage & storage, uint_t numberOfTargetProcesses );

/// \brief Distributes a second PrimitiveStorage equal to another one.
///
/// \param targetDistributionStorage this PrimitiveStorage specifies the desired distribution, it is not modified
/// \param storageToRedistribute this PrimitiveStorage is redistributed
MigrationMap_T copyDistribution( const PrimitiveStorage & targetDistributionStorage, PrimitiveStorage & storageToRedistribute );

/// \brief Same as copyDistribution but only creates mapping, does not perform any migration.
///
/// \param targetDistributionStorage this PrimitiveStorage specifies the desired distribution, it is not modified
/// \param storageToRedistribute this PrimitiveStorage shall be redistributed, but is not modified in this function
MigrationMap_T copyDistributionDry( const PrimitiveStorage & targetDistributionStorage, const PrimitiveStorage & storageToRedistribute );

}
}
}
