/*
* Copyright (c) 2017-2023 Dominik Thoennes.
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

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/debug/CheckFunctions.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::uint_c;

namespace hyteg {
void testNonDistributedSetupStorage()
{
   MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 2, 2, real_c( 1.0 ), real_c( 2.0 ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   PrimitiveStorage      primitiveStorage( setupStorage, 1 );

   std::vector<PrimitiveID> distributed;
   primitiveStorage.getNeighboringPrimitiveIDs(distributed);

   SetupPrimitiveStorage setupStorageRoot( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ), true );
   PrimitiveStorage      primitiveStorageRoot( setupStorageRoot, 1 );

   std::vector<PrimitiveID> rootOnly;
   primitiveStorageRoot.getNeighboringPrimitiveIDs(rootOnly);

   WALBERLA_CHECK_EQUAL(primitiveStorage.getPrimitiveIDs(), primitiveStorageRoot.getPrimitiveIDs())
   WALBERLA_CHECK_EQUAL(distributed, rootOnly)
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testNonDistributedSetupStorage();

   return EXIT_SUCCESS;
}