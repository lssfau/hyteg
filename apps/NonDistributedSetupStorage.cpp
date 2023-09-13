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

#include <unistd.h>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::uint_t;

namespace hyteg {
void runTest()
{
   //uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   std::shared_ptr< PrimitiveStorage > storage;
   {
      //const std::string meshFileName = "../../data/meshes/porous_fine.msh";
      //const std::string meshFileName = "../../data/meshes/bfs_126el.msh";
      const std::string meshFileName     = "../../data/meshes/tri_2el.msh";
      const std::string distributionFile = "../../output/PrimitiveStorageTestDistribution.csv";

      MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      //loadbalancing::greedy( setupStorage );

      WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

      //auto foo = setupStorage.getVerticesOnRank(1);
      WALBERLA_LOG_INFO_ON_ROOT( "Building PrimitiveStorage" );

      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }

   storage->checkConsistency();
   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_CRITICAL_SECTION_START
   WALBERLA_LOG_INFO( storageInfo );
   WALBERLA_CRITICAL_SECTION_END
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::runTest();
   return EXIT_SUCCESS;
}