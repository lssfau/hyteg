/*
* Copyright (c) 2017-2023 Nils Kohl.
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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

static void primitiveStorageParallelSetupWrite( const std::string& meshFile, uint_t numProcesses, const std::string& file )
{
   WALBERLA_ROOT_SECTION()
   {
      MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFile );
      SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
      loadbalancing::roundRobinVolume( setupStorage, numProcesses );

      WALBERLA_LOG_INFO( setupStorage );

      setupStorage.writeToFile( file, numProcesses );
   }
}

static void primitiveStorageParallelSetupRead( const std::string& meshFile, uint_t numProcesses, const std::string& file )
{
   // Build storage as usual.
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
   loadbalancing::roundRobinVolume( setupStorage, numProcesses );
   PrimitiveStorage storage( setupStorage );

   // Load storage from file
   PrimitiveStorage storageRead( file );

   auto info     = storage.getGlobalInfo();
   auto infoRead = storageRead.getGlobalInfo();

   WALBERLA_LOG_INFO_ON_ROOT( "Storage built during run time:" )
   WALBERLA_LOG_INFO_ON_ROOT( info );

   WALBERLA_LOG_INFO_ON_ROOT( "Storage read from file:" )
   WALBERLA_LOG_INFO_ON_ROOT( infoRead );

   auto ngf = storageRead.getNumberOfGlobalFaces();
   WALBERLA_CHECK_GREATER( ngf, 0 );

   WALBERLA_CHECK_EQUAL( info, infoRead );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const auto numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   hyteg::primitiveStorageParallelSetupWrite( "../../../data/meshes/bfs_126el.msh", numProcesses, "test_00.data" );
   hyteg::primitiveStorageParallelSetupWrite( "../../../data/meshes/3D/cube_24el.msh", numProcesses, "test_01.data" );

   hyteg::primitiveStorageParallelSetupRead( "../../../data/meshes/bfs_126el.msh", numProcesses, "test_00.data" );
   hyteg::primitiveStorageParallelSetupRead( "../../../data/meshes/3D/cube_24el.msh", numProcesses, "test_01.data" );

   return EXIT_SUCCESS;
}
