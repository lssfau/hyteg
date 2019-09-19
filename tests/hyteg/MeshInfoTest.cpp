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
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "core/Filesystem.h"

namespace hyteg {

static void testMeshInfo()
{
  const std::string meshFileDir = "../../data/meshes/";

  std::vector< std::string > gmshFiles;

  walberla::filesystem::recursive_directory_iterator dirIterator( meshFileDir );
  walberla::filesystem::recursive_directory_iterator dirIteratorEnd;

  for ( ; dirIterator != dirIteratorEnd; dirIterator++ )
  {
    std::string dirOrFilePath = dirIterator->path().string();
    if ( dirOrFilePath.compare( dirOrFilePath.size() - 4, 4, ".msh" ) == 0 )
    {
       gmshFiles.push_back( dirOrFilePath );
    }
  }

  for ( const auto & gmshFile : gmshFiles )
  {
    MeshInfo meshInfo = MeshInfo::fromGmshFile( gmshFile );
    SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
    loadbalancing::roundRobin( setupStorage );
    WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
    PrimitiveStorage storage( setupStorage );
  }
}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testMeshInfo();

   return EXIT_SUCCESS;
}
