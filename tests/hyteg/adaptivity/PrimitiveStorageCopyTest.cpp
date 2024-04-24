/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

static void testPrimitiveStorageCopy( const MeshInfo& mesh )
{
   const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   const uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   SetupPrimitiveStorage setupStorage( mesh, numProcesses );
   const auto            numGlobalPrimitives = setupStorage.getNumberOfPrimitives();

   loadbalancing::roundRobin( setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // create copy of the storage
   auto copiedStorage = storage->createCopy();

   WALBERLA_CHECK( storage->getPrimitiveIDs() == copiedStorage->getPrimitiveIDs() );

   loadbalancing::distributed::roundRobin( *copiedStorage );

   WALBERLA_CHECK( storage->getPrimitiveIDs() == copiedStorage->getPrimitiveIDs() );

   const auto migrationInfo = loadbalancing::distributed::roundRobin( *copiedStorage, 1 );

   if ( rank == 0 )
   {
      WALBERLA_CHECK_EQUAL( copiedStorage->getNumberOfLocalPrimitives(), numGlobalPrimitives );
   }
   else
   {
      WALBERLA_CHECK_EQUAL( copiedStorage->getNumberOfLocalPrimitives(), 0 );
   }

   WALBERLA_LOG_DEVEL("before copy dist")
   WALBERLA_MPI_BARRIER()

   WALBERLA_LOG_DEVEL( migrationInfo );

   loadbalancing::distributed::reverseDistribution( migrationInfo, *copiedStorage );

   WALBERLA_CHECK( storage->getPrimitiveIDs() == copiedStorage->getPrimitiveIDs() );


}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testPrimitiveStorageCopy( hyteg::MeshInfo::fromGmshFile( "../../meshes/3D/cube_24el.msh" ) );

   return EXIT_SUCCESS;
}
