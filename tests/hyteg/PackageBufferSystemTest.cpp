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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/communication/PackageBufferSystem.hpp"

namespace hyteg {

static void testAllToAll()
{
   WALBERLA_LOG_INFO_ON_ROOT( "All to all" )
   WALBERLA_MPI_BARRIER()

   const auto comm         = walberla::mpi::MPIManager::instance()->comm();
   const auto rank         = walberla::mpi::MPIManager::instance()->rank();
   const auto numProcesses = walberla::mpi::MPIManager::instance()->numProcesses();

   const auto numPackages = 3;

   PackageBufferSystem bs( comm, uint_c( numPackages * numProcesses ) );

   for ( uint_t targetProcess = 0; targetProcess < uint_c( numProcesses ); targetProcess++ )
   {
      for ( uint_t package = 0; package < numPackages; package++ )
      {
         auto& buffer = bs.getPackageSendBuffer( targetProcess );
         buffer << rank << package;
      }
   }

   bs.sendAll();

   std::vector< int > packageSum( uint_c( numProcesses ), 0 );

   while ( !bs.allPackagesReceived() )
   {
      auto   p = bs.getNextPackage();
      int    senderRank;
      uint_t content;
      p.buffer() >> senderRank;
      p.buffer() >> content;
      packageSum[uint_c( senderRank )] += int_c( content );
   }

   for ( auto ps : packageSum )
   {
      WALBERLA_CHECK_EQUAL( ps, ( ( numPackages - 1 ) * ( numPackages - 1 ) + ( numPackages - 1 ) ) / 2 );
   }
}

static void testAllToLower()
{
   WALBERLA_LOG_INFO_ON_ROOT( "All to lower" )
   WALBERLA_MPI_BARRIER()

   const auto comm         = walberla::mpi::MPIManager::instance()->comm();
   const auto rank         = walberla::mpi::MPIManager::instance()->rank();
   const auto numProcesses = walberla::mpi::MPIManager::instance()->numProcesses();

   const auto numPackages = 3;

   PackageBufferSystem bs( comm, uint_c( numPackages * ( numProcesses - rank ) ) );

   for ( uint_t targetProcess = 0; targetProcess < uint_c( rank + 1 ); targetProcess++ )
   {
      for ( uint_t package = 0; package < numPackages; package++ )
      {
         auto& buffer = bs.getPackageSendBuffer( targetProcess );
         buffer << rank << package;
      }
   }

   bs.sendAll();

   std::vector< int > packageSum( uint_c( numProcesses ), 0 );

   while ( !bs.allPackagesReceived() )
   {
      auto   p = bs.getNextPackage();
      int    senderRank;
      uint_t content;
      p.buffer() >> senderRank;
      p.buffer() >> content;
      WALBERLA_CHECK_GREATER_EQUAL( senderRank, rank );
      packageSum[uint_c( senderRank )] += int_c( content );
   }

   for ( int sr = 0; sr < numProcesses; sr++ )
   {
      WALBERLA_CHECK_EQUAL( packageSum[uint_c( sr )],
                            ( ( ( numPackages - 1 ) * ( numPackages - 1 ) + ( numPackages - 1 ) ) / 2 ) *
                                ( rank <= sr ? 1 : 0 ) );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testAllToAll();
   hyteg::testAllToLower();
   return EXIT_SUCCESS;
}
