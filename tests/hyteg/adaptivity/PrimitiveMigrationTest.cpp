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
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/indexing/Common.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::uint_t;

namespace hyteg {

struct Data
{
   uint_t primitiveIDInData = 345678;
};

class DataHandling : PrimitiveDataHandling< Data, Primitive >
{
 public:
   std::shared_ptr< Data > initialize( const Primitive* const primitive ) const
   {
      auto data               = std::make_shared< Data >();
      data->primitiveIDInData = uint_c( primitive->getID() );
      return data;
   }

   virtual void
       serialize( const Primitive* const primitive, const PrimitiveDataID< Data, Primitive >& id, SendBuffer& buffer ) const
   {
      auto data = primitive->getData( id );
      buffer << data->primitiveIDInData;
   }

   virtual void
       deserialize( const Primitive* const primitive, const PrimitiveDataID< Data, Primitive >& id, RecvBuffer& buffer ) const
   {
      auto data = primitive->getData( id );
      buffer >> data->primitiveIDInData;
   }
};

static void testPrimitiveMigration()
{
   uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   const std::string meshFileName = "../../data/meshes/3D/cube_24el.msh";

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   loadbalancing::roundRobin( setupStorage );

   //WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

   WALBERLA_LOG_INFO( "Building PrimitiveStorage" );

   std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

   writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition_before_migration" );

   PrimitiveDataID< Data, Primitive > dataID;
   auto                               dataHandling = std::make_shared< DataHandling >();
   storage->addPrimitiveData( dataID, dataHandling, "test data" );

   WALBERLA_MPI_SECTION()
   {
      std::vector< PrimitiveID > primitiveIDs;
      storage->getPrimitiveIDsGenerically< Primitive >( primitiveIDs );

      std::map< PrimitiveID, uint_t > migrationMap;
      uint_t                                  lel = 0;
      for( const auto& id : primitiveIDs )
      {
         uint_t targetRank = ++lel % numProcesses;
         // WALBERLA_LOG_INFO( "Migrating " << id.getID() << " to rank " << targetRank );
         migrationMap[id] = targetRank;
      }

      const auto numReceivingPrimitives = getNumReceivingPrimitives( migrationMap );
      const MigrationInfo migrationInfo( migrationMap, numReceivingPrimitives );
      WALBERLA_LOG_INFO( migrationInfo )
      storage->migratePrimitives( migrationInfo );

      PrimitiveStorage::PrimitiveMap primitives;
      storage->getPrimitives( primitives );
      for( const auto& it : primitives )
      {
         WALBERLA_CHECK_EQUAL( it.first, it.second->getData( dataID )->primitiveIDInData );
      }
   }

   writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition_after_migration" );
}

static void testPrimitiveMigrationMaps()
{
  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  const std::string meshFileName = "../../data/meshes/3D/cube_24el.msh";

  MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  loadbalancing::roundRobin( setupStorage );

  WALBERLA_LOG_INFO( "Building PrimitiveStorage" );

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

  writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition_before_migration" );

  PrimitiveDataID< LevelWiseMemory< std::map< indexing::IndexIncrement, uint_t > >, Cell > dataID;
  auto cellDataHandling =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< std::map< indexing::IndexIncrement, uint_t > >, Cell > >(2, 4 );
  storage->addCellData( dataID, cellDataHandling, "test data" );

  for ( auto it : storage->getCells() )
  {
     auto & cellData = it.second->getData( dataID )->getData( 2 );
     cellData[ indexing::IndexIncrement( {0, 2, 3} ) ] = 3;
     cellData[ indexing::IndexIncrement( {0, 2, 5} ) ] = 7;
  }

  for ( auto it : storage->getCells() )
  {
    auto cellData = it.second->getData( dataID )->getData( 2 );
    WALBERLA_LOG_INFO( "Cell " << it.first << "(0, 2, 3): " << cellData.at( indexing::IndexIncrement( {0, 2, 3} ) ) << ", (0, 2, 5): " << cellData.at( indexing::IndexIncrement( {0, 2, 5} ) )  )
  }

  WALBERLA_MPI_BARRIER()
  WALBERLA_LOG_INFO( "Migration..." )

  WALBERLA_MPI_SECTION()
  {
    std::vector< PrimitiveID > primitiveIDs;
    storage->getPrimitiveIDsGenerically< Primitive >( primitiveIDs );

    std::map< PrimitiveID, uint_t > migrationMap;
    uint_t                                  lel = 0;
    for( const auto& id : primitiveIDs )
    {
      uint_t targetRank = ++lel % numProcesses;
      // WALBERLA_LOG_INFO( "Migrating " << id.getID() << " to rank " << targetRank );
      migrationMap[id] = targetRank;
    }

    const auto numReceivingPrimitives = getNumReceivingPrimitives( migrationMap );
    const MigrationInfo migrationInfo( migrationMap, numReceivingPrimitives );
    WALBERLA_LOG_INFO( migrationInfo );
    storage->migratePrimitives( migrationInfo );

    for ( auto it : storage->getCells() )
    {
      auto cellData = it.second->getData( dataID )->getData( 2 );
      WALBERLA_CHECK_EQUAL( cellData.at( indexing::IndexIncrement( {0, 2, 3} ) ), 3 );
      WALBERLA_CHECK_EQUAL( cellData.at( indexing::IndexIncrement( {0, 2, 5} ) ), 7 );
      WALBERLA_LOG_INFO( "Cell " << it.first << "(0, 2, 3): " << cellData.at( indexing::IndexIncrement( {0, 2, 3} ) ) << ", (0, 2, 5): " << cellData.at( indexing::IndexIncrement( {0, 2, 5} ) )  )
    }
  }

  writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition_after_migration" );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::debug::enterTestMode();
   hyteg::testPrimitiveMigration();
   hyteg::testPrimitiveMigrationMaps();

   return EXIT_SUCCESS;
}
