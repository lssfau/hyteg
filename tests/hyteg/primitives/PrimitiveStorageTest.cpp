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
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

static void testPrimitiveStorage()
{
   uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   const std::string meshFileName = prependHyTeGMeshDir( "2D/porous_fine.msh" );
   // const std::string meshFileName = prependHyTeGMeshDir( "2D/bfs_126el.msh");
   // const std::string meshFileName = prependHyTeGMeshDir( "2D/tri_2el.msh");
   const std::string distributionFile = "../../output/PrimitiveStorageTestDistribution.csv";

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   WALBERLA_LOG_INFO_ON_ROOT( "LB start" );
   loadbalancing::greedy( setupStorage );
   WALBERLA_LOG_INFO_ON_ROOT( "LB end" );

   WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

   WALBERLA_LOG_INFO_ON_ROOT( "Building PrimitiveStorage" );

   std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

   WALBERLA_LOG_PROGRESS_ON_ROOT(
       "Checking that all primitives have been loadbalanced as expected, checking neighborhood on SetupStorage" );

   for ( auto it : setupStorage.getVertices() )
   {
      if ( setupStorage.getTargetRank( it.first ) == rank )
      {
         WALBERLA_CHECK( storage->vertexExistsLocally( it.first ) );
      }
      else
      {
         WALBERLA_CHECK( !storage->vertexExistsLocally( it.first ) );
      }

      WALBERLA_CHECK_EQUAL( it.second->getNumLowerDimNeighbors(), 0 );
      WALBERLA_CHECK_GREATER( it.second->getNumHigherDimNeighbors(), 0 );
   }

   for ( auto it : setupStorage.getEdges() )
   {
      if ( setupStorage.getTargetRank( it.first ) == rank )
      {
         WALBERLA_CHECK( storage->edgeExistsLocally( it.first ) );
      }
      else
      {
         WALBERLA_CHECK( !storage->edgeExistsLocally( it.first ) );
      }

      WALBERLA_CHECK_EQUAL( it.second->getNumLowerDimNeighbors(), 2 );
      WALBERLA_CHECK_GREATER( it.second->getNumHigherDimNeighbors(), 0 );
   }

   for ( auto it : setupStorage.getFaces() )
   {
      if ( setupStorage.getTargetRank( it.first ) == rank )
      {
         WALBERLA_CHECK( storage->faceExistsLocally( it.first ) );
      }
      else
      {
         WALBERLA_CHECK( !storage->faceExistsLocally( it.first ) );
      }

      WALBERLA_CHECK_EQUAL( it.second->getNumLowerDimNeighbors(), 3 );
      WALBERLA_CHECK_EQUAL( it.second->getNumHigherDimNeighbors(), 0 );
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Checking neighborhood on distributed storage" );
   for ( const auto& it : storage->getVertices() )
   {
      WALBERLA_CHECK_EQUAL( it.second->getNumLowerDimNeighbors(), 0 );
      WALBERLA_CHECK_GREATER( it.second->getNumHigherDimNeighbors(), 0 );
   }
   for ( const auto& it : storage->getEdges() )
   {
      WALBERLA_CHECK_EQUAL( it.second->getNumLowerDimNeighbors(), 2 );
      WALBERLA_CHECK_GREATER( it.second->getNumHigherDimNeighbors(), 0 );
   }
   for ( const auto& it : storage->getFaces() )
   {
      WALBERLA_CHECK_EQUAL( it.second->getNumLowerDimNeighbors(), 3 );
      WALBERLA_CHECK_EQUAL( it.second->getNumHigherDimNeighbors(), 0 );
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Testing generic getters" );

   std::vector< PrimitiveID > vertexIDs;
   std::vector< PrimitiveID > vertexIDsGeneric;
   storage->getVertexIDs( vertexIDs );
   storage->getPrimitiveIDsGenerically< Vertex >( vertexIDsGeneric );
   WALBERLA_CHECK_EQUAL( vertexIDs.size(), vertexIDsGeneric.size() );

   std::vector< PrimitiveID > edgeIDs;
   std::vector< PrimitiveID > edgeIDsGeneric;
   storage->getEdgeIDs( edgeIDs );
   storage->getPrimitiveIDsGenerically< Edge >( edgeIDsGeneric );
   WALBERLA_CHECK_EQUAL( edgeIDs.size(), edgeIDsGeneric.size() );

   std::vector< PrimitiveID > faceIDs;
   std::vector< PrimitiveID > faceIDsGeneric;
   storage->getFaceIDs( faceIDs );
   storage->getPrimitiveIDsGenerically< Face >( faceIDsGeneric );
   WALBERLA_CHECK_EQUAL( faceIDs.size(), faceIDsGeneric.size() );

   for ( const PrimitiveID& vertexID : vertexIDs )
   {
      WALBERLA_CHECK( storage->primitiveExistsLocallyGenerically< Primitive >( vertexID ) );
      WALBERLA_CHECK( storage->primitiveExistsLocallyGenerically< Vertex >( vertexID ) );
      WALBERLA_CHECK( !storage->primitiveExistsLocallyGenerically< Edge >( vertexID ) );
      WALBERLA_CHECK( !storage->primitiveExistsLocallyGenerically< Face >( vertexID ) );
      Vertex* vertex = storage->getPrimitiveGenerically< Vertex >( vertexID );
      WALBERLA_LOG_INFO( "" << vertex->getID() );
   }

   for ( const PrimitiveID& edgeID : edgeIDs )
   {
      WALBERLA_CHECK( storage->primitiveExistsLocallyGenerically< Primitive >( edgeID ) );
      WALBERLA_CHECK( !storage->primitiveExistsLocallyGenerically< Vertex >( edgeID ) );
      WALBERLA_CHECK( storage->primitiveExistsLocallyGenerically< Edge >( edgeID ) );
      WALBERLA_CHECK( !storage->primitiveExistsLocallyGenerically< Face >( edgeID ) );
      Edge* edge = storage->getPrimitiveGenerically< Edge >( edgeID );
      WALBERLA_LOG_INFO( "" << edge->getID() );
   }

   for ( const PrimitiveID& faceID : faceIDs )
   {
      WALBERLA_CHECK( storage->primitiveExistsLocallyGenerically< Primitive >( faceID ) );
      WALBERLA_CHECK( !storage->primitiveExistsLocallyGenerically< Vertex >( faceID ) );
      WALBERLA_CHECK( !storage->primitiveExistsLocallyGenerically< Edge >( faceID ) );
      WALBERLA_CHECK( storage->primitiveExistsLocallyGenerically< Face >( faceID ) );
      Face* face = storage->getPrimitiveGenerically< Face >( faceID );
      WALBERLA_LOG_INFO( "" << face->getID() );
   }

   writePrimitiveStorageDistributionCSV( storage, distributionFile );
   writeDomainPartitioningVTK( storage, "../../output/", "domain_decomposition" );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testPrimitiveStorage();

   return EXIT_SUCCESS;
}
