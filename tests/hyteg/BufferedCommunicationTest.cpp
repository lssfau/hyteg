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

#include <hyteg/communication/BufferedCommunication.hpp>
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/communication/PackInfo.hpp"


namespace hyteg {

struct VertexTestData
{
  PrimitiveID ownID;
  std::vector< PrimitiveID > edgeIDs;
};

struct EdgeTestData
{
  PrimitiveID ownID;
  std::vector< PrimitiveID > vertexIDs;
};

struct VertexTestDataHandling : OnlyInitializeDataHandling< VertexTestData, Vertex >
{
  virtual std::shared_ptr< VertexTestData > initialize( const Vertex * const primitive ) const
  {
    auto data = std::make_shared< VertexTestData >();
    data->ownID = primitive->getID();
    return data;
  }
};

struct EdgeTestDataHandling : OnlyInitializeDataHandling< EdgeTestData, Edge >
{
  virtual std::shared_ptr< EdgeTestData > initialize( const Edge * const primitive ) const
  {
    auto data = std::make_shared< EdgeTestData >();
    data->ownID = primitive->getID();
    return data;
  }
};


class TestPackInfo : public communication::PackInfo
{
public:

  TestPackInfo( PrimitiveDataID< VertexTestData, Vertex > & vertexDataID,
		PrimitiveDataID< EdgeTestData,   Edge >   & edgeDataID ) :
		  vertexDataID_( vertexDataID ), edgeDataID_( edgeDataID )
  {}

  virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
  {
    WALBERLA_UNUSED( receiver );
    VertexTestData * data = sender->getData( vertexDataID_ );
    buffer << data->ownID;
    // WALBERLA_LOG_INFO( "Packing | Vertex: " << sender->getID().getID() << ", Receiver: " << receiver.getID() << ", Data: " << data->someInt );
  }

  virtual void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
  {
    WALBERLA_UNUSED( sender );

    EdgeTestData * data = receiver->getData( edgeDataID_ );
    PrimitiveID vertexData;
    buffer >> vertexData;
    // WALBERLA_LOG_INFO( "Unpacking | Edge: " << receiver->getID().getID() << ", Data: " << vertexData );
    data->vertexIDs.push_back( vertexData );
  }

  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const
  {
    VertexTestData * vertexData = sender->getData( vertexDataID_ );
    EdgeTestData   * edgeData   = receiver->getData( edgeDataID_ );
    // WALBERLA_LOG_INFO( "Direct | Vertex: " << sender->getID().getID() << ", Edge: " << receiver->getID().getID() << ", Data: " << vertexData->someInt );
    edgeData->vertexIDs.push_back( vertexData->ownID );
  }


  virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
  {
    WALBERLA_UNUSED( receiver );
    EdgeTestData * data = sender->getData( edgeDataID_ );
    buffer << data->ownID;
  }

  virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
  {
    WALBERLA_UNUSED( sender );

    VertexTestData * data = receiver->getData( vertexDataID_ );
    PrimitiveID edgeData;
    buffer >> edgeData;
    data->edgeIDs.push_back( edgeData );
  }

  virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const
  {
    EdgeTestData   * edgeData   = sender->getData( edgeDataID_ );
    VertexTestData * vertexData = receiver->getData( vertexDataID_ );
    vertexData->edgeIDs.push_back( edgeData->ownID );
  }


  virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
  }


  virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
  }

private:

  PrimitiveDataID< VertexTestData, Vertex > vertexDataID_;
  PrimitiveDataID< EdgeTestData,   Edge >   edgeDataID_;

};

static void testBufferedCommunication()
{

  std::string meshFileName = "../../data/meshes/bfs_126el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  loadbalancing::roundRobin( setupStorage );

  // WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
  WALBERLA_MPI_BARRIER();

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

  communication::BufferedCommunicator communicator( storage );

  auto vertexTestDataHandling = std::make_shared< VertexTestDataHandling >();
  auto edgeTestDataHandling = std::make_shared< EdgeTestDataHandling >();

  PrimitiveDataID< VertexTestData, Vertex > vertexTestDataID;
  storage->addVertexData( vertexTestDataID, vertexTestDataHandling, "vertex data" );

  PrimitiveDataID< EdgeTestData, Edge > edgeTestDataID;
  storage->addEdgeData  ( edgeTestDataID, edgeTestDataHandling,   "edge data" );

  std::shared_ptr< TestPackInfo > testPackInfo( new TestPackInfo( vertexTestDataID, edgeTestDataID ) );

  for ( const auto & it : storage->getVertices() )
  {
    auto vertex = it.second;
    auto data = vertex->getData( vertexTestDataID );
    WALBERLA_CHECK_EQUAL( data->ownID, vertex->getID() );
  }

  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    auto data = edge->getData( edgeTestDataID );
    WALBERLA_CHECK_EQUAL( data->vertexIDs.size(), 0 );
  }

  communicator.addPackInfo( testPackInfo );
  communicator.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  communicator.enableTiming( timingTree );

  communicator.startCommunication< Vertex, Edge >();
  communicator.endCommunication< Vertex, Edge >();
#if 0
  communicator.startCommunication< Vertex, Edge >();
  communicator.endCommunication< Vertex, Edge >();
#endif

  communicator.startCommunication< Edge, Vertex >();

  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    auto data = edge->getData( edgeTestDataID );
    WALBERLA_CHECK_EQUAL( data->vertexIDs.size(), 2 );
    WALBERLA_CHECK_UNEQUAL( data->vertexIDs[0], data->vertexIDs[1], "Failing on Edge: " << it.first );

    for ( const auto & lowerDimNeighborID : edge->getLowerDimNeighbors() )
    {
      WALBERLA_CHECK( data->vertexIDs[0] == lowerDimNeighborID || data->vertexIDs[1] == lowerDimNeighborID, "Failing on Edge: " << lowerDimNeighborID );
    }

    // WALBERLA_LOG_INFO( "Edge " << edge->getID().getID() << " received: " << data->someInts[0] << ", " << data->someInts[1] );
  }

  communicator.endCommunication< Edge, Vertex >(); // checking interleaved communication and processing

  for ( const auto & it : storage->getVertices() )
  {
    auto vertex = it.second;
    auto data = vertex->getData( vertexTestDataID );
    WALBERLA_CHECK_GREATER( data->edgeIDs.size(), 0 );
    std::set< PrimitiveID > edgeIdsSet( data->edgeIDs.begin(), data->edgeIDs.end() );
    WALBERLA_CHECK_EQUAL( data->edgeIDs.size(), edgeIdsSet.size() );

    for ( const auto & higherDimNeighborID : vertex->getHigherDimNeighbors() )
    {
      WALBERLA_ABORT( "Here a primitive ID must be casted - check and fix test." );
      // TODO: comment in this line when fixed!
      // WALBERLA_CHECK_EQUAL( edgeIdsSet.count( uint_c( higherDimNeighborID.getID() ) ), 1 );
    }
  }

  WALBERLA_MPI_SECTION()
  {
    const walberla::WcTimingTree tt = timingTree->getReduced();
    WALBERLA_LOG_INFO_ON_ROOT( tt );
  }
  WALBERLA_NON_MPI_SECTION()
  {
    WALBERLA_LOG_INFO_ON_ROOT( *timingTree );
  }

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testBufferedCommunication();

   return EXIT_SUCCESS;
}
