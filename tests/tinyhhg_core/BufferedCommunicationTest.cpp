
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

struct VertexTestData
{
  uint_t someInt;
};

struct EdgeTestData
{
  std::vector< uint_t > someInts;
};

struct VertexTestDataHandling : OnlyInitializeDataHandling< VertexTestData, Vertex >
{
  virtual VertexTestData * initialize( const Vertex * const primitive ) const
  {
    VertexTestData * data = new VertexTestData;
    data->someInt = primitive->getID().getID();
    return data;
  }
};

struct EdgeTestDataHandling : OnlyInitializeDataHandling< EdgeTestData, Edge >
{
  virtual EdgeTestData * initialize( const Edge * const ) const
  {
    EdgeTestData * data = new EdgeTestData;
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

  virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer)
  {
    WALBERLA_UNUSED( receiver );

    VertexTestData * data = sender->getData( vertexDataID_ );
    buffer << data->someInt;
  }

  virtual void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer)
  {
    WALBERLA_UNUSED( sender );

    EdgeTestData * data = receiver->getData( edgeDataID_ );
    uint_t vertexData;
    buffer >> vertexData;
    data->someInts.push_back( vertexData );
  }

  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
  }


  virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
  }


  virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
  }


  virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer)
  {
    WALBERLA_UNUSED( sender   );
    WALBERLA_UNUSED( receiver );
    WALBERLA_UNUSED( buffer   );
  }

  virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver)
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

  uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  std::string meshFileName = "../../data/meshes/tri_2el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  RoundRobin loadbalancer;
  setupStorage.balanceLoad( loadbalancer, 0.0 );

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( rank, setupStorage ) );

  communication::BufferedCommunicator communicator( storage );

  VertexTestDataHandling vertexTestDataHandling;
  EdgeTestDataHandling   edgeTestDataHandling;

  auto vertexTestDataID = storage->addVertexData( vertexTestDataHandling, "vertex data" );
  auto edgeTestDataID   = storage->addEdgeData  ( edgeTestDataHandling,   "edge data" );

  std::shared_ptr< TestPackInfo > testPackInfo( new TestPackInfo( vertexTestDataID, edgeTestDataID ) );

  for ( auto it = storage->beginVertices(); it != storage->endVertices(); it++ )
  {
    Vertex * vertex = it->second;
    VertexTestData * data = vertex->getData( vertexTestDataID );
    WALBERLA_CHECK_EQUAL( data->someInt, vertex->getID().getID() );
  }

  for ( auto it = storage->beginEdges(); it != storage->endEdges(); it++ )
  {
    Edge * edge = it->second;
    EdgeTestData * data = edge->getData( edgeTestDataID );
    WALBERLA_CHECK_EQUAL( data->someInts.size(), 0 );
  }

  communicator.addPackInfo( testPackInfo );

  communicator.startCommunicationVertexToEdge();
  communicator.endCommunicationVertexToEdge();

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testBufferedCommunication();

   return EXIT_SUCCESS;
}
