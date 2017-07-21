
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

static void testPrimitiveStorage()
{
  uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  std::string meshFileName = "../../data/meshes/tri_2el.msh";

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  // uint_t balanceRank = 2;
  // AllBlocksOnOneRank loadbalancer( 2 );
  RoundRobin loadbalancer;
  setupStorage.balanceLoad( loadbalancer, 0.0 );

  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

  WALBERLA_LOG_INFO_ON_ROOT( "Building PrimitiveStorage" );

  PrimitiveStorage storage( rank, setupStorage );

  WALBERLA_LOG_PROGRESS_ON_ROOT( "Checking that all primitives have been loadbalanced as expected, checking neighborhood on SetupStorage" );

  for ( auto it = setupStorage.beginVertices(); it != setupStorage.endVertices(); it++ )
  {
    if ( setupStorage.getTargetRank( it->first ) == rank )
    {
      WALBERLA_CHECK( storage.vertexExistsLocally( it->first ) );
    }
    else
    {
      WALBERLA_CHECK( !storage.vertexExistsLocally( it->first ) );
    }

    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
    WALBERLA_CHECK_GREATER( it->second->getNumHigherDimNeighbors(), 0 );
  }

  for ( auto it = setupStorage.beginEdges(); it != setupStorage.endEdges(); it++ )
  {
    if ( setupStorage.getTargetRank( it->first ) == rank )
    {
      WALBERLA_CHECK( storage.edgeExistsLocally( it->first ) );
    }
    else
    {
      WALBERLA_CHECK( !storage.edgeExistsLocally( it->first ) );
    }

    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
    WALBERLA_CHECK_GREATER( it->second->getNumHigherDimNeighbors(), 0 );
  }

  for ( auto it = setupStorage.beginFaces(); it != setupStorage.endFaces(); it++ )
  {
    if ( setupStorage.getTargetRank( it->first ) == rank )
    {
      WALBERLA_CHECK( storage.faceExistsLocally( it->first ) );
    }
    else
    {
      WALBERLA_CHECK( !storage.faceExistsLocally( it->first ) );
    }

    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 3 );
    WALBERLA_CHECK_EQUAL( it->second->getNumHigherDimNeighbors(), 0 );
  }

  WALBERLA_LOG_PROGRESS_ON_ROOT( "Checking neighborhood on distributed storage" );
  for ( auto it = storage.beginVertices(); it != storage.endVertices(); it++ )
  {
	WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
	WALBERLA_CHECK_GREATER( it->second->getNumHigherDimNeighbors(), 0 );
  }
  for ( auto it = storage.beginEdges(); it != storage.endEdges(); it++ )
  {
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
    WALBERLA_CHECK_GREATER( it->second->getNumHigherDimNeighbors(), 0 );
  }
  for ( auto it = storage.beginFaces(); it != storage.endFaces(); it++ )
  {
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 3 );
    WALBERLA_CHECK_EQUAL( it->second->getNumHigherDimNeighbors(), 0 );
  }


  WALBERLA_LOG_PROGRESS_ON_ROOT( "Testing generic getters" );

  std::vector< PrimitiveID > vertexIDs;
  std::vector< PrimitiveID > vertexIDsGeneric;
  storage.getVertexIDs( vertexIDs );
  storage.getPrimitiveIDsGenerically< Vertex >( vertexIDsGeneric );
  WALBERLA_CHECK_EQUAL( vertexIDs.size(), vertexIDsGeneric.size() );

  std::vector< PrimitiveID > edgeIDs;
  std::vector< PrimitiveID > edgeIDsGeneric;
  storage.getEdgeIDs( edgeIDs );
  storage.getPrimitiveIDsGenerically< Edge >( edgeIDsGeneric );
  WALBERLA_CHECK_EQUAL( edgeIDs.size(), edgeIDsGeneric.size() );

  std::vector< PrimitiveID > faceIDs;
  std::vector< PrimitiveID > faceIDsGeneric;
  storage.getFaceIDs( faceIDs );
  storage.getPrimitiveIDsGenerically< Face >( faceIDsGeneric );
  WALBERLA_CHECK_EQUAL( faceIDs.size(), faceIDsGeneric.size() );

  for ( const PrimitiveID & vertexID : vertexIDs )
  {
    WALBERLA_CHECK(  storage.primitiveExistsLocallyGenerically< Primitive >( vertexID ) );
    WALBERLA_CHECK(  storage.primitiveExistsLocallyGenerically< Vertex >( vertexID ) );
    WALBERLA_CHECK( !storage.primitiveExistsLocallyGenerically< Edge >( vertexID ) );
    WALBERLA_CHECK( !storage.primitiveExistsLocallyGenerically< Face >( vertexID ) );
    Vertex * vertex = storage.getPrimitiveGenerically< Vertex >( vertexID );
    WALBERLA_LOG_INFO( "" << vertex->getID().getID() );
  }

  for ( const PrimitiveID & edgeID : edgeIDs )
  {
    WALBERLA_CHECK(  storage.primitiveExistsLocallyGenerically< Primitive >( edgeID ) );
    WALBERLA_CHECK( !storage.primitiveExistsLocallyGenerically< Vertex >( edgeID ) );
    WALBERLA_CHECK(  storage.primitiveExistsLocallyGenerically< Edge >( edgeID ) );
    WALBERLA_CHECK( !storage.primitiveExistsLocallyGenerically< Face >( edgeID ) );
    Edge * edge = storage.getPrimitiveGenerically< Edge >( edgeID );
    WALBERLA_LOG_INFO( "" << edge->getID().getID() );
  }

  for ( const PrimitiveID & faceID : faceIDs )
  {
    WALBERLA_CHECK(  storage.primitiveExistsLocallyGenerically< Primitive >( faceID ) );
    WALBERLA_CHECK( !storage.primitiveExistsLocallyGenerically< Vertex >( faceID ) );
    WALBERLA_CHECK( !storage.primitiveExistsLocallyGenerically< Edge >( faceID ) );
    WALBERLA_CHECK(  storage.primitiveExistsLocallyGenerically< Face >( faceID ) );
    Face * face = storage.getPrimitiveGenerically< Face >( faceID );
    WALBERLA_LOG_INFO( "" << face->getID().getID() );
  }


}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testPrimitiveStorage();

   return EXIT_SUCCESS;
}
