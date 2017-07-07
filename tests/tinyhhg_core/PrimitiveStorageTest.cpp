
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
    if ( it->second->getTargetRank() == rank )
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
    if ( it->second->getTargetRank() == rank )
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
    if ( it->second->getTargetRank() == rank )
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
