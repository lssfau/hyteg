
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

std::string getExampleMeshFileContent()
{
  std::string content =
  "$MeshFormat\n"
  "2.2 0 8\n"
  "$EndMeshFormat\n"
  "$Nodes\n"
  "4\n"
  "1 0 0 0\n"
  "2 1 0 0\n"
  "3 0 1 0\n"
  "4 1 1 0\n"
  "$EndNodes\n"
  "$Elements\n"
  "9\n"
  "1 15 2 1 1 1\n"
  "2 15 2 1 2 2\n"
  "3 15 2 1 3 3\n"
  "4 1 2 1 1 1 2\n"
  "5 1 2 1 2 2 4\n"
  "6 1 2 1 2 4 3\n"
  "7 1 2 1 3 3 1\n"
  "8 2 2 0 5 2 4 1\n"
  "9 2 2 0 5 1 4 3\n"
  "$EndElements\n";

  return content;
}

void writeTestMeshFile( const std::string & meshFileName )
{
  WALBERLA_ROOT_SECTION()
  {
    std::string meshFileContent = getExampleMeshFileContent();
    std::ofstream file( meshFileName );
    file << meshFileContent;
    file.close();
  }
}

static void testPrimitiveStorage()
{
  uint_t rank = walberla::mpi::MPIManager::instance()->rank();

  std::string meshFileName = "./tmpMeshFile.msh";
  writeTestMeshFile( meshFileName );

  WALBERLA_MPI_BARRIER();

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  // uint_t balanceRank = 2;
  // AllBlocksOnOneRank loadbalancer( 2 );
  RoundRobin loadbalancer;
  setupStorage.balanceLoad( loadbalancer, 0.0 );

  WALBERLA_LOG_INFO( setupStorage );

  WALBERLA_MPI_BARRIER();
  WALBERLA_LOG_INFO( "Building PrimitiveStorage" );

  PrimitiveStorage storage( rank, setupStorage );

  WALBERLA_LOG_PROGRESS( "Checking that all primitives have been loadbalanced as expected" );

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
