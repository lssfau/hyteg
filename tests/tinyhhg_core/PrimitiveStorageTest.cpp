
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
  std::string meshFileContent = getExampleMeshFileContent();
  std::ofstream file( meshFileName );
  file << meshFileContent;
  file.close();
}

static void testPrimitiveStorage()
{

  std::string meshFileName = "./tmpMeshFile.msh";
  writeTestMeshFile( meshFileName );

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo );

  uint_t balanceRank = 2;
  AllBlocksOnOneRank loadbalancer( 2 );
  setupStorage.balanceLoad( loadbalancer, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ), 0.0 );

  WALBERLA_LOG_INFO( setupStorage );

  PrimitiveStorage storage( balanceRank, setupStorage );

  WALBERLA_LOG_PROGRESS( "Checking that all primitives have been loadbalanced as expected" );

  for ( auto it = storage.beginVertices(); it != storage.endVertices(); it++ )
  {
    WALBERLA_LOG_PROGRESS( "Checking that all primitives have been loadbalanced as expected" );
    WALBERLA_CHECK_EQUAL( it->second->getRank(), balanceRank, "A vertex is not correctly balanced." );
  }

  for ( auto it = storage.beginEdges(); it != storage.endEdges(); it++ )
  {
    WALBERLA_CHECK_EQUAL( it->second->getRank(), balanceRank, "An edge is not correctly balanced." );
  }

  for ( auto it = storage.beginFaces(); it != storage.endFaces(); it++ )
  {
    WALBERLA_CHECK_EQUAL( it->second->getRank(), balanceRank, "A face is not correctly balanced." );
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
