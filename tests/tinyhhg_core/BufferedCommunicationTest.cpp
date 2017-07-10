
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

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
