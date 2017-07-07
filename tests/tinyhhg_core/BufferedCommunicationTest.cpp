
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

static void testBufferedCommunication()
{
  communication::BufferedCommunicator communicator;
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
