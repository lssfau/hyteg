
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFDataHandling.hpp"

namespace hhg {

static void testEdgeDoFFunctionMemorySize()
{
  WALBERLA_CHECK_EQUAL( EdgeDoFMacroVertexFunctionMemorySize( 2, 4 ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFMacroVertexFunctionMemorySize( 3, 4 ), 8 );

  WALBERLA_CHECK_EQUAL( EdgeDoFMacroEdgeFunctionMemorySize( 2, 1 ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFMacroEdgeFunctionMemorySize( 2, 2 ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFMacroEdgeFunctionMemorySize( 3, 1 ), 31 );
  WALBERLA_CHECK_EQUAL( EdgeDoFMacroEdgeFunctionMemorySize( 3, 2 ), 54 );

  WALBERLA_CHECK_EQUAL( EdgeDoFMacroFaceFunctionMemorySize( 2, 2 ),  30 );
  WALBERLA_CHECK_EQUAL( EdgeDoFMacroFaceFunctionMemorySize( 3, 2 ), 108 );
}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testEdgeDoFFunctionMemorySize();

   return EXIT_SUCCESS;
}
