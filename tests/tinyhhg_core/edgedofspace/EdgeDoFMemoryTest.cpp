
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFDataHandling.hpp"

namespace hhg {

static void testEdgeDoFFunctionMemorySize()
{
  WALBERLA_CHECK_EQUAL( edgeDoFMacroVertexFunctionMemorySize( 2, 4 ), 8 );
  WALBERLA_CHECK_EQUAL( edgeDoFMacroVertexFunctionMemorySize( 3, 4 ), 8 );

  WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 2, 1 ), 15 );
  WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 2, 2 ), 26 );
  WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 3, 1 ), 31 );
  WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 3, 2 ), 54 );

  WALBERLA_CHECK_EQUAL( edgeDoFMacroFaceFunctionMemorySize( 2, 2 ),  30 );
  WALBERLA_CHECK_EQUAL( edgeDoFMacroFaceFunctionMemorySize( 3, 2 ), 108 );
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
