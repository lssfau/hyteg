
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

namespace hhg {

static void testEdgeDoFFunctionMemorySize()
{
  using namespace edgedof;

  auto storage = PrimitiveStorage::createFromGmshFile( "../../data/meshes/quad_8el.msh" );

  for ( const auto & it : storage->getEdges() )
  {
    const Edge & edge = *it.second;

    if ( edge.getNumNeighborFaces() == 1 )
    {
      WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 2, edge ), 15 );
      WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 3, edge ), 31 );
    }
    else if ( edge.getNumNeighborFaces() == 2 )
    {
      WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 2, edge ), 26 );
      WALBERLA_CHECK_EQUAL( edgeDoFMacroEdgeFunctionMemorySize( 3, edge ), 54 );
    }
  }

  for ( const auto & it : storage->getFaces() )
  {
    const Face & face = *it.second;

    WALBERLA_CHECK_EQUAL( edgeDoFMacroFaceFunctionMemorySize( 2, face ),  30 );
    WALBERLA_CHECK_EQUAL( edgeDoFMacroFaceFunctionMemorySize( 3, face ), 108 );
  }

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
