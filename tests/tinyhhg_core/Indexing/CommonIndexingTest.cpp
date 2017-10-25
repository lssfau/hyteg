
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"

#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/VertexDoFIndexing.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

namespace hhg {

static void testCommonIndexing()
{
  using namespace indexing;

  WALBERLA_LOG_INFO_ON_ROOT( "Index P1      - face, level 3, (3, 3, center): " << VertexDoFOnMacroFaceIndexFromVertex< 3 >( 3, 3, stencilDirection::VERTEX_C ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Index EdgeDoF - face, level 3, (3, 3, center): " << EdgeDoFFaceIndexFromVertex< 3 >( 3, 3, stencilDirection::EDGE_HO_C ) );

  for ( const auto & it : VertexDoFFaceBorderIterator< 3 >( FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 1 ) )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it.col() << ", row = " << it.row() << " ( idx = " << VertexDoFOnMacroFaceIndexFromVertex< 3 >( it.col(), it.row(), stencilDirection::VERTEX_C ) << " ) " );
  }


}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testCommonIndexing();

   return EXIT_SUCCESS;
}
