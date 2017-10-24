
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

static void testCommonIndexing()
{
  using namespace indexing;

  WALBERLA_LOG_INFO_ON_ROOT( "Index P1      - face, level 3, (3, 3, center): " << P1FaceIndexFromVertex< 3 >( 3, 3, VERTEX_C ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Index EdgeDoF - face, level 3, (3, 3, center): " << EdgeDoFFaceIndexFromVertex< 3 >( 3, 3, EDGE_HO_C ) );

  for ( const auto & it : P1FaceBorderIterator< 3 >( Direction::DIAGONAL_BOTTOM_TO_TOP, 1 ) )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it[0] << ", row = " << it[1] << " ( idx = " << P1FaceIndexFromVertex< 3 >( it[0], it[1], VERTEX_C ) << " ) " );
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
