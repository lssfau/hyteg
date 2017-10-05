
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOnMacroFaceIndex.hpp"

namespace hhg {

static void testEdgeDoFsOnMacroFace()
{
  typedef stencilDirection sD;

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 0, 0, sD::EDGE_HO_E ),  0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 1, 0, sD::EDGE_HO_E ),  1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 0, sD::EDGE_HO_E ),  2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 1, sD::EDGE_HO_E ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 2, sD::EDGE_HO_E ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 0, 7, sD::EDGE_HO_E ), 35 );
  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing: EDGE_HO_E correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 1, 0, sD::EDGE_HO_W ),  0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 0, sD::EDGE_HO_W ),  1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 1, sD::EDGE_HO_W ),  9 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 2, sD::EDGE_HO_W ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 1, 7, sD::EDGE_HO_W ), 35 );
  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing: EDGE_HO_W correct!" );

#ifdef NDEBUG
  static_assert( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 0, 0, sD::EDGE_HO_E )  ==  0 , "EDGE_HO_E  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 1, sD::EDGE_HO_E )  == 10 , "EDGE_HO_E  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVertex< 3 >( 2, 2, sD::EDGE_HO_E )  == 17 , "EDGE_HO_E  cannot be statically computed by the compiler!" );
#endif

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testEdgeDoFsOnMacroFace();

   return EXIT_SUCCESS;
}
