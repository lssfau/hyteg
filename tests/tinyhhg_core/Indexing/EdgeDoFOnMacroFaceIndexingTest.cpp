
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOnMacroFaceIndex.hpp"

namespace hhg {

static void testEdgeDoFsOnMacroFace()
{
  typedef stencilDirection sD;

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 0, sD::EDGE_HO_E ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 0, sD::EDGE_HO_E ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 0, sD::EDGE_HO_E ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 0, sD::EDGE_HO_E ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 1, sD::EDGE_HO_E ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_HO_E ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_HO_E ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 2, sD::EDGE_HO_E ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_HO_E ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 3, sD::EDGE_HO_E ), 9 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 0, sD::EDGE_HO_E ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 0, sD::EDGE_HO_E ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 0, sD::EDGE_HO_E ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 0, sD::EDGE_HO_E ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 0, sD::EDGE_HO_E ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 0, sD::EDGE_HO_E ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 0, sD::EDGE_HO_E ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 0, sD::EDGE_HO_E ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 1, sD::EDGE_HO_E ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_HO_E ), 9 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_HO_E ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_HO_E ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_HO_E ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_HO_E ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_HO_E ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 2, sD::EDGE_HO_E ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_HO_E ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_HO_E ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_HO_E ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_HO_E ), 19 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_HO_E ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 3, sD::EDGE_HO_E ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_HO_E ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_HO_E ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_HO_E ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_HO_E ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 4, sD::EDGE_HO_E ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_HO_E ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_HO_E ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_HO_E ), 29 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_HO_E ), 30 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_HO_E ), 31 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_HO_E ), 32 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 6, sD::EDGE_HO_E ), 33 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_HO_E ), 34 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 7, sD::EDGE_HO_E ), 35 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing: EDGE_HO_E correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 0, sD::EDGE_HO_W ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 0, sD::EDGE_HO_W ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 0, sD::EDGE_HO_W ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 4, 0, sD::EDGE_HO_W ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_HO_W ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_HO_W ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 1, sD::EDGE_HO_W ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_HO_W ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 2, sD::EDGE_HO_W ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 3, sD::EDGE_HO_W ), 9 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 0, sD::EDGE_HO_W ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 0, sD::EDGE_HO_W ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 0, sD::EDGE_HO_W ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 0, sD::EDGE_HO_W ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 0, sD::EDGE_HO_W ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 0, sD::EDGE_HO_W ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 0, sD::EDGE_HO_W ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 8, 0, sD::EDGE_HO_W ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_HO_W ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_HO_W ), 9 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_HO_W ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_HO_W ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_HO_W ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_HO_W ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 1, sD::EDGE_HO_W ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_HO_W ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_HO_W ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_HO_W ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_HO_W ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_HO_W ), 19 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 2, sD::EDGE_HO_W ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_HO_W ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_HO_W ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_HO_W ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_HO_W ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 3, sD::EDGE_HO_W ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_HO_W ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_HO_W ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_HO_W ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 4, sD::EDGE_HO_W ), 29 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_HO_W ), 30 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_HO_W ), 31 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 5, sD::EDGE_HO_W ), 32 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_HO_W ), 33 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 6, sD::EDGE_HO_W ), 34 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 7, sD::EDGE_HO_W ), 35 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing: EDGE_HO_W correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 0, sD::EDGE_HO_NW ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 0, sD::EDGE_HO_NW ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 0, sD::EDGE_HO_NW ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_HO_NW ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_HO_NW ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_HO_NW ), 9 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 0, sD::EDGE_HO_NW ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 0, sD::EDGE_HO_NW ), 9 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 0, sD::EDGE_HO_NW ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 0, sD::EDGE_HO_NW ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 0, sD::EDGE_HO_NW ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 0, sD::EDGE_HO_NW ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 0, sD::EDGE_HO_NW ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_HO_NW ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_HO_NW ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_HO_NW ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_HO_NW ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_HO_NW ), 19 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_HO_NW ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_HO_NW ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_HO_NW ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_HO_NW ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_HO_NW ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_HO_NW ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_HO_NW ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_HO_NW ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_HO_NW ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_HO_NW ), 29 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_HO_NW ), 30 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_HO_NW ), 31 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_HO_NW ), 32 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_HO_NW ), 33 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_HO_NW ), 34 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_HO_NW ), 35 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing: EDGE_HO_NW correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 1, sD::EDGE_HO_SE ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_HO_SE ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_HO_SE ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 1, sD::EDGE_HO_SE ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 2, sD::EDGE_HO_SE ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_HO_SE ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 2, sD::EDGE_HO_SE ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 3, sD::EDGE_HO_SE ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 3, sD::EDGE_HO_SE ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 4, sD::EDGE_HO_SE ), 9 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 1, sD::EDGE_HO_SE ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_HO_SE ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_HO_SE ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_HO_SE ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_HO_SE ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_HO_SE ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_HO_SE ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 1, sD::EDGE_HO_SE ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 2, sD::EDGE_HO_SE ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_HO_SE ), 9 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_HO_SE ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_HO_SE ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_HO_SE ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_HO_SE ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 2, sD::EDGE_HO_SE ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 3, sD::EDGE_HO_SE ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_HO_SE ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_HO_SE ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_HO_SE ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_HO_SE ), 19 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 3, sD::EDGE_HO_SE ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 4, sD::EDGE_HO_SE ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_HO_SE ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_HO_SE ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_HO_SE ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 4, sD::EDGE_HO_SE ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_HO_SE ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_HO_SE ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_HO_SE ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 5, sD::EDGE_HO_SE ), 29 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 6, sD::EDGE_HO_SE ), 30 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_HO_SE ), 31 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 6, sD::EDGE_HO_SE ), 32 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 7, sD::EDGE_HO_SE ), 33 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 7, sD::EDGE_HO_SE ), 34 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 8, sD::EDGE_HO_SE ), 35 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing: EDGE_HO_SE correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 0, sD::EDGE_VE_N ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 0, sD::EDGE_VE_N ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 0, sD::EDGE_VE_N ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 0, sD::EDGE_VE_N ), 23 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 1, sD::EDGE_VE_N ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_VE_N ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_VE_N ), 26 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 2, sD::EDGE_VE_N ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_VE_N ), 28 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 3, sD::EDGE_VE_N ), 29 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing: EDGE_VE_N correct!" );

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
