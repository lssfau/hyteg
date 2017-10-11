
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

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_HO_E correct!" );

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

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_HO_W correct!" );

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

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_HO_NW correct!" );

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

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_HO_SE correct!" );

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

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 0, sD::EDGE_VE_N ), 72 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 0, sD::EDGE_VE_N ), 73 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 0, sD::EDGE_VE_N ), 74 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 0, sD::EDGE_VE_N ), 75 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 0, sD::EDGE_VE_N ), 76 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 0, sD::EDGE_VE_N ), 77 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 0, sD::EDGE_VE_N ), 78 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 0, sD::EDGE_VE_N ), 79 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 1, sD::EDGE_VE_N ), 80 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_VE_N ), 81 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_VE_N ), 82 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_VE_N ), 83 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_VE_N ), 84 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_VE_N ), 85 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_VE_N ), 86 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 2, sD::EDGE_VE_N ), 87 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_VE_N ), 88 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_VE_N ), 89 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_VE_N ), 90 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_VE_N ), 91 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_VE_N ), 92 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 3, sD::EDGE_VE_N ), 93 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_VE_N ), 94 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_VE_N ), 95 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_VE_N ), 96 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_VE_N ), 97 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 4, sD::EDGE_VE_N ), 98 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_VE_N ), 99 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_VE_N ), 100 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_VE_N ), 101 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_VE_N ), 102 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_VE_N ), 103 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_VE_N ), 104 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 6, sD::EDGE_VE_N ), 105 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_VE_N ), 106 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 7, sD::EDGE_VE_N ), 107 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_VE_N correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 1, sD::EDGE_VE_S ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_VE_S ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_VE_S ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 1, sD::EDGE_VE_S ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 2, sD::EDGE_VE_S ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_VE_S ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 2, sD::EDGE_VE_S ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 3, sD::EDGE_VE_S ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 3, sD::EDGE_VE_S ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 4, sD::EDGE_VE_S ), 29 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 1, sD::EDGE_VE_S ), 72 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_VE_S ), 73 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_VE_S ), 74 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_VE_S ), 75 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_VE_S ), 76 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_VE_S ), 77 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_VE_S ), 78 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 1, sD::EDGE_VE_S ), 79 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 2, sD::EDGE_VE_S ), 80 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_VE_S ), 81 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_VE_S ), 82 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_VE_S ), 83 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_VE_S ), 84 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_VE_S ), 85 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 2, sD::EDGE_VE_S ), 86 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 3, sD::EDGE_VE_S ), 87 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_VE_S ), 88 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_VE_S ), 89 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_VE_S ), 90 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_VE_S ), 91 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 3, sD::EDGE_VE_S ), 92 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 4, sD::EDGE_VE_S ), 93 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_VE_S ), 94 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_VE_S ), 95 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_VE_S ), 96 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 4, sD::EDGE_VE_S ), 97 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_VE_S ), 98 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_VE_S ), 99 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_VE_S ), 100 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 5, sD::EDGE_VE_S ), 101 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 6, sD::EDGE_VE_S ), 102 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_VE_S ), 103 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 6, sD::EDGE_VE_S ), 104 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 7, sD::EDGE_VE_S ), 105 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 7, sD::EDGE_VE_S ), 106 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 8, sD::EDGE_VE_S ), 107 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_VE_S correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 0, sD::EDGE_VE_NW ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 0, sD::EDGE_VE_NW ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 0, sD::EDGE_VE_NW ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 4, 0, sD::EDGE_VE_NW ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_VE_NW ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_VE_NW ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 1, sD::EDGE_VE_NW ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_VE_NW ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 2, sD::EDGE_VE_NW ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 3, sD::EDGE_VE_NW ), 29 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 0, sD::EDGE_VE_NW ), 72 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 0, sD::EDGE_VE_NW ), 73 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 0, sD::EDGE_VE_NW ), 74 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 0, sD::EDGE_VE_NW ), 75 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 0, sD::EDGE_VE_NW ), 76 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 0, sD::EDGE_VE_NW ), 77 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 0, sD::EDGE_VE_NW ), 78 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 8, 0, sD::EDGE_VE_NW ), 79 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_VE_NW ), 80 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_VE_NW ), 81 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_VE_NW ), 82 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_VE_NW ), 83 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_VE_NW ), 84 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_VE_NW ), 85 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 1, sD::EDGE_VE_NW ), 86 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_VE_NW ), 87 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_VE_NW ), 88 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_VE_NW ), 89 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_VE_NW ), 90 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_VE_NW ), 91 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 2, sD::EDGE_VE_NW ), 92 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_VE_NW ), 93 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_VE_NW ), 94 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_VE_NW ), 95 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_VE_NW ), 96 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 3, sD::EDGE_VE_NW ), 97 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_VE_NW ), 98 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_VE_NW ), 99 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_VE_NW ), 100 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 4, sD::EDGE_VE_NW ), 101 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_VE_NW ), 102 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_VE_NW ), 103 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 5, sD::EDGE_VE_NW ), 104 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_VE_NW ), 105 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 6, sD::EDGE_VE_NW ), 106 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 7, sD::EDGE_VE_NW ), 107 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_VE_NW correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 1, sD::EDGE_VE_SE ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_VE_SE ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_VE_SE ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 2, sD::EDGE_VE_SE ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_VE_SE ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 3, sD::EDGE_VE_SE ), 28 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 1, sD::EDGE_VE_SE ), 73 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_VE_SE ), 74 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_VE_SE ), 75 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_VE_SE ), 76 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_VE_SE ), 77 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_VE_SE ), 78 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_VE_SE ), 79 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 2, sD::EDGE_VE_SE ), 81 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_VE_SE ), 82 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_VE_SE ), 83 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_VE_SE ), 84 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_VE_SE ), 85 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_VE_SE ), 86 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 3, sD::EDGE_VE_SE ), 88 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_VE_SE ), 89 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_VE_SE ), 90 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_VE_SE ), 91 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_VE_SE ), 92 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 4, sD::EDGE_VE_SE ), 94 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_VE_SE ), 95 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_VE_SE ), 96 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_VE_SE ), 97 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_VE_SE ), 99 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_VE_SE ), 100 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_VE_SE ), 101 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 6, sD::EDGE_VE_SE ), 103 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_VE_SE ), 104 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 7, sD::EDGE_VE_SE ), 106 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_VE_SE correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 0, sD::EDGE_DI_NW ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 0, sD::EDGE_DI_NW ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 0, sD::EDGE_DI_NW ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 4, 0, sD::EDGE_DI_NW ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_DI_NW ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_DI_NW ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 1, sD::EDGE_DI_NW ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_DI_NW ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 2, sD::EDGE_DI_NW ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 3, sD::EDGE_DI_NW ), 19 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 0, sD::EDGE_DI_NW ), 36 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 0, sD::EDGE_DI_NW ), 37 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 0, sD::EDGE_DI_NW ), 38 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 0, sD::EDGE_DI_NW ), 39 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 0, sD::EDGE_DI_NW ), 40 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 0, sD::EDGE_DI_NW ), 41 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 0, sD::EDGE_DI_NW ), 42 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 8, 0, sD::EDGE_DI_NW ), 43 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_DI_NW ), 44 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_DI_NW ), 45 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_DI_NW ), 46 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_DI_NW ), 47 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_DI_NW ), 48 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_DI_NW ), 49 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 1, sD::EDGE_DI_NW ), 50 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_DI_NW ), 51 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_DI_NW ), 52 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_DI_NW ), 53 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_DI_NW ), 54 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_DI_NW ), 55 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 2, sD::EDGE_DI_NW ), 56 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_DI_NW ), 57 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_DI_NW ), 58 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_DI_NW ), 59 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_DI_NW ), 60 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 3, sD::EDGE_DI_NW ), 61 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_DI_NW ), 62 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_DI_NW ), 63 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_DI_NW ), 64 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 4, sD::EDGE_DI_NW ), 65 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_DI_NW ), 66 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_DI_NW ), 67 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 5, sD::EDGE_DI_NW ), 68 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_DI_NW ), 69 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 6, sD::EDGE_DI_NW ), 70 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 7, sD::EDGE_DI_NW ), 71 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_DI_NW correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 0, sD::EDGE_DI_NE ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 0, sD::EDGE_DI_NE ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 0, sD::EDGE_DI_NE ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 0, sD::EDGE_DI_NE ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 1, sD::EDGE_DI_NE ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_DI_NE ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_DI_NE ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 2, sD::EDGE_DI_NE ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_DI_NE ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 3, sD::EDGE_DI_NE ), 19 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 0, sD::EDGE_DI_NE ), 36 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 0, sD::EDGE_DI_NE ), 37 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 0, sD::EDGE_DI_NE ), 38 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 0, sD::EDGE_DI_NE ), 39 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 0, sD::EDGE_DI_NE ), 40 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 0, sD::EDGE_DI_NE ), 41 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 0, sD::EDGE_DI_NE ), 42 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 0, sD::EDGE_DI_NE ), 43 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 1, sD::EDGE_DI_NE ), 44 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_DI_NE ), 45 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_DI_NE ), 46 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_DI_NE ), 47 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_DI_NE ), 48 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_DI_NE ), 49 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_DI_NE ), 50 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 2, sD::EDGE_DI_NE ), 51 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_DI_NE ), 52 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_DI_NE ), 53 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_DI_NE ), 54 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_DI_NE ), 55 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_DI_NE ), 56 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 3, sD::EDGE_DI_NE ), 57 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_DI_NE ), 58 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_DI_NE ), 59 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_DI_NE ), 60 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_DI_NE ), 61 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 4, sD::EDGE_DI_NE ), 62 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_DI_NE ), 63 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_DI_NE ), 64 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_DI_NE ), 65 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_DI_NE ), 66 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_DI_NE ), 67 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_DI_NE ), 68 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 6, sD::EDGE_DI_NE ), 69 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_DI_NE ), 70 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 7, sD::EDGE_DI_NE ), 71 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_DI_NE correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 1, sD::EDGE_DI_SE ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_DI_SE ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_DI_SE ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 1, sD::EDGE_DI_SE ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 2, sD::EDGE_DI_SE ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_DI_SE ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 2, sD::EDGE_DI_SE ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 3, sD::EDGE_DI_SE ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 3, sD::EDGE_DI_SE ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 0, 4, sD::EDGE_DI_SE ), 19 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 1, sD::EDGE_DI_SE ), 36 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_DI_SE ), 37 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_DI_SE ), 38 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_DI_SE ), 39 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_DI_SE ), 40 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_DI_SE ), 41 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_DI_SE ), 42 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 1, sD::EDGE_DI_SE ), 43 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 2, sD::EDGE_DI_SE ), 44 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_DI_SE ), 45 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_DI_SE ), 46 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_DI_SE ), 47 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_DI_SE ), 48 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_DI_SE ), 49 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 2, sD::EDGE_DI_SE ), 50 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 3, sD::EDGE_DI_SE ), 51 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_DI_SE ), 52 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_DI_SE ), 53 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_DI_SE ), 54 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_DI_SE ), 55 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 3, sD::EDGE_DI_SE ), 56 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 4, sD::EDGE_DI_SE ), 57 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_DI_SE ), 58 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_DI_SE ), 59 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_DI_SE ), 60 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 4, sD::EDGE_DI_SE ), 61 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_DI_SE ), 62 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_DI_SE ), 63 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_DI_SE ), 64 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 5, sD::EDGE_DI_SE ), 65 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 6, sD::EDGE_DI_SE ), 66 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_DI_SE ), 67 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 6, sD::EDGE_DI_SE ), 68 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 7, sD::EDGE_DI_SE ), 69 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 7, sD::EDGE_DI_SE ), 70 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 8, sD::EDGE_DI_SE ), 71 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_DI_SE correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 1, sD::EDGE_DI_SW ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 1, sD::EDGE_DI_SW ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 3, 1, sD::EDGE_DI_SW ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 2, sD::EDGE_DI_SW ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 2, 2, sD::EDGE_DI_SW ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 2 >( 1, 3, sD::EDGE_DI_SW ), 17 );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 1, sD::EDGE_DI_SW ), 36 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_DI_SW ), 37 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 1, sD::EDGE_DI_SW ), 38 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 1, sD::EDGE_DI_SW ), 39 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 1, sD::EDGE_DI_SW ), 40 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 1, sD::EDGE_DI_SW ), 41 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 7, 1, sD::EDGE_DI_SW ), 42 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 2, sD::EDGE_DI_SW ), 44 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_DI_SW ), 45 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 2, sD::EDGE_DI_SW ), 46 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 2, sD::EDGE_DI_SW ), 47 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 2, sD::EDGE_DI_SW ), 48 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 6, 2, sD::EDGE_DI_SW ), 49 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 3, sD::EDGE_DI_SW ), 51 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 3, sD::EDGE_DI_SW ), 52 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 3, sD::EDGE_DI_SW ), 53 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 3, sD::EDGE_DI_SW ), 54 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 5, 3, sD::EDGE_DI_SW ), 55 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 4, sD::EDGE_DI_SW ), 57 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 4, sD::EDGE_DI_SW ), 58 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 4, sD::EDGE_DI_SW ), 59 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 4, 4, sD::EDGE_DI_SW ), 60 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 5, sD::EDGE_DI_SW ), 62 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 5, sD::EDGE_DI_SW ), 63 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 3, 5, sD::EDGE_DI_SW ), 64 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 6, sD::EDGE_DI_SW ), 66 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 6, sD::EDGE_DI_SW ), 67 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 1, 7, sD::EDGE_DI_SW ), 69 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertex): EDGE_DI_SW correct!" );

#ifdef NDEBUG
#ifndef _MSC_VER
  static_assert( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 0, sD::EDGE_HO_E )  ==   0 , "EDGE_HO_E  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 1, sD::EDGE_HO_E )  ==  10 , "EDGE_HO_E  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 2, 2, sD::EDGE_HO_E )  ==  17 , "EDGE_HO_E  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVertex < 3 >( 0, 5, sD::EDGE_VE_N )  == 102 , "EDGE_VE_N  cannot be statically computed by the compiler!" );
#endif
#endif

  /////////////////////////////////
  // access from horizontal edge //
  /////////////////////////////////

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 0, sD::EDGE_HO_C ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 0, sD::EDGE_HO_C ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 0, sD::EDGE_HO_C ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 3, 0, sD::EDGE_HO_C ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 1, sD::EDGE_HO_C ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 1, sD::EDGE_HO_C ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 1, sD::EDGE_HO_C ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 2, sD::EDGE_HO_C ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 2, sD::EDGE_HO_C ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 3, sD::EDGE_HO_C ), 9 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from horizontal edge): EDGE_HO_C correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 0, sD::EDGE_DI_N ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 0, sD::EDGE_DI_N ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 0, sD::EDGE_DI_N ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 3, 0, sD::EDGE_DI_N ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 1, sD::EDGE_DI_N ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 1, sD::EDGE_DI_N ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 1, sD::EDGE_DI_N ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 2, sD::EDGE_DI_N ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 2, sD::EDGE_DI_N ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 3, sD::EDGE_DI_N ), 19 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from horizontal edge): EDGE_DI_N correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 1, sD::EDGE_DI_S ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 1, sD::EDGE_DI_S ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 1, sD::EDGE_DI_S ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 2, sD::EDGE_DI_S ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 2, sD::EDGE_DI_S ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 3, sD::EDGE_DI_S ), 17 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from horizontal edge): EDGE_DI_S correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 0, sD::EDGE_VE_NW ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 0, sD::EDGE_VE_NW ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 0, sD::EDGE_VE_NW ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 3, 0, sD::EDGE_VE_NW ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 1, sD::EDGE_VE_NW ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 1, sD::EDGE_VE_NW ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 1, sD::EDGE_VE_NW ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 2, sD::EDGE_VE_NW ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 2, sD::EDGE_VE_NW ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 3, sD::EDGE_VE_NW ), 29 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from horizontal edge): EDGE_VE_NW correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 1, sD::EDGE_VE_SE ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 1, sD::EDGE_VE_SE ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 1, sD::EDGE_VE_SE ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 2, sD::EDGE_VE_SE ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 2, sD::EDGE_VE_SE ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 3, sD::EDGE_VE_SE ), 28 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from horizontal edge): EDGE_VE_SE correct!" );

#ifdef NDEBUG
#ifndef _MSC_VER
  static_assert( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 0, sD::EDGE_DI_N )  == 10, "EDGE_DI_N  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 1, 1, sD::EDGE_DI_S )  == 11, "EDGE_DI_S  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 2, 0, sD::EDGE_VE_NW ) == 22, "EDGE_VE_NW cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromHorizontalEdge < 2 >( 0, 3, sD::EDGE_VE_SE ) == 28, "EDGE_VE_SE cannot be statically computed by the compiler!" );
#endif
#endif

  ///////////////////////////////
  // access from vertical edge //
  ///////////////////////////////

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 0, sD::EDGE_VE_C ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 0, sD::EDGE_VE_C ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 0, sD::EDGE_VE_C ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 3, 0, sD::EDGE_VE_C ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 1, sD::EDGE_VE_C ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 1, sD::EDGE_VE_C ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 1, sD::EDGE_VE_C ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 2, sD::EDGE_VE_C ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 2, sD::EDGE_VE_C ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 3, sD::EDGE_VE_C ), 29 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertical edge): EDGE_VE_C correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 0, sD::EDGE_DI_E ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 0, sD::EDGE_DI_E ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 0, sD::EDGE_DI_E ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 3, 0, sD::EDGE_DI_E ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 1, sD::EDGE_DI_E ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 1, sD::EDGE_DI_E ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 1, sD::EDGE_DI_E ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 2, sD::EDGE_DI_E ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 2, sD::EDGE_DI_E ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 3, sD::EDGE_DI_E ), 19 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertical edge): EDGE_DI_E correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 0, sD::EDGE_DI_W ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 0, sD::EDGE_DI_W ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 3, 0, sD::EDGE_DI_W ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 1, sD::EDGE_DI_W ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 1, sD::EDGE_DI_W ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 2, sD::EDGE_DI_W ), 17 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertical edge): EDGE_DI_W correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 0, sD::EDGE_HO_SE ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 0, sD::EDGE_HO_SE ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 0, sD::EDGE_HO_SE ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 3, 0, sD::EDGE_HO_SE ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 1, sD::EDGE_HO_SE ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 1, sD::EDGE_HO_SE ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 1, sD::EDGE_HO_SE ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 2, sD::EDGE_HO_SE ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 2, sD::EDGE_HO_SE ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 3, sD::EDGE_HO_SE ), 9 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertical edge): EDGE_HO_SE correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 0, sD::EDGE_HO_NW ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 0, sD::EDGE_HO_NW ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 3, 0, sD::EDGE_HO_NW ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 1, sD::EDGE_HO_NW ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 1, sD::EDGE_HO_NW ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 2, sD::EDGE_HO_NW ), 9 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from vertical edge): EDGE_HO_NW correct!" );

#ifdef NDEBUG
#ifndef _MSC_VER
  static_assert( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 0, 0, sD::EDGE_DI_E )  == 10, "EDGE_DI_E  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 0, sD::EDGE_DI_W )  == 11, "EDGE_DI_W  cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 2, 0, sD::EDGE_HO_SE ) ==  2, "EDGE_HO_SE cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromVerticalEdge < 2 >( 1, 2, sD::EDGE_HO_NW ) ==  9, "EDGE_HO_NW cannot be statically computed by the compiler!" );
#endif
#endif

  ///////////////////////////////
  // access from diagonal edge //
  ///////////////////////////////

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 0, sD::EDGE_DI_C ), 10 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 0, sD::EDGE_DI_C ), 11 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 0, sD::EDGE_DI_C ), 12 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 3, 0, sD::EDGE_DI_C ), 13 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 1, sD::EDGE_DI_C ), 14 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 1, sD::EDGE_DI_C ), 15 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 1, sD::EDGE_DI_C ), 16 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 2, sD::EDGE_DI_C ), 17 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 2, sD::EDGE_DI_C ), 18 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 3, sD::EDGE_DI_C ), 19 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from diagonal edge): EDGE_DI_C correct!" );


  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 0, sD::EDGE_HO_N ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 0, sD::EDGE_HO_N ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 0, sD::EDGE_HO_N ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 1, sD::EDGE_HO_N ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 1, sD::EDGE_HO_N ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 2, sD::EDGE_HO_N ), 9 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from diagonal edge): EDGE_HO_N correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 0, sD::EDGE_HO_S ), 0 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 0, sD::EDGE_HO_S ), 1 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 0, sD::EDGE_HO_S ), 2 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 3, 0, sD::EDGE_HO_S ), 3 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 1, sD::EDGE_HO_S ), 4 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 1, sD::EDGE_HO_S ), 5 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 1, sD::EDGE_HO_S ), 6 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 2, sD::EDGE_HO_S ), 7 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 2, sD::EDGE_HO_S ), 8 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 3, sD::EDGE_HO_S ), 9 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from diagonal edge): EDGE_HO_S correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 0, sD::EDGE_VE_W ), 20 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 0, sD::EDGE_VE_W ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 0, sD::EDGE_VE_W ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 3, 0, sD::EDGE_VE_W ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 1, sD::EDGE_VE_W ), 24 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 1, sD::EDGE_VE_W ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 1, sD::EDGE_VE_W ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 2, sD::EDGE_VE_W ), 27 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 2, sD::EDGE_VE_W ), 28 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 3, sD::EDGE_VE_W ), 29 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from diagonal edge): EDGE_VE_W correct!" );

  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 0, sD::EDGE_VE_E ), 21 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 0, sD::EDGE_VE_E ), 22 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 0, sD::EDGE_VE_E ), 23 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 1, sD::EDGE_VE_E ), 25 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 1, sD::EDGE_VE_E ), 26 );
  WALBERLA_CHECK_EQUAL( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 2, sD::EDGE_VE_E ), 28 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro face indexing (from diagonal edge): EDGE_VE_E correct!" );

#ifdef NDEBUG
#ifndef _MSC_VER
  static_assert( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 0, sD::EDGE_HO_N ) ==  4, "EDGE_HO_N cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 2, 0, sD::EDGE_HO_S ) ==  2, "EDGE_HO_S cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 1, 1, sD::EDGE_VE_W ) == 25, "EDGE_VE_W cannot be statically computed by the compiler!" );
  static_assert( EdgeDoFOnMacroFace::indexFromDiagonalEdge < 2 >( 0, 1, sD::EDGE_VE_E ) == 25, "EDGE_VE_E cannot be statically computed by the compiler!" );
#endif
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
