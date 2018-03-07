
#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

static void testEdgeDoFsOnMacroEdge()
{
  typedef stencilDirection sD;

  ///////////////////////
  // index from vertex //
  ///////////////////////

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 0, sD::EDGE_HO_E ), 0 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_HO_E ), 1 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_HO_E ), 2 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_HO_E ), 3 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 0, sD::EDGE_HO_E ), 0 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_HO_E ), 1 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_HO_E ), 2 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_HO_E ), 3 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_HO_E ), 4 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_HO_E ), 5 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_HO_E ), 6 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_HO_E ), 7 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_HO_E correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_HO_W ), 0 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_HO_W ), 1 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_HO_W ), 2 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 4, sD::EDGE_HO_W ), 3 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_HO_W ), 0 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_HO_W ), 1 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_HO_W ), 2 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_HO_W ), 3 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_HO_W ), 4 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_HO_W ), 5 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_HO_W ), 6 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 8, sD::EDGE_HO_W ), 7 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_HO_W correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_HO_NW ), 15 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_HO_NW ), 16 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_HO_NW ), 17 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_HO_NW ), 31 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_HO_NW ), 32 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_HO_NW ), 33 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_HO_NW ), 34 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_HO_NW ), 35 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_HO_NW ), 36 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_HO_NW ), 37 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_HO_NW correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_HO_SE ), 4 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_HO_SE ), 5 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_HO_SE ), 6 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_HO_SE ), 8 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_HO_SE ), 9 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_HO_SE ), 10 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_HO_SE ), 11 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_HO_SE ), 12 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_HO_SE ), 13 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_HO_SE ), 14 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_HO_SE correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_VE_S ), 11 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_VE_S ), 12 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_VE_S ), 13 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 4, sD::EDGE_VE_S ), 14 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_VE_S ), 23 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_VE_S ), 24 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_VE_S ), 25 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_VE_S ), 26 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_VE_S ), 27 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_VE_S ), 28 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_VE_S ), 29 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 8, sD::EDGE_VE_S ), 30 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_VE_S correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 0, sD::EDGE_VE_N ), 22 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_VE_N ), 23 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_VE_N ), 24 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_VE_N ), 25 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 0, sD::EDGE_VE_N ), 46 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_VE_N ), 47 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_VE_N ), 48 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_VE_N ), 49 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_VE_N ), 50 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_VE_N ), 51 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_VE_N ), 52 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_VE_N ), 53 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_VE_N correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_VE_NW ), 22 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_VE_NW ), 23 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_VE_NW ), 24 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 4, sD::EDGE_VE_NW ), 25 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_VE_NW ), 46 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_VE_NW ), 47 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_VE_NW ), 48 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_VE_NW ), 49 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_VE_NW ), 50 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_VE_NW ), 51 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_VE_NW ), 52 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 8, sD::EDGE_VE_NW ), 53 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_VE_NW correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 0, sD::EDGE_VE_SE ), 11 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_VE_SE ), 12 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_VE_SE ), 13 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_VE_SE ), 14 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 0, sD::EDGE_VE_SE ), 23 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_VE_SE ), 24 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_VE_SE ), 25 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_VE_SE ), 26 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_VE_SE ), 27 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_VE_SE ), 28 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_VE_SE ), 29 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_VE_SE ), 30 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_VE_SE correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_DI_NW ), 18 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_DI_NW ), 19 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_DI_NW ), 20 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 4, sD::EDGE_DI_NW ), 21 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_DI_NW ), 38 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_DI_NW ), 39 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_DI_NW ), 40 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_DI_NW ), 41 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_DI_NW ), 42 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_DI_NW ), 43 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_DI_NW ), 44 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 8, sD::EDGE_DI_NW ), 45 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_DI_NW correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 0, sD::EDGE_DI_NE ), 18 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_DI_NE ), 19 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_DI_NE ), 20 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_DI_NE ), 21 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 0, sD::EDGE_DI_NE ), 38 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_DI_NE ), 39 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_DI_NE ), 40 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_DI_NE ), 41 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_DI_NE ), 42 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_DI_NE ), 43 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_DI_NE ), 44 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_DI_NE ), 45 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_DI_NE correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_DI_SW ), 7 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_DI_SW ), 8 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_DI_SW ), 9 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 4, sD::EDGE_DI_SW ), 10 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_DI_SW ), 15 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_DI_SW ), 16 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_DI_SW ), 17 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_DI_SW ), 18 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_DI_SW ), 19 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_DI_SW ), 20 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_DI_SW ), 21 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 8, sD::EDGE_DI_SW ), 22 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_DI_SW correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 0, sD::EDGE_DI_SE ), 7 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 1, sD::EDGE_DI_SE ), 8 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 2, sD::EDGE_DI_SE ), 9 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 2, 3, sD::EDGE_DI_SE ), 10 );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 0, sD::EDGE_DI_SE ), 15 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_DI_SE ), 16 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 2, sD::EDGE_DI_SE ), 17 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_DI_SE ), 18 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_DI_SE ), 19 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 5, sD::EDGE_DI_SE ), 20 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_DI_SE ), 21 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_DI_SE ), 22 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from vertex): EDGE_DI_SE correct!" );

#ifdef NDEBUG
  static_assert( edgedof::macroedge::indexFromVertex( 3, 3, sD::EDGE_HO_E)  ==  3 , "EDGE_HO_E  cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_HO_W)  ==  5 , "EDGE_HO_W  cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 8, sD::EDGE_VE_S)  == 30 , "EDGE_VE_S  cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_VE_N)  == 50 , "EDGE_VE_N  cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_DI_NW) == 43 , "EDGE_DI_NW cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 4, sD::EDGE_DI_NE) == 42 , "EDGE_DI_NE cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_DI_SW) == 21 , "EDGE_DI_SW cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 7, sD::EDGE_DI_SE) == 22 , "EDGE_DI_SE cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 8, sD::EDGE_HO_NW) == 38 , "EDGE_HO_NW cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_HO_SE) == 13 , "EDGE_HO_SE cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 1, sD::EDGE_VE_NW) == 46 , "EDGE_VE_NW cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromVertex( 3, 6, sD::EDGE_VE_SE) == 29 , "EDGE_VE_SE cannot be statically computed by the compiler!" );
#endif

  ////////////////////////////////
  // index from horizontal edge //
  ////////////////////////////////

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 0, sD::EDGE_HO_C ), 0 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 1, sD::EDGE_HO_C ), 1 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 2, sD::EDGE_HO_C ), 2 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 3, sD::EDGE_HO_C ), 3 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from horizontal edge): EDGE_HO_C correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 0, sD::EDGE_DI_N ), 18 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 1, sD::EDGE_DI_N ), 19 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 2, sD::EDGE_DI_N ), 20 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 3, sD::EDGE_DI_N ), 21 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from horizontal edge): EDGE_DI_N correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 0, sD::EDGE_DI_S ), 7 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 1, sD::EDGE_DI_S ), 8 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 2, sD::EDGE_DI_S ), 9 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 3, sD::EDGE_DI_S ), 10 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from horizontal edge): EDGE_DI_S correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 0, sD::EDGE_VE_NW ), 22 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 1, sD::EDGE_VE_NW ), 23 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 2, sD::EDGE_VE_NW ), 24 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 3, sD::EDGE_VE_NW ), 25 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from horizontal edge): EDGE_VE_NW correct!" );

  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 0, sD::EDGE_VE_SE ), 11 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 1, sD::EDGE_VE_SE ), 12 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 2, sD::EDGE_VE_SE ), 13 );
  WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexFromHorizontalEdge( 2, 3, sD::EDGE_VE_SE ), 14 );

  WALBERLA_LOG_INFO_ON_ROOT( "Edge DoFs on macro edge indexing (from horizontal edge): EDGE_VE_SE correct!" );

#ifdef NDEBUG
  static_assert( edgedof::macroedge::indexFromHorizontalEdge( 2, 0, sD::EDGE_DI_N)  == 18, "EDGE_DI_N  cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromHorizontalEdge( 2, 1, sD::EDGE_DI_S)  ==  8, "EDGE_DI_S  cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromHorizontalEdge( 2, 2, sD::EDGE_VE_NW) == 24, "EDGE_VE_NW cannot be statically computed by the compiler!" );
  static_assert( edgedof::macroedge::indexFromHorizontalEdge( 2, 3, sD::EDGE_VE_SE) == 14, "EDGE_VE_SE cannot be statically computed by the compiler!" );
#endif

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testEdgeDoFsOnMacroEdge();

   return EXIT_SUCCESS;
}
