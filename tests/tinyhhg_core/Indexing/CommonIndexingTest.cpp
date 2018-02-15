
#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp>
#include <iostream>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"
#include "core/logging/all.h"

#include "tinyhhg_core/levelinfo.hpp"

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/MacroCellIndexing.hpp"
#include "tinyhhg_core/indexing/Optimization.hpp"

namespace hhg {

static void testFaceBorderIterator( const std::vector< std::array< uint_t, 2 > > & expectedValues, const indexing::FaceBorderDirection & faceBorderDirection,
                                    const uint_t & width, const uint_t & offsetToCenter, const uint_t & offsetFromVertices )
{
  std::vector< std::array< uint_t, 2 > > iteratorResult;

  for ( const auto & it : indexing::FaceBorderIterator( width, faceBorderDirection, offsetToCenter, offsetFromVertices ) )
  {
    iteratorResult.push_back( {{ it.col(), it.row() }} );
  }

  WALBERLA_CHECK_EQUAL( iteratorResult.size(), expectedValues.size() );

  for ( uint_t i = 0; i < iteratorResult.size(); i++ )
  {
    WALBERLA_CHECK_EQUAL( iteratorResult[i][0], expectedValues[i][0] );
    WALBERLA_CHECK_EQUAL( iteratorResult[i][1], expectedValues[i][1] );
  }
}

static void testCellIterator( const std::vector< std::array< uint_t, 3 > > & expectedValues, const uint_t & width, const uint_t & offsetToCenter )
{
  std::vector< std::array< uint_t, 3 > > iteratorResult;

  for ( const auto & it : indexing::CellIterator( width, offsetToCenter ) )
  {
    iteratorResult.push_back( {{ it.col(), it.row(), it.dep() }} );
  }

  WALBERLA_CHECK_EQUAL( iteratorResult.size(), expectedValues.size() );

  for ( uint_t i = 0; i < iteratorResult.size(); i++ )
  {
    WALBERLA_CHECK_EQUAL( iteratorResult[i][0], expectedValues[i][0] );
    WALBERLA_CHECK_EQUAL( iteratorResult[i][1], expectedValues[i][1] );
    WALBERLA_CHECK_EQUAL( iteratorResult[i][2], expectedValues[i][2] );
  }
}

static void testCellBorderIterator( const std::vector< std::array< uint_t, 3 > > & expectedValues, const uint_t & width,
                                    const std::array< uint_t, 3 > & vertices, const uint_t & offsetToCenter = 0 )
{
  std::vector< std::array< uint_t, 3 > > iteratorResult;

  for ( const auto & it : indexing::CellBorderIterator( width, vertices, offsetToCenter ) )
  {
    iteratorResult.push_back( {{ it.col(), it.row(), it.dep() }} );
  }

  WALBERLA_CHECK_EQUAL( iteratorResult.size(), expectedValues.size() );

  for ( uint_t i = 0; i < iteratorResult.size(); i++ )
  {
    WALBERLA_CHECK_EQUAL( iteratorResult[i][0], expectedValues[i][0] );
    WALBERLA_CHECK_EQUAL( iteratorResult[i][1], expectedValues[i][1] );
    WALBERLA_CHECK_EQUAL( iteratorResult[i][2], expectedValues[i][2] );
  }
}
static void testCommonIndexing()
{
  using walberla::uint_t;

  uint_t testCol = 0;
  uint_t testRow = 0;

  // macro edge

  for ( const auto & it : indexing::EdgeIterator( 8, 0 ) )
  {
    WALBERLA_CHECK_EQUAL( testCol++, it.col() );
    WALBERLA_CHECK_EQUAL( 0,         it.row() );
  }
  WALBERLA_CHECK_EQUAL( testCol, 8 );

  testCol = 1;
  for ( const auto & it : indexing::EdgeIterator( 8, 1 ) )
  {
    WALBERLA_CHECK_EQUAL( testCol++, it.col() );
    WALBERLA_CHECK_EQUAL( 0,         it.row() );
  }
  WALBERLA_CHECK_EQUAL( testCol, 7 );

  // macro face

  testCol = 0;
  testRow = 0;
  for ( const auto & it : indexing::FaceIterator( 4, 0 ) )
  {
    WALBERLA_CHECK_EQUAL( testCol++, it.col() );
    WALBERLA_CHECK_EQUAL( testRow  , it.row() );

    if ( testCol + testRow == 4 )
    {
      testRow++;
      testCol = 0;
    }
  }
  WALBERLA_CHECK_EQUAL( testRow, 4 );
  WALBERLA_CHECK_EQUAL( testCol, 0 );

  testCol = 1;
  testRow = 1;
  for ( const auto & it : indexing::FaceIterator( 4, 1 ) )
  {
    testCol++;
    testRow++;
    WALBERLA_CHECK_EQUAL( it.col(), 1 );
    WALBERLA_CHECK_EQUAL( it.row(), 1 );
  }
  WALBERLA_CHECK_EQUAL( testRow, 2 );
  WALBERLA_CHECK_EQUAL( testCol, 2 );

  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{0, 0}}, {{1, 0}}, {{2, 0}}, {{3, 0}} }} ),           indexing::FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT, 4, 0, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{2, 0}}, {{1, 0}} }} ),                               indexing::FaceBorderDirection::BOTTOM_RIGHT_TO_LEFT, 4, 0, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{3, 1}}, {{2, 1}}, {{1, 1}} }} ),                     indexing::FaceBorderDirection::BOTTOM_RIGHT_TO_LEFT, 4, 1, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{2, 1}} }} ),                                         indexing::FaceBorderDirection::BOTTOM_RIGHT_TO_LEFT, 4, 1, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{4, 0}}, {{3, 1}}, {{2, 2}}, {{1, 3}}, {{0, 4}} }} ), indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 5, 0, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{3, 1}}, {{2, 2}}, {{1, 3}} }} ),                     indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 5, 0, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{3, 0}}, {{2, 1}}, {{1, 2}}, {{0, 3}} }} ),           indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 5, 1, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{1, 2}}, {{2, 1}} }} ),                               indexing::FaceBorderDirection::DIAGONAL_TOP_TO_BOTTOM, 5, 1, 1 );

  // macro cell

  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  1 >(),   1 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  2 >(),   4 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  3 >(),  10 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  4 >(),  20 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  5 >(),  35 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  6 >(),  56 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  7 >(),  84 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  8 >(), 120 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize<  9 >(), 165 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize< 10 >(), 220 );

  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex< 5 >( 0, 0, 0 ),  0 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex< 5 >( 2, 0, 0 ),  2 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex< 5 >( 1, 3, 0 ), 13 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex< 5 >( 1, 1, 1 ), 20 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex< 5 >( 1, 1, 2 ), 29 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex< 5 >( 0, 0, 4 ), 34 );

  testCellIterator( std::vector< std::array< uint_t, 3 > >( {{ {{ 0, 0, 0 }}, {{ 1, 0, 0 }}, {{ 2, 0, 0 }}, {{ 3, 0, 0 }}, {{ 0, 1, 0 }}, {{ 1, 1, 0 }},
                                                               {{ 2, 1, 0 }}, {{ 0, 2, 0 }}, {{ 1, 2, 0 }}, {{ 0, 3, 0 }}, {{ 0, 0, 1 }}, {{ 1, 0, 1 }},
                                                               {{ 2, 0, 1 }}, {{ 0, 1, 1 }}, {{ 1, 1, 1 }}, {{ 0, 2, 1 }}, {{ 0, 0, 2 }}, {{ 1, 0, 2 }},
                                                               {{ 0, 1, 2 }}, {{ 0, 0, 3 }} }} ),
                                                            4, 0 );

  testCellIterator( std::vector< std::array< uint_t, 3 > >( {{ {{ 1, 1, 1 }} }} ), 5, 1 );

  testCellIterator( std::vector< std::array< uint_t, 3 > >( {{ {{ 1, 1, 1 }}, {{ 2, 1, 1 }}, {{ 3, 1, 1 }}, {{ 4, 1, 1 }}, {{ 5, 1, 1 }}, {{ 1, 2, 1 }}, {{ 2, 2, 1 }}, {{ 3, 2, 1 }}, {{ 4, 2, 1 }},
                                                               {{ 1, 3, 1 }}, {{ 2, 3, 1 }}, {{ 3, 3, 1 }}, {{ 1, 4, 1 }}, {{ 2, 4, 1 }}, {{ 1, 5, 1 }},
                                                               {{ 1, 1, 2 }}, {{ 2, 1, 2 }}, {{ 3, 1, 2 }}, {{ 4, 1, 2 }}, {{ 1, 2, 2 }}, {{ 2, 2, 2 }}, {{ 3, 2, 2 }}, {{ 1, 3, 2 }}, {{ 2, 3, 2 }}, {{ 1, 4, 2 }},
                                                               {{ 1, 1, 3 }}, {{ 2, 1, 3 }}, {{ 3, 1, 3 }}, {{ 1, 2, 3 }}, {{ 2, 2, 3 }}, {{ 1, 3, 3 }},
                                                               {{ 1, 1, 4 }}, {{ 2, 1, 4 }}, {{ 1, 2, 4 }},
                                                               {{ 1, 1, 5 }} }} ), 9, 1 );

  testCellIterator( std::vector< std::array< uint_t, 3 > >( {{ {{ 2, 2, 2 }}, {{ 3, 2, 2 }}, {{ 4, 2, 2 }}, {{ 2, 3, 2 }}, {{ 3, 3, 2 }}, {{ 2, 4, 2 }},
                                                               {{ 2, 2, 3 }}, {{ 3, 2, 3 }}, {{ 2, 3, 3 }}, {{ 2, 2, 4 }} }} ), 9, 2 );

  // no offset to center

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 0, 0, 0 }}, {{ 1, 0, 0 }}, {{ 2, 0, 0 }}, {{ 3, 0, 0 }}, {{ 4, 0, 0 }}, {{ 0, 1, 0 }}, {{ 1, 1, 0 }}, {{ 2, 1, 0 }}, {{ 3, 1, 0 }}, {{ 0, 2, 0 }}, {{ 1, 2, 0 }}, {{ 2, 2, 0 }}, {{ 0, 3, 0 }}, {{ 1, 3, 0 }}, {{ 0, 4, 0 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 0, 1, 2}} ) );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 0, 0, 0 }}, {{ 0, 1, 0 }}, {{ 0, 2, 0 }}, {{ 0, 3, 0 }}, {{ 0, 4, 0 }}, {{ 1, 0, 0 }}, {{ 1, 1, 0 }}, {{ 1, 2, 0 }}, {{ 1, 3, 0 }}, {{ 2, 0, 0 }}, {{ 2, 1, 0 }}, {{ 2, 2, 0 }}, {{ 3, 0, 0 }}, {{ 3, 1, 0 }}, {{ 4, 0, 0 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 0, 2, 1}} ) );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 4, 0, 0 }}, {{ 3, 0, 0 }}, {{ 2, 0, 0 }}, {{ 1, 0, 0 }}, {{ 0, 0, 0 }}, {{ 3, 1, 0 }}, {{ 2, 1, 0 }}, {{ 1, 1, 0 }}, {{ 0, 1, 0 }}, {{ 2, 2, 0 }}, {{ 1, 2, 0 }}, {{ 0, 2, 0 }}, {{ 1, 3, 0 }}, {{ 0, 3, 0 }}, {{ 0, 4, 0 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 1, 0, 2}} ) );
  
  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 4, 0, 0 }}, {{ 3, 1, 0 }}, {{ 2, 2, 0 }}, {{ 1, 3, 0 }}, {{ 0, 4, 0 }}, {{ 3, 0, 1 }}, {{ 2, 1, 1 }}, {{ 1, 2, 1 }}, {{ 0, 3, 1 }}, {{ 2, 0, 2 }}, {{ 1, 1, 2 }}, {{ 0, 2, 2 }}, {{ 1, 0, 3 }}, {{ 0, 1, 3 }}, {{ 0, 0, 4 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 1, 2, 3 }} ) );
  
  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 0, 4, 0 }}, {{ 1, 3, 0 }}, {{ 2, 2, 0 }}, {{ 3, 1, 0 }}, {{ 4, 0, 0 }}, {{ 0, 3, 1 }}, {{ 1, 2, 1 }}, {{ 2, 1, 1 }}, {{ 3, 0, 1 }}, {{ 0, 2, 2 }}, {{ 1, 1, 2 }}, {{ 2, 0, 2 }}, {{ 0, 1, 3 }}, {{ 1, 0, 3 }}, {{ 0, 0, 4 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 2, 1, 3 }} ) );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 0, 0, 0 }}, {{ 0, 1, 0 }}, {{ 0, 2, 0 }}, {{ 0, 3, 0 }}, {{ 0, 4, 0 }}, {{ 0, 5, 0 }}, {{ 0, 6, 0 }}, {{ 0, 7, 0 }}, {{ 0, 8, 0 }}, {{ 1, 0, 0 }}, {{ 1, 1, 0 }}, {{ 1, 2, 0 }}, {{ 1, 3, 0 }}, {{ 1, 4, 0 }}, {{ 1, 5, 0 }}, {{ 1, 6, 0 }}, {{ 1, 7, 0 }}, {{ 2, 0, 0 }}, {{ 2, 1, 0 }}, {{ 2, 2, 0 }}, {{ 2, 3, 0 }}, {{ 2, 4, 0 }}, {{ 2, 5, 0 }}, {{ 2, 6, 0 }}, {{ 3, 0, 0 }}, {{ 3, 1, 0 }}, {{ 3, 2, 0 }}, {{ 3, 3, 0 }}, {{ 3, 4, 0 }}, {{ 3, 5, 0 }}, {{ 4, 0, 0 }}, {{ 4, 1, 0 }}, {{ 4, 2, 0 }}, {{ 4, 3, 0 }}, {{ 4, 4, 0 }}, {{ 5, 0, 0 }}, {{ 5, 1, 0 }}, {{ 5, 2, 0 }}, {{ 5, 3, 0 }}, {{ 6, 0, 0 }}, {{ 6, 1, 0 }}, {{ 6, 2, 0 }}, {{ 7, 0, 0 }}, {{ 7, 1, 0 }}, {{ 8, 0, 0 }}
  }} ), 9, std::array< uint_t, 3 >( {{ 0, 2, 1 }} ) );

  // with offset to center

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 0, 0, 1 }}, {{ 1, 0, 1 }}, {{ 2, 0, 1 }}, {{ 3, 0, 1 }}, {{ 0, 1, 1 }}, {{ 1, 1, 1 }}, {{ 2, 1, 1 }}, {{ 0, 2, 1 }}, {{ 1, 2, 1 }}, {{ 0, 3, 1 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 0, 1, 2}} ), 1 );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 0, 0, 1 }}, {{ 0, 1, 1 }}, {{ 0, 2, 1 }}, {{ 0, 3, 1 }}, {{ 1, 0, 1 }}, {{ 1, 1, 1 }}, {{ 1, 2, 1 }}, {{ 2, 0, 1 }}, {{ 2, 1, 1 }}, {{ 3, 0, 1 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 0, 2, 1}} ), 1 );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 3, 0, 1 }}, {{ 2, 0, 1 }}, {{ 1, 0, 1 }}, {{ 0, 0, 1 }}, {{ 2, 1, 1 }}, {{ 1, 1, 1 }}, {{ 0, 1, 1 }}, {{ 1, 2, 1 }}, {{ 0, 2, 1 }}, {{ 0, 3, 1 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 1, 0, 2}} ), 1 );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 3, 0, 0 }}, {{ 2, 1, 0 }}, {{ 1, 2, 0 }}, {{ 0, 3, 0 }}, {{ 2, 0, 1 }}, {{ 1, 1, 1 }}, {{ 0, 2, 1 }}, {{ 1, 0, 2 }}, {{ 0, 1, 2 }}, {{ 0, 0, 3 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 1, 2, 3 }} ), 1 );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 0, 3, 0 }}, {{ 1, 2, 0 }}, {{ 2, 1, 0 }}, {{ 3, 0, 0 }}, {{ 0, 2, 1 }}, {{ 1, 1, 1 }}, {{ 2, 0, 1 }}, {{ 0, 1, 2 }}, {{ 1, 0, 2 }}, {{ 0, 0, 3 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 2, 1, 3 }} ), 1 );

  testCellBorderIterator( std::vector< std::array< uint_t, 3 > >( {{
    {{ 1, 0, 0 }}, {{ 1, 0, 1 }}, {{ 1, 0, 2 }}, {{ 1, 0, 3 }}, {{ 1, 1, 0 }}, {{ 1, 1, 1 }}, {{ 1, 1, 2 }}, {{ 1, 2, 0 }}, {{ 1, 2, 1 }}, {{ 1, 3, 0 }}
  }} ), 5, std::array< uint_t, 3 >( {{ 0, 3, 2 }} ), 1 );
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
