/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <hyteg/edgedofspace/EdgeDoFIndexing.hpp>
#include <hyteg/p1functionspace/VertexDoFIndexing.hpp>
#include <iostream>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"
#include "core/logging/all.h"

#include "hyteg/Levelinfo.hpp"

#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/Optimization.hpp"

namespace hyteg {

static void testFaceBorderIterator( const std::vector< std::array< uint_t, 2 > > & expectedValues, const indexing::FaceBoundaryDirection & faceBorderDirection,
                                    const uint_t & width, const uint_t & offsetToCenter, const uint_t & offsetFromVertices )
{
  std::vector< std::array< idx_t, 2 > > iteratorResult;

  for ( const auto & it : indexing::FaceBoundaryIterator( width, faceBorderDirection, offsetToCenter, offsetFromVertices ) )
  {
    iteratorResult.push_back( {{ it.x(), it.y() }} );
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
  std::vector< std::array< idx_t, 3 > > iteratorResult;

  for ( const auto & it : indexing::CellIterator( width, offsetToCenter ) )
  {
    iteratorResult.push_back( {{ it.x(), it.y(), it.z() }} );
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
  std::vector< std::array< idx_t, 3 > > iteratorResult;

  for ( const auto & it : indexing::CellBoundaryIterator( width, vertices, offsetToCenter ) )
  {
    iteratorResult.push_back( {{ it.x(), it.y(), it.z() }} );
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
    WALBERLA_CHECK_EQUAL( testCol++, it.x() );
    WALBERLA_CHECK_EQUAL( 0,         it.y() );
  }
  WALBERLA_CHECK_EQUAL( testCol, 8 );

  testCol = 1;
  for ( const auto & it : indexing::EdgeIterator( 8, 1 ) )
  {
    WALBERLA_CHECK_EQUAL( testCol++, it.x() );
    WALBERLA_CHECK_EQUAL( 0,         it.y() );
  }
  WALBERLA_CHECK_EQUAL( testCol, 7 );

  // macro face

  testCol = 0;
  testRow = 0;
  for ( const auto & it : indexing::FaceIterator( 4, 0 ) )
  {
    WALBERLA_CHECK_EQUAL( testCol++, it.x() );
    WALBERLA_CHECK_EQUAL( testRow  , it.y() );

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
    WALBERLA_CHECK_EQUAL( it.x(), 1 );
    WALBERLA_CHECK_EQUAL( it.y(), 1 );
  }
  WALBERLA_CHECK_EQUAL( testRow, 2 );
  WALBERLA_CHECK_EQUAL( testCol, 2 );

  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{0, 0}}, {{1, 0}}, {{2, 0}}, {{3, 0}} }} ),           indexing::FaceBoundaryDirection::BOTTOM_LEFT_TO_RIGHT, 4, 0, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{2, 0}}, {{1, 0}} }} ),                               indexing::FaceBoundaryDirection::BOTTOM_RIGHT_TO_LEFT, 4, 0, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{2, 1}}, {{1, 1}}, {{0, 1}} }} ),                     indexing::FaceBoundaryDirection::BOTTOM_RIGHT_TO_LEFT, 4, 1, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{1, 1}} }} ),                                         indexing::FaceBoundaryDirection::BOTTOM_RIGHT_TO_LEFT, 4, 1, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{0, 0}}, {{0, 1}}, {{0, 2}}, {{0, 3}} }} ),           indexing::FaceBoundaryDirection::LEFT_BOTTOM_TO_TOP, 4, 0, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{1, 1}} }} ),                                         indexing::FaceBoundaryDirection::LEFT_BOTTOM_TO_TOP, 4, 1, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{1, 2}}, {{1, 1}}, {{1, 0}} }} ),                     indexing::FaceBoundaryDirection::LEFT_TOP_TO_BOTTOM, 4, 1, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{1, 1}} }} ),                                         indexing::FaceBoundaryDirection::LEFT_TOP_TO_BOTTOM, 4, 1, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{4, 0}}, {{3, 1}}, {{2, 2}}, {{1, 3}}, {{0, 4}} }} ), indexing::FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP, 5, 0, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{3, 1}}, {{2, 2}}, {{1, 3}} }} ),                     indexing::FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP, 5, 0, 1 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{3, 0}}, {{2, 1}}, {{1, 2}}, {{0, 3}} }} ),           indexing::FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP, 5, 1, 0 );
  testFaceBorderIterator( std::vector< std::array< uint_t, 2 > >( {{ {{1, 2}}, {{2, 1}} }} ),                               indexing::FaceBoundaryDirection::DIAGONAL_TOP_TO_BOTTOM, 5, 1, 1 );

  // macro cell

  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 1 ), 1 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 2 ), 4 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 3 ), 10 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 4 ), 20 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 5 ), 35 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 6 ), 56 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 7 ), 84 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 8 ), 120 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 9 ), 165 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellSize( 10 ), 220 );

  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex( 5, 0, 0, 0 ), 0 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex( 5, 2, 0, 0 ), 2 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex( 5, 1, 3, 0 ), 13 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex( 5, 1, 1, 1 ), 20 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex( 5, 1, 1, 2 ), 29 );
  WALBERLA_CHECK_EQUAL( indexing::macroCellIndex( 5, 0, 0, 4 ), 34 );

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

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testCommonIndexing();

   return EXIT_SUCCESS;
}
