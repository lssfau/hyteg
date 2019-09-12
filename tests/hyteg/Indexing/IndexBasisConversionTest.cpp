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

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"
#include "core/logging/all.h"

#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"

namespace hyteg {

using walberla::uint_t;

static void testBasisConversion( const std::array< uint_t, 4 > testBasis0,
                                 const std::array< uint_t, 4 > testBasis1,
                                 const uint_t & cellWidth,
                                 const uint_t & offsetToCenter )
{
  auto testIterator0 = indexing::CellBoundaryIterator( cellWidth, {{ testBasis0[0], testBasis0[1], testBasis0[2] }}, offsetToCenter );
  auto testIterator1 = indexing::CellBoundaryIterator( cellWidth, {{ testBasis1[0], testBasis1[1], testBasis1[2] }}, offsetToCenter );

  for( auto referenceIdx : indexing::CellBoundaryIterator( cellWidth, {{ 0, 1, 2 }}, offsetToCenter ) )
  {
    auto testIdx0 = *testIterator0;
    auto testIdx1 = *testIterator1;

    // conversion test0 -> reference
    auto idxConvertedToReference = indexing::basisConversion( testIdx0, {{ 0, 1, 2, 3 }}, testBasis0, cellWidth );
    WALBERLA_CHECK_EQUAL( referenceIdx, idxConvertedToReference );

    // conversion reference -> test1
    auto idxConvertedToSecondTestBasis = indexing::basisConversion( referenceIdx, testBasis1, {{ 0, 1, 2, 3 }}, cellWidth );
    WALBERLA_CHECK_EQUAL( testIdx1, idxConvertedToSecondTestBasis );

    testIterator0++;
    testIterator1++;
  }
}

}

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hyteg::testBasisConversion( {{1, 0, 2, 3}}, {{0, 1, 3, 2}}, 5, 0 );
  hyteg::testBasisConversion( {{1, 0, 2, 3}}, {{0, 1, 3, 2}}, 6, 0 );
  hyteg::testBasisConversion( {{1, 3, 2, 0}}, {{1, 3, 2, 0}}, 9, 1 );
  hyteg::testBasisConversion( {{0, 1, 2, 3}}, {{1, 3, 2, 0}}, 35, 1 );
  hyteg::testBasisConversion( {{3, 2, 1, 0}}, {{0, 1, 2, 3}}, 37, 2 );

  return EXIT_SUCCESS;
}