
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"
#include "core/logging/all.h"

#include "tinyhhg_core/indexing/DistanceCoordinateSystem.hpp"
#include "tinyhhg_core/indexing/MacroCellIndexing.hpp"

namespace hhg {

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
  hhg::testBasisConversion( {{1, 0, 2, 3}}, {{0, 1, 3, 2}}, 5, 0 );
  hhg::testBasisConversion( {{1, 0, 2, 3}}, {{0, 1, 3, 2}}, 6, 0 );
  hhg::testBasisConversion( {{1, 3, 2, 0}}, {{1, 3, 2, 0}}, 9, 1 );
  hhg::testBasisConversion( {{0, 1, 2, 3}}, {{1, 3, 2, 0}}, 35, 1 );
  hhg::testBasisConversion( {{3, 2, 1, 0}}, {{0, 1, 2, 3}}, 37, 2 );

  return EXIT_SUCCESS;
}