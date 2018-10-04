
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"
#include "core/logging/all.h"

#include "tinyhhg_core/indexing/DistanceCoordinateSystem.hpp"
#include "tinyhhg_core/indexing/MacroCellIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

namespace hhg {

using walberla::uint_t;

static void convertAndTest( edgedof::EdgeDoFOrientation src, uint_t v0, uint_t v1, uint_t v2, edgedof::EdgeDoFOrientation expectedResult )
{
  // vx are the local vertex IDs of the face in the macro-cell
  // TODO: call function here
  WALBERLA_CHECK_EQUAL( edgedof::convertEdgeDoFOrientation( src, v0, v1, v2 ), expectedResult, src << " " << v0 << " " << v1 << " " << v2 );
}

static void testEdgeOrientationConversion()
{
  for ( auto edgedofOrientation : edgedof::allEdgeDoFOrientations )
  {
    convertAndTest( edgedofOrientation, 0, 1, 2, edgedofOrientation );
  }

  for ( uint_t v0 = 0; v0 < 4; v0++ )
    for ( uint_t v1 = 0; v1 < 4; v1++ )
      for ( uint_t v2 = 0; v2 < 4; v2++ )
        convertAndTest( edgedof::EdgeDoFOrientation::XYZ, v0, v1, v2, edgedof::EdgeDoFOrientation::XYZ );

  convertAndTest( edgedof::EdgeDoFOrientation::Y, 1, 0, 2, edgedof::EdgeDoFOrientation::XY );
  convertAndTest( edgedof::EdgeDoFOrientation::Y, 1, 2, 0, edgedof::EdgeDoFOrientation::X );

  convertAndTest( edgedof::EdgeDoFOrientation::XY, 1, 2, 0, edgedof::EdgeDoFOrientation::Y );

  convertAndTest( edgedof::EdgeDoFOrientation::Z, 1, 2, 0, edgedof::EdgeDoFOrientation::XZ );
  convertAndTest( edgedof::EdgeDoFOrientation::Z, 1, 2, 3, edgedof::EdgeDoFOrientation::X );

  convertAndTest( edgedof::EdgeDoFOrientation::XZ, 1, 2, 3, edgedof::EdgeDoFOrientation::Y );
}

}

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  hhg::testEdgeOrientationConversion();

  return EXIT_SUCCESS;
}