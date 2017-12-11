
#include <iostream>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"
#include "core/logging/all.h"

#include "tinyhhg_core/levelinfo.hpp"

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/VertexDoFIndexing.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/indexing/Optimization.hpp"

namespace hhg {

static void testCommonIndexing()
{
  using namespace indexing;
  using walberla::uint_t;
  using walberla::real_t;

  WALBERLA_LOG_INFO_ON_ROOT( "Index P1      - face, level 3, (3, 3, center): " << vertexdof::macroface::indexFromVertex< 3 >( 3, 3, stencilDirection::VERTEX_C ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Index EdgeDoF - face, level 3, (3, 3, center): " << indexing::edgedof::macroface::indexFromVertex< 3 >( 3, 3, stencilDirection::EDGE_HO_E ) );

  for ( const auto & it : vertexdof::macroface::BorderIterator( 3, FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 1 ) )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it.col() << ", row = " << it.row() << " ( idx = " << vertexdof::macroface::indexFromVertex< 3 >( it.col(), it.row(), stencilDirection::VERTEX_C ) << " ) " );
  }

  WALBERLA_LOG_INFO_ON_ROOT( "FaceIterator:" )
  for ( const auto & it : FaceIterator( 9, 1 ) )
  {
    std::cout << "(" << it.col() << ", " << it.row() << ") -> " << macroFaceIndex< 9 >( it.col(), it.row() ) << "\n";

  }

  WALBERLA_LOG_INFO_ON_ROOT( "P1FaceIterator (inner face), accessing neighboring horizontal edges" );
  for ( const auto & it : vertexdof::macroface::Iterator( 3, 1 ) )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "Inner face, indexFromVertex (horizontal edge west) = " << indexing::edgedof::macroface::indexFromVertex< 3 >( it.col(), it.row(), stencilDirection::EDGE_HO_W ) );
  }

  const uint_t level = 3;
  const uint_t size = macroFaceSize< vertexdof::levelToWidth< level > >();

  real_t * a = new real_t[ size ];
  real_t * b = new real_t[ size ];

  walberla::WcTimingPool timer;


#if 0
  const uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      a[ macroFaceIndex< vertexdof::levelToWidth< level > >(j, i) ] = 0.0001;
      b[ macroFaceIndex< vertexdof::levelToWidth< level > >(j, i) ] = 0.0002;
    }
    --inner_rowsize;
  }

  timer[ "loop" ].start();

  real_t sp = 0.0;
  inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      sp += a[ macroFaceIndex< vertexdof::levelToWidth< level > >(j, i) ] * b[ macroFaceIndex< vertexdof::levelToWidth< level > >(j, i) ];
    }
    --inner_rowsize;
  }

#else

#if 0

  for (uint_t i = 0; i < optimization::unwrapNumRows< vertexdof::levelToWidth< level > >(); ++i) {
    for (uint_t j = 0; j < optimization::unwrapNumCols< vertexdof::levelToWidth< level > >(); ++j) {

      const uint_t actualRow = optimization::unwrapRow< vertexdof::levelToWidth< level > >(j, i);
      const uint_t actualCol = optimization::unwrapCol< vertexdof::levelToWidth< level > >(j, i);

      a[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ] = 0.0001;
      b[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ] = 0.0002;

    }
  }

  timer[ "loop" ].start();

  real_t sp = 0.0;

  for (uint_t i = 0; i < optimization::unwrapNumRows< vertexdof::levelToWidth< level > >(); ++i) {
    for (uint_t j = 0; j < optimization::unwrapNumCols< vertexdof::levelToWidth< level > >(); ++j) {

      const uint_t actualRow = optimization::unwrapRow< vertexdof::levelToWidth< level > >(j, i);
      const uint_t actualCol = optimization::unwrapCol< vertexdof::levelToWidth< level > >(j, i);

      sp += a[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ] * b[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ];
    }
  }

#else

  for ( const auto & it : vertexdof::macroface::Iterator( level, 1 ) )
  {
    a[ macroFaceIndex< vertexdof::levelToWidth< level > >(it.col(), it.row()) ] = 0.0001;
    b[ macroFaceIndex< vertexdof::levelToWidth< level > >(it.col(), it.row()) ] = 0.0002;
  }

  timer[ "loop" ].start();

  real_t sp = 0.0;

  for ( const auto & it : vertexdof::macroface::Iterator( level, 1 ) )
  {
    sp += a[ macroFaceIndex< vertexdof::levelToWidth< level > >(it.col(), it.row()) ] * b[ macroFaceIndex< vertexdof::levelToWidth< level > >(it.col(), it.row()) ];
  }

#endif

#endif

  timer[ "loop" ].end();

  std::cout << sp << std::endl;
  WALBERLA_LOG_INFO_ON_ROOT( "result = " << sp );
  WALBERLA_LOG_INFO_ON_ROOT( timer );

  delete[] a;
  delete[] b;



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
