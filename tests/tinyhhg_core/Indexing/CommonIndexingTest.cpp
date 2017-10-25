
#include <iostream>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"

#include "tinyhhg_core/tinyhhg.hpp"

#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/VertexDoFIndexing.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

namespace hhg {

static void testCommonIndexing()
{
  using namespace indexing;
  using walberla::uint_t;
  using walberla::real_t;

  WALBERLA_LOG_INFO_ON_ROOT( "Index P1      - face, level 3, (3, 3, center): " << vertexdof::macroface::indexFromVertex< 3 >( 3, 3, stencilDirection::VERTEX_C ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Index EdgeDoF - face, level 3, (3, 3, center): " << EdgeDoFFaceIndexFromVertex< 3 >( 3, 3, stencilDirection::EDGE_HO_C ) );

  for ( const auto & it : vertexdof::macroface::BorderIterator< 3 >( FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 1 ) )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it.col() << ", row = " << it.row() << " ( idx = " << vertexdof::macroface::indexFromVertex< 3 >( it.col(), it.row(), stencilDirection::VERTEX_C ) << " ) " );
  }

  const uint_t level = 3;
  const uint_t size = macroFaceSize< vertexdof::levelToWidth< level > >();

  WALBERLA_LOG_INFO_ON_ROOT( size );

  for ( uint_t i = 0; i < 2; i++ )
  {
    for ( uint_t j = 0; j < 5; j++ )
    {
      const uint_t actualRow = unwrapRow<4>(j, i);
      const uint_t actualCol = unwrapCol<4>(j, i);

      std::cout << macroFaceIndex<4>(actualCol, actualRow) << " ";
    }
    std::cout << "\n";
  }

  for ( uint_t i = 0; i < 3; i++ )
  {
    for ( uint_t j = 0; j < 5; j++ )
    {
      const uint_t actualRow = unwrapRow<5>(j, i);
      const uint_t actualCol = unwrapCol<5>(j, i);

      std::cout << macroFaceIndex<5>(actualCol, actualRow) << " ";
    }
    std::cout << "\n";
  }


  real_t * a = new real_t[ size ];
  real_t * b = new real_t[ size ];

  walberla::WcTimingPool timer;


#if 0
  const uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 0; i < rowsize; ++i) {
    for (uint_t j = 0; j < inner_rowsize; ++j) {
      a[ macroFaceIndex< levelToWidthVertexDoF< level > >(j, i) ] = 0.0001;
      b[ macroFaceIndex< levelToWidthVertexDoF< level > >(j, i) ] = 0.0002;
    }
    --inner_rowsize;
  }

  timer[ "loop" ].start();

  real_t sp = 0.0;
  inner_rowsize = rowsize;

  for (uint_t i = 0; i < rowsize; ++i) {
    for (uint_t j = 0; j < inner_rowsize; ++j) {
      sp += a[ macroFaceIndex< levelToWidthVertexDoF< level > >(j, i) ] * b[ macroFaceIndex< levelToWidthVertexDoF< level > >(i, j) ];
    }
    --inner_rowsize;
  }

#else

  for (uint_t i = 0; i < unwrapNumRows< vertexdof::levelToWidth< level > >(); ++i) {
    for (uint_t j = 0; j < unwrapNumCols< vertexdof::levelToWidth< level > >(); ++j) {

      const uint_t actualRow = unwrapRow< vertexdof::levelToWidth< level > >(j, i);
      const uint_t actualCol = unwrapCol< vertexdof::levelToWidth< level > >(j, i);

      a[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ] = 0.0001;
      b[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ] = 0.0002;

    }
  }

  timer[ "loop" ].start();

  real_t sp = 0.0;

  for (uint_t i = 0; i < unwrapNumRows< vertexdof::levelToWidth< level > >(); ++i) {
    for (uint_t j = 0; j < unwrapNumCols< vertexdof::levelToWidth< level > >(); ++j) {

      const uint_t actualRow = unwrapRow< vertexdof::levelToWidth< level > >(j, i);
      const uint_t actualCol = unwrapCol< vertexdof::levelToWidth< level > >(j, i);

      sp += a[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ] * b[ macroFaceIndex< vertexdof::levelToWidth< level > >(actualCol, actualRow) ];
    }
  }

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
