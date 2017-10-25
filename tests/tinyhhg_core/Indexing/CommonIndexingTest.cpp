
#include <iostream>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"

#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/VertexDoFIndexing.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

namespace hhg {

static void testCommonIndexing()
{
  using namespace indexing;
  using walberla::uint_t;

  WALBERLA_LOG_INFO_ON_ROOT( "Index P1      - face, level 3, (3, 3, center): " << VertexDoFOnMacroFaceIndexFromVertex< 3 >( 3, 3, stencilDirection::VERTEX_C ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Index EdgeDoF - face, level 3, (3, 3, center): " << EdgeDoFFaceIndexFromVertex< 3 >( 3, 3, stencilDirection::EDGE_HO_C ) );

  for ( const auto & it : VertexDoFFaceBorderIterator< 3 >( FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 1 ) )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it.col() << ", row = " << it.row() << " ( idx = " << VertexDoFOnMacroFaceIndexFromVertex< 3 >( it.col(), it.row(), stencilDirection::VERTEX_C ) << " ) " );
  }


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
