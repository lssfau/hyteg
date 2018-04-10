
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

#include "tinyhhg_core/p1functionspace/generated/p1_tet_diffusion.h"

namespace hhg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;

static void testVertexDoFStencilAssembly()
{
  typedef stencilDirection sd;

  const uint_t minLevel = 2;
  const uint_t maxLevel = 12;

  auto storage = PrimitiveStorage::createFromGmshFile( "../../data/meshes/3D/tet_1el.msh" );

  p1_tet_diffusion_cell_integral_0_otherwise ufcOperator;

  std::array< std::array< real_t, 15 >, maxLevel + 1 > stencilOnLevel;

  for ( uint_t level = minLevel; level <= maxLevel; level++ )
  {
    for ( const auto & cell : storage->getCells())
    {
      std::array< real_t, 15 > stencil = P1Elements::CellVertexDoF::assembleP1LocalStencil< p1_tet_diffusion_cell_integral_0_otherwise >( *cell.second, level, ufcOperator );
      stencilOnLevel[ level ] = stencil;
      real_t rowSum = real_c( 0 );

      for ( uint_t stencilIdx = 0; stencilIdx < 15; stencilIdx++ )
      {
        rowSum += stencil[stencilIdx];
        // WALBERLA_LOG_INFO( "Stencil entry on level " << level << ", idx " << stencilIdx << ": " << stencil[stencilIdx] );

        if ( level > minLevel )
        {
          // Checking if stencil weight scale with h
          WALBERLA_CHECK_FLOAT_EQUAL( 2.0 * stencil[ stencilIdx ], stencilOnLevel[ level - 1 ][ stencilIdx ] )
        }
      }

      // Checking system matrix row sum
      WALBERLA_CHECK_FLOAT_EQUAL( rowSum, 0.0 );

      // Checking stencil weight symmetry
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_W )], stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_E )] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_N )], stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_S )] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_NW )], stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_SE )] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_BC )], stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_FC )] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_BW )], stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_FE )] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_BS )], stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_FN )] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_BSW )], stencil[vertexdof::stencilIndexFromVertex( sd::VERTEX_FNE )] );
    }
  }
}

} // namespace hhg


int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testVertexDoFStencilAssembly();

  return EXIT_SUCCESS;
}
