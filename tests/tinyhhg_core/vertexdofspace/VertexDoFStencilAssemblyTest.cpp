
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

#include "tinyhhg_core/p1functionspace/generated/p1_tet_diffusion.h"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"

#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"

namespace hhg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;

static void testVertexDoFStencilAssembly()
{
  typedef stencilDirection sd;

  const uint_t minLevel = 2;
  const uint_t maxLevel = 12;

  auto storage = PrimitiveStorage::createFromGmshFile( "../../data/meshes/3D/regular_octahedron_8el.msh" );

  p1_tet_diffusion_cell_integral_0_otherwise ufcOperator;

  for ( const auto & it : storage->getVertices() )
  {
    const auto vertex = *it.second;

    if ( !storage->onBoundary( it.first ))
    {
      std::array< std::vector< real_t >, maxLevel + 1 > macroVertexStencilOnLevel;

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
        std::vector< real_t > stencil = P1Elements::P1Elements3D::assembleP1LocalStencil< p1_tet_diffusion_cell_integral_0_otherwise >( storage, vertex, indexing::Index( 0, 0, 0 ), level,
                                                                                                                                         ufcOperator );
        macroVertexStencilOnLevel[level] = stencil;
        real_t rowSum = real_c( 0 );

        for ( uint_t stencilIdx = 0; stencilIdx < stencil.size(); stencilIdx++ )
        {
          rowSum += stencil[stencilIdx];
          WALBERLA_LOG_INFO( "Stencil entry on level " << level << ", idx " << stencilIdx << ": " << stencil[stencilIdx] );

          if ( level > minLevel )
          {
            // Checking if stencil weight scale with h
            WALBERLA_CHECK_FLOAT_EQUAL( 2.0 * stencil[ stencilIdx ], macroVertexStencilOnLevel[ level - 1 ][ stencilIdx ] )
          }
        }

        // Checking system matrix row sum
        WALBERLA_LOG_DEVEL( "Vertex row sum: " << rowSum << "\n" );
        WALBERLA_CHECK_FLOAT_EQUAL( rowSum, 0.0 );
      }
    }
  }

  for ( const auto & it : storage->getEdges() )
  {
    const auto edge = *it.second;

    if ( !storage->onBoundary( it.first ) )
    {
      std::array< std::vector< real_t >, maxLevel + 1 > macroEdgeStencilOnLevel;

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
        std::vector< real_t > stencil = P1Elements::P1Elements3D::assembleP1LocalStencil< p1_tet_diffusion_cell_integral_0_otherwise >( storage, edge, indexing::Index( 1, 0, 0 ), level,
                                                                                                                                         ufcOperator );
        macroEdgeStencilOnLevel[level] = stencil;
        real_t rowSum = real_c( 0 );

        // test weights on edge
        for ( auto dir : std::vector< sd >( { sd::VERTEX_C, sd::VERTEX_W, sd::VERTEX_E } ))
        {
          const uint_t index = vertexdof::macroedge::stencilIndexOnEdge( dir );
          const real_t weight = stencil[index];
          WALBERLA_LOG_INFO( "Stencil entry on level " << level << ", on edge " << stencilDirectionToStr.at( dir ) << ": " << weight );
          rowSum += weight;
          if ( level > minLevel )
          {
            // Checking if stencil weight scale with h
            WALBERLA_CHECK_FLOAT_EQUAL( 2.0 * weight, macroEdgeStencilOnLevel[level - 1][index] )
          }
        }

        // test weights on neighbor faces
        for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
        {
          for ( auto dir : std::vector< sd >( { sd::VERTEX_W, sd::VERTEX_E } ))
          {
            const uint_t index = vertexdof::macroedge::stencilIndexOnNeighborFace( dir, neighborFace );
            const real_t weight = stencil[index];
            WALBERLA_LOG_INFO( "Stencil entry on level " << level << ", on face " << neighborFace << " " << stencilDirectionToStr.at( dir ) << ": " << weight );
            rowSum += weight;
            if ( level > minLevel )
            {
              // Checking if stencil weight scale with h
              WALBERLA_CHECK_FLOAT_EQUAL( 2.0 * weight, macroEdgeStencilOnLevel[level - 1][index] )
            }
          }
        }

        // test weights on neighbor cells
        for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
        {
          const uint_t index = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces());
          const real_t weight = stencil[index];
          WALBERLA_LOG_INFO( "Stencil entry on level " << level << ", on cell " << neighborCell << ": " << weight );
          rowSum += weight;
          if ( level > minLevel )
          {
            // Checking if stencil weight scale with h
            WALBERLA_CHECK_FLOAT_EQUAL( 2.0 * weight, macroEdgeStencilOnLevel[level - 1][index] )
          }
        }

        // Checking system matrix row sum
        WALBERLA_LOG_INFO( "Row sum = " << rowSum );
        WALBERLA_CHECK_FLOAT_EQUAL( rowSum, 0.0 );

        // Checking stencil weight symmetry
        WALBERLA_CHECK_FLOAT_EQUAL( stencil[vertexdof::macroedge::stencilIndexOnEdge( sd::VERTEX_W )], stencil[vertexdof::macroedge::stencilIndexOnEdge( sd::VERTEX_E )] );
        WALBERLA_LOG_INFO_ON_ROOT( "" );
      }
    }
  }

  for ( const auto & it : storage->getFaces() )
  {
    const auto face = *it.second;
    if ( face.getNumNeighborCells() == 2 )
    {
      std::array< vertexdof::macroface::StencilMap_T, maxLevel + 1 > macroFaceStencilOnLevel;

      for ( uint_t level = minLevel; level <= maxLevel; level++ ) {

        for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
        {
          auto neighborCell = storage->getCell( face.neighborCells().at( neighborCellID ) );
          auto vertexAssemblyIndexInCell =
            vertexdof::macroface::getIndexInNeighboringMacroCell( {1, 1, 0}, face, neighborCellID, *storage, level );
          macroFaceStencilOnLevel[level][neighborCellID] = P1Elements::P1Elements3D::assembleP1LocalStencilNew(
            storage, *neighborCell, vertexAssemblyIndexInCell, level, ufcOperator );
        }

        real_t rowSum = real_c( 0 );

        for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
        {
          for ( const auto itStencil : macroFaceStencilOnLevel[level][neighborCellID] )
          {
            const auto direction = itStencil.first;
            const auto weight = itStencil.second;
            rowSum += weight;
            WALBERLA_LOG_INFO( "Stencil entry on level " << level << ", direction " << direction << ": " << weight );

            if ( level > minLevel )
            {
              // Checking if stencil weight scale with h
              WALBERLA_CHECK_FLOAT_EQUAL( 2.0 * weight, macroFaceStencilOnLevel[level - 1][neighborCellID][direction] )
            }
          }
        }
        // Checking system matrix row sum
        WALBERLA_CHECK_FLOAT_EQUAL( rowSum, 0.0 );
      }
    }
  }

  for ( const auto & cell : storage->getCells() )
  {
    std::array< std::map< stencilDirection, real_t >, maxLevel + 1 > macroCellStencilOnLevel;

    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      std::map< stencilDirection, real_t > stencil = P1Elements::P1Elements3D::assembleP1LocalStencil< p1_tet_diffusion_cell_integral_0_otherwise >( storage, *cell.second, indexing::Index( 1, 1, 1 ), level, ufcOperator );
      macroCellStencilOnLevel[ level ] = stencil;
      real_t rowSum = real_c( 0 );

      for ( const auto itStencil : stencil )
      {
        const auto direction = itStencil.first;
        const auto weight = itStencil.second;

        rowSum += weight;
        // WALBERLA_LOG_INFO( "Stencil entry on level " << level << ", idx " << stencilIdx << ": " << stencil[stencilIdx] );

        if ( level > minLevel )
        {
          // Checking if stencil weight scale with h
          WALBERLA_CHECK_FLOAT_EQUAL( 2.0 * weight, macroCellStencilOnLevel[ level - 1 ][ direction ] )
        }
      }

      // Checking system matrix row sum
      WALBERLA_CHECK_FLOAT_EQUAL( rowSum, 0.0 );

      // Checking stencil weight symmetry
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[sd::VERTEX_W], stencil[sd::VERTEX_E] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[sd::VERTEX_N], stencil[sd::VERTEX_S] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[sd::VERTEX_NW], stencil[sd::VERTEX_SE] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[sd::VERTEX_TC], stencil[sd::VERTEX_BC] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[sd::VERTEX_TW], stencil[sd::VERTEX_BE] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[sd::VERTEX_TS], stencil[sd::VERTEX_BN] );
      WALBERLA_CHECK_FLOAT_EQUAL( stencil[sd::VERTEX_TSE], stencil[sd::VERTEX_BNW] );
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
