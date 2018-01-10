
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOperator.hpp"


namespace hhg {

static void testEdgeDoFToEdgeDoFOperator()
{
  const uint_t maxLevel = 4;

  MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto edge_dst = std::make_shared< EdgeDoFFunction< real_t > >( "edge_dst", storage, maxLevel, maxLevel );
  auto edge_src = std::make_shared< EdgeDoFFunction< real_t > >( "edge_src", storage, maxLevel, maxLevel );

  EdgeDoFOperator edgeToEdgeOperator( storage, maxLevel, maxLevel );

  // Test setup:
  // Writing different values to different kind of Edge DoF types (horizontal, vertical, diagonal)
  // Setting specific stencil weights other than 0.0.

  // Stencils

  const real_t macroEdgeHorizontalStencilValue = 1.1;
  const real_t macroEdgeVerticalStencilValue   = 1.2;
  const real_t macroEdgeDiagonalStencilValue   = 1.3;

  const real_t macroFaceHorizontalStencilValue = 1.4;
  const real_t macroFaceVerticalStencilValue   = 1.5;
  const real_t macroFaceDiagonalStencilValue   = 1.6;


  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    auto stencil = edge->getData( edgeToEdgeOperator.getEdgeStencilID_() )->getPointer( maxLevel );

    for ( const auto & stencilDir : indexing::edgedof::macroedge::neighborsOnEdgeFromHorizontalEdge )
    {
      stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeHorizontalStencilValue;
    }

    for ( const auto & stencilDir : indexing::edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }

    if ( edge->getNumNeighborFaces() == 2 )
    {
      for ( const auto & stencilDir : indexing::edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge )
      {
        if ( isDiagonalEdge( stencilDir ) )
        {
          stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeDiagonalStencilValue;
        }
        else if ( isHorizontalEdge( stencilDir ) )
        {
          stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeHorizontalStencilValue;
        }
        else if ( isVerticalEdge( stencilDir ) )
        {
          stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeVerticalStencilValue;
        }
        else
        {
          WALBERLA_ABORT( "invalid stencil direction" );
        }
      }
    }
  }

  for ( const auto & it : storage->getFaces() )
  {
    auto face = it.second;
    auto stencil = face->getData( edgeToEdgeOperator.getFaceStencilID_() )->getPointer( maxLevel );

    for ( const auto & stencilDir : indexing::edgedof::macroface::neighborsFromHorizontalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroFaceDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroFaceHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroFaceVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }

    for ( const auto & stencilDir : indexing::edgedof::macroface::neighborsFromDiagonalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromDiagonalEdge( stencilDir ) ] = macroFaceDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromDiagonalEdge( stencilDir ) ] = macroFaceHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromDiagonalEdge( stencilDir ) ] = macroFaceVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }

    for ( const auto & stencilDir : indexing::edgedof::macroface::neighborsFromVerticalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromVerticalEdge( stencilDir ) ] = macroFaceDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromVerticalEdge( stencilDir ) ] = macroFaceHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ indexing::edgedof::stencilIndexFromVerticalEdge( stencilDir ) ] = macroFaceVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }
  }

  // Interpolate src function
  const real_t edgeSrcValue = 0.5;

  std::function< real_t ( const Point3D & ) > initEdgeSrc = [ edgeSrcValue ]( const Point3D & ) -> real_t { return edgeSrcValue; };
  edge_src->interpolate( initEdgeSrc, maxLevel, DoFType::All );

  auto communicator = edge_src->getCommunicator( maxLevel );

  // Pull all halos
  communicator->communicate< Face, Edge >();
  communicator->communicate< Edge, Vertex >();
  communicator->communicate< Edge, Face >();

  edgeToEdgeOperator.apply( *edge_src, *edge_dst, maxLevel, DoFType::All, UpdateType::Replace );

  // Check macro edges
  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    auto edgeFunction = edge->getData( edge_dst->getEdgeDataID() );

    for ( const auto & idxIt : indexing::edgedof::macroedge::Iterator( maxLevel ) )
    {
      auto ptr = edgeFunction->getPointer( maxLevel );
      auto idx = indexing::edgedof::macroedge::indexFromHorizontalEdge< maxLevel >( idxIt.col(), stencilDirection::EDGE_HO_C );

      const real_t expectedValue = edgeSrcValue * ( macroEdgeHorizontalStencilValue + real_c( edge->getNumNeighborFaces() ) * ( macroEdgeDiagonalStencilValue + macroEdgeVerticalStencilValue ) );
      WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
    }
  }


  // Check macro faces
  for ( const auto & it : storage->getFaces() )
  {
    auto face = it.second;
    auto faceFunction = face->getData( edge_dst->getFaceDataID() );

    auto ptr = faceFunction->getPointer( maxLevel );

    for ( const auto & idxIt : indexing::edgedof::macroface::Iterator( maxLevel ) )
    {
      if ( idxIt.col() != 0 )
      {
        const auto idx = indexing::edgedof::macroface::indexFromVerticalEdge< maxLevel >( idxIt.col(), idxIt.row(), stencilDirection::EDGE_VE_C );

        const real_t expectedValue = edgeSrcValue * ( macroFaceVerticalStencilValue + 2.0 * ( macroFaceDiagonalStencilValue + macroFaceHorizontalStencilValue ) );
        WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
      }

      if ( idxIt.row() != 0 )
      {
        const auto idx = indexing::edgedof::macroface::indexFromHorizontalEdge< maxLevel >( idxIt.col(), idxIt.row(), stencilDirection::EDGE_HO_C );

        const real_t expectedValue = edgeSrcValue * ( macroFaceHorizontalStencilValue + 2.0 * ( macroFaceDiagonalStencilValue + macroFaceVerticalStencilValue ) );
        WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
      }

      if ( idxIt.row() + idxIt.col() < levelinfo::num_microedges_per_edge( maxLevel ) - 1 )
      {
        const auto idx = indexing::edgedof::macroface::indexFromDiagonalEdge< maxLevel >( idxIt.col(), idxIt.row(), stencilDirection::EDGE_DI_C );

        const real_t expectedValue = edgeSrcValue * ( macroFaceDiagonalStencilValue + 2.0 * ( macroFaceHorizontalStencilValue + macroFaceVerticalStencilValue ) );
        WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
      }
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
   hhg::testEdgeDoFToEdgeDoFOperator();

   return EXIT_SUCCESS;
}
