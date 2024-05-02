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
#include "core/timing/all.h"
#include "core/DataTypes.h"

#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;

namespace hyteg {

static void testEdgeDoFToVertexDoFOperator()
{
  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  MeshInfo mesh = MeshInfo::fromGmshFile( "../../meshes/quad_4el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto vertex_dst      = std::make_shared< P1Function< real_t > >     ( "vertex_dst", storage, minLevel, maxLevel );
  auto vertex_expected = std::make_shared< P1Function< real_t > >     ( "vertex_expected", storage, minLevel, maxLevel );
  auto edge_src        = std::make_shared< EdgeDoFFunction< real_t > >( "edge_src", storage, minLevel, maxLevel );

  GenericEdgeDoFToVertexDoFOperator edgeToVertexOperator( storage, minLevel, maxLevel );

  // Test setup:
  // Writing different values to different kind of Edge DoF types (horizontal, vertical, diagonal)
  // Setting specific stencil weights other than 0.0.

  // Stencils

  const real_t macroVertexStencilValue = 1.0;

  const real_t macroEdgeHorizontalStencilValue = 1.1;
  const real_t macroEdgeVerticalStencilValue   = 1.2;
  const real_t macroEdgeDiagonalStencilValue   = 1.3;

  const real_t macroFaceHorizontalStencilValue = 1.4;
  const real_t macroFaceVerticalStencilValue   = 1.5;
  const real_t macroFaceDiagonalStencilValue   = 1.6;

  for ( const auto & it : storage->getVertices() )
  {
    auto vertex = it.second;
    auto stencil = vertex->getData( edgeToVertexOperator.getVertexStencilID() );

    for ( uint_t stencilIdx = 0; stencilIdx < stencil->getSize( maxLevel ); stencilIdx++ )
    {
      stencil->getPointer( maxLevel )[ stencilIdx ] = macroVertexStencilValue;
    }
  }

  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    auto stencil = edge->getData( edgeToVertexOperator.getEdgeStencilID() )->getPointer( maxLevel );

    for ( const auto & stencilDir : edgedof::macroedge::neighborsOnEdgeFromVertex )
    {
      stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroEdgeHorizontalStencilValue;
    }

    for ( const auto & stencilDir : edgedof::macroedge::neighborsOnSouthFaceFromVertex )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroEdgeDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroEdgeHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroEdgeVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }

    if ( edge->getNumNeighborFaces() == 2 )
    {
      for ( const auto & stencilDir : edgedof::macroface::neighborsFromVertex )
      {
        if ( isDiagonalEdge( stencilDir ) )
        {
          stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroEdgeDiagonalStencilValue;
        }
        else if ( isHorizontalEdge( stencilDir ) )
        {
          stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroEdgeHorizontalStencilValue;
        }
        else if ( isVerticalEdge( stencilDir ) )
        {
          stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroEdgeVerticalStencilValue;
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
    auto stencil = face->getData( edgeToVertexOperator.getFaceStencilID() )->getPointer( maxLevel );

    for ( const auto & stencilDir : edgedof::macroface::neighborsFromVertex )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroFaceDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroFaceHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVertex( stencilDir ) ] = macroFaceVerticalStencilValue;
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

  edgeToVertexOperator.apply( *edge_src, *vertex_dst, maxLevel, DoFType::All, UpdateType::Replace );

  // Check macro vertices
  for ( const auto & it : storage->getVertices() )
  {
    auto vertex = it.second;
    auto vertexFunction = vertex->getData( vertex_dst->getVertexDataID() );
    WALBERLA_CHECK_FLOAT_EQUAL( vertexFunction->getPointer( maxLevel )[0], edgeSrcValue * macroVertexStencilValue * real_c( vertex->getNumNeighborEdges() + vertex->getNumNeighborFaces() ) );
  }

  // Check macro edges
  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    auto edgeFunction = edge->getData( vertex_dst->getEdgeDataID() );
    for ( const auto & idxIt : vertexdof::macroedge::Iterator( maxLevel, 1 ) )
    {
      auto ptr = edgeFunction->getPointer( maxLevel );
      auto idx = vertexdof::macroedge::index( maxLevel, idxIt.x());

      const real_t expectedValue = edgeSrcValue *
        ( 2.0 * real_c( edge->getNumNeighborFaces() ) * macroEdgeDiagonalStencilValue
        + ( 2.0 + real_c( edge->getNumNeighborFaces() ) ) * macroEdgeHorizontalStencilValue
        + 2.0 * real_c( edge->getNumNeighborFaces() ) * macroEdgeVerticalStencilValue );

      WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
    }
  }


  // Check macro faces
  for ( const auto & it : storage->getFaces() )
  {
    auto face = it.second;
    auto faceFunction = face->getData( vertex_dst->getFaceDataID() );
    for ( const auto & idxIt : vertexdof::macroface::Iterator( maxLevel, 1 ) )
    {
      auto ptr = faceFunction->getPointer( maxLevel );
      auto idx = vertexdof::macroface::index( maxLevel, idxIt.x(), idxIt.y());

      const real_t expectedValue = edgeSrcValue * 4.0 * ( macroFaceHorizontalStencilValue + macroFaceDiagonalStencilValue + macroFaceVerticalStencilValue );

      WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
    }
  }
}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testEdgeDoFToVertexDoFOperator();

   return EXIT_SUCCESS;
}
