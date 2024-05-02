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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constant_stencil_operator/EdgeDoFOperator.hpp"

using walberla::real_c;

namespace hyteg {

static void testEdgeDoFToEdgeDoFOperator()
{
  const uint_t maxLevel = 4;

  MeshInfo mesh = MeshInfo::fromGmshFile( "../../meshes/quad_4el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto edge_dst = std::make_shared< EdgeDoFFunction< real_t > >( "edge_dst", storage, maxLevel, maxLevel );
  auto edge_src = std::make_shared< EdgeDoFFunction< real_t > >( "edge_src", storage, maxLevel, maxLevel );

  GenericEdgeDoFOperator edgeToEdgeOperator( storage, maxLevel, maxLevel );

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

    auto stencil = edge->getData( edgeToEdgeOperator.getEdgeStencilID() )->getPointer( maxLevel );

    for ( const auto & stencilDir : edgedof::macroedge::neighborsOnEdgeFromHorizontalEdge )
    {
      stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeHorizontalStencilValue;
    }

    for ( const auto & stencilDir : edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }

    if ( edge->getNumNeighborFaces() == 2 )
    {
      for ( const auto & stencilDir : edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge )
      {
        if ( isDiagonalEdge( stencilDir ) )
        {
          stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeDiagonalStencilValue;
        }
        else if ( isHorizontalEdge( stencilDir ) )
        {
          stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeHorizontalStencilValue;
        }
        else if ( isVerticalEdge( stencilDir ) )
        {
          stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroEdgeVerticalStencilValue;
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
    auto stencil = face->getData( edgeToEdgeOperator.getFaceStencilID() )->getPointer( maxLevel );

    for ( const auto & stencilDir : edgedof::macroface::neighborsFromHorizontalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroFaceDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroFaceHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromHorizontalEdge( stencilDir ) ] = macroFaceVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }

    for ( const auto & stencilDir : edgedof::macroface::neighborsFromDiagonalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromDiagonalEdge( stencilDir ) ] = macroFaceDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromDiagonalEdge( stencilDir ) ] = macroFaceHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromDiagonalEdge( stencilDir ) ] = macroFaceVerticalStencilValue;
      }
      else
      {
        WALBERLA_ABORT( "invalid stencil direction" );
      }
    }

    for ( const auto & stencilDir : edgedof::macroface::neighborsFromVerticalEdge )
    {
      if ( isDiagonalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVerticalEdge( stencilDir ) ] = macroFaceDiagonalStencilValue;
      }
      else if ( isHorizontalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVerticalEdge( stencilDir ) ] = macroFaceHorizontalStencilValue;
      }
      else if ( isVerticalEdge( stencilDir ) )
      {
        stencil[ edgedof::stencilIndexFromVerticalEdge( stencilDir ) ] = macroFaceVerticalStencilValue;
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

  edgeToEdgeOperator.apply( *edge_src, *edge_dst, maxLevel, DoFType::All, UpdateType::Replace );

  // Check macro edges
  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    auto edgeFunction = edge->getData( edge_dst->getEdgeDataID() );

    for ( const auto & idxIt : edgedof::macroedge::Iterator( maxLevel ) )
    {
      auto ptr = edgeFunction->getPointer( maxLevel );
      auto idx = edgedof::macroedge::indexFromHorizontalEdge( maxLevel, idxIt.x(), stencilDirection::EDGE_HO_C );

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

    for ( const auto & idxIt : edgedof::macroface::Iterator( maxLevel ) )
    {
      if ( idxIt.x() != 0 )
      {
        const auto idx = edgedof::macroface::indexFromVerticalEdge( maxLevel, idxIt.x(), idxIt.y(), stencilDirection::EDGE_VE_C );

        const real_t expectedValue = edgeSrcValue * ( macroFaceVerticalStencilValue + 2.0 * ( macroFaceDiagonalStencilValue + macroFaceHorizontalStencilValue ) );
        WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
      }

      if ( idxIt.y() != 0 )
      {
        const auto idx = edgedof::macroface::indexFromHorizontalEdge( maxLevel, idxIt.x(), idxIt.y(), stencilDirection::EDGE_HO_C );

        const real_t expectedValue = edgeSrcValue * ( macroFaceHorizontalStencilValue + 2.0 * ( macroFaceDiagonalStencilValue + macroFaceVerticalStencilValue ) );
        WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
      }

      if ( idxIt.y() + idxIt.x() < idx_t( levelinfo::num_microedges_per_edge( maxLevel ) - 1 ) )
      {
        const auto idx = edgedof::macroface::indexFromDiagonalEdge( maxLevel, idxIt.x(), idxIt.y(), stencilDirection::EDGE_DI_C );

        const real_t expectedValue = edgeSrcValue * ( macroFaceDiagonalStencilValue + 2.0 * ( macroFaceHorizontalStencilValue + macroFaceVerticalStencilValue ) );
        WALBERLA_CHECK_FLOAT_EQUAL( ptr[ idx ], expectedValue );
      }
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
   hyteg::testEdgeDoFToEdgeDoFOperator();

   return EXIT_SUCCESS;
}
