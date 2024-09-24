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
#include "core/debug/all.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t level = 5;

   size_t mEperEdge = levelinfo::num_microedges_per_edge( level );

   EdgeDoFFunction< real_t >                                        x( "x", storage, level, level );
   std::vector< PrimitiveDataID< FunctionMemory< real_t >, Edge > > emptyEdgeIds;
   std::vector< PrimitiveDataID< FunctionMemory< real_t >, Face > > emptyFaceIds;

   for ( const auto& face : storage->getFaces() )
   {
      for ( size_t i = 0; i < levelinfo::num_microedges_per_face( level ); ++i )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( face.second->getData( x.getFaceDataID() )->getPointer( level )[i], 0.0 );
      }
   }
   for ( const auto& edge : storage->getEdges() )
   {
      for ( size_t i = 0; i < mEperEdge + edge.second->getNumHigherDimNeighbors() * ( mEperEdge * 2 + ( mEperEdge - 1 ) ); ++i )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( edge.second->getData( x.getEdgeDataID() )->getPointer( level )[i], 0.0 );
      }
   }

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > exact =
       []( const Point3D& xx, const std::vector< real_t >& ) { return 2 * xx[0] + xx[1]; };

   real_t value, xStepSize, yStepSize;
   for ( const auto& faceIter : storage->getFaces() )
   {
      auto face = faceIter.second;
      xStepSize =
          walberla::real_c( face->getCoordinates()[1][0] - face->getCoordinates()[0][0] ) / walberla::real_c( ( mEperEdge ) );
      yStepSize =
          walberla::real_c( face->getCoordinates()[2][1] - face->getCoordinates()[0][1] ) / walberla::real_c( ( mEperEdge ) );
      value = ( xStepSize / 2 + face->getCoordinates()[0][0] ) * 2 + face->getCoordinates()[0][0];

      edgedof::macroface::interpolate< real_t >( storage, level, *face, x.getFaceDataID(), emptyFaceIds, exact );
      for ( idx_t i = 0; i < idx_t( mEperEdge ); ++i )
      {
         for ( idx_t j = 0; j < idx_t( mEperEdge ) - i; ++j )
         {
            uint_t idx_ho = edgedof::macroface::indexFromHorizontalEdge( level, j, i, stencilDirection::EDGE_HO_C );
            uint_t idx_di = edgedof::macroface::indexFromDiagonalEdge( level, j, i, stencilDirection::EDGE_DI_C );
            uint_t idx_ve = edgedof::macroface::indexFromVerticalEdge( level, j, i, stencilDirection::EDGE_VE_C );
            if ( i == 0 )
            {
               WALBERLA_CHECK_FLOAT_EQUAL( face->getData( x.getFaceDataID() )->getPointer( level )[idx_ho],
                                           0.0,
                                           "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value );
            }
            else
            {
               WALBERLA_CHECK_FLOAT_EQUAL( face->getData( x.getFaceDataID() )->getPointer( level )[idx_ho],
                                           value,
                                           "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value );
            }
            if ( j == 0 )
            {
               WALBERLA_CHECK_FLOAT_EQUAL( face->getData( x.getFaceDataID() )->getPointer( level )[idx_ve],
                                           0.0,
                                           "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value );
            }
            else
            {
               WALBERLA_CHECK_FLOAT_EQUAL( face->getData( x.getFaceDataID() )->getPointer( level )[idx_ve],
                                           value + yStepSize / 2 - xStepSize,
                                           "i: " << i << " j: " << j << " idx: " << idx_ho << " value was "
                                                 << value + yStepSize / 2 - xStepSize );
            }
            if ( i + j == static_cast< idx_t >( levelinfo::num_microvertices_per_edge( level ) ) - 2 )
            {
               WALBERLA_CHECK_FLOAT_EQUAL( face->getData( x.getFaceDataID() )->getPointer( level )[idx_di],
                                           0.0,
                                           "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value );
            }
            else
            {
               WALBERLA_CHECK_FLOAT_EQUAL( face->getData( x.getFaceDataID() )->getPointer( level )[idx_di],
                                           value + yStepSize / 2,
                                           "i: " << i << " j: " << j << " idx: " << idx_ho << " value was "
                                                 << value + yStepSize / 2 );
            }
            value += 2 * xStepSize;
         }
         value = yStepSize * (real_t) ( i + 1 ) + xStepSize / 2 * 2;
      }
   }

   for ( const auto& edgeIter : storage->getEdges() )
   {
      auto edge = edgeIter.second;
      hyteg::edgedof::macroedge::interpolate< real_t >( storage, level, *edge, x.getEdgeDataID(), emptyEdgeIds, exact );
      value     = 2 * edge->getCoordinates()[0][0] + edge->getCoordinates()[0][1];
      xStepSize = edge->getDirection()[0] / walberla::real_c( ( mEperEdge ) );
      yStepSize = edge->getDirection()[1] / walberla::real_c( ( mEperEdge ) );
      value += xStepSize * 2 / 2;
      value += yStepSize / 2;
      for ( uint_t i = 0; i < mEperEdge; ++i )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             edge->getData( x.getEdgeDataID() )->getPointer( level )[i], value, "i: " << i << " edge: " << *edge );
         value += 2 * xStepSize;
         value += yStepSize;
      }
   }

   return EXIT_SUCCESS;
}
