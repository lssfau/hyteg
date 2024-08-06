/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr.
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
#include "core/debug/all.h"
#include "core/mpi/all.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;

using walberla::real_t;

template < uint_t Level >
void checkComm( std::string meshfile, bool bufferComm = false )
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( meshfile );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hyteg::EdgeDoFFunction< int > x( "x", storage, Level, Level );
   if ( bufferComm )
   {
      x.setLocalCommunicationMode( communication::BufferedCommunicator::BUFFERED_MPI );
   }

   int    num   = 0;
   uint_t check = 0;

   uint_t totalDoFs = hyteg::levelinfo::num_microedges_per_face( Level ) * storage->getNumberOfLocalFaces();

   for ( auto& edgeIt : storage->getEdges() )
   {
      hyteg::edgedof::macroedge::enumerate< int >( Level, *edgeIt.second, x.getEdgeDataID(), num );
      Edge& edge     = *edgeIt.second;
      int*  edgeData = edge.getData( x.getEdgeDataID() )->getPointer( Level );

      if ( edgeIt.second->getNumNeighborFaces() == 2 )
         totalDoFs -= levelinfo::num_microedges_per_edge( Level );

      for ( idx_t i = 0; i < idx_t( levelinfo::num_microedges_per_edge( Level ) ); ++i )
      {
         WALBERLA_CHECK_EQUAL( edgeData[i], check );
         WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_DI_S )],
                               0 );
         WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_VE_SE )],
                               0 );
         if ( i != 0 )
         {
            WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_SE )], 0 );
         }
         if ( edgeIt.second->getNumNeighborFaces() == 2 )
         {
            WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_DI_N )],
                                  0 );
            WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_VE_NW )],
                                  0 );
            if ( i != 0 )
            {
               WALBERLA_CHECK_EQUAL( edgeData[edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_NW )], 0 );
            }
         }
         check++;
      }
   }

   for ( auto& faceIt : storage->getFaces() )
   {
      hyteg::edgedof::macroface::enumerate< int >( Level, *faceIt.second, x.getFaceDataID(), num );
      size_t idxCounter    = 0;
      Face&  face          = *faceIt.second;
      int*   faceData      = face.getData( x.getFaceDataID() )->getPointer( Level );
      uint_t rowsize       = levelinfo::num_microedges_per_edge( Level );
      uint_t inner_rowsize = rowsize;
      for ( uint_t i = 0; i < rowsize; ++i )
      {
         for ( uint_t j = 0; j < inner_rowsize; ++j )
         {
            if ( i == 0 )
            {
               WALBERLA_CHECK_EQUAL( faceData[idxCounter], 0, "idxCounter was: " << idxCounter );
               ++idxCounter;
            }
            else
            {
               WALBERLA_CHECK_EQUAL( faceData[idxCounter], check, "idxCounter was: " << idxCounter );
               ++idxCounter;
               ++check;
            }
         }
         --inner_rowsize;
      }
      rowsize       = levelinfo::num_microedges_per_edge( Level );
      inner_rowsize = rowsize;
      for ( uint_t i = 0; i < rowsize; ++i )
      {
         for ( uint_t j = 0; j < inner_rowsize; ++j )
         {
            if ( ( j + i ) == ( levelinfo::num_microedges_per_edge( Level ) - 1 ) )
            {
               WALBERLA_CHECK_EQUAL( faceData[idxCounter], 0, "idxCounter was: " << idxCounter );
               ++idxCounter;
            }
            else
            {
               WALBERLA_CHECK_EQUAL( faceData[idxCounter], check, "idxCounter was: " << idxCounter );
               ++idxCounter;
               ++check;
            }
         }
         --inner_rowsize;
      }
      rowsize       = levelinfo::num_microedges_per_edge( Level );
      inner_rowsize = rowsize;
      for ( uint_t i = 0; i < rowsize; ++i )
      {
         for ( uint_t j = 0; j < inner_rowsize; ++j )
         {
            if ( j == 0 )
            {
               WALBERLA_CHECK_EQUAL( faceData[idxCounter], 0, "idxCounter was: " << idxCounter );
               ++idxCounter;
            }
            else
            {
               WALBERLA_CHECK_EQUAL( faceData[idxCounter], check, "idxCounter was: " << idxCounter );
               ++idxCounter;
               ++check;
            }
         }
         --inner_rowsize;
      }
   }

   //--check;
   WALBERLA_CHECK_EQUAL( check, totalDoFs );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::debug::enterTestMode();

   checkComm< 2 >( prependHyTeGMeshDir( "quad_2el.msh" ) );
   checkComm< 2 >( prependHyTeGMeshDir( "quad_2el.msh" ) );
   checkComm< 2 >( prependHyTeGMeshDir( "bfs_12el.msh" ) );

   checkComm< 3 >( prependHyTeGMeshDir( "tri_1el.msh" ) );
   checkComm< 3 >( prependHyTeGMeshDir( "quad_2el.msh" ) );
   checkComm< 3 >( prependHyTeGMeshDir( "bfs_12el.msh" ) );

   checkComm< 4 >( prependHyTeGMeshDir( "tri_1el.msh" ) );
   checkComm< 4 >( prependHyTeGMeshDir( "quad_2el.msh" ) );
   checkComm< 4 >( prependHyTeGMeshDir( "bfs_12el.msh" ) );
}
