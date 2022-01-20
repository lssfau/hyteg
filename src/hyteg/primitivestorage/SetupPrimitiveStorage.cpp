/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include <algorithm>
#include <iomanip>
#include <set>

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

SetupPrimitiveStorage::SetupPrimitiveStorage( const MeshInfo& meshInfo, const uint_t& numberOfProcesses )
: numberOfProcesses_( numberOfProcesses )
{
   WALBERLA_ASSERT_GREATER( numberOfProcesses_, 0, "Number of processes must be positive" );

   // since the MeshInfo IDs of the vertices do not necessarily
   // match the primitive IDs of the vertices in the SetupStorage, we need an assignment
   std::map< uint_t, PrimitiveID > meshVertexIDToPrimitiveID;

   // We cache the inserted primitives (edges, faces and cells) by filling
   // these maps with the surrounding vertexIDs as keys and the inserted
   // PrimitiveIDs as values.
   // This way we do not need to search for the neighboring lower level
   // primitives when building inner primitives.
   std::map< std::vector< PrimitiveID >, PrimitiveID > vertexIDsToEdgeIDs;
   std::map< std::vector< PrimitiveID >, PrimitiveID > vertexIDsToFaceIDs;
   auto findCachedPrimitiveID = []( const std::vector< PrimitiveID >&                          unsortedPrimitiveIDs,
                                    const std::map< std::vector< PrimitiveID >, PrimitiveID >& cache ) -> PrimitiveID {
      std::vector< PrimitiveID > sortedKey( unsortedPrimitiveIDs );
      std::sort( sortedKey.begin(), sortedKey.end() );
      WALBERLA_ASSERT_GREATER(
          cache.count( sortedKey ), 0, "Could not find primitive in cache during SetupStorage construction." );
      return cache.at( sortedKey );
   };

   // Adding vertices to storage
   const MeshInfo::VertexContainer vertices = meshInfo.getVertices();
   for ( const auto& it : vertices )
   {
      const MeshInfo::Vertex meshInfoVertex = it.second;

      PrimitiveID vertexID = generatePrimitiveID();

      meshVertexIDToPrimitiveID[meshInfoVertex.getID()] = vertexID;

      Point3D coordinates( meshInfoVertex.getCoordinates() );
      vertices_[vertexID.getID()] = std::make_shared< Vertex >( vertexID, coordinates );

      setMeshBoundaryFlag( vertexID.getID(), meshInfoVertex.getBoundaryFlag() );
   }

   // Adding edges to storage
   const MeshInfo::EdgeContainer edges = meshInfo.getEdges();
   for ( const auto& it : edges )
   {
      const MeshInfo::Edge meshInfoEdge = it.second;

      PrimitiveID edgeID = generatePrimitiveID();

      WALBERLA_ASSERT_EQUAL( meshInfoEdge.getVertices().size(), 2, "Edges are expected to have two vertices." );
      PrimitiveID vertexID0 = meshVertexIDToPrimitiveID[meshInfoEdge.getVertices().at( 0 )];
      PrimitiveID vertexID1 = meshVertexIDToPrimitiveID[meshInfoEdge.getVertices().at( 1 )];

      std::array< Point3D, 2 > coords;

      coords[0] = vertices_[vertexID0.getID()]->getCoordinates();
      coords[1] = vertices_[vertexID1.getID()]->getCoordinates();

      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID.getID() ), 0 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID0.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID1.getID() ), 1 );
      edges_[edgeID.getID()] = std::make_shared< Edge >( edgeID, vertexID0, vertexID1, coords );

      setMeshBoundaryFlag( edgeID.getID(), meshInfoEdge.getBoundaryFlag() );

      // Adding edge ID as neighbor to SetupVertices
      vertices_[vertexID0.getID()]->addEdge( edgeID );
      vertices_[vertexID1.getID()]->addEdge( edgeID );

      // Caching neighboring vertices
      std::vector< PrimitiveID > vertexIDs = {{vertexID0, vertexID1}};
      std::sort( vertexIDs.begin(), vertexIDs.end() );
      vertexIDsToEdgeIDs[vertexIDs] = edgeID;
   }

   // Adding faces to storage
   const MeshInfo::FaceContainer faces = meshInfo.getFaces();
   for ( const auto& it : faces )
   {
      const MeshInfo::Face meshInfoFace = it.second;

      PrimitiveID faceID = generatePrimitiveID();

      WALBERLA_ASSERT_EQUAL( meshInfoFace.getVertices().size(), 3, "Only supporting triangle faces." );
      PrimitiveID vertexID0 = meshVertexIDToPrimitiveID[meshInfoFace.getVertices().at( 0 )];
      PrimitiveID vertexID1 = meshVertexIDToPrimitiveID[meshInfoFace.getVertices().at( 1 )];
      PrimitiveID vertexID2 = meshVertexIDToPrimitiveID[meshInfoFace.getVertices().at( 2 )];

      WALBERLA_ASSERT_EQUAL( faces_.count( faceID.getID() ), 0 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID0.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID1.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID2.getID() ), 1 );

      PrimitiveID edgeID0 = findCachedPrimitiveID( {{vertexID0, vertexID1}}, vertexIDsToEdgeIDs );
      PrimitiveID edgeID1 = findCachedPrimitiveID( {{vertexID0, vertexID2}}, vertexIDsToEdgeIDs );
      PrimitiveID edgeID2 = findCachedPrimitiveID( {{vertexID1, vertexID2}}, vertexIDsToEdgeIDs );

      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID0.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID1.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID2.getID() ), 1 );

      // Edge Orientation
      std::array< int, 3 > edgeOrientation;

      PrimitiveID edge0Vertex0 = edges_[edgeID0.getID()]->getVertexID0();
      PrimitiveID edge0Vertex1 = edges_[edgeID0.getID()]->getVertexID1();
      PrimitiveID edge1Vertex0 = edges_[edgeID1.getID()]->getVertexID0();
      PrimitiveID edge1Vertex1 = edges_[edgeID1.getID()]->getVertexID1();
      PrimitiveID edge2Vertex0 = edges_[edgeID2.getID()]->getVertexID0();
      PrimitiveID edge2Vertex1 = edges_[edgeID2.getID()]->getVertexID1();

      WALBERLA_UNUSED( edge2Vertex1 );

      if ( edge0Vertex0 == vertexID0 )
      {
         WALBERLA_ASSERT( edge0Vertex1 == vertexID1 );
         edgeOrientation[0] = 1;
      }
      else
      {
         WALBERLA_ASSERT( edge0Vertex0 == vertexID1 );
         WALBERLA_ASSERT( edge0Vertex1 == vertexID0 );
         edgeOrientation[0] = -1;
      }

      if ( edge1Vertex0 == vertexID0 )
      {
         WALBERLA_ASSERT( edge1Vertex1 == vertexID2 );
         edgeOrientation[1] = 1;
      }
      else
      {
         WALBERLA_ASSERT( edge1Vertex0 == vertexID2 );
         WALBERLA_ASSERT( edge1Vertex1 == vertexID0 );
         edgeOrientation[1] = -1;
      }

      if ( edge2Vertex0 == vertexID1 )
      {
         WALBERLA_ASSERT( edge2Vertex1 == vertexID2 );
         edgeOrientation[2] = 1;
      }
      else
      {
         WALBERLA_ASSERT( edge2Vertex0 == vertexID2 );
         WALBERLA_ASSERT( edge2Vertex1 == vertexID1 );
         edgeOrientation[2] = -1;
      }

      // Corner coordinates
      std::array< Point3D, 3 >     coordinates;
      std::array< PrimitiveID, 3 > vertexIDs;
      std::vector< PrimitiveID >   verticesOnBoundary;
      std::vector< PrimitiveID >   edgesOnBoundary;

      if ( edgeOrientation[0] == 1 )
      {
         coordinates[0] = vertices_[edge0Vertex0.getID()]->getCoordinates();
         coordinates[1] = vertices_[edge0Vertex1.getID()]->getCoordinates();

         vertexIDs[0] = edge0Vertex0.getID();
         vertexIDs[1] = edge0Vertex1.getID();
      }
      else
      {
         coordinates[0] = vertices_[edge0Vertex1.getID()]->getCoordinates();
         coordinates[1] = vertices_[edge0Vertex0.getID()]->getCoordinates();

         vertexIDs[0] = edge0Vertex1.getID();
         vertexIDs[1] = edge0Vertex0.getID();
      }

      if ( edgeOrientation[1] == 1 )
      {
         coordinates[2] = vertices_[edge1Vertex1.getID()]->getCoordinates();

         vertexIDs[2] = edge1Vertex1.getID();
      }
      else
      {
         coordinates[2] = vertices_[edge1Vertex0.getID()]->getCoordinates();

         vertexIDs[2] = edge1Vertex0.getID();
      }

      faces_[faceID.getID()] =
          std::shared_ptr< Face >( new Face( faceID, vertexIDs, {{edgeID0, edgeID1, edgeID2}}, edgeOrientation, coordinates ) );

      setMeshBoundaryFlag( faceID.getID(), meshInfoFace.getBoundaryFlag() );

      // Adding face ID to vertices as neighbors
      vertices_[vertexIDs[0].getID()]->addFace( faceID );
      vertices_[vertexIDs[1].getID()]->addFace( faceID );
      vertices_[vertexIDs[2].getID()]->addFace( faceID );

      // Adding face ID to edges as neighbors
      edges_[edgeID0.getID()]->addFace( faceID );
      edges_[edgeID1.getID()]->addFace( faceID );
      edges_[edgeID2.getID()]->addFace( faceID );

      // Caching neighboring vertices
      std::vector< PrimitiveID > neighboringVertexIDs = {{vertexID0, vertexID1, vertexID2}};
      std::sort( neighboringVertexIDs.begin(), neighboringVertexIDs.end() );
      vertexIDsToFaceIDs[neighboringVertexIDs] = faceID;
   }

   // Adding cells to storage
   const MeshInfo::CellContainer cells = meshInfo.getCells();
   for ( const auto& it : cells )
   {
      const MeshInfo::Cell meshInfoCell = it.second;

      PrimitiveID cellID = generatePrimitiveID();

      WALBERLA_ASSERT_EQUAL( meshInfoCell.getVertices().size(), 4, "Only supporting tetrahedron cells." );

      PrimitiveID vertexID0 = meshVertexIDToPrimitiveID[meshInfoCell.getVertices().at( 0 )];
      PrimitiveID vertexID1 = meshVertexIDToPrimitiveID[meshInfoCell.getVertices().at( 1 )];
      PrimitiveID vertexID2 = meshVertexIDToPrimitiveID[meshInfoCell.getVertices().at( 2 )];
      PrimitiveID vertexID3 = meshVertexIDToPrimitiveID[meshInfoCell.getVertices().at( 3 )];

      PrimitiveID edgeID0 = findCachedPrimitiveID( {{vertexID0, vertexID1}}, vertexIDsToEdgeIDs );
      PrimitiveID edgeID1 = findCachedPrimitiveID( {{vertexID0, vertexID2}}, vertexIDsToEdgeIDs );
      PrimitiveID edgeID2 = findCachedPrimitiveID( {{vertexID1, vertexID2}}, vertexIDsToEdgeIDs );
      PrimitiveID edgeID3 = findCachedPrimitiveID( {{vertexID0, vertexID3}}, vertexIDsToEdgeIDs );
      PrimitiveID edgeID4 = findCachedPrimitiveID( {{vertexID1, vertexID3}}, vertexIDsToEdgeIDs );
      PrimitiveID edgeID5 = findCachedPrimitiveID( {{vertexID2, vertexID3}}, vertexIDsToEdgeIDs );

      PrimitiveID faceID0 = findCachedPrimitiveID( {{vertexID0, vertexID1, vertexID2}}, vertexIDsToFaceIDs );
      PrimitiveID faceID1 = findCachedPrimitiveID( {{vertexID0, vertexID1, vertexID3}}, vertexIDsToFaceIDs );
      PrimitiveID faceID2 = findCachedPrimitiveID( {{vertexID0, vertexID2, vertexID3}}, vertexIDsToFaceIDs );
      PrimitiveID faceID3 = findCachedPrimitiveID( {{vertexID1, vertexID2, vertexID3}}, vertexIDsToFaceIDs );

      WALBERLA_ASSERT_EQUAL( cells_.count( cellID.getID() ), 0 );

      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID0.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID1.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID2.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID3.getID() ), 1 );

      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID0.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID1.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID2.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID3.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID4.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( edges_.count( edgeID5.getID() ), 1 );

      WALBERLA_ASSERT_EQUAL( faces_.count( faceID0.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( faces_.count( faceID1.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( faces_.count( faceID2.getID() ), 1 );
      WALBERLA_ASSERT_EQUAL( faces_.count( faceID3.getID() ), 1 );

      std::vector< PrimitiveID > cellVertices = {{vertexID0, vertexID1, vertexID2, vertexID3}};
      std::vector< PrimitiveID > cellEdges    = {{edgeID0, edgeID1, edgeID2, edgeID3, edgeID4, edgeID5}};
      std::vector< PrimitiveID > cellFaces    = {{faceID0, faceID1, faceID2, faceID3}};

      for ( const auto& id : cellVertices )
      {
         WALBERLA_ASSERT( vertexExists( id ) );
         vertices_[id.getID()]->addCell( cellID );
      }
      for ( const auto& id : cellEdges )
      {
         WALBERLA_ASSERT( edgeExists( id ) );
         edges_[id.getID()]->addCell( cellID );
      }
      for ( const auto& id : cellFaces )
      {
         WALBERLA_ASSERT( faceExists( id ) );
         faces_[id.getID()]->addCell( cellID );
      }

      std::array< Point3D, 4 > cellCoordinates = {{vertices_.at( vertexID0.getID() )->getCoordinates(),
                                                   vertices_.at( vertexID1.getID() )->getCoordinates(),
                                                   vertices_.at( vertexID2.getID() )->getCoordinates(),
                                                   vertices_.at( vertexID3.getID() )->getCoordinates()}};

      std::array< std::map< uint_t, uint_t >, 6 > edgeLocalVertexToCellLocalVertexMaps;

      // edgeLocalVertexToCellLocalVertexMaps[ cellLocalEdgeID ][ edgeLocalVertexID ] = cellLocalVertexID;

      edgeLocalVertexToCellLocalVertexMaps[0][edges_.at( edgeID0.getID() )->vertex_index( vertexID0 )] = 0;
      edgeLocalVertexToCellLocalVertexMaps[0][edges_.at( edgeID0.getID() )->vertex_index( vertexID1 )] = 1;

      edgeLocalVertexToCellLocalVertexMaps[1][edges_.at( edgeID1.getID() )->vertex_index( vertexID0 )] = 0;
      edgeLocalVertexToCellLocalVertexMaps[1][edges_.at( edgeID1.getID() )->vertex_index( vertexID2 )] = 2;

      edgeLocalVertexToCellLocalVertexMaps[2][edges_.at( edgeID2.getID() )->vertex_index( vertexID1 )] = 1;
      edgeLocalVertexToCellLocalVertexMaps[2][edges_.at( edgeID2.getID() )->vertex_index( vertexID2 )] = 2;

      edgeLocalVertexToCellLocalVertexMaps[3][edges_.at( edgeID3.getID() )->vertex_index( vertexID0 )] = 0;
      edgeLocalVertexToCellLocalVertexMaps[3][edges_.at( edgeID3.getID() )->vertex_index( vertexID3 )] = 3;

      edgeLocalVertexToCellLocalVertexMaps[4][edges_.at( edgeID4.getID() )->vertex_index( vertexID1 )] = 1;
      edgeLocalVertexToCellLocalVertexMaps[4][edges_.at( edgeID4.getID() )->vertex_index( vertexID3 )] = 3;

      edgeLocalVertexToCellLocalVertexMaps[5][edges_.at( edgeID5.getID() )->vertex_index( vertexID2 )] = 2;
      edgeLocalVertexToCellLocalVertexMaps[5][edges_.at( edgeID5.getID() )->vertex_index( vertexID3 )] = 3;

      std::array< std::map< uint_t, uint_t >, 4 > faceLocalVertexToCellLocalVertexMaps;

      // faceLocalVertexToCellLocalVertexMaps[ cellLocalFaceID ][ faceLocalVertexID ] = cellLocalVertexID;

      faceLocalVertexToCellLocalVertexMaps[0][faces_.at( faceID0.getID() )->vertex_index( vertexID0 )] = 0;
      faceLocalVertexToCellLocalVertexMaps[0][faces_.at( faceID0.getID() )->vertex_index( vertexID1 )] = 1;
      faceLocalVertexToCellLocalVertexMaps[0][faces_.at( faceID0.getID() )->vertex_index( vertexID2 )] = 2;

      faceLocalVertexToCellLocalVertexMaps[1][faces_.at( faceID1.getID() )->vertex_index( vertexID0 )] = 0;
      faceLocalVertexToCellLocalVertexMaps[1][faces_.at( faceID1.getID() )->vertex_index( vertexID1 )] = 1;
      faceLocalVertexToCellLocalVertexMaps[1][faces_.at( faceID1.getID() )->vertex_index( vertexID3 )] = 3;

      faceLocalVertexToCellLocalVertexMaps[2][faces_.at( faceID2.getID() )->vertex_index( vertexID0 )] = 0;
      faceLocalVertexToCellLocalVertexMaps[2][faces_.at( faceID2.getID() )->vertex_index( vertexID2 )] = 2;
      faceLocalVertexToCellLocalVertexMaps[2][faces_.at( faceID2.getID() )->vertex_index( vertexID3 )] = 3;

      faceLocalVertexToCellLocalVertexMaps[3][faces_.at( faceID3.getID() )->vertex_index( vertexID1 )] = 1;
      faceLocalVertexToCellLocalVertexMaps[3][faces_.at( faceID3.getID() )->vertex_index( vertexID2 )] = 2;
      faceLocalVertexToCellLocalVertexMaps[3][faces_.at( faceID3.getID() )->vertex_index( vertexID3 )] = 3;

      cells_[cellID.getID()] = std::make_shared< Cell >( cellID,
                                                         cellVertices,
                                                         cellEdges,
                                                         cellFaces,
                                                         cellCoordinates,
                                                         edgeLocalVertexToCellLocalVertexMaps,
                                                         faceLocalVertexToCellLocalVertexMaps );

      setMeshBoundaryFlag( cellID.getID(), meshInfoCell.getBoundaryFlag() );
   }

   // add indirect neighbor faces
   for ( auto& it : faces_ )
   {
      auto faceID = it.first;
      auto face   = it.second;

      std::set< PrimitiveID > indirectNeighborsSet;

      for ( const auto& vertexID : face->neighborVertices() )
      {
         auto vertex = getVertex( vertexID );
         for ( const auto& neighborFaceID : vertex->neighborFaces() )
         {
            if ( neighborFaceID != faceID )
            {
               indirectNeighborsSet.insert( neighborFaceID );
            }
         }
      }

      face->indirectNeighborFaceIDs_.clear();
      face->indirectNeighborFaceIDs_.insert(
          face->indirectNeighborFaceIDs_.begin(), indirectNeighborsSet.begin(), indirectNeighborsSet.end() );
   }

   // add indirect neighbor cells
   for ( auto& it : cells_ )
   {
      auto cellID = it.first;
      auto cell   = it.second;

      std::set< PrimitiveID > indirectNeighborsSet;

      for ( const auto& vertexID : cell->neighborVertices() )
      {
         auto vertex = getVertex( vertexID );
         for ( const auto& neighborCellID : vertex->neighborCells() )
         {
            if ( neighborCellID != cellID )
            {
               indirectNeighborsSet.insert( neighborCellID );
            }
         }
      }

      cell->indirectNeighborCellIDs_.clear();
      cell->indirectNeighborCellIDs_.insert(
          cell->indirectNeighborCellIDs_.begin(), indirectNeighborsSet.begin(), indirectNeighborsSet.end() );
   }

   loadbalancing::roundRobin( *this );
}

SetupPrimitiveStorage::SetupPrimitiveStorage( const VertexMap& vertices,
                                              const EdgeMap&   edges,
                                              const FaceMap&   faces,
                                              const CellMap&   cells,
                                              const uint_t&    numberOfProcesses )
: numberOfProcesses_( numberOfProcesses )
, vertices_( vertices )
, edges_( edges )
, faces_( faces )
, cells_( cells )
{
   loadbalancing::roundRobin( *this );
}

Primitive* SetupPrimitiveStorage::getPrimitive( const PrimitiveID& id )
{
   if ( vertexExists( id ) )
   {
      return getVertex( id );
   }
   if ( edgeExists( id ) )
   {
      return getEdge( id );
   }
   if ( faceExists( id ) )
   {
      return getFace( id );
   }
   if ( cellExists( id ) )
   {
      return getCell( id );
   }
   return nullptr;
}

const Primitive* SetupPrimitiveStorage::getPrimitive( const PrimitiveID& id ) const
{
   if ( vertexExists( id ) )
   {
      return getVertex( id );
   }
   if ( edgeExists( id ) )
   {
      return getEdge( id );
   }
   if ( faceExists( id ) )
   {
      return getFace( id );
   }
   if ( cellExists( id ) )
   {
      return getCell( id );
   }
   return nullptr;
}

uint_t SetupPrimitiveStorage::getNumCellsOnRank( uint_t rank ) const
{
   uint_t n = 0;
   for ( const auto& it : cells_ )
   {
      if ( getTargetRank( it.first ) == rank )
      {
         n++;
      }
   }
   return n;
}

uint_t SetupPrimitiveStorage::getNumFacesOnRank( uint_t rank ) const
{
   uint_t n = 0;
   for ( const auto& it : faces_ )
   {
      if ( getTargetRank( it.first ) == rank )
      {
         n++;
      }
   }
   return n;
}

uint_t SetupPrimitiveStorage::getNumEdgesOnRank( uint_t rank ) const
{
   uint_t n = 0;
   for ( const auto& it : edges_ )
   {
      if ( getTargetRank( it.first ) == rank )
      {
         n++;
      }
   }
   return n;
}

uint_t SetupPrimitiveStorage::getNumVerticesOnRank( uint_t rank ) const
{
   uint_t n = 0;
   for ( const auto& it : vertices_ )
   {
      if ( getTargetRank( it.first ) == rank )
      {
         n++;
      }
   }
   return n;
}

void SetupPrimitiveStorage::assembleRankToSetupPrimitivesMap( RankToSetupPrimitivesMap& rankToSetupPrimitivesMap ) const
{
   rankToSetupPrimitivesMap.clear();

   PrimitiveMap setupPrimitives;
   getSetupPrimitives( setupPrimitives );
   for ( uint_t rank = 0; rank < numberOfProcesses_; rank++ )
   {
      rankToSetupPrimitivesMap[rank] = std::vector< PrimitiveID::IDType >();
      for ( auto setupPrimitive : setupPrimitives )
      {
         if ( rank == getTargetRank( setupPrimitive.first ) )
         {
            rankToSetupPrimitivesMap[rank].push_back( setupPrimitive.first );
         }
      }
   }

   WALBERLA_ASSERT_LESS_EQUAL( rankToSetupPrimitivesMap.size(), numberOfProcesses_ );
}

uint_t SetupPrimitiveStorage::getNumberOfEmptyProcesses() const
{
   uint_t                   numberOfEmptyProcesses = 0;
   RankToSetupPrimitivesMap rankToSetupPrimitivesMap;
   assembleRankToSetupPrimitivesMap( rankToSetupPrimitivesMap );
   for ( auto const& rankToSetupPrimitives : rankToSetupPrimitivesMap )
   {
      if ( rankToSetupPrimitives.second.size() == 0 )
      {
         numberOfEmptyProcesses++;
      }
   }
   return numberOfEmptyProcesses;
}

uint_t SetupPrimitiveStorage::getMinPrimitivesPerRank() const
{
   uint_t                   minNumberOfPrimitives = std::numeric_limits< uint_t >::max();
   RankToSetupPrimitivesMap rankToSetupPrimitivesMap;
   assembleRankToSetupPrimitivesMap( rankToSetupPrimitivesMap );
   for ( auto const& rankToSetupPrimitives : rankToSetupPrimitivesMap )
   {
      minNumberOfPrimitives = std::min( rankToSetupPrimitives.second.size(), minNumberOfPrimitives );
   }
   return minNumberOfPrimitives;
}

uint_t SetupPrimitiveStorage::getMaxPrimitivesPerRank() const
{
   uint_t                   maxNumberOfPrimitives = 0;
   RankToSetupPrimitivesMap rankToSetupPrimitivesMap;
   assembleRankToSetupPrimitivesMap( rankToSetupPrimitivesMap );
   for ( auto const& rankToSetupPrimitives : rankToSetupPrimitivesMap )
   {
      maxNumberOfPrimitives = std::max( rankToSetupPrimitives.second.size(), maxNumberOfPrimitives );
   }
   return maxNumberOfPrimitives;
}

real_t SetupPrimitiveStorage::getAvgPrimitivesPerRank() const
{
   return real_t( getNumberOfPrimitives() ) / real_t( numberOfProcesses_ );
}

uint_t SetupPrimitiveStorage::getNumVerticesOnBoundary() const
{
   return uint_c( std::count_if(
       vertices_.begin(), vertices_.end(), [this]( const std::pair< PrimitiveID::IDType, std::shared_ptr< Vertex > >& mapEntry ) {
          return onBoundary( PrimitiveID( mapEntry.first ), true );
       } ) );
}

uint_t SetupPrimitiveStorage::getNumEdgesOnBoundary() const
{
   return uint_c( std::count_if(
       edges_.begin(), edges_.end(), [this]( const std::pair< PrimitiveID::IDType, std::shared_ptr< Edge > >& mapEntry ) {
          return onBoundary( PrimitiveID( mapEntry.first ), true );
       } ) );
}

uint_t SetupPrimitiveStorage::getNumFacesOnBoundary() const
{
   return uint_c( std::count_if(
       faces_.begin(), faces_.end(), [this]( const std::pair< PrimitiveID::IDType, std::shared_ptr< Face > >& mapEntry ) {
          return onBoundary( PrimitiveID( mapEntry.first ), true );
       } ) );
}

uint_t SetupPrimitiveStorage::getNumCellsOnBoundary() const
{
   return uint_c( std::count_if(
       cells_.begin(), cells_.end(), [this]( const std::pair< PrimitiveID::IDType, std::shared_ptr< Cell > >& mapEntry ) {
          return onBoundary( PrimitiveID( mapEntry.first ), true );
       } ) );
}

void SetupPrimitiveStorage::toStream( std::ostream& os, bool verbose ) const
{
   os << "SetupPrimitiveStorage:\n";

   os << " - Processes (overall): " << std::setw( 10 ) << numberOfProcesses_ << "\n";
   os << " - Processes (empty)  : " << std::setw( 10 ) << getNumberOfEmptyProcesses() << "\n";

   os << " - Number of...\n"
      << "   +  Vertices: " << std::setw( 10 ) << vertices_.size() << " | " << getNumVerticesOnBoundary() << "/"
      << vertices_.size() << " on boundary\n"
      << "   +     Edges: " << std::setw( 10 ) << edges_.size() << " | " << getNumEdgesOnBoundary() << "/" << edges_.size()
      << " on boundary\n"
      << "   +     Faces: " << std::setw( 10 ) << faces_.size() << " | " << getNumFacesOnBoundary() << "/" << faces_.size()
      << " on boundary\n"
      << "   +     Cells: " << std::setw( 10 ) << cells_.size() << " | " << getNumCellsOnBoundary() << "/" << cells_.size()
      << " on boundary\n";

   os << " - Primitives per process...\n"
      << "   +      min: " << std::setw( 10 ) << getMinPrimitivesPerRank() << "\n"
      << "   +      max: " << std::setw( 10 ) << getMaxPrimitivesPerRank() << "\n"
      << "   +      avg: " << std::setw( 10 ) << getAvgPrimitivesPerRank() << "";

   if ( verbose )
   {
      os << "\n";
      os << "Vertices:   ID | Target Rank | Position  | Neighbor Edges \n"
         << "---------------------------------------------------------\n";
      for ( auto it = vertices_.begin(); it != vertices_.end(); it++ )
      {
         Point3D coordinates = it->second->getCoordinates();
         os << "          " << std::setw( 4 ) << it->first << " | " << std::setw( 11 ) << getTargetRank( it->first ) << " | "
            << coordinates << " | ";
         for ( const auto& neighborEdgeID : it->second->getHigherDimNeighbors() )
         {
            os << neighborEdgeID.getID() << " ";
         }
         os << "\n";
      }
      os << "\n";

      os << "Edges:      ID | Target Rank | VertexID_0 | VertexID_1 | mesh boundary flag   | Neighbor Faces \n"
         << "----------------------------------------------------------------------------------------------\n";
      for ( auto it = edges_.begin(); it != edges_.end(); it++ )
      {
         os << "          " << std::setw( 4 ) << it->first << " | " << std::setw( 11 ) << getTargetRank( it->first ) << " | "
            << std::setw( 10 ) << it->second->getVertexID0().getID() << " | " << std::setw( 10 )
            << it->second->getVertexID1().getID() << " | " << std::setw( 20 ) << it->second->getMeshBoundaryFlag() << " | ";
         for ( const auto& neighborFaceID : it->second->getHigherDimNeighbors() )
         {
            os << neighborFaceID.getID() << " ";
         }
         os << "\n";
      }
      os << "\n";

      os << "Faces:      ID | Target Rank | EdgeID_0 | EdgeID_1 | EdgeID_2\n"
         << "-------------------------------------------------------------\n";
      for ( auto it = faces_.begin(); it != faces_.end(); it++ )
      {
         os << "          " << std::setw( 4 ) << it->first << " | " << std::setw( 11 ) << getTargetRank( it->first ) << " | "
            << std::setw( 8 ) << it->second->getEdgeID0().getID() << " | " << std::setw( 8 ) << it->second->getEdgeID1().getID()
            << " | " << std::setw( 8 ) << it->second->getEdgeID2().getID() << "\n";
      }
      os << "\n";

      if ( cells_.size() > 0 )
      {
         os << "Cells:      ID | Target Rank | FaceID_0 | FaceID_1 | FaceID_2 | FaceID_3\n"
            << "------------------------------------------------------------------------\n";
         for ( auto it = cells_.begin(); it != cells_.end(); it++ )
         {
            os << "          " << std::setw( 4 ) << it->first << " | " << std::setw( 11 ) << getTargetRank( it->first ) << " | "
               << std::setw( 8 ) << it->second->neighborFaces()[0].getID() << " | " << std::setw( 8 )
               << it->second->neighborFaces()[1].getID() << " | " << std::setw( 8 ) << it->second->neighborFaces()[2].getID()
               << " | " << std::setw( 8 ) << it->second->neighborFaces()[3].getID() << "\n";
         }
      }
   }
}

void SetupPrimitiveStorage::getSetupPrimitives( PrimitiveMap& setupPrimitiveMap ) const
{
   setupPrimitiveMap.clear();

   setupPrimitiveMap.insert( vertices_.begin(), vertices_.end() );
   setupPrimitiveMap.insert( edges_.begin(), edges_.end() );
   setupPrimitiveMap.insert( faces_.begin(), faces_.end() );
   setupPrimitiveMap.insert( cells_.begin(), cells_.end() );

   WALBERLA_ASSERT_EQUAL( setupPrimitiveMap.size(), vertices_.size() + edges_.size() + faces_.size() + cells_.size() );
}

PrimitiveID SetupPrimitiveStorage::generatePrimitiveID() const
{
   PrimitiveID newID( getNumberOfPrimitives() );
   WALBERLA_ASSERT( !primitiveExists( newID ) );
   return newID;
}

void SetupPrimitiveStorage::setMeshBoundaryFlagsOnBoundary( const uint_t& meshBoundaryFlagOnBoundary,
                                                            const uint_t& meshBoundaryFlagInner,
                                                            const bool&   highestDimensionAlwaysInner )
{
   PrimitiveMap primitives;
   getSetupPrimitives( primitives );
   for ( auto& primitive : primitives )
   {
      primitive.second->meshBoundaryFlag_ =
          onBoundary( primitive.first, highestDimensionAlwaysInner ) ? meshBoundaryFlagOnBoundary : meshBoundaryFlagInner;
   }
}

void SetupPrimitiveStorage::setMeshBoundaryFlagsInner( const uint_t& meshBoundaryFlagInner,
                                                       const bool&   highestDimensionAlwaysInner )
{
   PrimitiveMap primitives;
   getSetupPrimitives( primitives );
   for ( auto& primitive : primitives )
   {
      primitive.second->meshBoundaryFlag_ = onBoundary( primitive.first, highestDimensionAlwaysInner ) ?
                                                primitive.second->getMeshBoundaryFlag() :
                                                meshBoundaryFlagInner;
   }
}

void SetupPrimitiveStorage::setMeshBoundaryFlagsByVertexLocation( const uint_t& meshBoundaryFlag,
                                                                  const std::function< bool( const Point3D& x ) >& onBoundary,
                                                                  const bool&                                      allVertices )
{
   auto cond = [allVertices, onBoundary]( const std::vector< Point3D >& coordinates ) {
      if ( allVertices )
         return std::all_of( coordinates.begin(), coordinates.end(), onBoundary );
      else
         return std::any_of( coordinates.begin(), coordinates.end(), onBoundary );
   };

   for ( const auto& p : vertices_ )
   {
      if ( cond( {p.second->getCoordinates()} ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }

   for ( const auto& p : edges_ )
   {
      if ( cond( std::vector< Point3D >( p.second->getCoordinates().begin(), p.second->getCoordinates().end() ) ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }

   for ( const auto& p : faces_ )
   {
      if ( cond( std::vector< Point3D >( p.second->getCoordinates().begin(), p.second->getCoordinates().end() ) ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }

   for ( const auto& p : cells_ )
   {
      if ( cond( std::vector< Point3D >( p.second->getCoordinates().begin(), p.second->getCoordinates().end() ) ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }
}

void SetupPrimitiveStorage::setMeshBoundaryFlagsByCentroidLocation( const uint_t& meshBoundaryFlag,
                                                                    const std::function< bool( const Point3D& x ) >& onBoundary )
{
   auto centroid = []( const std::vector< Point3D >& coordinates ) -> Point3D {
      Point3D c( {0, 0, 0} );
      for ( const auto& p : coordinates )
         c += p;
      c *= 1. / static_cast< double >( coordinates.size() );
      return c;
   };

   for ( const auto& p : vertices_ )
   {
      const auto c = centroid( {p.second->getCoordinates()} );
      if ( onBoundary( c ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }

   for ( const auto& p : edges_ )
   {
      const auto c = centroid( std::vector< Point3D >( p.second->getCoordinates().begin(), p.second->getCoordinates().end() ) );
      if ( onBoundary( c ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }

   for ( const auto& p : faces_ )
   {
      const auto c = centroid( std::vector< Point3D >( p.second->getCoordinates().begin(), p.second->getCoordinates().end() ) );
      if ( onBoundary( c ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }

   for ( const auto& p : cells_ )
   {
      const auto c = centroid( std::vector< Point3D >( p.second->getCoordinates().begin(), p.second->getCoordinates().end() ) );
      if ( onBoundary( c ) )
         setMeshBoundaryFlag( p.first, meshBoundaryFlag );
   }
}

bool SetupPrimitiveStorage::onBoundary( const PrimitiveID& primitiveID, const bool& highestDimensionAlwaysInner ) const
{
   WALBERLA_ASSERT( primitiveExists( primitiveID ) );

   if ( getNumberOfCells() == 0 )
   {
      // 2D
      if ( highestDimensionAlwaysInner && faceExists( primitiveID ) )
      {
         return false;
      }
      if ( edgeExists( primitiveID ) )
      {
         const auto edge = getEdge( primitiveID );
         WALBERLA_ASSERT_GREATER( edge->getNumNeighborFaces(), 0 );
         WALBERLA_ASSERT_LESS_EQUAL( edge->getNumNeighborFaces(), 2 );
         return edge->getNumNeighborFaces() == 1;
      }
      else
      {
         const auto                 primitive = getPrimitive( primitiveID );
         std::vector< PrimitiveID > neighborEdges;
         primitive->getNeighborEdges( neighborEdges );
         for ( auto it : neighborEdges )
         {
            if ( onBoundary( it ) )
            {
               return true;
            }
         }
         return false;
      }
   }
   else
   {
      // 3D
      if ( highestDimensionAlwaysInner && cellExists( primitiveID ) )
      {
         return false;
      }
      if ( faceExists( primitiveID ) )
      {
         const auto face = getFace( primitiveID );
         WALBERLA_ASSERT_GREATER( face->getNumNeighborCells(), 0 );
         WALBERLA_ASSERT_LESS_EQUAL( face->getNumNeighborCells(), 2 );
         return face->getNumNeighborCells() == 1;
      }
      else
      {
         const auto                 primitive = getPrimitive( primitiveID );
         std::vector< PrimitiveID > neighborFaces;
         primitive->getNeighborFaces( neighborFaces );
         for ( auto it : neighborFaces )
         {
            if ( onBoundary( it ) )
            {
               return true;
            }
         }
         return false;
      }
   }
}

} // namespace hyteg
