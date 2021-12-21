/*
 * Copyright (c) 2021 Benjamin Mann
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

#pragma once
#include "mesh.hpp"

// todo fix communication bug

namespace hyteg {
namespace adaptiveRefinement {

template < class K_Simplex >
uint_t K_Mesh< K_Simplex >::updateSetupStorage( const EdgeData& edges,
                                                const FaceData& faces,
                                                const CellData& cells,
                                                const uint_t&   n_processes )
{
   SetupPrimitiveStorage::VertexMap vertices_sps;
   SetupPrimitiveStorage::EdgeMap   edges_sps;
   SetupPrimitiveStorage::FaceMap   faces_sps;
   SetupPrimitiveStorage::CellMap   cells_sps;

   // give each primitive a running id
   uint_t id = 0;

   WALBERLA_LOG_INFO( ">>>>>>> add vertices <<<<<<<<<<<<<<\n\n" );
   //****** add vertices to storage ******

   for ( const auto& vtx : _vertices )
   {
      PrimitiveID vtxID( id );

      // add new vertex
      auto primitive   = std::make_shared< Vertex >( vtxID, vtx );
      vertices_sps[id] = primitive;

      // add properties
      primitive->meshBoundaryFlag_ = _vertexBoundaryFlag[id];
      primitive->geometryMap_      = _geometryMap[_vertexGeometryMap[id]];

      ++id;
   }

   WALBERLA_LOG_INFO( ">>>>>>>> add edges <<<<<<<<<<<<<<\n\n" );
   //****** add edges to storage ******

   // identify edges with their vertex IDs
   std::map< Idx< 2 >, PrimitiveID > vertexIDsToEdgeID;

   for ( size_t e = 0; e < edges.size(); ++e )
   {
      constexpr uint_t K = 1;

      PrimitiveID edgeID( id );

      // simplexData
      auto v             = edges.get_vertices( e );
      auto geometryMapID = edges.getGeometryMap( e );
      auto boundaryFlag  = edges.getBoundaryFlag( e );

      // extract data from simplex
      // if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
      // {
      //    edge->setPrimitiveID( edgeID );
      //    v             = edge->get_vertices();
      //    coords        = edge->coordinates( vertices );
      //    geometryMapID = edge->getGeometryMap();
      //    // boundaryFlag = edge->getBoundaryFlag() // todo
      // }
      // walberla::mpi::broadcastObject( v );
      // walberla::mpi::broadcastObject( coords );
      // walberla::mpi::broadcastObject( geometryMapID );
      // walberla::mpi::broadcastObject( boundaryFlag );

      // vertex coordinates and IDs
      std::array< Point3D, K + 1 >     coords;
      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( size_t i = 0; i < K + 1; ++i )
      {
         coords[i]    = _vertices[v[i]];
         vertexIDs[i] = PrimitiveID( v[i] );
      }

      // add new edge
      auto primitive = std::make_shared< Edge >( edgeID, vertexIDs[0], vertexIDs[1], coords );
      edges_sps[id]  = primitive;

      // add properties
      primitive->meshBoundaryFlag_ = boundaryFlag;
      primitive->geometryMap_      = _geometryMap[geometryMapID];

      // Adding edge ID as neighbor to SetupVertices
      for ( const auto& vertexID : vertexIDs )
      {
         vertices_sps[vertexID.getID()]->addEdge( edgeID );
      }

      // Caching neighboring vertices
      vertexIDsToEdgeID[v] = edgeID;

      ++id;
   }

   WALBERLA_LOG_INFO( ">>>>>>>> add faces <<<<<<<<<<<<<<\n\n" );
   //****** add faces to storage ******

   // write back ID of first face
   auto face0 = id;

   // identify faces with their vertex IDs
   std::map< Idx< 3 >, PrimitiveID > vertexIDsToFaceID;

   for ( size_t f = 0; f < faces.size(); ++f )
   {
      constexpr uint_t K = 2;

      PrimitiveID faceID( id );

      // simplexData
      auto v             = faces.get_vertices( f );
      auto geometryMapID = faces.getGeometryMap( f );
      auto boundaryFlag  = faces.getBoundaryFlag( f );

      // // extract data from simplex
      // if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
      // {
      //    face->setPrimitiveID( faceID );
      //    v             = face->get_vertices();
      //    coords        = face->coordinates( _vertices );
      //    geometryMapID = face->getGeometryMap();
      //    // boundaryFlag = face->getBoundaryFlag() // todo
      // }
      // walberla::mpi::broadcastObject( v );
      // walberla::mpi::broadcastObject( coords );
      // walberla::mpi::broadcastObject( geometryMapID );
      // walberla::mpi::broadcastObject( boundaryFlag );

      // vertex coordinates and IDs
      std::array< Point3D, K + 1 >     coords;
      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( size_t i = 0; i < K + 1; ++i )
      {
         coords[i]    = _vertices[v[i]];
         vertexIDs[i] = PrimitiveID( v[i] );
      }

      // ordering of edges
      constexpr std::array< std::array< uint_t, 2 >, K + 1 > edgeOrder{ { { 0ul, 1ul }, { 0ul, 2ul }, { 1ul, 2ul } } };

      // edge IDs
      std::array< PrimitiveID, K + 1 > edgeIDs;
      for ( uint_t i = 0; i < K + 1; ++i )
      {
         edgeIDs[i] = vertexIDsToEdgeID[{ v[edgeOrder[i][0]], v[edgeOrder[i][1]] }];
      }

      // edgeIDs[0] = vertexIDsToEdgeID[{ v[0], v[1] }];
      // edgeIDs[1] = vertexIDsToEdgeID[{ v[0], v[2] }];
      // edgeIDs[2] = vertexIDsToEdgeID[{ v[1], v[2] }];

      // edge orientation
      std::array< int, K + 1 > edgeOrientation;
      for ( uint_t i = 0; i < K + 1; ++i )
      {
         std::vector< PrimitiveID > edgeVertices;
         edges_sps[edgeIDs[i].getID()]->getNeighborVertices( edgeVertices );

         if ( edgeVertices[0].getID() == v[edgeOrder[i][0]] )
            edgeOrientation[i] = 1;
         else
            edgeOrientation[i] = -1;
      }

      // add new face
      auto primitive = std::make_shared< Face >( faceID, vertexIDs, edgeIDs, edgeOrientation, coords );
      faces_sps[id]  = primitive;

      // add properties
      primitive->meshBoundaryFlag_ = boundaryFlag;
      primitive->geometryMap_      = _geometryMap[geometryMapID];

      // Adding face ID to vertices as neighbors
      for ( const auto& vertexID : vertexIDs )
      {
         vertices_sps[vertexID.getID()]->addFace( faceID );
      }
      // Adding face ID to edges as neighbors
      for ( const auto& edgeID : edgeIDs )
      {
         edges_sps[edgeID.getID()]->addFace( faceID );
      }

      // Caching neighboring vertices
      vertexIDsToFaceID[v] = faceID;

      ++id;
   }

   WALBERLA_LOG_INFO( ">>>>>>>> add cells <<<<<<<<<<<<<<\n\n" );
   //****** add cells to storage ******

   // write back ID of first cell
   auto cell0 = id;

   for ( uint_t c = 0; c < cells.size(); ++c )
   {
      constexpr uint_t K = 3;

      PrimitiveID cellID( id );

      // simplexData
      auto v             = cells.get_vertices( c );
      auto geometryMapID = cells.getGeometryMap( c );
      auto boundaryFlag  = cells.getBoundaryFlag( c );

      // if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
      // {
      //    cell->setPrimitiveID( cellID );
      //    v             = cell->get_vertices();
      //    coords        = cell->coordinates( _vertices );
      //    geometryMapID = cell->getGeometryMap();
      //    // boundaryFlag = cell->getBoundaryFlag(); // todo
      // }
      // walberla::mpi::broadcastObject( v );
      // walberla::mpi::broadcastObject( coords );
      // walberla::mpi::broadcastObject( geometryMapID );
      // walberla::mpi::broadcastObject( boundaryFlag );

      // vertex IDs
      // vertex coordinates and IDs
      std::array< Point3D, K + 1 > coords;
      std::vector< PrimitiveID >   vertexIDs( K + 1 );
      for ( size_t i = 0; i < K + 1; ++i )
      {
         coords[i]    = _vertices[v[i]];
         vertexIDs[i] = PrimitiveID( v[i] );
      }

      // ordering of edges
      constexpr std::array< std::array< uint_t, 2 >, 6 > edgeOrder{
          { { 0ul, 1ul }, { 0ul, 2ul }, { 1ul, 2ul }, { 0ul, 3ul }, { 1ul, 3ul }, { 2ul, 3ul } } };

      // edge IDs
      std::vector< PrimitiveID > edgeIDs( 6 );
      for ( uint_t i = 0; i < 6; ++i )
      {
         edgeIDs[i] = vertexIDsToEdgeID[{ v[edgeOrder[i][0]], v[edgeOrder[i][1]] }];
      }
      // edgeIDs[0] = vertexIDsToEdgeID[{ v[0], v[1] }];
      // edgeIDs[1] = vertexIDsToEdgeID[{ v[0], v[2] }];
      // edgeIDs[2] = vertexIDsToEdgeID[{ v[1], v[2] }];
      // edgeIDs[3] = vertexIDsToEdgeID[{ v[0], v[3] }];
      // edgeIDs[4] = vertexIDsToEdgeID[{ v[1], v[3] }];
      // edgeIDs[5] = vertexIDsToEdgeID[{ v[2], v[3] }];

      // ordering of faces
      constexpr std::array< std::array< uint_t, 3 >, K + 1 > faceOrder{
          { { 0ul, 1ul, 2ul }, { 0ul, 1ul, 3ul }, { 0ul, 2ul, 3ul }, { 1ul, 2ul, 3ul } } };

      // face IDs
      std::vector< PrimitiveID > faceIDs( K + 1 );
      for ( uint_t i = 0; i < K + 1; ++i )
      {
         faceIDs[i] = vertexIDsToFaceID[{ v[faceOrder[i][0]], v[faceOrder[i][1]], v[faceOrder[i][2]] }];
      }
      // faceIDs[0] = vertexIDsToFaceID[{ v[0], v[1], v[2] }];
      // faceIDs[1] = vertexIDsToFaceID[{ v[0], v[1], v[3] }];
      // faceIDs[2] = vertexIDsToFaceID[{ v[0], v[2], v[3] }];
      // faceIDs[3] = vertexIDsToFaceID[{ v[1], v[2], v[3] }];

      std::array< std::map< uint_t, uint_t >, 6 > edgeLocalVertexToCellLocalVertexMaps;

      // edgeLocalVertexToCellLocalVertexMaps[ cellLocalEdgeID ][ edgeLocalVertexID ] = cellLocalVertexID;

      for ( uint_t i = 0; i < 6; ++i )
      {
         for ( const auto& j : edgeOrder[i] )
         {
            edgeLocalVertexToCellLocalVertexMaps[i][edges_sps.at( edgeIDs[i].getID() )->vertex_index( vertexIDs[j] )] = j;
         }
      }

      // edgeLocalVertexToCellLocalVertexMaps[0][edges_sps.at( edgeIDs[0].getID() )->vertex_index( vertexIDs[0] )] = 0;
      // edgeLocalVertexToCellLocalVertexMaps[0][edges_sps.at( edgeIDs[0].getID() )->vertex_index( vertexIDs[1] )] = 1;

      // edgeLocalVertexToCellLocalVertexMaps[1][edges_sps.at( edgeIDs[1].getID() )->vertex_index( vertexIDs[0] )] = 0;
      // edgeLocalVertexToCellLocalVertexMaps[1][edges_sps.at( edgeIDs[1].getID() )->vertex_index( vertexIDs[2] )] = 2;

      // edgeLocalVertexToCellLocalVertexMaps[2][edges_sps.at( edgeIDs[2].getID() )->vertex_index( vertexIDs[1] )] = 1;
      // edgeLocalVertexToCellLocalVertexMaps[2][edges_sps.at( edgeIDs[2].getID() )->vertex_index( vertexIDs[2] )] = 2;

      // edgeLocalVertexToCellLocalVertexMaps[3][edges_sps.at( edgeIDs[3].getID() )->vertex_index( vertexIDs[0] )] = 0;
      // edgeLocalVertexToCellLocalVertexMaps[3][edges_sps.at( edgeIDs[3].getID() )->vertex_index( vertexIDs[3] )] = 3;

      // edgeLocalVertexToCellLocalVertexMaps[4][edges_sps.at( edgeIDs[4].getID() )->vertex_index( vertexIDs[1] )] = 1;
      // edgeLocalVertexToCellLocalVertexMaps[4][edges_sps.at( edgeIDs[4].getID() )->vertex_index( vertexIDs[3] )] = 3;

      // edgeLocalVertexToCellLocalVertexMaps[5][edges_sps.at( edgeIDs[5].getID() )->vertex_index( vertexIDs[2] )] = 2;
      // edgeLocalVertexToCellLocalVertexMaps[5][edges_sps.at( edgeIDs[5].getID() )->vertex_index( vertexIDs[3] )] = 3;

      std::array< std::map< uint_t, uint_t >, 4 > faceLocalVertexToCellLocalVertexMaps;

      // faceLocalVertexToCellLocalVertexMaps[ cellLocalFaceID ][ faceLocalVertexID ] = cellLocalVertexID;

      for ( uint_t i = 0; i < K + 1; ++i )
      {
         for ( auto& j : faceOrder[i] )
         {
            faceLocalVertexToCellLocalVertexMaps[i][faces_sps.at( faceIDs[i].getID() )->vertex_index( vertexIDs[j] )] = j;
         }
      }

      // faceLocalVertexToCellLocalVertexMaps[0][faces_sps.at( faceIDs[0].getID() )->vertex_index( vertexIDs[0] )] = 0;
      // faceLocalVertexToCellLocalVertexMaps[0][faces_sps.at( faceIDs[0].getID() )->vertex_index( vertexIDs[1] )] = 1;
      // faceLocalVertexToCellLocalVertexMaps[0][faces_sps.at( faceIDs[0].getID() )->vertex_index( vertexIDs[2] )] = 2;

      // faceLocalVertexToCellLocalVertexMaps[1][faces_sps.at( faceIDs[1].getID() )->vertex_index( vertexIDs[0] )] = 0;
      // faceLocalVertexToCellLocalVertexMaps[1][faces_sps.at( faceIDs[1].getID() )->vertex_index( vertexIDs[1] )] = 1;
      // faceLocalVertexToCellLocalVertexMaps[1][faces_sps.at( faceIDs[1].getID() )->vertex_index( vertexIDs[3] )] = 3;

      // faceLocalVertexToCellLocalVertexMaps[2][faces_sps.at( faceIDs[2].getID() )->vertex_index( vertexIDs[0] )] = 0;
      // faceLocalVertexToCellLocalVertexMaps[2][faces_sps.at( faceIDs[2].getID() )->vertex_index( vertexIDs[2] )] = 2;
      // faceLocalVertexToCellLocalVertexMaps[2][faces_sps.at( faceIDs[2].getID() )->vertex_index( vertexIDs[3] )] = 3;

      // faceLocalVertexToCellLocalVertexMaps[3][faces_sps.at( faceIDs[3].getID() )->vertex_index( vertexIDs[1] )] = 1;
      // faceLocalVertexToCellLocalVertexMaps[3][faces_sps.at( faceIDs[3].getID() )->vertex_index( vertexIDs[2] )] = 2;
      // faceLocalVertexToCellLocalVertexMaps[3][faces_sps.at( faceIDs[3].getID() )->vertex_index( vertexIDs[3] )] = 3;

      // add new cell
      auto primitive = std::make_shared< Cell >( cellID,
                                                 vertexIDs,
                                                 edgeIDs,
                                                 faceIDs,
                                                 coords,
                                                 edgeLocalVertexToCellLocalVertexMaps,
                                                 faceLocalVertexToCellLocalVertexMaps );

      cells_sps[id] = primitive;

      // add properties
      primitive->meshBoundaryFlag_ = boundaryFlag;
      primitive->geometryMap_      = _geometryMap[geometryMapID];

      // Adding cell ID to vertices as neighbors
      for ( const auto& vertexID : vertexIDs )
      {
         vertices_sps[vertexID.getID()]->addCell( cellID );
      }
      // Adding cell ID to edges as neighbors
      for ( const auto& edgeID : edgeIDs )
      {
         edges_sps[edgeID.getID()]->addCell( cellID );
      }
      // Adding cell ID to faces as neighbors
      for ( const auto& faceID : faceIDs )
      {
         faces_sps[faceID.getID()]->addCell( cellID );
      }

      ++id;
   }

   WALBERLA_LOG_INFO( ">>>>>>>> add indirect nbrs <<<<<<<<<<<<<<\n\n" );
   //****** add indirect neighbor faces ******

   for ( const auto& [faceID, face] : faces_sps )
   {
      std::set< PrimitiveID > indirectNeighborsSet;

      for ( const auto& vertexID : face->neighborVertices() )
      {
         auto vertex = vertices_sps[vertexID.getID()];
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

   //****** add indirect neighbor cells ******

   for ( const auto& [cellID, cell] : cells_sps )
   {
      std::set< PrimitiveID > indirectNeighborsSet;

      for ( const auto& vertexID : cell->neighborVertices() )
      {
         auto vertex = vertices_sps[vertexID.getID()];
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

   WALBERLA_LOG_INFO( ">>>>>>>> setupStorage <<<<<<<<<<<<<<\n\n" );
   //****** construct new setupStorage ******
   _setupStorage = SetupPrimitiveStorage( vertices_sps, edges_sps, faces_sps, cells_sps, n_processes );

   WALBERLA_LOG_INFO( ">>>>>>>> return <<<<<<<<<<<<<<\n\n" );
   return ( cells.size() == 0 ) ? face0 : cell0;
}

} // namespace adaptiveRefinement
} // namespace hyteg