/*
 * Copyright (c) 2021-2022 Benjamin Mann
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

#include "mesh.hpp"

#include <algorithm>
#include <assert.h>
#include <core/mpi/Broadcast.h>
#include <iostream>
#include <map>
#include <numeric>
#include <type_traits>

// #include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "refine_cell.hpp"
#include "simplexFactory.hpp"

namespace hyteg {
namespace adaptiveRefinement {

template < class K_Simplex >
K_Mesh< K_Simplex >::K_Mesh( const SetupPrimitiveStorage& setupStorage )
: _n_processes( setupStorage.getNumberOfProcesses() )
{
   // copy geometrymaps for each primitive in initial setupStorage
   SetupPrimitiveStorage::PrimitiveMap setupPrimitives;
   setupStorage.getSetupPrimitives( setupPrimitives );
   for ( auto& [id, primitive] : setupPrimitives )
   {
      _geometryMap[id] = primitive->getGeometryMap();
   }
   _geometryMap[size_t( -1 )] = nullptr; // used for uninitialized values

   // internal data structures are only required on rank_0
   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      // extract vertices
      _n_vertices = setupStorage.getVertices().size();
      _vertices.resize( _n_vertices );
      _vertexBoundaryFlag.resize( _n_vertices );
      _vertexGeometryMap.resize( _n_vertices );

      // [0,1,...,n-1]
      std::vector< uint_t > vtxIndices( _n_vertices );
      // convert PrimitiveID::ID to Mesh::vertexID
      std::map< PrimitiveID::IDType, uint_t > conversion;

      // initialize vertices
      uint_t idx = 0;
      for ( auto& [id, vtx] : setupStorage.getVertices() )
      {
         // extract coordinates of vertex
         _vertices[idx] = vtx->getCoordinates();
         // extract geometrymap of vertex
         _vertexGeometryMap[idx] = id;
         // extract boundaryFlag of vertex
         _vertexBoundaryFlag[idx] = vtx->getMeshBoundaryFlag();
         // prepare element setup
         conversion[id]  = idx;
         vtxIndices[idx] = idx;
         ++idx;
      }

      // initialize all edges, faces and cells

      SimplexFactory fac( nullptr, vtxIndices );
      // simplex factory does not store cells
      std::set< std::shared_ptr< Simplex3 > > cells;
      std::vector< PrimitiveID >              v;

      for ( auto& [id, edge] : setupStorage.getEdges() )
      {
         edge->getNeighborVertices( v );
         auto myEdge = fac.make_edge( conversion[v[0].getID()], conversion[v[1].getID()] );
         myEdge->setPrimitiveID( id );
         myEdge->setGeometryMap( id );
         myEdge->setBoundaryFlag( edge->getMeshBoundaryFlag() );
      }

      for ( auto& [id, face] : setupStorage.getFaces() )
      {
         face->getNeighborVertices( v );
         auto myFace = fac.make_face( conversion[v[0].getID()], conversion[v[1].getID()], conversion[v[2].getID()] );
         myFace->setPrimitiveID( id );
         myFace->setGeometryMap( id );
         myFace->setBoundaryFlag( face->getMeshBoundaryFlag() );
      }

      for ( auto& [id, cell] : setupStorage.getCells() )
      {
         cell->getNeighborVertices( v );
         auto myCell = fac.make_cell(
             conversion[v[0].getID()], conversion[v[1].getID()], conversion[v[2].getID()], conversion[v[3].getID()] );
         myCell->setPrimitiveID( id );
         myCell->setGeometryMap( id );
         myCell->setBoundaryFlag( cell->getMeshBoundaryFlag() );
         cells.insert( myCell );
      }

      // insert volume elements into _T
      if constexpr ( std::is_same_v< K_Simplex, Simplex2 > )
      {
         if ( !cells.empty() )
         {
            WALBERLA_ABORT( "Adaptive 2D mesh requires SetupPrimitiveStorage without any cells!" );
         }
         for ( auto& p : fac.faces() )
         {
            _T.insert( p.second );
         }
      }
      if constexpr ( std::is_same_v< K_Simplex, Simplex3 > )
      {
         if ( cells.empty() )
         {
            WALBERLA_ABORT( "Adaptive 3D mesh requires SetupPrimitiveStorage containing at least one cell!" );
         }

         _T = cells;
      }

      _n_elements = _T.size();
   }

   // broadcast required values to all processes
   walberla::mpi::broadcastObject( _n_elements );
   walberla::mpi::broadcastObject( _n_vertices );
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::refineRG( const std::vector< PrimitiveID >& elements_to_refine, uint_t n_el_max )
{
   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      // pessimistic estimate! in most cases the actual growth will be significantly smaller
      const uint_t est_growth_factor = ( K_Simplex::DIM == 2 ) ? 18 : 138;

      uint_t predict = _T.size();

      /* green elements must not be refined any further to
         prevent mesh degeneration
      */
      remove_green_edges();

      // get elements corresponding to given IDs
      auto R = init_R( elements_to_refine );

      // unprocessed elements
      std::set< std::shared_ptr< K_Simplex > > U = _T;
      // refined elements
      std::set< std::shared_ptr< K_Simplex > > refined;

      /* successively apply recursive red-refinement to parts of R
         until predicted n_el exceeds n_el_max
      */
      auto prt      = R.begin();
      auto prt_size = R.size();
      while ( prt_size > 0 && prt != R.end() && predict < n_el_max )
      {
         // move to next chunk
         auto done = prt;
         prt_size  = ( n_el_max - predict ) / est_growth_factor;
         prt       = done + std::min( int64_t( prt_size ), int64_t( R.end() - done ) );

         std::set< std::shared_ptr< K_Simplex > > R_prt( done, prt );
         R_prt.erase( nullptr );

         /* recursively apply red refinement for elements
         that otherwise would be subject to multiple
         green refinement steps later on
         */
         while ( !R_prt.empty() )
         {
            refined.merge( refine_red( R_prt, U ) );
            R_prt = find_elements_for_red_refinement( U );
         }

         // predict number of elements after required green step
         predict = U.size() + refined.size() + predict_n_el_green( U );
      }

      // apply green refinement
      refined.merge( refine_green( U ) );

      // update current configuration
      _T = U;
      _T.merge( refined );
      _n_vertices = _vertices.size();
      _n_elements = _T.size();
   }

   walberla::mpi::broadcastObject( _n_vertices );
   walberla::mpi::broadcastObject( _n_elements );
}

template < class K_Simplex >
std::shared_ptr< PrimitiveStorage > K_Mesh< K_Simplex >::make_storage()
{
   // extract connectivity, geometry and boundary data and add PrimitiveIDs
   std::vector< EdgeData > edges;
   std::vector< FaceData > faces;
   std::vector< CellData > cells;
   extract_data( edges, faces, cells );

   // broadcast data to all processes
   walberla::mpi::broadcastObject( _vertices );
   walberla::mpi::broadcastObject( _vertexGeometryMap );
   walberla::mpi::broadcastObject( _vertexBoundaryFlag );
   walberla::mpi::broadcastObject( edges );
   walberla::mpi::broadcastObject( faces );
   walberla::mpi::broadcastObject( cells );

   // apply loadbalancing
   std::vector< uint_t > vertices_targetRank( _vertices.size() );
   loadbalancing( vertices_targetRank, edges, faces, cells, _n_processes );

   // todo use this method instead
   return make_localPrimitives( vertices_targetRank, edges, faces, cells );

   // // create local primitives and associated setupStorage
   // auto setupStorage = make_localSetupStorage( vertices_targetRank, edges, faces, cells );

   // // copy targetrank to setupStorage
   // for ( uint_t id = 0; id < vertices_targetRank.size(); ++id )
   // {
   //    setupStorage->setTargetRank( PrimitiveID( id ), vertices_targetRank[id] );
   // }
   // for ( auto& p : edges )
   // {
   //    setupStorage->setTargetRank( p.getPrimitiveID(), p.getTargetRank() );
   // }
   // for ( auto& p : faces )
   // {
   //    setupStorage->setTargetRank( p.getPrimitiveID(), p.getTargetRank() );
   // }
   // for ( auto& p : cells )
   // {
   //    setupStorage->setTargetRank( p.getPrimitiveID(), p.getTargetRank() );
   // }
   // return std::make_shared< PrimitiveStorage >( *setupStorage );
}

template < class K_Simplex >
std::shared_ptr< SetupPrimitiveStorage >
    K_Mesh< K_Simplex >::make_localSetupStorage( const std::vector< uint_t >&   vertices_targetRank,
                                                 const std::vector< EdgeData >& edges,
                                                 const std::vector< FaceData >& faces,
                                                 const std::vector< CellData >& cells ) const
{
   // todo only create local and halo primitives

   SetupPrimitiveStorage::VertexMap vertices_sps;
   SetupPrimitiveStorage::EdgeMap   edges_sps;
   SetupPrimitiveStorage::FaceMap   faces_sps;
   SetupPrimitiveStorage::CellMap   cells_sps;

   // identify edges with their vertex IDs
   std::map< Idx< 2 >, PrimitiveID > vertexIDsToEdgeID;
   // identify faces with their vertex IDs
   std::map< Idx< 3 >, PrimitiveID > vertexIDsToFaceID;

   // auto rank = uint_t( walberla::mpi::MPIManager::instance()->rank() );

   // convenient acces to SimplexData
   auto get_edge = [&]( uint_t id ) -> const EdgeData& {
      const uint_t id0 = edges[0].getPrimitiveID().getID();
      return edges[id - id0];
   };
   auto get_face = [&]( uint_t id ) -> const FaceData& {
      const uint_t id0 = faces[0].getPrimitiveID().getID();
      return faces[id - id0];
   };
   // auto get_cell = [&]( uint_t id ) -> const CellData& {
   //    const uint_t id0 = cells[0].getPrimitiveID().getID();
   //    return cells[id - id0];
   // };

   auto add_vertex = [&]( uint_t id ) {
      // add new vertex
      vertices_sps[id] = std::make_shared< Vertex >( PrimitiveID( id ), _vertices[id] );
      // add properties
      vertices_sps[id]->meshBoundaryFlag_ = _vertexBoundaryFlag[id];
      vertices_sps[id]->geometryMap_      = _geometryMap.at( _vertexGeometryMap[id] );
   };

   auto add_edge = [&]( const EdgeData& edge ) {
      constexpr uint_t K = 1;

      // vertex coordinates and IDs
      auto v      = edge.get_vertices();
      auto coords = edge.get_coordinates( _vertices );

      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         vertexIDs[i] = PrimitiveID( v[i] );
      }

      // add new edge
      auto id       = edge.getPrimitiveID().getID();
      edges_sps[id] = std::make_shared< Edge >( edge.getPrimitiveID(), vertexIDs[0], vertexIDs[1], coords );

      // add properties
      edges_sps[id]->meshBoundaryFlag_ = edge.getBoundaryFlag();
      edges_sps[id]->geometryMap_      = _geometryMap.at( edge.getGeometryMap() );

      // Adding edge ID as neighbor to SetupVertices
      for ( const auto& vertexID : vertexIDs )
      {
         if ( vertices_sps.count( vertexID.getID() ) == 0 )
         {
            add_vertex( vertexID.getID() );
         }
         vertices_sps[vertexID.getID()]->addEdge( edge.getPrimitiveID() );
      }

      // Caching neighboring vertices
      vertexIDsToEdgeID[v] = edge.getPrimitiveID();
   };

   auto add_face = [&]( const FaceData& face ) {
      constexpr uint_t K = 2;

      // vertex coordinates and IDs
      auto v      = face.get_vertices();
      auto coords = face.get_coordinates( _vertices );

      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         vertexIDs[i] = PrimitiveID( v[i] );
      }

      // ordering of edges
      constexpr std::array< std::array< uint_t, 2 >, K + 1 > edgeOrder{ { { 0ul, 1ul }, { 0ul, 2ul }, { 1ul, 2ul } } };

      // edge IDs
      std::array< PrimitiveID, K + 1 > edgeIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         edgeIDs[i] = vertexIDsToEdgeID[{ v[edgeOrder[i][0]], v[edgeOrder[i][1]] }];
      }

      // edge orientation
      std::array< int, K + 1 > edgeOrientation;
      for ( uint_t i = 0; i <= K; ++i )
      {
         std::vector< PrimitiveID > edgeVertices;
         edges_sps[edgeIDs[i].getID()]->getNeighborVertices( edgeVertices );

         if ( edgeVertices[0].getID() == v[edgeOrder[i][0]] )
            edgeOrientation[i] = 1;
         else
            edgeOrientation[i] = -1;
      }

      // add new face
      auto id       = face.getPrimitiveID().getID();
      faces_sps[id] = std::make_shared< Face >( face.getPrimitiveID(), vertexIDs, edgeIDs, edgeOrientation, coords );

      // add properties
      faces_sps[id]->meshBoundaryFlag_ = face.getBoundaryFlag();
      faces_sps[id]->geometryMap_      = _geometryMap.at( face.getGeometryMap() );

      // Adding face ID to edges as neighbors
      for ( const auto& edgeID : edgeIDs )
      {
         if ( edges_sps.count( edgeID.getID() ) == 0 )
         {
            add_edge( get_edge( edgeID.getID() ) );
         }
         edges_sps[edgeID.getID()]->addFace( face.getPrimitiveID() );
      }
      // Adding face ID to vertices as neighbors
      for ( const auto& vertexID : vertexIDs )
      {
         // note: if vertex didn't exist before, it has been added along with the edge
         vertices_sps[vertexID.getID()]->addFace( face.getPrimitiveID() );
      }

      // Caching neighboring vertices
      vertexIDsToFaceID[v] = face.getPrimitiveID();
   };

   auto add_cell = [&]( const CellData& cell ) {
      constexpr uint_t K = 3;

      // vertex coordinates and IDs
      auto v      = cell.get_vertices();
      auto coords = cell.get_coordinates( _vertices );

      std::vector< PrimitiveID > vertexIDs( K + 1 );
      for ( uint_t i = 0; i <= K; ++i )
      {
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

      // ordering of faces
      constexpr std::array< std::array< uint_t, 3 >, K + 1 > faceOrder{
          { { 0ul, 1ul, 2ul }, { 0ul, 1ul, 3ul }, { 0ul, 2ul, 3ul }, { 1ul, 2ul, 3ul } } };

      // face IDs
      std::vector< PrimitiveID > faceIDs( K + 1 );
      for ( uint_t i = 0; i <= K; ++i )
      {
         faceIDs[i] = vertexIDsToFaceID[{ v[faceOrder[i][0]], v[faceOrder[i][1]], v[faceOrder[i][2]] }];
      }

      std::array< std::map< uint_t, uint_t >, 6 > edgeLocalVertexToCellLocalVertexMaps;

      // edgeLocalVertexToCellLocalVertexMaps[ cellLocalEdgeID ][ edgeLocalVertexID ] = cellLocalVertexID;

      for ( uint_t i = 0; i < 6; ++i )
      {
         for ( const auto& j : edgeOrder[i] )
         {
            edgeLocalVertexToCellLocalVertexMaps[i][edges_sps.at( edgeIDs[i].getID() )->vertex_index( vertexIDs[j] )] = j;
         }
      }

      std::array< std::map< uint_t, uint_t >, 4 > faceLocalVertexToCellLocalVertexMaps;

      // faceLocalVertexToCellLocalVertexMaps[ cellLocalFaceID ][ faceLocalVertexID ] = cellLocalVertexID;

      for ( uint_t i = 0; i < K + 1; ++i )
      {
         for ( auto& j : faceOrder[i] )
         {
            faceLocalVertexToCellLocalVertexMaps[i][faces_sps.at( faceIDs[i].getID() )->vertex_index( vertexIDs[j] )] = j;
         }
      }

      // add new cell
      auto id       = cell.getPrimitiveID().getID();
      cells_sps[id] = std::make_shared< Cell >( cell.getPrimitiveID(),
                                                vertexIDs,
                                                edgeIDs,
                                                faceIDs,
                                                coords,
                                                edgeLocalVertexToCellLocalVertexMaps,
                                                faceLocalVertexToCellLocalVertexMaps );

      // add properties
      cells_sps[id]->meshBoundaryFlag_ = cell.getBoundaryFlag();
      cells_sps[id]->geometryMap_      = _geometryMap.at( cell.getGeometryMap() );

      // Adding cell ID to faces as neighbors
      for ( const auto& faceID : faceIDs )
      {
         if ( faces_sps.count( faceID.getID() ) == 0 )
         {
            add_face( get_face( faceID.getID() ) );
         }
         faces_sps[faceID.getID()]->addCell( cell.getPrimitiveID() );
      }
      // Adding cell ID to edges as neighbors
      for ( const auto& edgeID : edgeIDs )
      {
         edges_sps[edgeID.getID()]->addCell( cell.getPrimitiveID() );
      }
      // Adding cell ID to vertices as neighbors
      for ( const auto& vertexID : vertexIDs )
      {
         vertices_sps[vertexID.getID()]->addCell( cell.getPrimitiveID() );
      }
   };

   //****** add vertices to storage ******

   for ( uint_t id = 0; id < _vertices.size(); ++id )
   {
      add_vertex( id );
   }

   //****** add edges to storage ******

   for ( auto& edge : edges )
   {
      add_edge( edge );
   }

   //****** add faces to storage ******

   for ( auto& face : faces )
   {
      add_face( face );
   }

   //****** add cells to storage ******

   for ( auto& cell : cells )
   {
      add_cell( cell );
   }

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

   //****** construct new SetupPrimitiveStorage ******
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( vertices_sps, edges_sps, faces_sps, cells_sps, _n_processes );
   return setupStorage;
}

template < class K_Simplex >
std::shared_ptr< PrimitiveStorage > K_Mesh< K_Simplex >::make_localPrimitives( const std::vector< uint_t >& vertices_targetRank,
                                                                               std::vector< EdgeData >&     edges,
                                                                               std::vector< FaceData >&     faces,
                                                                               std::vector< CellData >&     cells ) const
{
   // todo FIX THIS METHOD
   auto rank = uint_t( walberla::mpi::MPIManager::instance()->rank() );

   // ****** create maps M: primitiveID -> simplexData ******

   // acces EdgeData via PrimitiveID
   auto get_edge = [&]( uint_t id ) -> EdgeData* {
      const uint_t id0 = edges[0].getPrimitiveID().getID();
      return &( edges[id - id0] );
   };
   // acces FaceData via PrimitiveID
   auto get_face = [&]( uint_t id ) -> FaceData* {
      const uint_t id0 = faces[0].getPrimitiveID().getID();
      return &( faces[id - id0] );
   };
   // acces CellData via PrimitiveID
   // auto get_cell = [&]( uint_t id ) -> CellData* {
   //    const uint_t id0 = cells[0].getPrimitiveID().getID();
   //    return &( cells[id - id0] );
   // };

   // ****** create maps M: vertexIds -> primitiveID ******

   // identify edges with their vertex IDs
   std::map< Idx< 2 >, PrimitiveID > vertexIDsToEdgeID;
   for ( auto& edge : edges )
   {
      vertexIDsToEdgeID[edge.get_vertices()] = edge.getPrimitiveID();
   }
   // identify faces with their vertex IDs
   std::map< Idx< 3 >, PrimitiveID > vertexIDsToFaceID;
   for ( auto& face : faces )
   {
      vertexIDsToFaceID[face.get_vertices()] = face.getPrimitiveID();
   }
   // identify cells with their vertex IDs
   std::map< Idx< 4 >, PrimitiveID > vertexIDsToCellID;
   for ( auto& cell : cells )
   {
      vertexIDsToCellID[cell.get_vertices()] = cell.getPrimitiveID();
   }

   // ****** find simplices required for halos ******

   std::vector< bool >   vtxOnHalo( _vertices.size(), false );
   hyteg::MigrationMap_T nbrRanks;
   // for ( uint_t i = 0; i < _vertices.size(); ++i )
   // {
   //    nbrRanks[vtx] = vertices_targetRank[vtx];
   // }
   // for ( auto& edge : edges )
   // {
   //    nbrRanks[edge.getPrimitiveID().getID()] = edge.getTargetRank();
   // }
   // for ( auto& face : faces )
   // {
   //    nbrRanks[face.getPrimitiveID().getID()] = face.getTargetRank();
   // }
   // for ( auto& cell : cells )
   // {
   //    nbrRanks[cell.getPrimitiveID().getID()] = cell.getTargetRank();
   // }

   for ( auto& edge : edges )
   {
      auto& v = edge.get_vertices();

      for ( auto& vtx : v )
      {
         if ( edge.getTargetRank() == rank && vertices_targetRank[vtx] != rank )
         {
            vtxOnHalo[vtx] = true;
            nbrRanks[vtx]  = vertices_targetRank[vtx];
         }
         if ( edge.getTargetRank() != rank && vertices_targetRank[vtx] == rank )
         {
            edge.setHaloInfo( true );
            nbrRanks[edge.getPrimitiveID().getID()] = edge.getTargetRank();
         }
      }
   }

   for ( auto& face : faces )
   {
      auto& v = face.get_vertices();

      for ( auto& vtx : v )
      {
         if ( face.getTargetRank() == rank && vertices_targetRank[vtx] != rank )
         {
            vtxOnHalo[vtx] = true;
            nbrRanks[vtx]  = vertices_targetRank[vtx];
         }
         if ( face.getTargetRank() != rank && vertices_targetRank[vtx] == rank )
         {
            face.setHaloInfo( true );
            nbrRanks[face.getPrimitiveID().getID()] = face.getTargetRank();
         }
      }

      for ( uint_t i = 0; i < 3; ++i )
      {
         Idx< 2 >  edgeVtxs{ v[i], v[( i + 1 ) % 3] };
         EdgeData* edge = get_edge( vertexIDsToEdgeID[edgeVtxs].getID() );

         if ( face.getTargetRank() == rank && edge->getTargetRank() != rank )
         {
            edge->setHaloInfo( true );
            nbrRanks[edge->getPrimitiveID().getID()] = edge->getTargetRank();
         }
         if ( face.getTargetRank() != rank && edge->getTargetRank() == rank )
         {
            face.setHaloInfo( true );
            nbrRanks[face.getPrimitiveID().getID()] = face.getTargetRank();
         }
      }
   }

   for ( auto& cell : cells )
   {
      auto& v = cell.get_vertices();

      for ( auto& vtx : v )
      {
         if ( cell.getTargetRank() == rank && vertices_targetRank[vtx] != rank )
         {
            vtxOnHalo[vtx] = true;
            nbrRanks[vtx]  = vertices_targetRank[vtx];
         }
         if ( cell.getTargetRank() != rank && vertices_targetRank[vtx] == rank )
         {
            cell.setHaloInfo( true );
            nbrRanks[cell.getPrimitiveID().getID()] = cell.getTargetRank();
         }
      }

      for ( uint_t i = 0; i < 4; ++i )
      {
         for ( uint_t j = i + 1; j < 4; ++j )
         {
            Idx< 2 >  edgeVtxs{ v[i], v[j] };
            EdgeData* edge = get_edge( vertexIDsToEdgeID[edgeVtxs].getID() );

            if ( cell.getTargetRank() == rank && edge->getTargetRank() != rank )
            {
               edge->setHaloInfo( true );
               nbrRanks[edge->getPrimitiveID().getID()] = edge->getTargetRank();
            }
            if ( cell.getTargetRank() != rank && edge->getTargetRank() == rank )
            {
               cell.setHaloInfo( true );
               nbrRanks[cell.getPrimitiveID().getID()] = cell.getTargetRank();
            }
         }
      }

      for ( uint_t j = 0; j < 4; ++j )
      {
         Idx< 3 >  faceVtxs{ v[j], v[( j + 1 ) % 4], v[( j + 2 ) % 4] };
         FaceData* face = get_face( vertexIDsToFaceID[faceVtxs].getID() );

         if ( cell.getTargetRank() == rank && face->getTargetRank() != rank )
         {
            face->setHaloInfo( true );
            nbrRanks[face->getPrimitiveID().getID()] = face->getTargetRank();
         }
         if ( cell.getTargetRank() != rank && face->getTargetRank() == rank )
         {
            cell.setHaloInfo( true );
            nbrRanks[cell.getPrimitiveID().getID()] = cell.getTargetRank();
         }
      }
   }

   // ****** create primitives ******

   PrimitiveStorage::VertexMap vtxs_ps;
   PrimitiveStorage::EdgeMap   edges_ps;
   PrimitiveStorage::FaceMap   faces_ps;
   PrimitiveStorage::CellMap   cells_ps;

   PrimitiveStorage::VertexMap nbrVtxs_ps;
   PrimitiveStorage::EdgeMap   nbrEdges_ps;
   PrimitiveStorage::FaceMap   nbrFaces_ps;
   PrimitiveStorage::CellMap   nbrCells_ps;

   auto add_vertex = [&]( PrimitiveStorage::VertexMap& map, uint_t id ) {
      // add new vertex
      map[id] = std::make_shared< Vertex >( PrimitiveID( id ), _vertices[id] );
      // add properties
      map[id]->meshBoundaryFlag_ = _vertexBoundaryFlag[id];
      map[id]->geometryMap_      = _geometryMap.at( _vertexGeometryMap[id] );
   };

   // create local and halo vertices
   for ( uint_t i = 0; i < _vertices.size(); ++i )
   {
      if ( vertices_targetRank[i] == rank )
      {
         add_vertex( vtxs_ps, i );
      }
      else if ( vtxOnHalo[i] )
      {
         add_vertex( nbrVtxs_ps, i );
      }
   }

   auto add_edge = [&]( PrimitiveStorage::EdgeMap& map, const EdgeData& edge ) {
      constexpr uint_t K = 1;

      // vertex coordinates and IDs
      auto v      = edge.get_vertices();
      auto coords = edge.get_coordinates( _vertices );

      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         vertexIDs[i] = PrimitiveID( v[i] );
      }

      // add new edge
      auto id = edge.getPrimitiveID().getID();
      map[id] = std::make_shared< Edge >( edge.getPrimitiveID(), vertexIDs[0], vertexIDs[1], coords );

      // add properties
      map[id]->meshBoundaryFlag_ = edge.getBoundaryFlag();
      map[id]->geometryMap_      = _geometryMap.at( edge.getGeometryMap() );
   };

   // create local and halo edges
   for ( auto& edge : edges )
   {
      if ( edge.getTargetRank() == rank )
      {
         add_edge( edges_ps, edge );
      }
      else if ( edge.getHaloInfo() )
      {
         add_edge( nbrEdges_ps, edge );
      }
   }

   auto add_face = [&]( PrimitiveStorage::FaceMap& map, const FaceData& face ) {
      constexpr uint_t K = 2;

      // vertex coordinates and IDs
      auto v      = face.get_vertices();
      auto coords = face.get_coordinates( _vertices );

      std::array< PrimitiveID, K + 1 > vertexIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         vertexIDs[i] = PrimitiveID( v[i] );
      }

      // ordering of edges
      constexpr std::array< std::array< uint_t, 2 >, K + 1 > edgeOrder{ { { 0ul, 1ul }, { 0ul, 2ul }, { 1ul, 2ul } } };

      // edge IDs
      std::array< PrimitiveID, K + 1 > edgeIDs;
      for ( uint_t i = 0; i <= K; ++i )
      {
         edgeIDs[i] = vertexIDsToEdgeID[{ v[edgeOrder[i][0]], v[edgeOrder[i][1]] }];
      }

      // edge orientation
      std::array< int, K + 1 > edgeOrientation;
      for ( uint_t i = 0; i <= K; ++i )
      {
         EdgeData* edge     = get_edge( edgeIDs[i].getID() );
         auto      edgeVtx0 = edge->get_vertices()[0];
         // auto edgeVtx0 = get_edge( edgeIDs[i] )->get_vertices()[0];
         edgeOrientation[i] = ( edgeVtx0 == v[edgeOrder[i][0]] ) ? 1 : -1;
      }

      // add new face
      auto id = face.getPrimitiveID().getID();
      map[id] = std::make_shared< Face >( face.getPrimitiveID(), vertexIDs, edgeIDs, edgeOrientation, coords );

      // add properties
      map[id]->meshBoundaryFlag_ = face.getBoundaryFlag();
      map[id]->geometryMap_      = _geometryMap.at( face.getGeometryMap() );
   };

   // create local and halo faces
   for ( auto& face : faces )
   {
      if ( face.getTargetRank() == rank )
      {
         add_face( faces_ps, face );
      }
      else if ( face.getHaloInfo() )
      {
         add_face( nbrFaces_ps, face );
      }
   }

   auto add_cell = [&]( PrimitiveStorage::CellMap& map, const CellData& cell ) {
      constexpr uint_t K = 3;

      // vertex coordinates and IDs
      auto v      = cell.get_vertices();
      auto coords = cell.get_coordinates( _vertices );

      std::vector< PrimitiveID > vertexIDs( K + 1 );
      for ( uint_t i = 0; i <= K; ++i )
      {
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

      // ordering of faces
      constexpr std::array< std::array< uint_t, 3 >, K + 1 > faceOrder{
          { { 0ul, 1ul, 2ul }, { 0ul, 1ul, 3ul }, { 0ul, 2ul, 3ul }, { 1ul, 2ul, 3ul } } };

      // face IDs
      std::vector< PrimitiveID > faceIDs( K + 1 );
      for ( uint_t i = 0; i <= K; ++i )
      {
         faceIDs[i] = vertexIDsToFaceID[{ v[faceOrder[i][0]], v[faceOrder[i][1]], v[faceOrder[i][2]] }];
      }

      std::array< std::map< uint_t, uint_t >, 6 > edgeLocalVertexToCellLocalVertexMaps;

      // edgeLocalVertexToCellLocalVertexMaps[ cellLocalEdgeID ][ edgeLocalVertexID ] = cellLocalVertexID;

      for ( uint_t i = 0; i < 6; ++i )
      {
         for ( const auto& j : edgeOrder[i] )
         {
            EdgeData* edge         = get_edge( edgeIDs[i].getID() );
            auto&     edgeVtxs     = edge->get_vertices();
            auto      vtxIdxOnEdge = std::find( edgeVtxs.begin(), edgeVtxs.end(), vertexIDs[j].getID() ) - edgeVtxs.begin();

            WALBERLA_ASSERT_LESS( vtxIdxOnEdge, edgeVtxs.size() );

            edgeLocalVertexToCellLocalVertexMaps[i][vtxIdxOnEdge] = j;
         }
      }

      std::array< std::map< uint_t, uint_t >, 4 > faceLocalVertexToCellLocalVertexMaps;

      // faceLocalVertexToCellLocalVertexMaps[ cellLocalFaceID ][ faceLocalVertexID ] = cellLocalVertexID;

      for ( uint_t i = 0; i < K + 1; ++i )
      {
         for ( auto& j : faceOrder[i] )
         {
            FaceData* face         = get_face( faceIDs[i].getID() );
            auto&     faceVtxs     = face->get_vertices();
            auto      vtxIdxOnFace = std::find( faceVtxs.begin(), faceVtxs.end(), vertexIDs[j].getID() ) - faceVtxs.begin();

            WALBERLA_ASSERT_LESS( vtxIdxOnFace, faceVtxs.size() );

            faceLocalVertexToCellLocalVertexMaps[i][vtxIdxOnFace] = j;
         }
      }

      // add new cell
      auto id = cell.getPrimitiveID().getID();
      map[id] = std::make_shared< Cell >( cell.getPrimitiveID(),
                                          vertexIDs,
                                          edgeIDs,
                                          faceIDs,
                                          coords,
                                          edgeLocalVertexToCellLocalVertexMaps,
                                          faceLocalVertexToCellLocalVertexMaps );

      // add properties
      map[id]->meshBoundaryFlag_ = cell.getBoundaryFlag();
      map[id]->geometryMap_      = _geometryMap.at( cell.getGeometryMap() );
   };

   // create local and halo cells
   for ( auto& cell : cells )
   {
      if ( cell.getTargetRank() == rank )
      {
         add_cell( cells_ps, cell );
      }
      else if ( cell.getHaloInfo() )
      {
         add_cell( nbrCells_ps, cell );
      }
   }

   // ****** add neighborhood information to primitives ******

   // add neighbor edges to vertices
   for ( auto& edge : edges )
   {
      auto& v = edge.get_vertices();

      for ( auto& vtx : v )
      {
         if ( vertices_targetRank[vtx] == rank )
         {
            vtxs_ps[vtx]->addEdge( edge.getPrimitiveID() );

            WALBERLA_ASSERT( edge.getHaloInfo() || edge.getTargetRank() == rank );
         }
         if ( vtxOnHalo[vtx] )
         {
            nbrVtxs_ps[vtx]->addEdge( edge.getPrimitiveID() );
         }
      }
   }

   // add neighbor faces to vertices and edges
   for ( auto& face : faces )
   {
      auto& v = face.get_vertices();

      for ( auto& vtx : v )
      {
         if ( vertices_targetRank[vtx] == rank )
         {
            vtxs_ps[vtx]->addFace( face.getPrimitiveID() );
            WALBERLA_ASSERT( face.getHaloInfo() || face.getTargetRank() == rank );
         }
         if ( vtxOnHalo[vtx] )
         {
            nbrVtxs_ps[vtx]->addFace( face.getPrimitiveID() );
         }
      }

      for ( uint_t i = 0; i < 3; ++i )
      {
         Idx< 2 >  edgeVtxs{ v[i], v[( i + 1 ) % 3] };
         EdgeData* edge = get_edge( vertexIDsToEdgeID[edgeVtxs].getID() );

         if ( edge->getTargetRank() == rank )
         {
            edges_ps[edge->getPrimitiveID().getID()]->addFace( face.getPrimitiveID() );
            WALBERLA_ASSERT( face.getHaloInfo() || face.getTargetRank() == rank );
         }
         if ( edge->getHaloInfo() )
         {
            nbrEdges_ps[edge->getPrimitiveID().getID()]->addFace( face.getPrimitiveID() );
         }
      }
   }

   // add neighbor cells to vertices, edges and faces
   for ( auto& cell : cells )
   {
      auto& v = cell.get_vertices();

      for ( auto& vtx : v )
      {
         if ( vertices_targetRank[vtx] == rank )
         {
            vtxs_ps[vtx]->addCell( cell.getPrimitiveID() );
            WALBERLA_ASSERT( cell.getHaloInfo() || cell.getTargetRank() == rank );
         }
         if ( vtxOnHalo[vtx] )
         {
            nbrVtxs_ps[vtx]->addCell( cell.getPrimitiveID() );
         }
      }

      for ( uint_t i = 0; i < 4; ++i )
      {
         for ( uint_t j = i + 1; j < 4; ++j )
         {
            Idx< 2 >  edgeVtxs{ v[i], v[j] };
            EdgeData* edge = get_edge( vertexIDsToEdgeID[edgeVtxs].getID() );

            if ( edge->getTargetRank() == rank )
            {
               edges_ps[edge->getPrimitiveID().getID()]->addCell( cell.getPrimitiveID() );
               WALBERLA_ASSERT( cell.getHaloInfo() || cell.getTargetRank() == rank );
            }
            if ( edge->getHaloInfo() )
            {
               nbrEdges_ps[edge->getPrimitiveID().getID()]->addCell( cell.getPrimitiveID() );
            }
         }
      }

      for ( uint_t j = 0; j < 4; ++j )
      {
         Idx< 3 >  faceVtxs{ v[j], v[( j + 1 ) % 4], v[( j + 2 ) % 4] };
         FaceData* face = get_face( vertexIDsToFaceID[faceVtxs].getID() );

         if ( face->getTargetRank() == rank )
         {
            faces_ps[face->getPrimitiveID().getID()]->addCell( cell.getPrimitiveID() );
            WALBERLA_ASSERT( cell.getHaloInfo() || cell.getTargetRank() == rank );
         }
         if ( face->getHaloInfo() )
         {
            nbrFaces_ps[face->getPrimitiveID().getID()]->addCell( cell.getPrimitiveID() );
         }
      }
   }

   // ****** create PrimitiveStorage ******

   return std::make_shared< PrimitiveStorage >(
       vtxs_ps, edges_ps, faces_ps, cells_ps, nbrVtxs_ps, nbrEdges_ps, nbrFaces_ps, nbrCells_ps, nbrRanks, cells.size() > 0 );
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::extract_data( std::vector< EdgeData >& edgeData,
                                        std::vector< FaceData >& faceData,
                                        std::vector< CellData >& cellData ) const
{
   std::set< std::shared_ptr< Simplex1 > > edges;
   std::set< std::shared_ptr< Simplex2 > > faces;
   std::set< std::shared_ptr< Simplex3 > > cells;

   // collect cells
   if constexpr ( std::is_same_v< K_Simplex, Simplex3 > )
   {
      cells = _T;
   }
   // collect faces
   if constexpr ( std::is_same_v< K_Simplex, Simplex2 > )
   {
      faces = _T;
   }
   for ( auto& cell : cells )
   {
      for ( auto& face : cell->get_faces() )
      {
         faces.insert( face );
      }
   }
   // collect edges
   for ( auto& face : faces )
   {
      for ( auto& edge : face->get_edges() )
      {
         edges.insert( edge );
      }
   }

   // collect data and add PrimitiveIDs
   uint_t id = _n_vertices;
   // collect edgedata
   for ( auto& edge : edges )
   {
      edge->setPrimitiveID( id );
      edgeData.push_back( EdgeData( edge.get() ) );
      ++id;
   }
   // collect facedata
   for ( auto& face : faces )
   {
      face->setPrimitiveID( id );
      faceData.push_back( FaceData( face.get() ) );
      ++id;
   }
   // collect celldata
   for ( auto& cell : cells )
   {
      cell->setPrimitiveID( id );
      cellData.push_back( CellData( cell.get() ) );
      ++id;
   }
}

template < class K_Simplex >
inline std::vector< std::shared_ptr< K_Simplex > >
    K_Mesh< K_Simplex >::init_R( const std::vector< PrimitiveID >& primitiveIDs ) const
{
   auto to_add = [&]( std::shared_ptr< K_Simplex > el, uint_t i ) -> bool {
      if ( el->has_children() )
      {
         // add element if it has a child with id=id
         for ( auto& child : el->get_children() )
         {
            if ( child->getPrimitiveID() == primitiveIDs[i] )
            {
               return true;
            }
         }
         return false;
      }
      else
      {
         // add element if id=id
         return ( el->getPrimitiveID() == primitiveIDs[i] );
      }
   };

   std::vector< std::shared_ptr< K_Simplex > > R( primitiveIDs.size() );
   for ( auto& el : _T )
   {
      for ( uint_t i = 0; i < R.size(); ++i )
      {
         if ( to_add( el, i ) )
         {
            R[i] = el;
            break;
         }
      }
   }

   return R;
}

template < class K_Simplex >
std::set< std::shared_ptr< K_Simplex > > K_Mesh< K_Simplex >::refine_red( const std::set< std::shared_ptr< K_Simplex > >& R,
                                                                          std::set< std::shared_ptr< K_Simplex > >&       U )
{
   std::set< std::shared_ptr< K_Simplex > > refined;

   for ( auto& el : R )
   {
      // mark el as processed
      if ( U.erase( el ) == 0 )
      {
         // for el ∉ U: don't try to refine
         continue;
      }

      // remove green edges
      bool check_subelements = el->kill_children();

      // apply regular refinement to el
      std::set< std::shared_ptr< K_Simplex > > subelements;
      if constexpr ( std::is_same_v< K_Simplex, Simplex2 > )
      {
         subelements = refine_face_red( _vertices, _vertexGeometryMap, _vertexBoundaryFlag, el );
      }
      if constexpr ( std::is_same_v< K_Simplex, Simplex3 > )
      {
         subelements = refine_cell_red( _vertices, _vertexGeometryMap, _vertexBoundaryFlag, el );
      }

      // if this red step only replaced a green step, subelements must be checked again
      if ( check_subelements )
      {
         U.merge( subelements );
      }
      else
      {
         refined.merge( subelements );
      }
   }

   return refined;
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::remove_green_edges()
{
   auto T_cpy = _T;

   for ( auto& el : T_cpy )
   {
      if ( el->has_green_edge() )
      {
         _T.erase( el );
         _T.insert( el->get_parent() );
      }
   }
}

template <>
std::set< std::shared_ptr< Simplex2 > >
    K_Mesh< Simplex2 >::find_elements_for_red_refinement( const std::set< std::shared_ptr< Simplex2 > >& U )
{
   std::set< std::shared_ptr< Simplex2 > > R;

   for ( auto& el : U )
   {
      if ( el->vertices_on_edges() > 1 )
      {
         R.insert( el );
      }
   }

   return R;
}

template <>
std::set< std::shared_ptr< Simplex3 > >
    K_Mesh< Simplex3 >::find_elements_for_red_refinement( const std::set< std::shared_ptr< Simplex3 > >& U )
{
   std::set< std::shared_ptr< Simplex3 > > R;

   for ( auto& el : U )
   {
      uint_t n_red = 0;

      for ( auto& face : el->get_faces() )
      {
         if ( face->vertices_on_edges() > 1 )
         {
            if ( face->get_children().size() == 2 )
            {
               // remove green edge from face
               face->kill_children();
            }

            if ( !face->has_children() )
            {
               // apply red refinement to face
               refine_face_red( _vertices, _vertexGeometryMap, _vertexBoundaryFlag, face );
            }

            ++n_red;
         }
      }

      // if more than one face has been red-refined, mark cell for red refinement
      if ( n_red > 1 )
      {
         R.insert( el );
      }
   }

   return R;
}

template < class K_Simplex >
uint_t K_Mesh< K_Simplex >::predict_n_el_green( const std::set< std::shared_ptr< K_Simplex > >& U ) const
{
   uint_t n_el_green = 0;

   for ( auto& el : U )
   {
      uint_t new_vertices = el->vertices_on_edges();

      n_el_green += ( new_vertices == 2 ) ? 3 : new_vertices;
   }

   return n_el_green;
}

template <>
std::set< std::shared_ptr< Simplex2 > > K_Mesh< Simplex2 >::refine_green( std::set< std::shared_ptr< Simplex2 > >& U )
{
   std::set< std::shared_ptr< Simplex2 > > refined;
   std::set< std::shared_ptr< Simplex2 > > U_cpy = U;

   for ( auto& el : U_cpy )
   {
      // count number of new vertices on the edges of el
      uint_t new_vertices = el->vertices_on_edges();

      if ( new_vertices > 0 )
      {
         WALBERLA_ASSERT( !el->has_green_edge() );
         WALBERLA_ASSERT( new_vertices == 1 );

         /* if green refinement had been applied to the same element before,
            nothing has to be done
         */
         if ( el->has_children() )
         {
            for ( auto& child : el->get_children() )
            {
               refined.insert( child );
            }
         }
         else
         {
            refined.merge( refine_face_green( el ) );
         }

         // mark el as processed
         U.erase( el );
      }
   }

   return refined;
}

template <>
std::set< std::shared_ptr< Simplex3 > > K_Mesh< Simplex3 >::refine_green( std::set< std::shared_ptr< Simplex3 > >& U )
{
   std::set< std::shared_ptr< Simplex3 > > refined;
   std::set< std::shared_ptr< Simplex3 > > U_cpy = U;

   auto keepChildren = [&]( std::shared_ptr< Simplex3 > el ) {
      for ( auto& child : el->get_children() )
      {
         refined.insert( child );
      }
   };

   for ( auto& el : U_cpy )
   {
      uint_t new_vertices = el->vertices_on_edges();

      switch ( new_vertices )
      {
      case 0:
         continue;
         break;

      case 1:
         if ( el->has_children() )
         {
            WALBERLA_ASSERT( el->get_children().size() == 2 );
            keepChildren( el );
         }
         else
         {
            refined.merge( refine_cell_green_1( el ) );
         }
         break;

      case 2:
         if ( el->has_children() && el->get_children().size() == 4 )
         {
            keepChildren( el );
         }
         else
         {
            WALBERLA_ASSERT( el->get_children().size() == 0 || el->get_children().size() == 2 );
            el->kill_children();
            refined.merge( refine_cell_green_2( el ) );
         }
         break;

      case 3:
         if ( el->has_children() && el->get_children().size() == 4 )
         {
            keepChildren( el );
         }
         else
         {
            WALBERLA_ASSERT( el->get_children().size() == 0 || el->get_children().size() == 2 );
            el->kill_children();
            refined.merge( refine_cell_green_3( el ) );
         }
         break;

      default:
         WALBERLA_ASSERT( new_vertices <= 3 );
         break;
      }

      // mark el as processed
      if ( new_vertices > 0 )
      {
         U.erase( el );
      }
   }

   return refined;
}

template < class K_Simplex >
std::pair< real_t, real_t > K_Mesh< K_Simplex >::min_max_angle() const
{
   std::pair< real_t, real_t > mm{ 10, 0 };

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         auto mm_el = el->min_max_angle( _vertices );

         mm.first  = std::min( mm.first, mm_el.first );
         mm.second = std::max( mm.second, mm_el.second );
      }
   }

   walberla::mpi::broadcastObject( mm );

   return mm;
}

template < class K_Simplex >
real_t K_Mesh< K_Simplex >::volume() const
{
   real_t v_tot = 0;

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         v_tot += el->volume( _vertices );
      }
   }

   walberla::mpi::broadcastObject( v_tot );

   return v_tot;
}

template class K_Mesh< Simplex2 >;
template class K_Mesh< Simplex3 >;

} // namespace adaptiveRefinement
} // namespace hyteg
