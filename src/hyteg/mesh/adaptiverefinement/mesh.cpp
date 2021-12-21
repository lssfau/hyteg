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

#include "mesh.hpp"

#include <algorithm>
#include <assert.h>
#include <core/mpi/Broadcast.h>
#include <iostream>
#include <map>
#include <type_traits>

#include "convertToSetupStorage.hpp"
#include "refine_cell.hpp"
#include "simplexFactory.hpp"

namespace hyteg {
namespace adaptiveRefinement {

template < class K_Simplex >
K_Mesh< K_Simplex >::K_Mesh( const SetupPrimitiveStorage& setupStorage )
: _setupStorage( setupStorage )
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

      // todo handle boundary flag
      for ( auto& [id, edge] : setupStorage.getEdges() )
      {
         edge->getNeighborVertices( v );
         fac.useGeometryMap( id );
         auto myEdge = fac.make_edge( conversion[v[0].getID()], conversion[v[1].getID()] );
         myEdge->setPrimitiveID( id );
      }

      for ( auto& [id, face] : setupStorage.getFaces() )
      {
         face->getNeighborVertices( v );
         fac.useGeometryMap( id );
         auto myFace = fac.make_face( conversion[v[0].getID()], conversion[v[1].getID()], conversion[v[2].getID()] );
         myFace->setPrimitiveID( id );
      }

      for ( auto& [id, cell] : setupStorage.getCells() )
      {
         cell->getNeighborVertices( v );
         fac.useGeometryMap( id );
         auto myCell = fac.make_cell(
             conversion[v[0].getID()], conversion[v[1].getID()], conversion[v[2].getID()], conversion[v[3].getID()] );
         myCell->setPrimitiveID( id );
         cells.insert( myCell );
      }

      // insert volume elements into _T
      if constexpr ( std::is_same_v< K_Simplex, Simplex2 > )
      {
         if ( !cells.empty() )
         {
            WALBERLA_ABORT( "Adaptive 2D mesh requires MeshInfo without any cells!" );
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
            WALBERLA_ABORT( "Adaptive 3D mesh requires MeshInfo containing at least one cell!" );
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
SetupPrimitiveStorage& K_Mesh< K_Simplex >::refineRG( const std::vector< PrimitiveID >& elements_to_refine )
{
   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      auto R = init_R( elements_to_refine );
      // remove green edges
      remove_green_edges( R );

      /* recursively apply red refinement for elements
        that otherwise would be subject to multiple
        green refinement steps
      */
      std::set< std::shared_ptr< K_Simplex > > U = _T;
      std::set< std::shared_ptr< K_Simplex > > refined;
      while ( !R.empty() )
      {
         refined.merge( refine_red( R, U ) );

         R = find_elements_for_red_refinement( U );
      }

      // apply green refinement
      refined.merge( refine_green( U ) );

      // update current configuration
      _T = U;
      _T.merge( refined );
      _n_vertices = _vertices.size();
      _n_elements = _T.size();
   }

   WALBERLA_LOG_INFO( ">>>>>>>> refinement done <<<<<<<<<<<<<<\n\n" );

   // extract connectivity, geometry and boundary data
   EdgeData edges;
   FaceData faces;
   CellData cells;
   extract_data( edges, faces, cells );

   WALBERLA_LOG_INFO( ">>>>>>>> extraction done <<<<<<<<<<<<<<\n\n" );

   // broadcast data to all processes
   walberla::mpi::broadcastObject( _n_elements );
   walberla::mpi::broadcastObject( _n_vertices );
   walberla::mpi::broadcastObject( _vertices );
   walberla::mpi::broadcastObject( _vertexGeometryMap );
   walberla::mpi::broadcastObject( _vertexBoundaryFlag );
   edges.broadcast();
   faces.broadcast();
   cells.broadcast();

   WALBERLA_LOG_INFO( ">>>>>>>> broadcast done <<<<<<<<<<<<<<\n\n" );

   // update storage
   auto id = updateSetupStorage( edges, faces, cells, _setupStorage.getNumberOfProcesses() );

   // todo clear vertex data on rank!=0

   WALBERLA_LOG_INFO( ">>>>>>>> storage done <<<<<<<<<<<<<<\n\n" );

   // insert PrimitiveIDs to volume elements
   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         el->setPrimitiveID( id );
         ++id;
      }
   }

   return _setupStorage;
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::extract_data( EdgeData& edgeData, FaceData& faceData, CellData& cellData ) const
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

   // collect celldata
   for ( auto& cell : cells )
   {
      cellData.add( cell.get() );
   }
   // collect facedata
   for ( auto& face : faces )
   {
      faceData.add( face.get() );
   }
   // collect edgedata
   for ( auto& edge : edges )
   {
      edgeData.add( edge.get() );
   }
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::collect_elements( std::set< std::shared_ptr< Simplex3 > >& cells,
                                            std::set< std::shared_ptr< Simplex2 > >& faces,
                                            std::set< std::shared_ptr< Simplex1 > >& edges ) const
{
   cells.clear();
   faces.clear();
   edges.clear();

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
}

template <>
inline hyteg::MeshInfo K_Mesh< Simplex2 >::export_meshInfo()
{
   std::vector< std::array< uint_t, 3 > > faces;

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         auto& v = el->get_vertices();
         faces.push_back( { uint_t( v[0] ), uint_t( v[1] ), uint_t( v[2] ) } );
      }
   }

   walberla::mpi::broadcastObject( faces );
   walberla::mpi::broadcastObject( _vertices );

   return MeshInfo::fromFaceData( _vertices, faces );
}

template <>
inline hyteg::MeshInfo K_Mesh< Simplex3 >::export_meshInfo()
{
   std::vector< std::array< uint_t, 4 > > cells;

   if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   {
      for ( auto& el : _T )
      {
         auto& v = el->get_vertices();
         cells.push_back( { uint_t( v[0] ), uint_t( v[1] ), uint_t( v[2] ), uint_t( v[3] ) } );
      }
   }

   walberla::mpi::broadcastObject( cells );
   walberla::mpi::broadcastObject( _vertices );

   return MeshInfo::fromCellData( _vertices, cells );
}

template < class K_Simplex >
inline std::set< std::shared_ptr< K_Simplex > >
    K_Mesh< K_Simplex >::init_R( const std::vector< PrimitiveID >& primitiveIDs ) const
{
   std::set< std::shared_ptr< K_Simplex > > R;

   auto barycenter_R = compute_barycenters( primitiveIDs );

   for ( auto& el : _T )
   {
      for ( const auto& bc : barycenter_R )
      {
         // ||bc_el - bc||/||bc||
         auto err = ( el->barycenter( _vertices ) - bc ).norm() / bc.norm();
         if ( err < 1e-14 )
         {
            R.insert( el );
            break;
         }
      }
   }

   return R;
}

template <>
inline std::vector< Point3D > K_Mesh< Simplex2 >::compute_barycenters( const std::vector< PrimitiveID >& primitiveIDs ) const
{
   std::vector< Point3D > barycenter_R( primitiveIDs.size() );

   for ( uint_t i = 0; i < primitiveIDs.size(); ++i )
   {
      auto face = _setupStorage.getFace( primitiveIDs[i] );
      if ( face == nullptr )
      {
         WALBERLA_ABORT( "All primitiveIDs for refineRG() must correspond to a face in Mesh::setupStorage()!" );
      }
      barycenter_R[i] = Simplex2::barycenter( face->getCoordinates() );
   }

   return barycenter_R;
}

template <>
inline std::vector< Point3D > K_Mesh< Simplex3 >::compute_barycenters( const std::vector< PrimitiveID >& primitiveIDs ) const
{
   std::vector< Point3D > barycenter_R( primitiveIDs.size() );

   for ( uint_t i = 0; i < primitiveIDs.size(); ++i )
   {
      auto cell = _setupStorage.getCell( primitiveIDs[i] );
      if ( cell == nullptr )
      {
         WALBERLA_ABORT( "All primitiveIDs for refineRG() must correspond to a cell in Mesh::setupStorage()!" );
      }
      barycenter_R[i] = Simplex3::barycenter( cell->getCoordinates() );
   }

   return barycenter_R;
}

template < class K_Simplex >
std::set< std::shared_ptr< K_Simplex > > K_Mesh< K_Simplex >::refine_red( const std::set< std::shared_ptr< K_Simplex > >& R,
                                                                          std::set< std::shared_ptr< K_Simplex > >&       U )
{
   std::set< std::shared_ptr< K_Simplex > > refined;

   for ( auto& el : R )
   {
      // remove green edges
      bool check_subelements = el->kill_children();

      auto subelements = refine_element_red( el );

      // mark el as processed
      U.erase( el );

      // mark subelements as unprocessed if necessary
      if ( check_subelements )
      {
         U.merge( subelements );
      }

      refined.merge( subelements );
   }

   return refined;
}

template < class K_Simplex >
void K_Mesh< K_Simplex >::remove_green_edges( std::set< std::shared_ptr< K_Simplex > >& R )
{
   auto T_cpy = _T;

   for ( auto& el : T_cpy )
   {
      if ( el->has_green_edge() )
      {
         _T.erase( el );
         _T.insert( el->get_parent() );

         if ( R.erase( el ) )
         {
            R.insert( el->get_parent() );
         }
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
               refine_face_red( _vertices, _vertexGeometryMap, face );
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
      U.erase( el );
   }

   return refined;
}

template <>
std::set< std::shared_ptr< Simplex3 > > K_Mesh< Simplex3 >::refine_element_red( std::shared_ptr< Simplex3 > element )
{
   return refine_cell_red( _vertices, _vertexGeometryMap, element );
}

template <>
std::set< std::shared_ptr< Simplex2 > > K_Mesh< Simplex2 >::refine_element_red( std::shared_ptr< Simplex2 > element )
{
   return refine_face_red( _vertices, _vertexGeometryMap, element );
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

// SimplexData

template < class K_Simplex >
template < size_t J >
inline void K_Mesh< K_Simplex >::SimplexData< J >::broadcast()
{
   walberla::mpi::broadcastObject( this->vertices );
   walberla::mpi::broadcastObject( this->geometryMap );
   walberla::mpi::broadcastObject( this->boundaryFlag );
}

} // namespace adaptiveRefinement
} // namespace hyteg
