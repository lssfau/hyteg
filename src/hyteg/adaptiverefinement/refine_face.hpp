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

#include <set>

#include "refine_edge.hpp"

namespace hyteg {
namespace adaptiveRefinement {

/* apply red refinement to face and add required vertices to vertices
      2
      |\
      | \
      |  \
   [4]|   \ [3]
      |    \
      | (2) \
      |      \
      |       \
      |        \
      |   [7]   \
      5|----------4
      |\        | \
      | \       |  \
      |  \  (3) |   \
      |   \[8]  |[6] \ [2]
   [5]|    \    |     \
      | (0) \   | (1)  \
      |      \  |       \
      |       \ |        \
      0---------3---------1
         [0]        [1]
   @param vertices      global coordinates of all vertices in the mesh
   @param vtxData       geometrymap, boundaryflag etc. of all vertices in the mesh
   @param face          subject to refinement
   @return sub-elements
*/
inline std::set< std::shared_ptr< Simplex2 > >
    refine_face_red( std::vector< Point3D >& vertices, std::vector< VertexData >& vtxData, std::shared_ptr< Simplex2 > face )
{
   // === prepare sets of vertices and edges for face->split() ===

   // collect old vertices ======================

   std::vector< uint_t > ref_vertices( 6 );
   for ( uint_t i = 0; i < 3; ++i )
   {
      ref_vertices[i] = face->get_vertices()[i];
   }

   // bisect edges ======================

   std::vector< std::shared_ptr< Simplex1 > > ref_edges( 9 );

   for ( uint_t i = 0; i < 3; ++i )
   {
      auto edge = face->get_edges()[i];

      if ( edge->has_children() )
      {
         ref_vertices[3 + i] = edge->get_midpoint_idx();
      }
      else
      {
         ref_vertices[3 + i] = bisect_edge( vertices, vtxData, edge );
      }

      auto child_edges     = edge->get_children_sorted( { ref_vertices[i], ref_vertices[( i + 1 ) % 3] } );
      ref_edges[2 * i + 0] = child_edges[0];
      ref_edges[2 * i + 1] = child_edges[1];
   }

   // add inner edges ======================

   auto childIDs = face->getPrimitiveID().createChildren();

   for ( uint_t i = 0; i < 3; ++i )
   {
      uint_t i0        = 3 + i;
      uint_t i1        = 3 + ( ( i + 1 ) % 3 );
      ref_edges[6 + i] = std::make_shared< Simplex1 >( ref_vertices[i0],
                                                       ref_vertices[i1],
                                                       nullptr,
                                                       RED,
                                                       face->getGeometryMap(),
                                                       face->getBoundaryFlag(),
                                                       face->getTargetRank() );
      ref_edges[6 + i]->setPrimitiveID( childIDs[i] );
   }

   // === split face ===

   // face_0
   std::array< uint_t, 3 >                      v0{ ref_vertices[0], ref_vertices[3], ref_vertices[5] };
   std::array< std::shared_ptr< Simplex1 >, 3 > e0{ ref_edges[0], ref_edges[8], ref_edges[5] };
   face->add_child( std::make_shared< Simplex2 >( v0, e0, face ) );
   // face_1
   std::array< uint_t, 3 >                      v1{ ref_vertices[3], ref_vertices[1], ref_vertices[4] };
   std::array< std::shared_ptr< Simplex1 >, 3 > e1{ ref_edges[1], ref_edges[2], ref_edges[6] };
   face->add_child( std::make_shared< Simplex2 >( v1, e1, face ) );
   // face_2
   std::array< uint_t, 3 >                      v2{ ref_vertices[5], ref_vertices[4], ref_vertices[2] };
   std::array< std::shared_ptr< Simplex1 >, 3 > e2{ ref_edges[7], ref_edges[3], ref_edges[4] };
   face->add_child( std::make_shared< Simplex2 >( v2, e2, face ) );
   // face_3
   std::array< uint_t, 3 >                      v3{ ref_vertices[4], ref_vertices[5], ref_vertices[3] };
   std::array< std::shared_ptr< Simplex1 >, 3 > e3{ ref_edges[7], ref_edges[8], ref_edges[6] };
   face->add_child( std::make_shared< Simplex2 >( v3, e3, face ) );

   auto& children = face->get_children();
   for ( uint_t i = 0; i < children.size(); ++i )
   {
      children[i]->setPrimitiveID( childIDs[i + 3] );
   }

   return std::set< std::shared_ptr< Simplex2 > >( face->get_children().begin(), face->get_children().end() );
}

/* apply green refinement to face
      0
         |\
         | \
         |  \ [3]
         |   \
         |    \
      [0]| (0) \3
         |    / \
         |   /   \
         |  /[4]  \ [2]
         | /       \
         |/   (1)   \
         1-----------2
               [1]

      @return sub-elements
   */
inline std::set< std::shared_ptr< Simplex2 > > refine_face_green( std::shared_ptr< Simplex2 > face )
{
   // === prepare sets of vertices and edges for face->split() ===

   // get index of split edge ======================
   int e_split = -1;
   for ( uint_t i = 0; i < 3; ++i )
   {
      if ( face->get_edges()[i]->has_children() )
      {
         e_split = int( i );
         break;
      }
   }
   WALBERLA_ASSERT( e_split >= 0 );
   auto edge = face->get_edges()[uint_t( e_split )];

   // index offset from the standard configuration (e_split = 2)
   const uint_t offset = uint_t( 1 + e_split );

   // get vertex indices ======================

   std::vector< uint_t > ref_vertices( 4 );

   for ( uint_t i = 0; i < 3; ++i )
   {
      ref_vertices[i] = face->get_vertices()[( i + offset ) % 3];
   }
   ref_vertices[3] = edge->get_midpoint_idx();

   // get edges ======================

   std::vector< std::shared_ptr< Simplex1 > > ref_edges( 5 );

   // old edges
   ref_edges[0] = face->get_edges()[( 0 + offset ) % 3];
   ref_edges[1] = face->get_edges()[( 1 + offset ) % 3];

   // child edges of split edge
   auto child_edges = edge->get_children_sorted( { ref_vertices[2], ref_vertices[0] } );
   ref_edges[2]     = child_edges[0];
   ref_edges[3]     = child_edges[1];

   // add new edge
   ref_edges[4] = std::make_shared< Simplex1 >(
       ref_vertices[1], ref_vertices[3], nullptr, GREEN, face->getGeometryMap(), face->getBoundaryFlag(), face->getTargetRank() );

   // === split face ===

   // face_0
   std::array< uint_t, 3 >                      v0{ ref_vertices[0], ref_vertices[1], ref_vertices[3] };
   std::array< std::shared_ptr< Simplex1 >, 3 > e0{ ref_edges[0], ref_edges[4], ref_edges[3] };
   face->add_child( std::make_shared< Simplex2 >( v0, e0, face ) );
   // face_1
   std::array< uint_t, 3 >                      v1{ ref_vertices[1], ref_vertices[2], ref_vertices[3] };
   std::array< std::shared_ptr< Simplex1 >, 3 > e1{ ref_edges[1], ref_edges[2], ref_edges[4] };
   face->add_child( std::make_shared< Simplex2 >( v1, e1, face ) );

   // add IDs to children
   auto childIDs = face->getPrimitiveID().createChildren();
   ref_edges[4]->setPrimitiveID( childIDs[0] );
   auto& children = face->get_children();
   for ( uint_t i = 0; i < children.size(); ++i )
   {
      children[i]->setPrimitiveID( childIDs[i + 3] );
   }

   return std::set< std::shared_ptr< Simplex2 > >( face->get_children().begin(), face->get_children().end() );
}

} // namespace adaptiveRefinement
} // namespace hyteg
