/*
 * Copyright (c) 2021-2024 Benjamin Mann
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

#include "refine_face.hpp"
#include "simplexFactory.hpp"

namespace hyteg {
namespace adaptiveRefinement {

/* apply red refinement to cell and add required vertices
   @param coords     global coordinates of all vertices in the mesh
   @param vtxs       geometrymap, boundaryflag etc. of all vertices in the mesh
   @param cell          subject to refinement
   @return sub-elements
*/
inline std::set< std::shared_ptr< Simplex3 > >
    refine_cell_red( EnumeratedList< Point3D >& coords, EnumeratedList< VertexData >& vtxs, std::shared_ptr< Simplex3 > cell )
{
   // === split faces ===
   for ( auto& face : cell->get_faces() )
   {
      if ( !face->has_children() )
      {
         refine_face_red( coords, vtxs, face );
      }
   }

   // === collect vertex indices ===
   std::vector< uint_t > ref_vertices( 10 ); // [0 1 2 3 01 02 03 12 13 23]
   // map edge-vertices to ref-idx of midpoint, e.g., 02->5
   auto refIdx = []( uint_t i, uint_t j ) {
      if ( i * j == 0 )
         return 3 + i + j;
      else
         return 4 + i + j;
   };
   const auto _01 = refIdx( 0, 1 );
   const auto _02 = refIdx( 0, 2 );
   const auto _03 = refIdx( 0, 3 );
   const auto _12 = refIdx( 1, 2 );
   const auto _13 = refIdx( 1, 3 );
   const auto _23 = refIdx( 2, 3 );

   // old vertices
   for ( uint_t i = 0; i < 4; ++i )
   {
      ref_vertices[i] = cell->get_vertices()[i];
   }
   // new vertices
   for ( uint_t i = 0; i < 4; ++i )
   {
      for ( uint_t j = i + 1; j < 4; ++j )
      {
         ref_vertices[refIdx( i, j )] = cell->get_Edge( ref_vertices[i], ref_vertices[j] )->get_midpoint_idx();
         ;
      }
   }

   SimplexFactory fac( cell, ref_vertices );

   // === collect edges ===

   // edge-children
   for ( uint_t i = 0; i < 4; ++i )
   {
      for ( uint_t j = i + 1; j < 4; ++j )
      {
         auto ij         = refIdx( i, j );
         auto parentEdge = cell->get_Edge( ref_vertices[i], ref_vertices[j] );                      // [i,j]
         auto children   = parentEdge->get_children_sorted( { ref_vertices[i], ref_vertices[j] } ); // {[i,ij], [ij,j]}
         fac.add_edge( i, ij, children[0] );
         fac.add_edge( j, ij, children[1] );
      }
   }

   // edges on face interior
   for ( uint_t i = 0; i < 4; ++i )
   {
      for ( uint_t j = i + 1; j < 4; ++j )
      {
         for ( uint_t k = j + 1; k < 4; ++k )
         {
            auto ij         = refIdx( i, j );
            auto jk         = refIdx( j, k );
            auto ki         = refIdx( k, i );
            auto parentFace = cell->get_Face( ref_vertices[i], ref_vertices[j], ref_vertices[k] ); // [i,j,k]
            auto innerChild = parentFace->get_children()[3];                                       // [ij, jk, ki]
            fac.add_edge( ij, jk, innerChild->get_Edge( ref_vertices[ij], ref_vertices[jk] ) );
            fac.add_edge( jk, ki, innerChild->get_Edge( ref_vertices[jk], ref_vertices[ki] ) );
            fac.add_edge( ki, ij, innerChild->get_Edge( ref_vertices[ki], ref_vertices[ij] ) );
         }
      }
   }

   // === new edge in cell interior ===
   // following Bey
   auto newEdge = fac.make_edge( _02, _13 );

   // === collect faces ===

   for ( uint_t i = 0; i < 4; ++i )
   {
      for ( uint_t j = i + 1; j < 4; ++j )
      {
         for ( uint_t k = j + 1; k < 4; ++k )
         {
            auto ij         = refIdx( i, j );
            auto jk         = refIdx( j, k );
            auto ki         = refIdx( k, i );
            auto parentFace = cell->get_Face( ref_vertices[i], ref_vertices[j], ref_vertices[k] ); // [i,j,k]
            auto childFaces = parentFace->get_children_sorted( { ref_vertices[i], ref_vertices[j], ref_vertices[k] } );
            fac.add_face( i, ij, ki, childFaces[0] );
            fac.add_face( j, ij, jk, childFaces[1] );
            fac.add_face( k, ki, jk, childFaces[2] );
            fac.add_face( ij, jk, ki, childFaces[3] );
         }
      }
   }

   // === new inner faces ===

   // inner faces "cutting off" the parent cells vertices
   std::array< std::shared_ptr< Simplex2 >, 8 > newFaces;
   newFaces[0] = fac.make_face( _01, _02, _03 );
   newFaces[1] = fac.make_face( _01, _12, _13 );
   newFaces[2] = fac.make_face( _02, _12, _23 );
   newFaces[3] = fac.make_face( _03, _13, _23 );

   // inner faces connecting the new inner edge with the remaining edge midpoints
   std::vector< uint_t > edgeMidpts{ _01, _03, _12, _23 };
   for ( uint_t i = 0; i < edgeMidpts.size(); ++i )
   {
      newFaces[i + 4] = fac.make_face( edgeMidpts[i], _02, _13 );
   }

   // === new cells ===

   // following Bey:

   cell->add_child( fac.make_cell( 0, _01, _02, _03 ) );
   cell->add_child( fac.make_cell( _01, 1, _12, _13 ) );
   cell->add_child( fac.make_cell( _02, _12, 2, _23 ) );
   cell->add_child( fac.make_cell( _03, _13, _23, 3 ) );

   cell->add_child( fac.make_cell( _01, _02, _03, _13 ) );
   cell->add_child( fac.make_cell( _01, _02, _12, _13 ) );
   cell->add_child( fac.make_cell( _02, _03, _13, _23 ) );
   cell->add_child( fac.make_cell( _02, _12, _13, _23 ) );

   // === child IDs ===

   auto childIDs = cell->getPrimitiveID().createChildren();
   newEdge->setPrimitiveID( childIDs[0] );
   for ( uint_t i = 0; i < newFaces.size(); ++i )
   {
      newFaces[i]->setPrimitiveID( childIDs[i + 1] );
   }
   auto& children = cell->get_children();
   for ( uint_t i = 0; i < children.size(); ++i )
   {
      children[i]->setPrimitiveID( childIDs[i + 1 + 8] );
   }

   return std::set< std::shared_ptr< Simplex3 > >( cell->get_children().begin(), cell->get_children().end() );
}

/* apply green refinement to cell with 1 new vertex
      @return sub-elements
   */
inline std::set< std::shared_ptr< Simplex3 > > refine_cell_green_1( std::shared_ptr< Simplex3 > cell )
{
   auto faces = cell->refined_faces();
   WALBERLA_ASSERT( faces[GREEN].size() == 2 && faces[RED].size() == 0 );

   // === split faces ===

   for ( auto& face : faces[GREEN] )
   {
      if ( !face->has_children() )
      {
         refine_face_green( face );
      }
   }

   // === collect vertex indices ===
   std::vector< uint_t > ref_vertices( 5 );

   // vertices at split edge
   std::shared_ptr< Simplex1 > split_edge;

   for ( auto& edge : cell->get_edges() )
   {
      if ( edge->has_children() )
      {
         ref_vertices[4] = edge->get_midpoint_idx();
         split_edge      = edge;
         break;
      }
   }
   WALBERLA_ASSERT( split_edge != nullptr );

   ref_vertices[1] = split_edge->get_vertices()[0];
   ref_vertices[2] = split_edge->get_vertices()[1];

   auto is_unbranded = [&]( uint_t vtx ) {
      return ( vtx != ref_vertices[1] && vtx != ref_vertices[2] && vtx != ref_vertices[4] );
   };
   std::vector< uint_t > unbranded;
   std::copy_if( cell->get_vertices().begin(), cell->get_vertices().end(), std::back_inserter( unbranded ), is_unbranded );

   ref_vertices[0] = unbranded[0];
   ref_vertices[3] = unbranded[1];

   SimplexFactory fac( cell, ref_vertices );

   // === collect faces ===

   auto face0        = cell->get_Face( ref_vertices[0], ref_vertices[1], ref_vertices[2] );
   auto child_faces0 = face0->get_children_sorted( { ref_vertices[1], ref_vertices[2], ref_vertices[0] } );

   auto face2        = cell->get_Face( ref_vertices[3], ref_vertices[1], ref_vertices[2] );
   auto child_faces2 = face2->get_children_sorted( { ref_vertices[1], ref_vertices[2], ref_vertices[3] } );

   fac.add_face( 0, 1, 3, cell->get_Face( ref_vertices[0], ref_vertices[1], ref_vertices[3] ) );
   fac.add_face( 0, 2, 3, cell->get_Face( ref_vertices[0], ref_vertices[2], ref_vertices[3] ) );
   fac.add_face( 0, 1, 4, child_faces0[0] );
   fac.add_face( 0, 2, 4, child_faces0[1] );
   fac.add_face( 1, 3, 4, child_faces2[0] );
   fac.add_face( 2, 3, 4, child_faces2[1] );

   // === collect edges ===

   fac.add_edge( 0, 1, cell->get_Edge( ref_vertices[0], ref_vertices[1] ) );
   fac.add_edge( 0, 2, cell->get_Edge( ref_vertices[0], ref_vertices[2] ) );
   fac.add_edge( 0, 3, cell->get_Edge( ref_vertices[0], ref_vertices[3] ) );
   fac.add_edge( 1, 3, cell->get_Edge( ref_vertices[1], ref_vertices[3] ) );
   fac.add_edge( 2, 3, cell->get_Edge( ref_vertices[2], ref_vertices[3] ) );
   fac.add_edge( 0, 4, child_faces0[0]->get_Edge( ref_vertices[0], ref_vertices[4] ) );
   fac.add_edge( 1, 4, child_faces0[0]->get_Edge( ref_vertices[1], ref_vertices[4] ) );
   fac.add_edge( 2, 4, child_faces0[1]->get_Edge( ref_vertices[2], ref_vertices[4] ) );
   fac.add_edge( 3, 4, child_faces2[0]->get_Edge( ref_vertices[3], ref_vertices[4] ) );

   // === new interior face ===

   auto newFace = fac.make_face( 0, 3, 4 );

   // === split cell ===

   cell->add_child( fac.make_cell( 0, 1, 4, 3 ) );
   cell->add_child( fac.make_cell( 0, 2, 4, 3 ) );

   // === child IDs ===

   auto childIDs = cell->getPrimitiveID().createChildren();
   newFace->setPrimitiveID( childIDs[1] );
   auto& children = cell->get_children();
   for ( uint_t i = 0; i < children.size(); ++i )
   {
      children[i]->setPrimitiveID( childIDs[i + 1 + 8 + 8] ); // we require different ids for different types of refinement
   }

   return std::set< std::shared_ptr< Simplex3 > >( cell->get_children().begin(), cell->get_children().end() );
}

/* apply green refinement to cell with 2 new vertices
      @return sub-elements
   */
inline std::set< std::shared_ptr< Simplex3 > > refine_cell_green_2( std::shared_ptr< Simplex3 > cell )
{
   auto faces = cell->refined_faces();
   WALBERLA_ASSERT( faces[GREEN].size() == 4 && faces[RED].size() == 0 );

   // === split faces ===

   for ( auto& face : faces[GREEN] )
   {
      if ( !face->has_children() )
      {
         refine_face_green( face );
      }
   }

   // === collect vertex indices ===
   std::vector< uint_t > ref_vertices( 6 );

   // vertices at split edges
   std::array< std::shared_ptr< Simplex1 >, 2 > split_edges;

   uint_t split_idx = 0;
   for ( auto& edge : cell->get_edges() )
   {
      if ( split_idx >= 2 )
         break;

      if ( edge->has_children() )
      {
         ref_vertices[4 + split_idx] = edge->get_midpoint_idx();
         split_edges[split_idx]      = edge;
         ++split_idx;
      }
   }
   WALBERLA_ASSERT( split_edges[0] && split_edges[1] );

   ref_vertices[1] = split_edges[0]->get_vertices()[0];
   ref_vertices[2] = split_edges[0]->get_vertices()[1];

   ref_vertices[0] = split_edges[1]->get_vertices()[0];
   ref_vertices[3] = split_edges[1]->get_vertices()[1];

   SimplexFactory fac( cell, ref_vertices );

   // === collect faces ===

   auto face0        = cell->get_Face( ref_vertices[0], ref_vertices[1], ref_vertices[2] );
   auto child_faces0 = face0->get_children_sorted( { ref_vertices[1], ref_vertices[2], ref_vertices[0] } );

   auto face1        = cell->get_Face( ref_vertices[0], ref_vertices[1], ref_vertices[3] );
   auto child_faces1 = face1->get_children_sorted( { ref_vertices[0], ref_vertices[3], ref_vertices[1] } );

   auto face2        = cell->get_Face( ref_vertices[3], ref_vertices[1], ref_vertices[2] );
   auto child_faces2 = face2->get_children_sorted( { ref_vertices[1], ref_vertices[2], ref_vertices[3] } );

   auto face3        = cell->get_Face( ref_vertices[0], ref_vertices[2], ref_vertices[3] );
   auto child_faces3 = face3->get_children_sorted( { ref_vertices[0], ref_vertices[3], ref_vertices[2] } );

   fac.add_face( 0, 1, 4, child_faces0[0] );
   fac.add_face( 0, 2, 4, child_faces0[1] );
   fac.add_face( 0, 1, 5, child_faces1[0] );
   fac.add_face( 1, 3, 5, child_faces1[1] );
   fac.add_face( 1, 3, 4, child_faces2[0] );
   fac.add_face( 2, 3, 4, child_faces2[1] );
   fac.add_face( 0, 2, 5, child_faces3[0] );
   fac.add_face( 2, 3, 5, child_faces3[1] );

   // === collect edges ===

   fac.add_edge( 0, 1, cell->get_Edge( ref_vertices[0], ref_vertices[1] ) );
   fac.add_edge( 0, 2, cell->get_Edge( ref_vertices[0], ref_vertices[2] ) );
   fac.add_edge( 1, 3, cell->get_Edge( ref_vertices[1], ref_vertices[3] ) );
   fac.add_edge( 2, 3, cell->get_Edge( ref_vertices[2], ref_vertices[3] ) );
   fac.add_edge( 0, 4, child_faces0[0]->get_Edge( ref_vertices[0], ref_vertices[4] ) );
   fac.add_edge( 1, 4, child_faces0[0]->get_Edge( ref_vertices[1], ref_vertices[4] ) );
   fac.add_edge( 2, 4, child_faces0[1]->get_Edge( ref_vertices[2], ref_vertices[4] ) );
   fac.add_edge( 3, 4, child_faces2[0]->get_Edge( ref_vertices[3], ref_vertices[4] ) );
   fac.add_edge( 0, 5, child_faces1[0]->get_Edge( ref_vertices[0], ref_vertices[5] ) );
   fac.add_edge( 1, 5, child_faces1[0]->get_Edge( ref_vertices[1], ref_vertices[5] ) );
   fac.add_edge( 2, 5, child_faces3[0]->get_Edge( ref_vertices[2], ref_vertices[5] ) );
   fac.add_edge( 3, 5, child_faces3[1]->get_Edge( ref_vertices[3], ref_vertices[5] ) );

   // new edge in cell interior
   auto newEdge = fac.make_edge( 4, 5, GREEN );

   // === new faces in cell interior

   std::array< std::shared_ptr< Simplex2 >, 4 > newFaces;
   for ( uint_t k = 0; k < 4; ++k )
   {
      newFaces[k] = fac.make_face( k, 4, 5 );
   }

   // === split cell ===

   cell->add_child( fac.make_cell( 0, 1, 4, 5 ) );
   cell->add_child( fac.make_cell( 0, 2, 4, 5 ) );
   cell->add_child( fac.make_cell( 1, 3, 4, 5 ) );
   cell->add_child( fac.make_cell( 2, 3, 4, 5 ) );

   // === child IDs ===

   auto childIDs = cell->getPrimitiveID().createChildren();
   newEdge->setPrimitiveID( childIDs[0] );
   for ( uint_t i = 0; i < newFaces.size(); ++i )
   {
      newFaces[i]->setPrimitiveID( childIDs[i + 1] );
   }
   auto& children = cell->get_children();
   for ( uint_t i = 0; i < children.size(); ++i )
   {
      children[i]->setPrimitiveID( childIDs[i + 1 + 8 + 8 + 2] ); // we require different ids for different types of refinement
   }

   return std::set< std::shared_ptr< Simplex3 > >( cell->get_children().begin(), cell->get_children().end() );
}

/* apply green refinement to cell with 3 new vertices
      @return sub-elements
   */
inline std::set< std::shared_ptr< Simplex3 > > refine_cell_green_3( std::shared_ptr< Simplex3 > cell )
{
   auto faces = cell->refined_faces();
   WALBERLA_ASSERT( faces[GREEN].size() == 3 && faces[RED].size() == 1 );

   // === split faces ===

   for ( auto& face : faces[GREEN] )
   {
      if ( !face->has_children() )
      {
         refine_face_green( face );
      }
   }

   // === collect vertex indices ===
   std::vector< uint_t > ref_vertices( 7 );

   // vertices at red face
   auto red_vtxs  = faces[RED][0]->get_vertices();
   auto red_edges = faces[RED][0]->get_edges_sorted( red_vtxs );

   for ( uint_t k = 0; k < 3; ++k )
   {
      ref_vertices[k]     = red_vtxs[k];
      ref_vertices[k + 4] = red_edges[k]->get_midpoint_idx();
   }

   // top vertex
   auto old_vtxs = cell->get_vertices();

   for ( uint_t k = 0; k < 4; ++k )
   {
      bool found = true;

      for ( uint_t bottom = 0; bottom < 3; ++bottom )
      {
         if ( old_vtxs[k] == ref_vertices[bottom] )
         {
            found = false;
         }
      }

      if ( found )
      {
         ref_vertices[3] = old_vtxs[k];
         break;
      }
   }

   SimplexFactory fac( cell, ref_vertices );

   // === collect faces ===

   auto face0        = cell->get_Face( ref_vertices[0], ref_vertices[1], ref_vertices[2] );
   auto child_faces0 = face0->get_children_sorted( { ref_vertices[0], ref_vertices[1], ref_vertices[2] } );

   auto face1        = cell->get_Face( ref_vertices[0], ref_vertices[1], ref_vertices[3] );
   auto child_faces1 = face1->get_children_sorted( { ref_vertices[0], ref_vertices[1], ref_vertices[3] } );

   auto face2        = cell->get_Face( ref_vertices[1], ref_vertices[2], ref_vertices[3] );
   auto child_faces2 = face2->get_children_sorted( { ref_vertices[1], ref_vertices[2], ref_vertices[3] } );

   auto face3        = cell->get_Face( ref_vertices[2], ref_vertices[0], ref_vertices[3] );
   auto child_faces3 = face3->get_children_sorted( { ref_vertices[2], ref_vertices[0], ref_vertices[3] } );

   fac.add_face( 0, 4, 6, child_faces0[0] );
   fac.add_face( 1, 4, 5, child_faces0[1] );
   fac.add_face( 2, 5, 6, child_faces0[2] );
   fac.add_face( 4, 5, 6, child_faces0[3] );

   fac.add_face( 0, 3, 4, child_faces1[0] );
   fac.add_face( 1, 3, 4, child_faces1[1] );

   fac.add_face( 1, 3, 5, child_faces2[0] );
   fac.add_face( 2, 3, 5, child_faces2[1] );

   fac.add_face( 2, 3, 6, child_faces3[0] );
   fac.add_face( 0, 3, 6, child_faces3[1] );

   // === collect edges ===

   // old edges
   fac.add_edge( 0, 3, cell->get_Edge( ref_vertices[0], ref_vertices[3] ) );
   fac.add_edge( 1, 3, cell->get_Edge( ref_vertices[1], ref_vertices[3] ) );
   fac.add_edge( 2, 3, cell->get_Edge( ref_vertices[2], ref_vertices[3] ) );
   // red edges
   fac.add_edge( 0, 4, child_faces0[0]->get_Edge( ref_vertices[0], ref_vertices[4] ) );
   fac.add_edge( 0, 6, child_faces0[0]->get_Edge( ref_vertices[0], ref_vertices[6] ) );
   fac.add_edge( 1, 4, child_faces0[1]->get_Edge( ref_vertices[1], ref_vertices[4] ) );
   fac.add_edge( 1, 5, child_faces0[1]->get_Edge( ref_vertices[1], ref_vertices[5] ) );
   fac.add_edge( 2, 5, child_faces0[2]->get_Edge( ref_vertices[2], ref_vertices[5] ) );
   fac.add_edge( 2, 6, child_faces0[2]->get_Edge( ref_vertices[2], ref_vertices[6] ) );
   fac.add_edge( 4, 5, child_faces0[3]->get_Edge( ref_vertices[4], ref_vertices[5] ) );
   fac.add_edge( 5, 6, child_faces0[3]->get_Edge( ref_vertices[5], ref_vertices[6] ) );
   fac.add_edge( 4, 6, child_faces0[3]->get_Edge( ref_vertices[6], ref_vertices[4] ) );
   // green edges
   fac.add_edge( 3, 4, child_faces1[0]->get_Edge( ref_vertices[3], ref_vertices[4] ) );
   fac.add_edge( 3, 5, child_faces2[0]->get_Edge( ref_vertices[3], ref_vertices[5] ) );
   fac.add_edge( 3, 6, child_faces3[0]->get_Edge( ref_vertices[3], ref_vertices[6] ) );

   // === new faces in cell interior

   std::array< std::shared_ptr< Simplex2 >, 3 > newFaces;

   newFaces[0] = fac.make_face( 3, 4, 5 );
   newFaces[1] = fac.make_face( 3, 5, 6 );
   newFaces[2] = fac.make_face( 3, 4, 6 );

   // === split cell ===

   cell->add_child( fac.make_cell( 0, 4, 6, 3 ) );
   cell->add_child( fac.make_cell( 1, 4, 5, 3 ) );
   cell->add_child( fac.make_cell( 2, 5, 6, 3 ) );
   cell->add_child( fac.make_cell( 4, 5, 6, 3 ) );

   // === child IDs ===

   auto childIDs = cell->getPrimitiveID().createChildren();
   for ( uint_t i = 0; i < newFaces.size(); ++i )
   {
      newFaces[i]->setPrimitiveID( childIDs[i + 1] );
   }
   auto& children = cell->get_children();
   for ( uint_t i = 0; i < children.size(); ++i )
   {
      children[i]->setPrimitiveID(
          childIDs[i + 1 + 8 + 8 + 2 + 4] ); // we require different ids for different types of refinement
   }

   return std::set< std::shared_ptr< Simplex3 > >( cell->get_children().begin(), cell->get_children().end() );
}

} // namespace adaptiveRefinement
} // namespace hyteg
