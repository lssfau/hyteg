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

#include "simplex.hpp"

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>

using walberla::real_c;

namespace hyteg {
namespace adaptiveRefinement {

// ============= K-Simplex =================
template < uint_t K, class K_Simplex >
const std::vector< std::shared_ptr< K_Simplex > >
    Simplex< K, K_Simplex >::get_children_sorted( const std::array< uint_t, K + 1 >& vertices ) const
{
   WALBERLA_ASSERT( K < 3 ); // only implemented for 1 and 2 dimensions

   WALBERLA_ASSERT( has_children() );

   std::vector< std::shared_ptr< K_Simplex > > sorted( K + 1 );

   for ( uint_t i = 0; i <= K; ++i )
   {
      for ( auto& child : _children )
      {
         WALBERLA_ASSERT( child != nullptr );
         if ( child->has_vertex( vertices[uint_t( i )] ) )
         {
            sorted[i] = child;
            break;
         }
      }

      WALBERLA_ASSERT( sorted[i] != nullptr );
   }

   for ( uint_t i = K + 1; i < _children.size(); ++i )
   {
      sorted.push_back( _children[i] );
   }

   return sorted;
}

template < uint_t K, class K_Simplex >
void Simplex< K, K_Simplex >::add_child( const std::shared_ptr< K_Simplex >& child )
{
   _children.push_back( child );
}

template < uint_t K, class K_Simplex >
bool Simplex< K, K_Simplex >::kill_children()
{
   if ( _children.empty() )
   {
      return false;
   }
   else
   {
      _children.clear();
      return true;
   }
}

template < uint_t K, class K_Simplex >
bool Simplex< K, K_Simplex >::has_grandkids() const
{
   for ( auto& child : _children )
   {
      if ( child->has_children() )
      {
         return true;
      }
   }
   return false;
}

template < uint_t K, class K_Simplex >
bool Simplex< K, K_Simplex >::has_vertex( uint_t idx ) const
{
   for ( auto& vtx : _vertices )
   {
      if ( vtx == idx )
         return true;
   }
   return false;
}

template < uint_t K, class K_Simplex >
inline Point3D Simplex< K, K_Simplex >::barycenter( const std::array< Point3D, K + 1 >& vertices )
{
   Point3D bc( 0, 0, 0 );
   for ( auto& vtx : vertices )
   {
      bc += vtx;
   }
   bc *= ( real_c( 1.0 ) / ( K + 1 ) );
   return bc;
}

template < uint_t K, class K_Simplex >
inline real_t Simplex< K, K_Simplex >::radius( const std::array< Point3D, K + 1 >& vertices )
{
   Point3D bc = barycenter( vertices );
   real_t  r  = real_c( 0.0 );
   for ( auto& vtx : vertices )
   {
      r = std::max( r, ( vtx - bc ).norm() );
   }
   return r;
}

template < uint_t K, class K_Simplex >
inline real_t Simplex< K, K_Simplex >::volume( const std::array< Point3D, K + 1 >& vertices )
{
   auto x = vertices[1] - vertices[0];

   if constexpr ( K == 1 )
   {
      // length of edge
      return x.norm();
   }
   else
   {
      auto y = vertices[2] - vertices[0];

      if constexpr ( K == 2 )
      {
         // ||X x Y|| / 2
         return crossProduct( x, y ).norm() / real_c( 2.0 );
      }
      if constexpr ( K == 3 )
      {
         auto z = vertices[3] - vertices[0];
         // |(X x Y)/2 * Z/3|
         return std::abs( crossProduct( x, y ).dot( z ) ) / real_c( 6.0 );
      }
   }

   return real_c( 0.0 );
}

template < uint_t K, class K_Simplex >
inline std::array< Point3D, K + 1 > Simplex< K, K_Simplex >::coordinates( const EnumeratedList< Point3D >& nodes ) const
{
   std::array< Point3D, K + 1 > vertices;
   for ( uint_t i = 0; i < K + 1; ++i )
   {
      vertices[i] = nodes[_vertices[i]];
   }
   return vertices;
}

template < uint_t K, class K_Simplex >
Point3D Simplex< K, K_Simplex >::barycenter( const EnumeratedList< Point3D >& nodes ) const
{
   return barycenter( this->coordinates( nodes ) );
}

template < uint_t K, class K_Simplex >
real_t Simplex< K, K_Simplex >::radius( const EnumeratedList< Point3D >& nodes ) const
{
   return radius( this->coordinates( nodes ) );
}

template < uint_t K, class K_Simplex >
real_t Simplex< K, K_Simplex >::volume( const EnumeratedList< Point3D >& nodes ) const
{
   return volume( this->coordinates( nodes ) );
}

template class Simplex< 1, Simplex1 >;
template class Simplex< 2, Simplex2 >;
template class Simplex< 3, Simplex3 >;

// ============= 1-Simplex (edge) =================

uint_t Simplex1::inner_vertices() const
{
   if ( has_children() )
   {
      uint_t n = 1;
      for ( auto& child : _children )
      {
         n += child->inner_vertices();
      }
      return n;
   }

   return 0;
}

// ============= 2-Simplex (face) =================

Simplex2::Simplex2( const std::array< uint_t, 3 >&                      vertices,
                    const std::array< std::shared_ptr< Simplex1 >, 3 >& edges,
                    std::shared_ptr< Simplex2 >                         parent,
                    const PrimitiveID&                                  geometryMap,
                    uint_t                                              boundaryFlag,
                    uint_t                                              targetRank )
: Simplex< 2, Simplex2 >( vertices, parent, geometryMap, boundaryFlag, targetRank )
, _edges( edges )
{
   for ( uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_ASSERT( _edges[i] != nullptr );

      WALBERLA_ASSERT( _edges[i]->has_vertex( _vertices[i] ) );
      WALBERLA_ASSERT( _edges[i]->has_vertex( _vertices[( i + 1 ) % 3] ) );
   }
}

uint_t Simplex2::vertices_on_edges() const
{
   uint_t n = 0;

   for ( auto& edge : _edges )
   {
      n += edge->inner_vertices();
   }

   return n;
}

bool Simplex2::has_green_edge() const
{
   for ( auto& edge : _edges )
   {
      if ( edge->is_green() )
      {
         return true;
      }
   }
   return false;
}

std::array< std::shared_ptr< Simplex1 >, 3 > Simplex2::get_edges_sorted( const std::array< uint_t, 3 >& vertices ) const
{
   std::array< std::shared_ptr< Simplex1 >, 3 > sorted;

   for ( uint_t i = 0; i < 3; ++i )
   {
      for ( auto& edge : _edges )
      {
         WALBERLA_ASSERT( edge != nullptr );
         if ( edge->has_vertex( vertices[i] ) && edge->has_vertex( vertices[( i + 1 ) % 3] ) )
         {
            sorted[i] = edge;
         }
      }

      WALBERLA_ASSERT( sorted[i] != nullptr );
   }

   return sorted;
}

std::shared_ptr< Simplex1 > Simplex2::get_Edge( uint_t a, uint_t b ) const
{
   for ( auto& edge : _edges )
   {
      WALBERLA_ASSERT( edge != nullptr );
      if ( edge->has_vertex( a ) && edge->has_vertex( b ) )
      {
         return edge;
      }
   }

   return nullptr;
}

std::pair< real_t, real_t > Simplex2::min_max_angle( const EnumeratedList< Point3D >& nodes ) const
{
   std::pair< real_t, real_t > mm{ 10, 0 };

   for ( uint_t i = 0; i < 3; ++i )
   {
      uint_t j = ( i + 1 ) % 3;
      uint_t k = ( i + 2 ) % 3;

      auto vi = nodes[_vertices[i]];
      auto vj = nodes[_vertices[j]];
      auto vk = nodes[_vertices[k]];

      auto v1      = vj - vi;
      auto v2      = vk - vi;
      auto angle_i = std::acos( v1.dot( v2 ) / ( v1.norm() * v2.norm() ) );

      mm.first  = std::min( mm.first, angle_i );
      mm.second = std::max( mm.second, angle_i );
   }

   return mm;
}

// ============= 3-Simplex (cell) =================

Simplex3::Simplex3( const std::array< uint_t, 4 >&                      vertices,
                    const std::array< std::shared_ptr< Simplex1 >, 6 >& edges,
                    const std::array< std::shared_ptr< Simplex2 >, 4 >& faces,
                    std::shared_ptr< Simplex3 >                         parent,
                    const PrimitiveID&                                  geometryMap,
                    uint_t                                              boundaryFlag,
                    uint_t                                              targetRank )
: Simplex< 3, Simplex3 >( vertices, parent, geometryMap, boundaryFlag, targetRank )
, _edges( edges )
, _faces( faces )
{
   // check ordering of edges and faces
   // should be: E = [01 02 03 12 13 23]
   //            F = [012 013 023 123]
   uint_t e = 0, f = 0;
   // e and f are unused in Release build
   WALBERLA_UNUSED( e );
   WALBERLA_UNUSED( f );
   for ( uint_t i = 0; i < 4; ++i )
   {
      for ( uint_t j = i + 1; j < 4; ++j )
      {
         WALBERLA_ASSERT( _edges[e] != nullptr );
         WALBERLA_ASSERT( _edges[e]->has_vertex( _vertices[i] ) );
         WALBERLA_ASSERT( _edges[e]->has_vertex( _vertices[j] ) );
         ++e;

         for ( uint_t k = j + 1; k < 4; ++k )
         {
            WALBERLA_ASSERT( _faces[f] != nullptr );
            WALBERLA_ASSERT( _faces[f]->has_vertex( _vertices[i] ) );
            WALBERLA_ASSERT( _faces[f]->has_vertex( _vertices[j] ) );
            WALBERLA_ASSERT( _faces[f]->has_vertex( _vertices[k] ) );
            ++f;
         }
      }
   }
}

bool Simplex3::has_green_edge() const
{
   for ( auto& edge : _edges )
   {
      if ( edge->is_green() )
      {
         return true;
      }
   }
   return false;
}

std::map< Color, std::vector< std::shared_ptr< Simplex2 > > > Simplex3::refined_faces() const
{
   std::map< Color, std::vector< std::shared_ptr< Simplex2 > > > map;
   map[GREEN] = std::vector< std::shared_ptr< Simplex2 > >();
   map[RED]   = std::vector< std::shared_ptr< Simplex2 > >();

   for ( auto& face : _faces )
   {
      if ( face->vertices_on_edges() == 1 )
      {
         map[GREEN].push_back( face );
      }
      else if ( face->vertices_on_edges() > 1 )
      {
         map[RED].push_back( face );
      }
   }

   return map;
}

uint_t Simplex3::vertices_on_edges() const
{
   uint_t n = 0;

   for ( auto& edge : _edges )
   {
      n += edge->inner_vertices();
   }

   return n;
}

std::pair< real_t, real_t > Simplex3::min_max_angle( const EnumeratedList< Point3D >& nodes ) const
{
   std::pair< real_t, real_t > mm{ 10, 0 };

   for ( auto& el : _faces )
   {
      auto mm_el = el->min_max_angle( nodes );

      mm.first  = std::min( mm.first, mm_el.first );
      mm.second = std::max( mm.second, mm_el.second );
   }

   return mm;
}

std::shared_ptr< Simplex1 > Simplex3::get_Edge( uint_t a, uint_t b ) const
{
   for ( auto& edge : _edges )
   {
      WALBERLA_ASSERT( edge != nullptr );
      if ( edge->has_vertex( a ) && edge->has_vertex( b ) )
      {
         return edge;
      }
   }

   return nullptr;
}

std::shared_ptr< Simplex2 > Simplex3::get_Face( uint_t a, uint_t b, uint_t c ) const
{
   for ( auto& face : _faces )
   {
      WALBERLA_ASSERT( face != nullptr );
      if ( face->has_vertex( a ) && face->has_vertex( b ) && face->has_vertex( c ) )
      {
         return face;
      }
   }

   return nullptr;
}

} // namespace adaptiveRefinement
} // namespace hyteg