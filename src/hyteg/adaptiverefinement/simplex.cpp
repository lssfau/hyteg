#include "simplex.hpp"

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
namespace hyteg {
namespace adaptiveRefinement {

// ============= K-Simplex =================
template < int K, class K_Simplex >
const std::vector< std::shared_ptr< K_Simplex > >
    Simplex< K, K_Simplex >::get_children_sorted( const std::array< int, K + 1 >& vertices ) const
{
   assert( K < 3 ); // only implemented for 1 and 2 dimensions

   assert( has_children() );

   std::vector< std::shared_ptr< K_Simplex > > sorted( K + 1 );

   for ( int i = 0; i <= K; ++i )
   {
      for ( auto& child : _children )
      {
         if ( child->has_vertex( vertices[i] ) )
         {
            sorted[i] = child;
            break;
         }
      }

      assert( sorted[i] );
   }

   for ( int i = K + 1; i < int( _children.size() ); ++i )
   {
      sorted.push_back( _children[i] );
   }

   return sorted;
}

template < int K, class K_Simplex >
void Simplex< K, K_Simplex >::add_child( const std::shared_ptr< K_Simplex >& child )
{
   _children.push_back( child );
}

template < int K, class K_Simplex >
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

template < int K, class K_Simplex >
bool Simplex< K, K_Simplex >::has_vertex( int idx ) const
{
   for ( auto& vtx : _vertices )
   {
      if ( vtx == idx )
         return true;
   }
   return false;
}

template < int K, class K_Simplex >
Point3D Simplex< K, K_Simplex >::barycenter( const std::vector< Point3D >& nodes ) const
{
   Point3D bc({ 0, 0, 0 });
   for ( auto& vtxID : _vertices )
   {
      bc += nodes[vtxID];
   }
   bc *= ( 1.0 / ( K + 1 ) );
   return bc;
}

template class Simplex< 1, Simplex1 >;
template class Simplex< 2, Simplex2 >;
template class Simplex< 3, Simplex3 >;

// ============= 1-Simplex (edge) =================

int Simplex1::inner_vertices() const
{
   if ( has_children() )
   {
      int n = 1;
      for ( auto& child : _children )
      {
         n += child->inner_vertices();
      }
      return n;
   }

   return 0;
}

// ============= 2-Simplex (face) =================

Simplex2::Simplex2( const std::array< int, 3 >&                         vertices,
                    const std::array< std::shared_ptr< Simplex1 >, 3 >& edges,
                    std::shared_ptr< Simplex2 >                         parent )
: Simplex< 2, Simplex2 >( vertices, parent )
, _edges( edges )
{
   for ( int i = 0; i < 3; ++i )
   {
      assert( _edges[i]->has_vertex( _vertices[i] ) );
      assert( _edges[i]->has_vertex( _vertices[( i + 1 ) % 3] ) );
   }
}

int Simplex2::vertices_on_edges() const
{
   int n = 0;

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

std::array< std::shared_ptr< Simplex1 >, 3 > Simplex2::get_edges_sorted( const std::array< int, 3 >& vertices ) const
{
   std::array< std::shared_ptr< Simplex1 >, 3 > sorted;

   for ( int i = 0; i < 3; ++i )
   {
      for ( auto& edge : _edges )
      {
         if ( edge->has_vertex( vertices[i] ) and edge->has_vertex( vertices[( i + 1 ) % 3] ) )
         {
            sorted[i] = edge;
         }
      }

      assert( sorted[i] );
   }

   return sorted;
}

template <>
double Simplex< 2, Simplex2 >::volume( const std::vector< Point3D >& nodes ) const
{
   auto& a = nodes[_vertices[0]];
   auto& b = nodes[_vertices[1]];
   auto& c = nodes[_vertices[2]];

   auto x = c - a;
   auto y = b - a;

   // ||X x Y|| / 2
   return crossProduct(x, y).norm() / 2;
}

std::shared_ptr< Simplex1 > Simplex2::get_Edge( int a, int b ) const
{
   for ( auto& edge : _edges )
   {
      if ( edge->has_vertex( a ) and edge->has_vertex( b ) )
      {
         return edge;
      }
   }

   return nullptr;
}

std::pair< real_t, real_t > Simplex2::min_max_angle( const std::vector< Point3D >& nodes ) const
{
   std::pair< real_t, real_t > mm{ 10, 0 };

   for ( int i = 0; i < 3; ++i )
   {
      int j = ( i + 1 ) % 3;
      int k = ( i + 2 ) % 3;

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

Simplex3::Simplex3( const std::array< int, 4 >&                         vertices,
                    const std::array< std::shared_ptr< Simplex1 >, 6 >& edges,
                    const std::array< std::shared_ptr< Simplex2 >, 4 >& faces,
                    std::shared_ptr< Simplex3 >                         parent )
: Simplex< 3, Simplex3 >( vertices, parent )
, _edges( edges )
, _faces( faces )
{
   // check edges
   for ( int i = 0; i < 3; ++i )
   {
      // edge 0--2
      assert( _edges[i]->has_vertex( _vertices[i] ) );
      assert( _edges[i]->has_vertex( _vertices[( i + 1 ) % 3] ) );
      // edge 3--5
      assert( _edges[i + 3]->has_vertex( _vertices[i] ) );
      assert( _edges[i + 3]->has_vertex( _vertices[3] ) );
   }

   // check faces
   // face 0
   for ( int j = 0; j < 3; ++j )
   {
      assert( _faces[0]->has_vertex( _vertices[j] ) );
   }
   // face 1--3
   for ( int i = 1; i < 4; ++i )
   {
      assert( _faces[i]->has_vertex( _vertices[i - 1] ) );
      assert( _faces[i]->has_vertex( _vertices[i % 3] ) );
      assert( _faces[i]->has_vertex( _vertices[3] ) );
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

int Simplex3::vertices_on_edges() const
{
   int n = 0;

   for ( auto& edge : _edges )
   {
      n += edge->inner_vertices();
   }

   return n;
}

template <>
double Simplex< 3, Simplex3 >::volume( const std::vector< Point3D >& nodes ) const
{
   auto& a = nodes[_vertices[0]];
   auto& b = nodes[_vertices[1]];
   auto& c = nodes[_vertices[2]];
   auto& d = nodes[_vertices[3]];

   auto x = a - d;
   auto y = b - d;
   auto z = c - d;

   // |X * (Y x Z)| / 6
   return std::abs( x.dot( crossProduct(y, z) ) ) / 6;
}

std::pair< real_t, real_t > Simplex3::min_max_angle( const std::vector< Point3D >& nodes ) const
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

std::shared_ptr< Simplex1 > Simplex3::get_Edge( int a, int b ) const
{
   for ( auto& edge : _edges )
   {
      if ( edge->has_vertex( a ) and edge->has_vertex( b ) )
      {
         return edge;
      }
   }

   return nullptr;
}

std::shared_ptr< Simplex2 > Simplex3::get_Face( int a, int b, int c ) const
{
   for ( auto& face : _faces )
   {
      if ( face->has_vertex( a ) and face->has_vertex( b ) and face->has_vertex( c ) )
      {
         return face;
      }
   }

   return nullptr;
}

} // namespace adaptiveRefinement
} // namespace hyteg