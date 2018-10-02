#pragma once

#include "core/DataTypes.h"

namespace hhg
{
namespace levelinfo
{

using walberla::uint_t;
using walberla::uint_c;

constexpr inline uint_t num_microvertices_per_vertex(uint_t /*level*/)
{
  return 1;
}

constexpr inline uint_t num_microvertices_per_edge(uint_t level)
{
  return ( (1u << level) + 1u);
  //return (uint_t) std::pow(2, level) + 1;
}

constexpr inline uint_t num_microvertices_per_edge_from_width( uint_t width )
{
  return width;
}

constexpr inline uint_t num_microedges_per_edge(uint_t level)
{
  return num_microvertices_per_edge(level) - 1;
}

constexpr inline uint_t num_microedges_per_edge_from_width( const uint_t & width)
{
  return num_microvertices_per_edge_from_width( width ) - 1;
}

constexpr inline uint_t num_microedges_per_face_from_width( const uint_t & width)
{
  return 3 * ( ( ( num_microedges_per_edge_from_width( width ) + 1 ) * num_microedges_per_edge_from_width( width ) ) / 2 );
}

constexpr inline uint_t num_microedges_per_face(uint_t level)
{
  return 3 * ( ( ( num_microedges_per_edge( level ) + 1 ) * num_microedges_per_edge( level ) ) / 2 );
}

constexpr inline uint_t num_microvertices_per_face_from_width( const uint_t & width )
{
  return ((width * (width + 1u)) >> 1u);
}

constexpr inline uint_t num_microvertices_per_face(uint_t level)
{
  const uint_t width = ( 1u << level ) + 1u;
  return num_microvertices_per_face_from_width( width );
}

constexpr inline uint_t num_microfaces_per_face(uint_t level)
{
  //pow(4, level)
  return (1u << (2u*level));
}

constexpr inline uint_t num_microvertices_per_cell_from_width( const uint_t & width )
{
  return ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
}

constexpr inline uint_t num_microvertices_per_cell( const uint_t & level )
{
  const uint_t width = ( 1u << level ) + 1;
  return ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
}

constexpr inline uint_t num_microcells_per_cell( const uint_t & level )
{
  // num_microcells_per_cell = 8 ^ level
  return (1 << (3 * level));
}

constexpr inline uint_t num_microcells_per_cell_from_width( const uint_t & width )
{
  const uint_t whiteUp   = num_microvertices_per_cell_from_width( width - 1 );
  const uint_t whiteDown = num_microvertices_per_cell_from_width( width - 3 );
  const uint_t others    = num_microvertices_per_cell_from_width( width - 2 );
  return whiteUp + whiteDown + 4 * others;
}


}

}

