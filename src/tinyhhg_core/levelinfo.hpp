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

constexpr inline uint_t num_microedges_per_edge(uint_t level)
{
  return num_microvertices_per_edge(level) - 1;
}

constexpr inline uint_t num_microedges_per_face(uint_t level)
{
  return 3 * ( ( ( num_microedges_per_edge( level ) + 1 ) * num_microedges_per_edge( level ) ) / 2 );
}

constexpr inline uint_t num_microvertices_per_face(uint_t level)
{
  //(pow(2,level) + 1) * (pow(2,level) + 2) / 2
  return ((((1u << level) + 1u) * ((1u << level) + 2u)) >> 1u);
}

constexpr inline uint_t num_microfaces_per_face(uint_t level)
{
  //pow(4, level)
  return (1u << (2u*level));
}

constexpr inline uint_t num_microvertices_per_cell( const uint_t & level )
{
  const uint_t width = ( 1u << level ) + 1;
  return ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
}

}
}

