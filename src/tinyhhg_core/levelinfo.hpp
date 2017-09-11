#ifndef LEVELINFO_HPP
#define LEVELINFO_HPP

#include <cmath>
#include "core/DataTypes.h"

namespace hhg
{
namespace levelinfo
{

using walberla::uint_t;

constexpr inline uint_t num_microvertices_per_vertex(uint_t /*level*/)
{
  return 1;
}

constexpr inline uint_t num_microvertices_per_edge(uint_t level)
{
  return (uint_t) std::pow(2, level) + 1;
}

constexpr inline uint_t num_microedges_per_edge(uint_t level)
{
  return num_microvertices_per_edge(level) - 1;
}

constexpr inline uint_t num_microvertices_per_face(uint_t level)
{
  return (uint_t) ((std::pow(2, level)+1) * (std::pow(2, level-1) + 1));
}

constexpr inline uint_t num_microfaces_per_face(uint_t level)
{
  return (uint_t) (std::pow(4, level));
}

}
}

#endif /* LEVELINFO_HPP */
