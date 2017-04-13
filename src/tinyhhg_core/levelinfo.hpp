#ifndef LEVELINFO_HPP
#define LEVELINFO_HPP

#include <cmath>

namespace hhg
{
namespace levelinfo
{

inline size_t num_microvertices_per_vertex(size_t level)
{
  return 1;
}

inline size_t num_microvertices_per_edge(size_t level)
{
  return std::pow(2, level) + 1;
}

inline size_t num_microvertices_per_face(size_t level)
{
  return (std::pow(2, level)+1) * (std::pow(2, level-1) + 1);
}

inline size_t num_microfaces_per_face(size_t level)
{
  return std::pow(4, level);
}

}
}

#endif /* LEVELINFO_HPP */