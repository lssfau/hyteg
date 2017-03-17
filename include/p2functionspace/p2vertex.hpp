#ifndef P2VERTEX_HPP
#define P2VERTEX_HPP

#include <levelinfo.hpp>
#include <comm.hpp>

#include <fmt/format.h>

namespace hhg
{

/// P2Vertex namespace for P2 macro-vertex kernels
/// 
/// [Vertex dof, Edges: [Edge dof, Vertex ghost dof], Faces: [Edge ghost dof]]
namespace P2Vertex
{

inline void allocate(Vertex& vertex, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  vertex.data.push_back(std::vector<double*>());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    size_t num_deps = 3 * vertex.edges.size();
    size_t total_n_dofs = levelinfo::num_microvertices_per_vertex(level) + num_deps;
    double* new_data = new double[total_n_dofs];
    memset(new_data, 0, total_n_dofs * sizeof(double));
    vertex.data[memory_id].push_back(new_data);
  }
}

inline void free(Vertex& vertex, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    delete[] vertex.data[memory_id][level - minLevel];
  }
}

inline void interpolate(Vertex& vertex, size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level)
{
  vertex.data[memory_id][level-2][0] = expr(vertex.coords);

  size_t offset = 1;
  for (Edge* edge : vertex.edges)
  {
    Point3D x(vertex.coords);
    double orientation = 1.0;

    if (edge->vertex_index(vertex) == 1)
    {
      orientation = -1.0;
    }

    Point3D direction(orientation * edge->direction);
    direction /= levelinfo::num_microedges_per_edge(level);

    vertex.data[memory_id][level-2][offset] = expr(x + 0.5 * direction);
    ++offset;

    // vertex.data[memory_id][level-2][offset] = expr(x + direction);
    ++offset;
  }
}

}
}

#endif /* P2VERTEX_HPP */