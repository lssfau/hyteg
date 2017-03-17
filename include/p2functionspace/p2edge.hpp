#ifndef P2EDGE_HPP
#define P2EDGE_HPP

#include <levelinfo.hpp>
#include <comm.hpp>

namespace hhg
{
namespace P2Edge
{

inline void allocate(Edge& edge, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  edge.data.push_back(std::vector<double*>());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    size_t n_dofs_per_edge = 2 * levelinfo::num_microvertices_per_edge(level) - 1;
    size_t n_dofs_per_edge_nbr_1 = n_dofs_per_edge - 1;
    size_t n_dofs_per_edge_nbr_2 = n_dofs_per_edge - 2;
    size_t num_deps = edge.faces.size();
    size_t total_n_dofs = n_dofs_per_edge + num_deps * (n_dofs_per_edge_nbr_1 + n_dofs_per_edge_nbr_2);
    double* new_data = new double[total_n_dofs];
    memset(new_data, 0, total_n_dofs * sizeof(double));
    edge.data[memory_id].push_back(new_data);
  }
}

inline void free(Edge& edge, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    delete[] edge.data[memory_id][level - minLevel];
  }
}

inline void interpolate(Edge& edge, size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);
  Point3D x = edge.v0->coords;
  Point3D dx = edge.direction / (double) (rowsize - 1);
  x += dx;

  for (size_t i = 2; i < rowsize-2; ++i)
  {
    edge.data[memory_id][level-2][i] = expr(x);
    x += dx;
  }

  Point3D v(edge.faces[0]->get_vertex_opposite_to_edge(edge)->coords);
  Point3D dirv((v - edge.v0->coords) / ((double) rowsize - 1));

  x = edge.v0->coords + dirv;

  size_t offset = rowsize;

  for (size_t i = 1; i < rowsize-1-1; ++i)
  {
    edge.data[memory_id][level-2][offset+i] = expr(x);
    x += dx;
  }

  offset += rowsize-1 + rowsize-2;

  if (edge.faces.size() == 2)
  {
    v = edge.faces[1]->get_vertex_opposite_to_edge(edge)->coords;
    dirv = (v - edge.v0->coords) / ((double) rowsize - 1);

    for (size_t i = 1; i < rowsize-1-1; ++i)
    {
      edge.data[memory_id][level-2][offset+i] = expr(x);
      x += dx;
    }
  }

}

}
}

#endif /* P2EDGE_HPP */