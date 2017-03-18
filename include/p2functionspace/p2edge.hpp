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

inline void interpolate(Edge& edge, size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level, bool boundaryDofs)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);
  Point3D x = edge.v0->coords;
  Point3D dx = edge.direction / (double) (rowsize - 1);
  x += dx;

  if (boundaryDofs)
  {
    for (size_t i = 1; i < rowsize-1; ++i)
    {
      edge.data[memory_id][level-2][i] = expr(x);
      x += dx;
    }
  }

  Point3D v(edge.faces[0]->get_vertex_opposite_to_edge(edge)->coords);
  Point3D dirv((v - edge.v0->coords) / ((double) rowsize - 1));

  x = edge.v0->coords + dx + dirv;

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

    x = edge.v0->coords + dx + dirv;

    for (size_t i = 1; i < rowsize-1-1; ++i)
    {
      edge.data[memory_id][level-2][offset+i] = expr(x);
      x += dx;
    }
  }

}

inline void pull_vertices(Edge& edge, size_t memory_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);

  if (edge.v0->rank == hhg::Comm::get().rk)
  {
    // local information
    if (edge.rank == hhg::Comm::get().rk)
    {
      edge.data[memory_id][level-2][0] = edge.v0->data[memory_id][level-2][0];
    }
    else
    {
      MPI_Send(&edge.v0->data[memory_id][level-2][0], 1, MPI_DOUBLE, edge.rank, 0, MPI_COMM_WORLD);
    }
  }
  else if (edge.rank == hhg::Comm::get().rk)
  {
    MPI_Recv(&edge.data[memory_id][level-2][0], 1, MPI_DOUBLE, edge.v0->rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (edge.v1->rank == hhg::Comm::get().rk)
  {
    // local information
    if (edge.rank == hhg::Comm::get().rk)
    {
      edge.data[memory_id][level-2][rowsize-1] = edge.v1->data[memory_id][level-2][0];
    }
    else
    {
      MPI_Send(&edge.v1->data[memory_id][level-2][0], 1, MPI_DOUBLE, edge.rank, 0, MPI_COMM_WORLD);
    }
  }
  else if (edge.rank == hhg::Comm::get().rk)
  {
    MPI_Recv(&edge.data[memory_id][level-2][rowsize-1], 1, MPI_DOUBLE, edge.v1->rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

inline void print(Edge & edge, size_t memory_id, size_t level) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);
  int midpos = 0;
  int face1pos = rowsize;
  int face1ghostpos = face1pos + rowsize - 1;
  int face2pos = face1ghostpos + rowsize - 2;
  int face2ghostpos = face2pos + rowsize - 1;
  auto mid = edge.data[memory_id][level - 2];
  auto face1 = edge.data[memory_id][level - 2];
  double * face2;
  ///This vector is used as a dummy in the case that only one edges exists;
  auto eights = std::vector<double>(rowsize,NAN);
  if(edge.faces.size() == 2) {
    face2 = edge.data[memory_id][level - 2];
  } else {
    face2 = eights.data();
    face2pos = 0;
    face2ghostpos = 0;
  }

  for (size_t i = 0; i < rowsize - 2; ++i) {
    fmt::print("{:<8} {:<8} {:<8} {:<8} {:<8}\n", face1[face1ghostpos++],
               face1[face1pos++], mid[midpos++],
               face2[face2pos++], face2[face2ghostpos++]);
  }
  fmt::print("{:<8} {:<8} {:<8} {:<8} {:<8}\n","",face1[face1pos],mid[midpos++], face2[face2pos],"");
  fmt::print("{:<8} {:<8} {:<8} {:<8} {:<8}\n","","",mid[midpos],"","");
}
}
}



#endif /* P2EDGE_HPP */
