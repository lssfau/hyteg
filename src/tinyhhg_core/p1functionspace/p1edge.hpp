#ifndef P1EDGE_HPP
#define P1EDGE_HPP

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/comm.hpp"

namespace hhg
{
namespace P1Edge
{

inline void allocate(Edge& edge, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  edge.data.push_back(std::vector<double*>());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    size_t n_dofs_per_edge = levelinfo::num_microvertices_per_edge(level);
    size_t n_dofs_per_edge_nbr = n_dofs_per_edge - 1;
    size_t num_deps = edge.faces.size();
    size_t total_n_dofs = n_dofs_per_edge + num_deps * n_dofs_per_edge_nbr;
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
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D x = edge.v0->coords;
  Point3D dx = edge.direction / (double) (rowsize - 1);
  x += dx;

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    edge.data[memory_id][level-2][i] = expr(x);
    x += dx;
  }
}

inline void pull_vertices(Edge& edge, size_t memory_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

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

inline void assign(Edge& edge, const std::vector<double>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    double tmp = scalars[0] * edge.data[src_ids[0]][level-2][i];

    for (size_t k = 1; k < src_ids.size(); ++k)
    {
      tmp += scalars[k] * edge.data[src_ids[k]][level-2][i];
    }

    edge.data[dst_id][level-2][i] = tmp;
  }
}

inline void add(Edge& edge, const std::vector<double>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    double tmp = 0.0;

    for (size_t k = 0; k < src_ids.size(); ++k)
    {
      tmp += scalars[k] * edge.data[src_ids[k]][level-2][i];
    }

    edge.data[dst_id][level-2][i] += tmp;
  }
}

inline double dot(Edge& edge, size_t lhs_id, size_t rhs_id, size_t level)
{
  double sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    sp += edge.data[lhs_id][level-2][i] * edge.data[rhs_id][level-2][i];
  }

  return sp;
}

inline void apply(Edge& edge, size_t opr_id, size_t src_id, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  double* opr_data = edge.opr_data[opr_id][level-2];
  double* src = edge.data[src_id][level-2];
  double* dst = edge.data[dst_id][level-2];

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    dst[i] = opr_data[2] * src[i-1] + opr_data[3] * src[i] + opr_data[4] * src[i+1];
    dst[i] += opr_data[0] * src[rowsize + i - 1] + opr_data[1] * src[rowsize + i];

    if (edge.faces.size() == 2)
    {
      dst[i] += opr_data[5] * src[rowsize + rowsize - 1 + i - 1] + opr_data[6] * src[rowsize + rowsize - 1 + i];
    }
  }
}

inline void smooth_gs(Edge& edge, size_t opr_id, size_t dst_id, size_t rhs_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  double* opr_data = edge.opr_data[opr_id][level-2];
  double* dst = edge.data[dst_id][level-2];
  double* rhs = edge.data[rhs_id][level-2];

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    dst[i] = rhs[i] - opr_data[2] * dst[i-1] - opr_data[4] * dst[i+1];
    dst[i] -= opr_data[0] * dst[rowsize + i - 1] + opr_data[1] * dst[rowsize + i];

    if (edge.faces.size() == 2)
    {
      dst[i] -= opr_data[5] * dst[rowsize + rowsize - 1 + i - 1] + opr_data[6] * dst[rowsize + rowsize - 1 + i];
    }

    dst[i] /= opr_data[3];
  }
}

inline void pull_halos(Edge& edge, size_t memory_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t rowsize_halo = rowsize - 1;

  size_t offset = rowsize;
  int rk = hhg::Comm::get().rk;

  auto pull = [rowsize, rowsize_halo, level](Edge& edge, double* edge_data, Face* face, double* face_data)
  {
    if (&edge == face->edges[0])
    {
      if (face->edge_orientation[0] == 1)
      {
        for (size_t i = 0; i < rowsize_halo; ++i)
        {
          edge_data[i] = face_data[rowsize + i];
        }
      }
      else
      {
        for (size_t i = 0; i < rowsize_halo; ++i)
        {
          edge_data[i] = face_data[rowsize + rowsize_halo - 1 - i];
        }
      }
    }
    else if (&edge == face->edges[1])
    {
      if (face->edge_orientation[1] == 1)
      {
        size_t idx = rowsize - 2;
        for (size_t i = 0; i < rowsize_halo; ++i)
        {
          edge_data[i] = face_data[idx];
          idx += rowsize - 1 - i;
        }
      }
      else
      {
        size_t idx = levelinfo::num_microvertices_per_face(level) - 3;
        for (size_t i = 0; i < rowsize_halo; ++i)
        {
          edge_data[i] = face_data[idx];
          idx -= i + 2;
        }
      }
    }
    else
    {
      if (face->edge_orientation[2] == 1)
      {
        size_t idx = levelinfo::num_microvertices_per_face(level) - 2;
        for (size_t i = 0; i < rowsize_halo; ++i)
        {
          edge_data[i] = face_data[idx];
          idx -= i+3;
        }
      }
      else
      {
        size_t idx = 1;
        for (size_t i = 0; i < rowsize_halo; ++i)
        {
          edge_data[i] = face_data[idx];
          idx += rowsize-i;
        }
      }
    }
  };
  
  for (Face* face : edge.faces)
  {
    if (edge.rank == rk)
    {
      double* edge_data = edge.data[memory_id][level-2];

      if (face->rank == rk)
      {
        double* face_data = face->data[memory_id][level-2];
        pull(edge, &edge_data[offset], face, face_data);
        offset += rowsize_halo;
      }
      else
      {
        MPI_Recv(&edge_data[offset], rowsize_halo, MPI_DOUBLE, face->rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        offset += rowsize_halo;
      }
    }
    else if (face->rank == rk)
    {
      double* face_data = face->data[memory_id][level-2];
      double* tmp = new double[rowsize_halo];
      pull(edge, tmp, face, face_data);
      MPI_Send(tmp, rowsize_halo, MPI_DOUBLE, edge.rank, 0, MPI_COMM_WORLD);
      delete[] tmp;
    }
  }
}

inline void prolongate(Edge& edge, size_t memory_id, size_t level)
{
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level);
  size_t i_fine = 1;

  double* edge_data_f = edge.data[memory_id][level-2+1];
  double* edge_data_c = edge.data[memory_id][level-2];

  for (size_t i_coarse = 0; i_coarse < rowsize_coarse-1; ++i_coarse)
  {
    edge_data_f[i_fine] = 0.5 * (edge_data_c[i_coarse] + edge_data_c[i_coarse+1]);
    edge_data_f[i_fine+1] = edge_data_c[i_coarse+1];
    i_fine += 2;
  }
}

inline void restrict(Edge& edge, size_t memory_id, size_t level)
{
  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(level);
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level-1);

  double* edge_data_f = edge.data[memory_id][level-2];
  double* edge_data_c = edge.data[memory_id][level-2-1];

  size_t i_fine = 2;
  size_t i_off = 1;

  for (size_t i_coarse = 1; i_coarse < rowsize_coarse-1; ++i_coarse)
  {
    // mid edge
    edge_data_c[i_coarse] = 0.5 * edge_data_f[i_fine - 1] + edge_data_f[i_fine] + 0.5 * edge_data_f[i_fine + 1];

    for (size_t off_edge = 0; off_edge < edge.faces.size(); ++off_edge)
    {
      edge_data_c[i_coarse] += 0.5 * edge_data_f[rowsize_fine + off_edge * (rowsize_fine-1) + i_off] + 0.5 * edge_data_f[rowsize_fine + off_edge * (rowsize_fine-1) + i_off + 1];
    }

    i_fine += 2;
    i_off += 2;
  }
}

}
}

#endif /* P1EDGE_HPP */