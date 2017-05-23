#ifndef P1FACE_HPP
#define P1FACE_HPP

#include "tinyhhg_core/levelinfo.hpp"

namespace hhg
{
namespace P1Face
{

enum Dir
{
  S  = 0,
  SE = 1,
  W  = 2,
  C  = 3,
  E  = 4,
  NW = 5,
  N  = 6
};

const Dir neighbors_with_center[] = {S, SE, W, C, E, NW, N};
const Dir neighbors[] = {S, SE, W, E, NW, N};

template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  size_t h = levelinfo::num_microvertices_per_edge(Level);
  size_t n = h * (h + 1) / 2;
  size_t center = (n - (h-row)*(h-row+1)/2) + col;
  switch (dir) {
    case C:
      return center;
    case N:
      return center + h - row;
    case E:
      return center + 1;
    case S:
      return center - h - 1 + row;
    case W:
      return center - 1;
    case SE:
      return center - h + row;
    case NW:
      return center + h - row - 1;
  }
  return 0;
}


//FIXME this can be removed after we are in waberla namespace
using namespace walberla::mpistubs;

inline void allocate(Face& face, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  face.data.push_back(std::vector<walberla::real_t*>());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    size_t total_n_dofs = levelinfo::num_microvertices_per_face(level);
    walberla::real_t* new_data = new walberla::real_t[total_n_dofs];
    memset(new_data, 0, total_n_dofs * sizeof(walberla::real_t));
    face.data[memory_id].push_back(new_data);
  }
}

inline void free(Face& face, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    delete[] face.data[memory_id][level - minLevel];
  }
}

template<size_t Level>
inline void interpolate(Face& face, size_t memory_id, std::function<walberla::real_t(const hhg::Point3D&)>& expr)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;

  if (face.edge_orientation[0] == 1)
  {
    x0 = face.edges[0]->v0->coords;
  }
  else
  {
    x0 = face.edges[0]->v1->coords;
  }

  Point3D d0 = face.edge_orientation[0] * face.edges[0]->direction / (walberla::real_c(rowsize-1));
  Point3D d2 = -face.edge_orientation[2] * face.edges[2]->direction / (walberla::real_c(rowsize-1));

  size_t mr_c = 1 + rowsize;
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize-3; ++i)
  {
    x = x0;
    x += (i+1) * d2 + d0;

    for (size_t j = 0; j < inner_rowsize-3; ++j)
    {
      face.data[memory_id][Level-2][mr_c] = expr(x);
      x += d0;
      mr_c += 1;
    }

    mr_c += 2;
    inner_rowsize -= 1;
  }
}

template<size_t Level>
inline void pull_edges(Face& face, size_t memory_id)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  walberla::real_t* edge_data_0 = NULL;
  walberla::real_t* edge_data_1 = NULL;
  walberla::real_t* edge_data_2 = NULL;

  MPI_Request req0;
  MPI_Request req1;
  MPI_Request req2;

  int rk = walberla::mpi::MPIManager::instance()->rank();

  if (face.edges[0]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_0 = face.edges[0]->data[memory_id][Level-2];
    }
    else
    {
      MPI_Send(&face.edges[0]->data[memory_id][Level-2][0], rowsize, walberla::MPITrait< walberla::real_t >::type(), face.rank, face.edges[0]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_0 = new walberla::real_t[rowsize];
    MPI_Irecv(edge_data_0, rowsize, walberla::MPITrait< walberla::real_t >::type(), face.edges[0]->rank, face.edges[0]->id, MPI_COMM_WORLD, &req0);
  }

  if (face.edges[1]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_1 = face.edges[1]->data[memory_id][Level-2];
    }
    else
    {
      MPI_Send(&face.edges[1]->data[memory_id][Level-2][0], rowsize, walberla::MPITrait< walberla::real_t >::type(), face.rank, face.edges[1]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_1 = new walberla::real_t[rowsize];
    MPI_Irecv(edge_data_1, rowsize, walberla::MPITrait< walberla::real_t >::type(), face.edges[1]->rank, face.edges[1]->id, MPI_COMM_WORLD, &req1);
  }

  if (face.edges[2]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_2 = face.edges[2]->data[memory_id][Level-2];
    }
    else
    {
      MPI_Send(&face.edges[2]->data[memory_id][Level-2][0], rowsize, walberla::MPITrait< walberla::real_t >::type(), face.rank, face.edges[2]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_2 = new walberla::real_t[rowsize];
    MPI_Irecv(edge_data_2, rowsize, walberla::MPITrait< walberla::real_t >::type(), face.edges[2]->rank, face.edges[2]->id, MPI_COMM_WORLD, &req2);
  }

  if (face.rank == rk)
  {
    walberla::real_t* face_data = face.data[memory_id][Level-2];

    if (face.edges[0]->rank != rk)
    {
      MPI_Wait(&req0,MPI_STATUS_IGNORE);
    }

    if (face.edges[1]->rank != rk)
    {
      MPI_Wait(&req1,MPI_STATUS_IGNORE);
    }

    if (face.edges[2]->rank != rk)
    {
      MPI_Wait(&req2,MPI_STATUS_IGNORE);
    }

    // edge 0
    if (face.edge_orientation[0] == 1)
    {
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[i] = edge_data_0[i];
      }
    }
    else
    {
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[i] = edge_data_0[rowsize - 1 - i];
      }
    }

    if (face.edges[0]->rank != rk)
    {
      delete[] edge_data_0;
    }

    // edge 1
    if (face.edge_orientation[1] == 1)
    {
      size_t idx = rowsize - 1;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_1[i];
        idx += rowsize - 1 - i;
      }
    }
    else
    {
      size_t idx = levelinfo::num_microvertices_per_face(Level) - 1;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_1[i];
        idx -= i + 1;
      }
    }

    if (face.edges[1]->rank != rk)
    {
      delete[] edge_data_1;
    }

    // edge 2
    if (face.edge_orientation[2] == 1)
    {
      size_t idx = levelinfo::num_microvertices_per_face(Level) - 1;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_2[i];
        idx -= i+2;
      }
    }
    else
    {
      size_t idx = 0;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_2[i];
        idx += rowsize-i;
      }
    }

    if (face.edges[2]->rank != rk)
    {
      delete[] edge_data_2;
    }
  }
}

template<size_t Level>
inline void assign(Face& face, const std::vector<walberla::real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {

      walberla::real_t tmp = face.data[src_ids[0]][Level-2][index<Level>(i, j, C)];

      for (size_t k = 1; k < src_ids.size(); ++k)
      {
        tmp += scalars[k] * face.data[src_ids[k]][Level-2][index<Level>(i, j, C)];
      }

      face.data[dst_id][Level-2][index<Level>(i, j, C)] = tmp;

    }
    --inner_rowsize;
  }
}

template<size_t Level>
inline void add(Face& face, const std::vector<walberla::real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {

      walberla::real_t tmp = 0.0;

      for (size_t k = 0; k < src_ids.size(); ++k)
      {
        tmp += scalars[k] * face.data[src_ids[k]][Level-2][index<Level>(i, j, C)];
      }

      face.data[dst_id][Level-2][index<Level>(i, j, C)] += tmp;

    }
    --inner_rowsize;
  }
}

template<size_t Level>
inline walberla::real_t dot(Face& face, size_t lhs_id, size_t rhs_id)
{
  walberla::real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      sp += face.data[lhs_id][Level-2][index<Level>(i, j, C)] * face.data[rhs_id][Level-2][index<Level>(i, j, C)];
    }
    --inner_rowsize;
  }

  return sp;
}

template<size_t Level>
inline void apply(Face& face, size_t opr_id, size_t src_id, size_t dst_id, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  walberla::real_t* opr_data = face.opr_data[opr_id][Level-2];
  walberla::real_t* src = face.data[src_id][Level-2];
  walberla::real_t* dst = face.data[dst_id][Level-2];

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      walberla::real_t tmp = 0.0;

      for (auto neighbor : neighbors_with_center)
      {
        tmp += opr_data[neighbor] * src[index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[index<Level>(i, j, C)] = tmp;
      }
      else {
        dst[index<Level>(i, j, C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

template<size_t Level>
inline void smooth_gs(Face& face, size_t opr_id, size_t dst_id, size_t rhs_id)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  walberla::real_t* opr_data = face.opr_data[opr_id][Level-2];
  walberla::real_t* dst = face.data[dst_id][Level-2];
  walberla::real_t* rhs = face.data[rhs_id][Level-2];

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      walberla::real_t tmp = rhs[index<Level>(i, j, C)];

      for (auto neighbor : neighbors)
      {
        tmp -= opr_data[neighbor] * dst[index<Level>(i, j, neighbor)];
      }

      dst[index<Level>(i, j, C)] = tmp / opr_data[C];
    }
    --inner_rowsize;
  }
}

template<size_t Level>
inline void prolongate(Face& face, size_t memory_id)
{
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(Level);
  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(Level+1);

  walberla::real_t* face_data_f = face.data[memory_id][Level-2+1];
  walberla::real_t* face_data_c = face.data[memory_id][Level-2];

  size_t mr_c = 1;
  size_t mr_f = rowsize_fine + 2;

  size_t i_rowsize_coarse = rowsize_coarse;

  for (size_t i_coarse = 0; i_coarse < rowsize_coarse-2; ++i_coarse)
  {
    for (size_t j_coarse = 0; j_coarse < i_rowsize_coarse-3; ++j_coarse)
    {
      face_data_f[mr_f] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse]);
      face_data_f[mr_f-1] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse - 1]);
      face_data_f[mr_f + rowsize_fine - 1 - 1] = 0.5 * (face_data_c[mr_c + i_rowsize_coarse] + face_data_c[mr_c + i_rowsize_coarse - 1]);

      face_data_f[mr_f + rowsize_fine - 1] = face_data_c[mr_c + i_rowsize_coarse];

      mr_c += 1;
      mr_f += 2;
    }

    face_data_f[mr_f] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse]);
    face_data_f[mr_f-1] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse - 1]);
    face_data_f[mr_f + rowsize_fine - 1 - 1] = 0.5 * (face_data_c[mr_c + i_rowsize_coarse] + face_data_c[mr_c + i_rowsize_coarse - 1]);

    mr_c += 3;
    mr_f += rowsize_fine - 1 + 3;

    rowsize_fine -= 2;
    i_rowsize_coarse -= 1;
  }
}

template<size_t Level>
inline void restrict(Face& face, size_t memory_id)
{
  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(Level);
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(Level-1);

  walberla::real_t* face_data_f = face.data[memory_id][Level-2];
  walberla::real_t* face_data_c = face.data[memory_id][Level-2-1];

  size_t mr_c = 1 + rowsize_coarse;

  size_t br_f = rowsize_fine + 2;
  size_t mr_f = br_f + rowsize_fine - 1;
  size_t tr_f = mr_f + rowsize_fine - 2;

  size_t i_rowsize_coarse = rowsize_coarse;

  for (size_t i = 0; i < rowsize_coarse - 3; ++i)
  {
    for (size_t j = 0; j < i_rowsize_coarse - 3; ++j)
    {
      face_data_c[mr_c] = 0.5 * (face_data_f[br_f] + face_data_f[br_f+1]);
      face_data_c[mr_c] += 0.5 * face_data_f[mr_f-1] + face_data_f[mr_f] + 0.5 * face_data_f[mr_f+1];
      face_data_c[mr_c] += 0.5 * (face_data_f[tr_f-1] + face_data_f[tr_f]);

      mr_c += 1;

      br_f += 2;
      mr_f += 2;
      tr_f += 2;
    }

    br_f += (rowsize_fine-2) + 4;
    mr_f += (rowsize_fine-3) + 3;
    tr_f += (rowsize_fine-4) + 2;

    mr_c += 2;

    rowsize_fine -= 2;
    i_rowsize_coarse -= 1;
  }
}

/// Checks if a given index is a the boundary of the face
/// \param index The index which should be checked
/// \param length Size of the triangle in the first dimension
bool is_boundary(size_t index, size_t length)
{
  if(index < length) return true;
  while(index >= length){
    index -= length;
    length--;
  }
  return(index == 0 || index == (length -1));
}

template<size_t Level>
inline void printmatrix(Face& face, size_t opr_id, size_t src_id)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  walberla::real_t* opr_data = face.opr_data[opr_id][Level-2];
  walberla::real_t* src = face.data[src_id][Level-2];
  size_t br = 1;
  size_t mr = 1 + rowsize ;
  size_t tr = mr + (rowsize - 1);

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[br], opr_data[0]);
      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[br+1], opr_data[1]);
      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[mr-1], opr_data[2]);
      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[mr], opr_data[3]);
      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[mr+1], opr_data[4]);
      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[tr-1], opr_data[5]);
      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[tr], opr_data[6]);

      br += 1;
      mr += 1;
      tr += 1;
    }

    br += 3;
    mr += 2;
    tr += 1;
    --inner_rowsize;
  }
}

}// namespace P1Face
}// namespace hhg


#endif /* P1FACE_HPP */
