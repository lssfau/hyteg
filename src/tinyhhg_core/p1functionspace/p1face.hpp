#ifndef P1FACE_HPP
#define P1FACE_HPP

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg
{
namespace P1Face
{
//FIXME this can be removed after we are in waberla namespace
using namespace walberla::mpistubs;

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


inline void allocate(Face& face, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  face.memory.push_back(new FaceP1Memory());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    getFaceP1Memory(face, memory_id)->addlevel(level);
  }
}

inline void free(Face& face, size_t memory_id)
{
  delete face.memory[memory_id];
  face.memory[memory_id] = nullptr;
}

inline void interpolate(Face& face, size_t memory_id, std::function<real_t(const hhg::Point3D&)>& expr, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
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
      getFaceP1Memory(face, memory_id)->data[level][mr_c] = expr(x);
      x += d0;
      mr_c += 1;
    }

    mr_c += 2;
    inner_rowsize -= 1;
  }
}

inline void pull_edges(Face& face, size_t memory_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  real_t* edge_data_0 = NULL;
  real_t* edge_data_1 = NULL;
  real_t* edge_data_2 = NULL;

  MPI_Request req0;
  MPI_Request req1;
  MPI_Request req2;

  int rk = walberla::mpi::MPIManager::instance()->rank();

  if (face.edges[0]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_0 = getEdgeP1Memory(*face.edges[0], memory_id)->data[level];
    }
    else
    {
      MPI_Send(&getEdgeP1Memory(*face.edges[0], memory_id)->data[level][0], rowsize, walberla::MPITrait< real_t >::type(), face.rank, face.edges[0]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_0 = new real_t[rowsize];
    MPI_Irecv(edge_data_0, rowsize, walberla::MPITrait< real_t >::type(), face.edges[0]->rank, face.edges[0]->id, MPI_COMM_WORLD, &req0);
  }

  if (face.edges[1]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_1 = getEdgeP1Memory(*face.edges[1], memory_id)->data[level];
    }
    else
    {
      MPI_Send(&getEdgeP1Memory(*face.edges[1], memory_id)->data[level][0], rowsize, walberla::MPITrait< real_t >::type(), face.rank, face.edges[1]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_1 = new real_t[rowsize];
    MPI_Irecv(edge_data_1, rowsize, walberla::MPITrait< real_t >::type(), face.edges[1]->rank, face.edges[1]->id, MPI_COMM_WORLD, &req1);
  }

  if (face.edges[2]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_2 = getEdgeP1Memory(*face.edges[2], memory_id)->data[level];
    }
    else
    {
      MPI_Send(&getEdgeP1Memory(*face.edges[2], memory_id)->data[level][0], rowsize, walberla::MPITrait< real_t >::type(), face.rank, face.edges[2]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_2 = new real_t[rowsize];
    MPI_Irecv(edge_data_2, rowsize, walberla::MPITrait< real_t >::type(), face.edges[2]->rank, face.edges[2]->id, MPI_COMM_WORLD, &req2);
  }

  if (face.rank == rk)
  {
    real_t* face_data = getFaceP1Memory(face, memory_id)->data[level];

    if (face.edges[0]->rank != rk)
    {
      MPI_Wait(&req0, MPI_STATUS_IGNORE);
    }

    if (face.edges[1]->rank != rk)
    {
      MPI_Wait(&req1, MPI_STATUS_IGNORE);
    }

    if (face.edges[2]->rank != rk)
    {
      MPI_Wait(&req2, MPI_STATUS_IGNORE);
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
      size_t idx = levelinfo::num_microvertices_per_face(level) - 1;
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
      size_t idx = levelinfo::num_microvertices_per_face(level) - 1;
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

inline void assign(Face& face, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      real_t tmp = scalars[0] * getFaceP1Memory(face, src_ids[0])->data[level][mr];

      for (size_t k = 1; k < src_ids.size(); ++k)
      {
        tmp += scalars[k] * getFaceP1Memory(face, src_ids[k])->data[level][mr];
      }
      getFaceP1Memory(face, dst_id)->data[level][mr] = tmp;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

inline void add(Face& face, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      real_t tmp = 0.0;

      for (size_t k = 0; k < src_ids.size(); ++k)
      {
        tmp += scalars[k] * getFaceP1Memory(face, src_ids[k])->data[level][mr];
      }

      getFaceP1Memory(face, dst_id)->data[level][mr] += tmp;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

inline real_t dot(Face& face, size_t lhs_id, size_t rhs_id, size_t level)
{
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      sp += getFaceP1Memory(face, lhs_id)->data[level][mr] * getFaceP1Memory(face, rhs_id)->data[level][mr];
      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }

  return sp;
}

template<size_t Level>
inline void apply_tmpl(Face& face, size_t opr_id, size_t src_id, size_t dst_id, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  real_t* opr_data = getFaceStencilMemory(face, opr_id)->data[Level];
  real_t* src = getFaceP1Memory(face, src_id)->data[Level];
  real_t* dst = getFaceP1Memory(face, dst_id)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      tmp = opr_data[C] * src[index<Level>(i, j, C)];

      for (auto neighbor : neighbors)
      {
        tmp += opr_data[neighbor] * src[index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[index<Level>(i, j, C)] = tmp;
      } else if (update == Add) {
        dst[index<Level>(i, j, C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)

template<size_t Level>
inline void smooth_gs_tmpl(Face& face, size_t opr_id, size_t dst_id, size_t rhs_id)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  real_t* opr_data = getFaceStencilMemory(face, opr_id)->data[Level];
  real_t* dst = getFaceP1Memory(face, dst_id)->data[Level];
  real_t* rhs = getFaceP1Memory(face, rhs_id)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      tmp = rhs[index<Level>(i, j, C)];

      for (auto neighbor : neighbors)
      {
        tmp -= opr_data[neighbor] * dst[index<Level>(i, j, neighbor)];
      }

      dst[index<Level>(i, j, C)] = tmp / opr_data[C];
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, smooth_gs_tmpl, smooth_gs)

inline void prolongate(Face& face, size_t memory_id, size_t level)
{
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level);
  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(level+1);

  real_t* face_data_f = getFaceP1Memory(face, memory_id)->data[level+1];
  real_t* face_data_c = getFaceP1Memory(face, memory_id)->data[level];

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

inline void restrict(Face& face, size_t memory_id, size_t level)
{
  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(level);
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level-1);

  real_t* face_data_f = getFaceP1Memory(face, memory_id)->data[level];
  real_t* face_data_c = getFaceP1Memory(face, memory_id)->data[level-1];

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

inline void printmatrix(Face& face, size_t opr_id, size_t src_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  real_t* opr_data = getFaceStencilMemory(face, opr_id)->data[level];
  real_t* src = getFaceP1Memory(face, src_id)->data[level];
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
