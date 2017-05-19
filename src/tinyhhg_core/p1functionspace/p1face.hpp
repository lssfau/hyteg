#ifndef P1FACE_HPP
#define P1FACE_HPP

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg
{
namespace P1Face
{

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

inline void interpolate(Face& face, size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level)
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

  Point3D d0 = face.edge_orientation[0] * face.edges[0]->direction / (rowsize-1);
  Point3D d2 = -face.edge_orientation[2] * face.edges[2]->direction / (rowsize-1);

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
  //WALBERLA_LOG_DEVEL("Started Pull Edges in p1face");
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  double* edge_data_0 = NULL;
  double* edge_data_1 = NULL;
  double* edge_data_2 = NULL;

  MPI_Request req0;
  MPI_Request req1;
  MPI_Request req2;

  int rk = walberla::mpi::MPIManager::instance()->rank();

  /*if (face.edges[0]->memory[memory_id]->type != P1)
    WALBERLA_LOG_WARNING("IN p1face: memory had not the right type!")
  if (face.edges[1]->memory[memory_id]->type != P1)
    WALBERLA_LOG_WARNING("IN p1face: memory had not the right type!")
  if (face.edges[1]->memory[memory_id]->type != P1)
    WALBERLA_LOG_WARNING("IN p1face: memory had not the right type!")*/
  if (face.edges[0]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_0 = getEdgeP1Memory(*face.edges[0], memory_id)->data[level];
    }
    else
    {
      //WALBERLA_LOG_DEVEL("Sending edge 0");
      MPI_Send(&getEdgeP1Memory(*face.edges[0], memory_id)->data[level][0], rowsize, MPI_DOUBLE, face.rank, face.edges[0]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    //WALBERLA_LOG_DEVEL("Receiving edge 0");
    edge_data_0 = new double[rowsize];
    MPI_Irecv(edge_data_0, rowsize, MPI_DOUBLE, face.edges[0]->rank, face.edges[0]->id, MPI_COMM_WORLD, &req0);
  }

  if (face.edges[1]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_1 = getEdgeP1Memory(*face.edges[1], memory_id)->data[level];
    }
    else
    {
      //WALBERLA_LOG_DEVEL("Sending edge 1");
      MPI_Send(&getEdgeP1Memory(*face.edges[1], memory_id)->data[level][0], rowsize, MPI_DOUBLE, face.rank, face.edges[1]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_1 = new double[rowsize];
    //WALBERLA_LOG_DEVEL("Receiving edge 1");
    MPI_Irecv(edge_data_1, rowsize, MPI_DOUBLE, face.edges[1]->rank, face.edges[1]->id, MPI_COMM_WORLD, &req1);
  }

  if (face.edges[2]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_2 = getEdgeP1Memory(*face.edges[2], memory_id)->data[level];
    }
    else
    {
      //WALBERLA_LOG_DEVEL("Sending edge 2");
      MPI_Send(&getEdgeP1Memory(*face.edges[2], memory_id)->data[level][0], rowsize, MPI_DOUBLE, face.rank, face.edges[2]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_2 = new double[rowsize];
    //WALBERLA_LOG_DEVEL("Receiving edge 2");
    MPI_Irecv(edge_data_2, rowsize, MPI_DOUBLE, face.edges[2]->rank, face.edges[2]->id, MPI_COMM_WORLD, &req2);
  }

  if (face.rank == rk)
  {
    double* face_data = getFaceP1Memory(face, memory_id)->data[level];

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
  //WALBERLA_LOG_DEVEL("Finished Pull Edges in p1face");
}

inline void assign(Face& face, const std::vector<double>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      double tmp = scalars[0] * getFaceP1Memory(face, src_ids[0])->data[level][mr];

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

inline void add(Face& face, const std::vector<double>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      double tmp = 0.0;

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

inline double dot(Face& face, size_t lhs_id, size_t rhs_id, size_t level)
{
  double sp = 0.0;
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

inline void apply(Face& face, size_t opr_id, size_t src_id, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  double* opr_data = getFaceStencilMemory(face, opr_id)->data[level];
  double* src = getFaceP1Memory(face, src_id)->data[level];
  double* dst = getFaceP1Memory(face, dst_id)->data[level];

  size_t br = 1;
  size_t mr = 1 + rowsize ;
  size_t tr = mr + (rowsize - 1);

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      dst[mr] = opr_data[0] * src[br] + opr_data[1] * src[br+1]
                + opr_data[2] * src[mr-1] + opr_data[3] * src[mr] + opr_data[4] * src[mr+1]
                + opr_data[5] * src[tr-1] + opr_data[6] * src[tr];
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

inline void smooth_gs(Face& face, size_t opr_id, size_t dst_id, size_t rhs_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  double* opr_data = getFaceStencilMemory(face, opr_id)->data[level];
  double* dst = getFaceP1Memory(face, dst_id)->data[level];
  double* rhs = getFaceP1Memory(face, rhs_id)->data[level];

  size_t br = 1;
  size_t mr = 1 + rowsize ;
  size_t tr = mr + (rowsize - 1);

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      dst[mr] = (rhs[mr] - opr_data[0] * dst[br] - opr_data[1] * dst[br+1]
                - opr_data[2] * dst[mr-1] - opr_data[4] * dst[mr+1]
                - opr_data[5] * dst[tr-1] - opr_data[6] * dst[tr]) / opr_data[3];
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

inline void prolongate(Face& face, size_t memory_id, size_t level)
{
  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level);
  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(level+1);

  double* face_data_f = getFaceP1Memory(face, memory_id)->data[level+1];
  double* face_data_c = getFaceP1Memory(face, memory_id)->data[level];

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

  double* face_data_f = getFaceP1Memory(face, memory_id)->data[level];
  double* face_data_c = getFaceP1Memory(face, memory_id)->data[level-1];

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

}// namespace P1Face
}// namespace hhg


#endif /* P1FACE_HPP */
