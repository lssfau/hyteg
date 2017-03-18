#ifndef P2FACE_HPP
#define P2FACE_HPP

#include <levelinfo.hpp>
#include <comm.hpp>

namespace hhg
{
namespace P2Face
{

inline void allocate(Face& face, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  face.data.push_back(std::vector<double*>());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    size_t total_n_dofs = levelinfo::num_microvertices_per_face(level) + levelinfo::num_microedges_per_face(level);
    double* new_data = new double[total_n_dofs];
    memset(new_data, 0, total_n_dofs * sizeof(double));
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

inline void interpolate(Face& face, size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);;
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

  size_t mr_c = 2 + rowsize + rowsize - 1;
  size_t inner_rowsize = rowsize-1;

  for (size_t i = 0; i < rowsize-5; ++i)
  {
    x = x0;
    x += (i+1) * d2 *2 + d0 * 2;

    for (size_t j = 0; j < inner_rowsize-5; ++j)
    {
      face.data[memory_id][level-2][mr_c] = expr(x);
      x += d0;
      mr_c += 1;
    }

    mr_c += 4;
    inner_rowsize -= 1;
  }
}
inline void print(Face & face, size_t memory_id, size_t level) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);
  auto facedata = face.data[memory_id][level - 2];
  int pos = 0;

  for (size_t i = 0; i < rowsize; ++i) {

    for (size_t j = 0; j < rowsize - i; ++j) {
      fmt::print("{:<8}", facedata[pos++]);
    }
    fmt::print("\n");
  }
}

inline void pull_edges(Face& face, size_t memory_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);
  double* edge_data_0 = NULL;
  double* edge_data_1 = NULL;
  // double* edge_data_2 = NULL;

  MPI_Request req0;
  MPI_Request req1;
  // MPI_Request req2;

  int rk = hhg::Comm::get().rk;
  
  if (face.edges[0]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_0 = face.edges[0]->data[memory_id][level-2];
    }
    else
    {
      MPI_Send(&face.edges[0]->data[memory_id][level-2][0], rowsize, MPI_DOUBLE, face.rank, face.edges[0]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_0 = new double[rowsize];
    MPI_Irecv(edge_data_0, rowsize, MPI_DOUBLE, face.edges[0]->rank, face.edges[0]->id, MPI_COMM_WORLD, &req0);
  }

  if (face.edges[1]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_1 = face.edges[1]->data[memory_id][level-2];
    }
    else
    {
      MPI_Send(&face.edges[1]->data[memory_id][level-2][0], rowsize, MPI_DOUBLE, face.rank, face.edges[1]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_1 = new double[rowsize];
    MPI_Irecv(edge_data_1, rowsize, MPI_DOUBLE, face.edges[1]->rank, face.edges[1]->id, MPI_COMM_WORLD, &req1);
  }

  // if (face.edges[2]->rank == rk)
  // {
  //   if (face.rank == rk)
  //   {
  //     edge_data_2 = face.edges[2]->data[memory_id][level-2];
  //   }
  //   else
  //   {
  //     MPI_Send(&face.edges[2]->data[memory_id][level-2][0], rowsize, MPI_DOUBLE, face.rank, face.edges[2]->id, MPI_COMM_WORLD);
  //   }
  // }
  // else if (face.rank == rk)
  // {
  //   edge_data_2 = new double[rowsize];
  //   MPI_Irecv(edge_data_2, rowsize, MPI_DOUBLE, face.edges[2]->rank, face.edges[2]->id, MPI_COMM_WORLD, &req2);
  // }

  if (face.rank == rk)
  {
    double* face_data = face.data[memory_id][level-2];

    if (face.edges[0]->rank != rk)
    {
      MPI_Wait(&req0, MPI_STATUS_IGNORE);
    }

    if (face.edges[1]->rank != rk)
    {
      MPI_Wait(&req1, MPI_STATUS_IGNORE);
    }

    // if (face.edges[2]->rank != rk)
    // {
    //   MPI_Wait(&req2, MPI_STATUS_IGNORE);
    // }

    // edge 0
    if (face.edge_orientation[0] == 1)
    {
      for (size_t i = 0; i < rowsize + rowsize - 1; ++i)
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

      for (size_t i = 0; i < rowsize-1; ++i)
      {
        face_data[i + rowsize] = edge_data_0[(rowsize + rowsize - 1) - 1 - i];
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

        fmt::print("idx  = {}\n", idx);
      }

      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_1[i];
        idx += rowsize - 1 - i;

        fmt::print("idx2 = {}\n", idx);
      }
    }
    // else
    // {
    //   size_t idx = levelinfo::num_microvertices_per_face(level) - 1;
    //   for (size_t i = 0; i < rowsize; ++i)
    //   {
    //     face_data[idx] = edge_data_1[i];
    //     idx -= i + 1;
    //   }
    // }

    // if (face.edges[1]->rank != rk)
    // {
    //   delete[] edge_data_1;
    // }

    // // edge 2
    // if (face.edge_orientation[2] == 1)
    // {
    //   size_t idx = levelinfo::num_microvertices_per_face(level) - 1;
    //   for (size_t i = 0; i < rowsize; ++i)
    //   {
    //     face_data[idx] = edge_data_2[i];
    //     idx -= i+2;
    //   }
    // }
    // else
    // {
    //   size_t idx = 0;
    //   for (size_t i = 0; i < rowsize; ++i)
    //   {
    //     face_data[idx] = edge_data_2[i];
    //     idx += rowsize-i;
    //   }
    // }

    // if (face.edges[2]->rank != rk)
    // {
    //   delete[] edge_data_2;
    // }
  }
}

}
}

#endif /* P2FACE_HPP */