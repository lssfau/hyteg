#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "p1bubbleedgecomm.hpp"
#include "p1bubblefacecomm.hpp"

namespace hhg
{
namespace P1BubbleEdge
{
//FIXME this can be removed after we moved into walberla namespace
using namespace walberla::mpistubs;

inline void allocate(Edge& edge, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  edge.memory.push_back(new EdgeP1BubbleFunctionMemory());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    P1Bubble::getEdgeFunctionMemory(edge, memory_id)->addlevel(level, edge.faces.size());
  }
}

inline void free(Edge& edge, size_t memory_id)
{
  delete edge.memory[memory_id];
  edge.memory[memory_id] = nullptr;
}

inline void interpolate(Edge& edge, size_t memory_id, std::function<real_t(const hhg::Point3D&)>& expr, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D x = edge.v0->coords;
  Point3D dx = edge.direction / (real_t) (rowsize - 1);
  x += dx;

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level][i] = expr(x);
    x += dx;
  }
}

inline void pull_vertices(Edge& edge, size_t memory_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  if (edge.v0->rank == walberla::mpi::MPIManager::instance()->rank())
  {
    // local information
    if (edge.rank == walberla::mpi::MPIManager::instance()->rank())
    {
      P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level][0] = P1Bubble::getVertexFunctionMemory(*edge.v0, memory_id)->data[level][0];
    }
    else
    {
      MPI_Send(&P1Bubble::getVertexFunctionMemory(*edge.v0, memory_id)->data[level][0], 1, walberla::MPITrait< real_t >::type(), edge.rank, 0, MPI_COMM_WORLD);
    }
  }
  else if (edge.rank == walberla::mpi::MPIManager::instance()->rank())
  {
    MPI_Recv(&P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level][0], 1, walberla::MPITrait< real_t >::type(), edge.v0->rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (edge.v1->rank == walberla::mpi::MPIManager::instance()->rank())
  {
    // local information
    if (edge.rank == walberla::mpi::MPIManager::instance()->rank())
    {
      P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level][rowsize-1] = P1Bubble::getVertexFunctionMemory(*edge.v1, memory_id)->data[level][0];
    }
    else
    {
      MPI_Send(&P1Bubble::getVertexFunctionMemory(*edge.v1, memory_id)->data[level][0], 1, walberla::MPITrait< real_t >::type(), edge.rank, 0, MPI_COMM_WORLD);
    }
  }
  else if (edge.rank == walberla::mpi::MPIManager::instance()->rank())
  {
    MPI_Recv(&P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level][rowsize-1], 1, walberla::MPITrait< real_t >::type(), edge.v1->rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

inline void assign(Edge& edge, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    real_t tmp = scalars[0] * P1Bubble::getEdgeFunctionMemory(edge, src_ids[0])->data[level][i];

    for (size_t k = 1; k < src_ids.size(); ++k)
    {
      tmp += scalars[k] * P1Bubble::getEdgeFunctionMemory(edge, src_ids[k])->data[level][i];
    }

    P1Bubble::getEdgeFunctionMemory(edge, dst_id)->data[level][i] = tmp;
  }
}

inline void add(Edge& edge, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    real_t tmp = 0.0;

    for (size_t k = 0; k < src_ids.size(); ++k)
    {
      tmp += scalars[k] * P1Bubble::getEdgeFunctionMemory(edge, src_ids[k])->data[level][i];
    }

    P1Bubble::getEdgeFunctionMemory(edge, dst_id)->data[level][i] += tmp;
  }
}

inline real_t dot(Edge& edge, size_t lhs_id, size_t rhs_id, size_t level)
{
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    sp += P1Bubble::getEdgeFunctionMemory(edge, lhs_id)->data[level][i] * P1Bubble::getEdgeFunctionMemory(edge, rhs_id)->data[level][i];
  }

  return sp;
}

template<size_t Level>
inline void apply_tmpl(Edge& edge, size_t opr_id, size_t src_id, size_t dst_id, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto& edge_vertex_stencil = P1Bubble::getEdgeStencilMemory(edge, opr_id)->data[Level];
  auto& src = P1Bubble::getEdgeFunctionMemory(edge, src_id)->data[Level];
  auto& dst = P1Bubble::getEdgeFunctionMemory(edge, dst_id)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    tmp = edge_vertex_stencil[EdgeCoordsVertex::VERTEX_C] * src[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)];

    for (auto neighbor : EdgeCoordsVertex::neighbors_edge)
    {
      tmp += edge_vertex_stencil[neighbor] * src[EdgeCoordsVertex::index<Level>(i, neighbor)];
    }

    for (auto neighbor : EdgeCoordsVertex::neighbors_south)
    {
      tmp += edge_vertex_stencil[neighbor] * src[EdgeCoordsVertex::index<Level>(i, neighbor)];
    }

    if (edge.faces.size() == 2)
    {
      for (auto neighbor : EdgeCoordsVertex::neighbors_north)
      {
        tmp += edge_vertex_stencil[neighbor] * src[EdgeCoordsVertex::index<Level>(i, neighbor)];
      }
    }

    if (update == Replace) {
      dst[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)] = tmp;
    } else if (update == Add) {
      dst[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)] += tmp;
    }
  }
}

SPECIALIZE(void, apply_tmpl, apply)

//inline void smooth_gs(Edge& edge, size_t opr_id, size_t dst_id, size_t rhs_id, size_t level)
//{
//  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
//
//  real_t* opr_data = getEdgeStencilMemory(edge, opr_id)->data[level];
//  real_t* dst = getEdgeP1Memory(edge, dst_id)->data[level];
//  real_t* rhs = getEdgeP1Memory(edge, rhs_id)->data[level];
//
//  for (size_t i = 1; i < rowsize-1; ++i)
//  {
//    dst[i] = rhs[i] - opr_data[2] * dst[i-1] - opr_data[4] * dst[i+1];
//    dst[i] -= opr_data[0] * dst[rowsize + i - 1] + opr_data[1] * dst[rowsize + i];
//
//    if (edge.faces.size() == 2)
//    {
//      dst[i] -= opr_data[5] * dst[rowsize + rowsize - 1 + i - 1] + opr_data[6] * dst[rowsize + rowsize - 1 + i];
//    }
//
//    dst[i] /= opr_data[3];
//  }
//}

inline void pull_halos(Edge& edge, size_t memory_id, size_t level)
{
  //TODO this is WIP only works with one mpi rank!
  walberla::mpi::SendBuffer sb;
  uint_t numberOfFaces = edge.faces.size();
  for(uint_t i = 0; i < numberOfFaces; ++i){
    hhg::P1BubbleFace::packData(level,*edge.faces[i],memory_id,sb,edge);
  }
  walberla::mpi::RecvBuffer rb(sb);
  for(uint_t i = 0; i < numberOfFaces; ++i){
    unpackFaceData(level,edge,memory_id,rb,*edge.faces[i]);
  }
}

//inline void prolongate(Edge& edge, size_t memory_id, size_t level)
//{
//  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level);
//  size_t i_fine = 1;
//
//  real_t* edge_data_f = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level+1];
//  real_t* edge_data_c = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level];
//
//  for (size_t i_coarse = 0; i_coarse < rowsize_coarse-1; ++i_coarse)
//  {
//    edge_data_f[i_fine] = 0.5 * (edge_data_c[i_coarse] + edge_data_c[i_coarse+1]);
//    edge_data_f[i_fine+1] = edge_data_c[i_coarse+1];
//    i_fine += 2;
//  }
//}

//inline void restrict(Edge& edge, size_t memory_id, size_t level)
//{
//  size_t rowsize_fine = levelinfo::num_microvertices_per_edge(level);
//  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(level-1);
//
//  real_t* edge_data_f = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level];
//  real_t* edge_data_c = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level-1];
//
//  size_t i_fine = 2;
//  size_t i_off = 1;
//
//  for (size_t i_coarse = 1; i_coarse < rowsize_coarse-1; ++i_coarse)
//  {
//    // mid edge
//    edge_data_c[i_coarse] = 0.5 * edge_data_f[i_fine - 1] + edge_data_f[i_fine] + 0.5 * edge_data_f[i_fine + 1];
//
//    for (size_t off_edge = 0; off_edge < edge.faces.size(); ++off_edge)
//    {
//      edge_data_c[i_coarse] += 0.5 * edge_data_f[rowsize_fine + off_edge * (rowsize_fine-1) + i_off] + 0.5 * edge_data_f[rowsize_fine + off_edge * (rowsize_fine-1) + i_off + 1];
//    }
//
//    i_fine += 2;
//    i_off += 2;
//  }
//}

//inline void printmatrix(Edge& edge, size_t opr_id, size_t src_id, size_t level)
//{
//  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
//
//  real_t* opr_data = P1Bubble::getEdgeStencilMemory(edge, opr_id)->data[level];
//  real_t* src = P1Bubble::getEdgeFunctionMemory(edge, src_id)->data[level];
//
//  for (size_t i = 1; i < rowsize-1; ++i)
//  {
//
//    fmt::printf("%d\t%d\t%e\n", (size_t)src[i], (size_t)src[i-1], opr_data[2]);
//    fmt::printf("%d\t%d\t%e\n", (size_t)src[i], (size_t)src[i], opr_data[3]);
//    fmt::printf("%d\t%d\t%e\n", (size_t)src[i], (size_t)src[i+1], opr_data[4]);
//
//    fmt::printf("%d\t%d\t%e\n", (size_t)src[i], (size_t)src[rowsize + i - 1], opr_data[0]);
//    fmt::printf("%d\t%d\t%e\n", (size_t)src[i], (size_t)src[rowsize + i], opr_data[1]);
//
//    if (edge.faces.size() == 2)
//    {
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[i], (size_t)src[rowsize + rowsize - 1 + i - 1], opr_data[5]);
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[i], (size_t)src[rowsize + rowsize - 1 + i], opr_data[6]);
//    }
//  }
//}

}
}
