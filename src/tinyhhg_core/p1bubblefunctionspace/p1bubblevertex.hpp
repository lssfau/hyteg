#pragma once

#include <fmt/format.h>

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "p1bubblevertexcomm.hpp"
#include "p1bubbleedgecomm.hpp"

#include <core/mpi/MPIWrapper.h>

namespace hhg
{

/// P1Vertex namespace for P1 macro-vertex kernels
namespace P1BubbleVertex
{
//FIXME this can be removed after me moved into walberla namespace
using namespace walberla::mpistubs;

/// Allocate memory for P1 macro-vertex including halos
/// \param vertex Reference to Vertex the allocated memory will belong to
/// \param memory_id Index of the \ref data
/// \param minLevel The minimum level allocated
/// \param maxLevel The maximum level allocated
///
/// This function allocates (1 + number of adjacent edges of \p vertex) real_ts for each level within the range of \p minLevel and \p maxLevel.
/// The pointers to these arrays are saved in the Vertex' data vector at index \p memory_id.
inline void allocate(Vertex& vertex, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  vertex.memory.push_back(new VertexP1BubbleFunctionMemory());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    P1Bubble::getVertexFunctionMemory(vertex, memory_id)->addlevel(level, vertex.edges.size() + vertex.faces.size());
  }
}

inline void free(Vertex& vertex, size_t memory_id)
{
  delete vertex.memory[memory_id];
  vertex.memory[memory_id] = nullptr;
}

inline void interpolate(Vertex& vertex, size_t memory_id, std::function<real_t(const hhg::Point3D&)>& expr, size_t level)
{
  P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level][0] = expr(vertex.coords);
}

inline void assign(Vertex& vertex, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level, DoFType flag)
{
  real_t tmp;

  if (testFlag(vertex.type, flag)) {
    tmp = scalars[0] * P1Bubble::getVertexFunctionMemory(vertex, src_ids[0])->data[level][0];

    for (size_t i = 1; i < src_ids.size(); ++i) {
      tmp += scalars[i] * P1Bubble::getVertexFunctionMemory(vertex, src_ids[i])->data[level][0];
    }

    P1Bubble::getVertexFunctionMemory(vertex, dst_id)->data[level][0] = tmp;
  }

  if (testFlag(hhg::Inner, flag)) {
    size_t offset = 1 + vertex.edges.size();

    for (size_t f = 0; f < vertex.faces.size(); ++f) {
      tmp = scalars[0] * P1Bubble::getVertexFunctionMemory(vertex, src_ids[0])->data[level][offset];

      for (size_t i = 1; i < src_ids.size(); ++i) {
        tmp += scalars[i] * P1Bubble::getVertexFunctionMemory(vertex, src_ids[i])->data[level][offset];
      }

      P1Bubble::getVertexFunctionMemory(vertex, dst_id)->data[level][offset] = tmp;
      ++offset;
    }
  }
}

inline void add(Vertex& vertex, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level, DoFType flag)
{
  real_t tmp;

  if (testFlag(vertex.type, flag)) {
    tmp = 0.0;

    for (size_t i = 0; i < src_ids.size(); ++i) {
      tmp += scalars[i] * P1Bubble::getVertexFunctionMemory(vertex, src_ids[i])->data[level][0];
    }

    P1Bubble::getVertexFunctionMemory(vertex, dst_id)->data[level][0] += tmp;
  }

  if (testFlag(hhg::Inner, flag)) {
    size_t offset = 1 + vertex.edges.size();

    for (size_t f = 0; f < vertex.faces.size(); ++f) {
      real_t tmp = 0.0;

      for (size_t i = 0; i < src_ids.size(); ++i) {
        tmp += scalars[i] * P1Bubble::getVertexFunctionMemory(vertex, src_ids[i])->data[level][offset];
      }

      P1Bubble::getVertexFunctionMemory(vertex, dst_id)->data[level][offset] += tmp;
      ++offset;
    }
  }
}

inline real_t dot(Vertex& vertex, size_t lhs_id, size_t rhs_id, size_t level, DoFType flag)
{
  auto& lhs_data = P1Bubble::getVertexFunctionMemory(vertex, lhs_id)->data[level];
  auto& rhs_data = P1Bubble::getVertexFunctionMemory(vertex, rhs_id)->data[level];

  real_t sp = 0.0;

  if (testFlag(vertex.type, flag)) {
    sp += lhs_data[0] * rhs_data[0];
  }

  if (testFlag(hhg::Inner, flag)) {
    size_t offset = 1 + vertex.edges.size();

    for (size_t f = 0; f < vertex.faces.size(); ++f) {
      sp += lhs_data[offset] * rhs_data[offset];
      ++offset;
    }
  }

  return sp;
}

inline void apply(Vertex& vertex, size_t opr_id, size_t src_id, size_t dst_id, size_t level, UpdateType update, DoFType flag)
{
  auto& stencil_stack = P1Bubble::getVertexStencilMemory(vertex, opr_id)->data[level];
  auto& src = P1Bubble::getVertexFunctionMemory(vertex, src_id)->data[level];
  auto& dst = P1Bubble::getVertexFunctionMemory(vertex, dst_id)->data[level];

  if (testFlag(vertex.type, flag)) {
    // apply first stencil to vertex dof
    auto &opr_data = stencil_stack[0];

    if (update == Replace) {
      dst[0] = opr_data[0] * src[0];
    } else if (update == Add) {
      dst[0] += opr_data[0] * src[0];
    }

    for (size_t i = 0; i < vertex.edges.size() + vertex.faces.size(); ++i) {
      dst[0] += opr_data[i + 1] * src[i + 1];
    }
  }

  if (testFlag(hhg::Inner, flag)) {
    // apply remaining stencils to adjacent cell dofs
    size_t offset = 1 + vertex.edges.size();

    for (size_t f = 0; f < vertex.faces.size(); ++f) {
      auto &opr_data = stencil_stack[f + 1];
      if (update == Replace) {
        dst[offset] = opr_data[0] * src[offset];
      } else if (update == Add) {
        dst[offset] += opr_data[0] * src[offset];
      }

      Face *face = vertex.faces[f];
      size_t v_i = face->vertex_index(vertex);

      dst[offset] += opr_data[v_i + 1] * src[0];

      std::vector<Edge *> adj_edges = face->adjacent_edges(vertex);

      // iterate over adjacent edges
      for (Edge *edge : adj_edges) {
        size_t edge_idx = vertex.edge_index(*edge) + 1;
        Vertex *vertex_j = edge->get_opposite_vertex(vertex);

        size_t v_j = face->vertex_index(*vertex_j);
        dst[offset] += opr_data[v_j + 1] * src[edge_idx];
      }
    }
  }
}

//inline void smooth_gs(Vertex& vertex, size_t opr_id, size_t f_id, size_t rhs_id, size_t level)
//{
//  real_t* opr_data = getVertexStencilMemory(vertex, opr_id)->data[level];
//  real_t* dst = getVertexP1Memory(vertex, f_id)->data[level];
//  real_t* rhs = getVertexP1Memory(vertex, rhs_id)->data[level];
//
//  dst[0] = rhs[0];
//
//  for (size_t i = 0; i < vertex.edges.size(); ++i)
//  {
//    dst[0] -= opr_data[i+1] * dst[i+1];
//  }
//
//  dst[0] /= opr_data[0];
//}

/// Function testing the perfomance and features of the walberla Buffersystem instead of manual mpi comm
//inline void pull_halos_bs(Vertex& vertex, size_t memory_id, size_t level, walberla::mpi::BufferSystem bs) {
//  size_t i = 1;
//  int rk = walberla::mpi::MPIManager::instance()->rank();
//
//  for (Edge* edge: vertex.edges) {
//    if(vertex.rank != rk && edge->rank == rk)
//    {
//      if (edge->vertex_index(vertex) == 0)
//      {
//        bs.sendBuffer(vertex.rank) << edge->data[memory_id][level - 2][1];
//      }
//      else
//      {
//        bs.sendBuffer(vertex.rank) << edge->data[memory_id][level - 2][levelinfo::num_microvertices_per_edge(level) - 2];
//      }
//    }
//    else if (vertex.rank == rk)
//    {
//      if (edge->rank == rk)
//      {
//        if(edge->vertex_index(vertex) == 0)
//        {
//          vertex.data[memory_id][level-2][i] = edge->data[memory_id][level-2][1];
//        }
//        else
//        {
//          vertex.data[memory_id][level-2][i] = edge->data[memory_id][level-2][levelinfo::num_microvertices_per_edge(level) - 2];
//        }
//      }
//      else
//      {
//        MPI_Recv(&vertex.data[memory_id][level-2][i], 1, walberla::MPITrait< real_t >::type(), edge->rank, i, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//      }
//    }
//  }
//}

inline void pull_halos(Vertex& vertex, size_t memory_id, size_t level)
{
  walberla::mpi::SendBuffer sb;
  for(hhg::Edge* edge : vertex.edges){
    hhg::P1BubbleEdge::packDataforVertex(*edge,memory_id,sb,level,vertex);
  }
  walberla::mpi::RecvBuffer rb(sb);
  for(hhg::Edge* edge : vertex.edges) {
    unpackEdgeData(level,vertex,memory_id,rb,*edge);
  }
}

inline void print(Vertex & vertex, size_t memory_id, size_t level) {
  auto &vertexData = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level];
  fmt::print("{:*^80}\n",fmt::format(" Vertex with id: {} ", vertex.id));
  fmt::print("Center dof: {}\n" , vertexData[0]);
  fmt::print("{:<8} {:<8} {:<8} {:<8}\n","vertex","cell","localID","globalID");
  for(size_t i = 0; i < vertex.edges.size(); ++i){
    fmt::print("{:<8} {:<8} {:<8} {:<8}\n",vertexData[i + 1],vertexData[vertex.edges.size() + i],i,vertex.edges[i]->id);
  }
  fmt::print("{:*^80}\n","*");
}


//inline void prolongate(Vertex& vertex, size_t memory_id, size_t level)
//{
//  P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level+1][0] = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level][0];
//}

//inline void restrict(Vertex& vertex, size_t memory_id, size_t level)
//{
//  real_t* vertex_data_f = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level];
//  real_t* vertex_data_c = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level-1];
//
//  vertex_data_c[0] = vertex_data_f[0];
//
//  size_t i = 1;
//  for (Edge* edge : vertex.edges)
//  {
//    vertex_data_c[0] += 0.5 * vertex_data_f[i];
//    i += 1;
//  }
//}

inline void printFunctionMemory(Vertex& vertex, uint_t id, uint_t level)
{
  auto& vertexData = P1Bubble::getVertexFunctionMemory(vertex, id)->data[level];

  std::cout <<  std::string(10,'*');
  std::cout << " Vertex ID: " << vertex.getID().getID();
  std::cout << " Center: " << vertexData[0];
  std::cout << " Memory ID: "<< id;
  std::cout <<  " " << std::string(10,'*') << std::endl;
  std::cout << "Edge ID: |" << " Vertex " << std::endl;
  for(uint_t i = 0; i < vertex.edges.size(); ++i){
    std::cout << std::left << std::setw(9) << vertex.edges[i]->getID().getID() << "|"  << vertexData[1+i] << std::endl;
  }
  std::cout << "Face ID: |" << " Cell " << std::endl;
  for(uint_t i = 0; i < vertex.faces.size(); ++i){
    std::cout << std::left << std::setw(9) << vertex.faces[i]->getID().getID() << "|"  << vertexData[1+vertex.edges.size()+i] << std::endl;
  }
  std::cout <<  std::string(100,'*') << std::endl;
}

}
}