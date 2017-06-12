#pragma once

#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "p1bubblefaceindex.hpp"
#include "core/mpi/RecvBuffer.h"

namespace hhg
{
namespace P1BubbleFace
{

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

/*!
 * Packs data from \ref face for \ref edge into \ref sendBuffer.
 * @param face Face
 * @param sendBuffer Buffer to pack into
 * @param memory_id Memory id of the data
 * @param level Multigrid level
 * @param edge corresponding Edge
 */
template<size_t Level>
inline void packEdgeHaloData_tmpl(Face &face, uint_t memory_id, walberla::mpi::RecvBuffer &sendBuffer, const Edge &edge){

}
SPECIALIZE(void, packEdgeHaloData_tmpl, packEdgeHaloData)

/*!
 * Unpacks data from \ref edge into the \ref face from \ref recvBuffer.
 * Only unpacks data owned by edge but not halo data
 * @param face Face to unpack to
 * @param recvBuffer Buffer to unpack from
 * @param memory_id Memory id of the data
 * @param level Multigrid level
 * @param edge corresponding Edge
 */
template<size_t Level>
inline void unpackEdgeData_tmpl(Face &face, uint_t memory_id, walberla::mpi::RecvBuffer &recvBuffer, const Edge &edge){
  real_t* face_data = P1Bubble::getFaceFunctionMemory(face, memory_id)->data[Level];
  uint_t edge_index = face.edge_index(edge);
  int edge_orientation = face.edge_orientation[edge_index];
  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(Level);
  uint_t j;
  switch(edge_index){
    case 0:
      for(uint_t i = 0; i < v_perEdge; ++i ) {
        edge_orientation == 1 ? j = i : j = v_perEdge-1 - i; // checks if memory is in reverse order
        recvBuffer >> face_data[CoordsVertex::index<Level>(0,j,CoordsVertex::VERTEX_C)];
      }
      break;
    case 1:
      for(uint_t i = 0; i < v_perEdge; ++i ) {
        edge_orientation == 1 ? j = i : j = v_perEdge - i; // checks if memory is in reverse order
        uint_t idx = CoordsVertex::index<Level>(j,v_perEdge-1-j,CoordsVertex::VERTEX_C);
        recvBuffer >> face_data[idx];
      }
      break;
    case 2:
      for(uint_t i = 0; i < v_perEdge; ++i ) {
        edge_orientation == 1 ? j = v_perEdge-1 - i : j = i; // checks if memory is in reverse order
        recvBuffer >> face_data[CoordsVertex::index<Level>(j,0,CoordsVertex::VERTEX_C)];
      }
      break;
    default:
      WALBERLA_LOG_WARNING("unpackData: " << edge << " is not Part of " << face);
      break;
  }
}
SPECIALIZE(void, unpackEdgeData_tmpl, unpackEdgeData)
}
}