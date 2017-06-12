#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "p1bubbleedge.hpp"

namespace hhg {
namespace P1BubbleEdge {

using walberla::uint_t;

/*!
 * Packs data on the edge into the \ref sendBuffer.
 * Only packs data owned by edge but not halo data
 * @param edge Edge
 * @param sendBuffer Buffer to pack into
 * @param memory_id Memory id of the data
 * @param level Multigrid level
 */
inline void packData(Edge &edge, uint_t memory_id, walberla::mpi::SendBuffer &sendBuffer, uint_t level) {
  real_t* edge_data = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level];
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  //TODO change to index function
  for (uint_t i = 0; i < rowsize; ++i) {
    sendBuffer << edge_data[i];
  }
}

}// namespace P1BubbleEdge
}// namespace hhg
