#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"

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
 * @param reverse if set to TRUE the data will be packed onto the buffer in reverse order
 */
void packData( Edge& edge, uint_t memory_id, walberla::mpi::SendBuffer & sendBuffer, uint_t level, bool reverse){
  real_t* edge_data = getEdgeP1BubbleMemory(edge, memory_id)->data[level];
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  if(reverse) {
    for (uin_t i = rowsize-1; i >= 0; --i) {
      sendBuffer << edge_data[i];
    }
  } else {
    for (uin_t i = 0; i < rowsize; ++i) {
      sendBuffer << edge_data[i];
    }
  }
}

}// namespace P1BubbleEdge
}// namespace hhg
