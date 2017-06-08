#pragma once

#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"

namespace hhg {
namespace P1BubbleFace {

using walberla::real_t;
using walberla::uint_t;

/*!
 * Unpacks data into the \ref face from \ref recvBuffer.
 * Only packs data owned by edge but not halo data
 * @param face Face to unpack to
 * @param recvBuffer Buffer to unpack from
 * @param memory_id Memory id of the data
 * @param level Multigrid level
 * @param edge Index of the corresponding edge
 */
void packData( Face& face, uint_t memory_id, walberla::mpi::RecvBuffer & recvBuffer, uint_t level, const Edge& edge){

}

}
}