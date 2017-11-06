#pragma once


#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMemory.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

namespace hhg {
namespace edgedof {
namespace macroedge {

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Edge &edge,
                          const PrimitiveDataID < FunctionMemory< ValueType >, Edge> &dstId,
                          uint_t& num)
{
  ValueType *dst = edge.getData(dstId)->getPointer(Level);

  for(uint_t i = 0 ; i < levelinfo::num_microedges_per_edge( Level ) ; ++i){
    dst[hhg::indexing::edgedof::macroedge::horizontalIndex< Level >(i)] = num;
    ++num;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate );

} ///namespace macroedge
} ///namespace edgedof
} ///namespace hhg