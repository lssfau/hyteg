#pragma once
#include "DGEdgeIndex.hpp"

namespace hhg {
namespace DGEdge {

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Edge &edge,
                      const PrimitiveDataID < FunctionMemory< ValueType >, Edge> &dstId,
                      uint_t& num)
{
  ValueType* dst = edge.getData(dstId)->getPointer( Level );
  std::vector< stencilDirection > dirs;
  dirs.push_back(stencilDirection::CELL_GRAY_SE);
  if(edge.getNumHigherDimNeighbors() == 2){
    dirs.push_back(stencilDirection::CELL_GRAY_NE);
  }
  for(auto dir : dirs) {
    for (uint_t i = 1; i < (hhg::levelinfo::num_microvertices_per_edge( Level ) - 2); ++i) {
      dst[indexDGFaceFromVertex< Level >(i, dir)] = num;
      ++num;
    }
  }
}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )


}//namespace DGEdge
}//namespace hhg