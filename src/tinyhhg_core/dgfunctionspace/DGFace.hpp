#pragma once
#include "DGFaceIndex.hpp"

namespace hhg {
namespace DGFace {

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Face &face,
                          const PrimitiveDataID < FunctionMemory< ValueType >, Face> &dstId,
                          uint_t& num)
{
  ValueType *dst = face.getData(dstId)->getPointer(Level);
  //the outermost gray cells belong to the edge
  const uint_t grayRowsize = levelinfo::num_microvertices_per_edge(Level) - 1;
  uint_t inner_rowsize = grayRowsize - 2;
  for(uint_t col = 1; col < grayRowsize - 2; ++col){
    for(uint_t row = 1; row < inner_rowsize; ++row){
      dst[indexDGFaceFromVertex< Level >(col,row,stencilDirection::CELL_GRAY_NE)] = num;
      ++num;
    }
    --inner_rowsize;
  }

  const uint_t blueRowsize = levelinfo::num_microvertices_per_edge(Level) - 2;
  inner_rowsize = blueRowsize;
  for(uint_t col = 0; col < blueRowsize; ++col){
    for(uint_t row = 0; row < inner_rowsize; ++row){
      dst[indexDGFaceFromGrayDGface< Level >(col,row,stencilDirection::CELL_BLUE_E)] = num;
      ++num;
    }
    --inner_rowsize;
  }

}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )

}//namespace DGFace
}//namespace hhg