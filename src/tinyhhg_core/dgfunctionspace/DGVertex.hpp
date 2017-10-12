#pragma once

namespace hhg {
namespace DGVertex {

template< typename ValueType >
inline void enumerate(Vertex &vertex, const PrimitiveDataID< FunctionMemory< ValueType >, Vertex> &dstId, uint_t level, uint_t& num){
  auto dst = vertex.getData(dstId)->getPointer( level );
  //for each adjacent edge there are two DoF where the first one is owned by the vertex
  for( uint_t i = 0; i < vertex.getNumHigherDimNeighbors(); ++i){
    dst[i * 2] = static_cast< ValueType >(num++);
  }

}





}// namespace DGVertex
}//namespace hhg