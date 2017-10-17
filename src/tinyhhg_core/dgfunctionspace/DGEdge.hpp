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


template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Edge &edge,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge>& edgeMemoryId,
                            std::function<ValueType(const hhg::Point3D &)> &expr,
                            const std::shared_ptr< PrimitiveStorage >& storage ) {

  auto edgeMemory = edge.getData(edgeMemoryId)->getPointer( Level );


  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  Point3D x;
  Point3D x0 = edge.getCoordinates()[0];
  Point3D dx = edge.getDirection()/(walberla::real_c(rowsize - 1));

  Point3D face0d0 = (storage.get()->getFace(edge.neighborFaces()[0].getID())->getCoordinates()[1] -
                     storage.get()->getFace(edge.neighborFaces()[0].getID())->getCoordinates()[0])/
                    (walberla::real_c(rowsize - 1));
  Point3D face0d2 = (storage.get()->getFace(edge.neighborFaces()[0].getID())->getCoordinates()[1] -
                     storage.get()->getFace(edge.neighborFaces()[0].getID())->getCoordinates()[0])/
                    (walberla::real_c(rowsize - 1));

  Point3D face1d0 = (storage.get()->getFace(edge.neighborFaces()[1].getID())->getCoordinates()[1] -
                     storage.get()->getFace(edge.neighborFaces()[1].getID())->getCoordinates()[0])/
                    (walberla::real_c(rowsize - 1));
  Point3D face1d2 = (storage.get()->getFace(edge.neighborFaces()[1].getID())->getCoordinates()[1] -
                     storage.get()->getFace(edge.neighborFaces()[1].getID())->getCoordinates()[0])/
                    (walberla::real_c(rowsize - 1));

  // gray south cells
  x = x0 + 1.0/3.0 * (face0d0 + face0d2);
  for (size_t i = 1; i < rowsize - 3; ++i) {
    edgeMemory[DGEdge::indexDGFaceFromVertex< Level >(i,stencilDirection::CELL_GRAY_SE)] = expr(x);
    x += dx;
  }

  // gray north cells
  x = x0 + 1.0/3.0 * (face1d0 + face1d2);
  for (size_t i = 1; i < rowsize - 3; ++i) {
    edgeMemory[DGEdge::indexDGFaceFromVertex< Level >(i,stencilDirection::CELL_GRAY_NE)] = expr(x);
    x += dx;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, interpolateTmpl, interpolate)


template< typename ValueType, uint_t Level >
inline void assignTmpl(Edge &edge,
                       const std::vector<ValueType>& scalars,
                       const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &srcIds,
                       const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto dst = edge.getData(dstId)->getPointer(Level);

  // gray south cells
  for (size_t i = 1; i < rowsize - 3; ++i) {
    uint_t cellIndex = DGEdge::indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE);
    ValueType tmp = scalars[0]*edge.getData(srcIds[0])->getPointer( Level )[cellIndex];
    for (uint_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( Level )[cellIndex];
    }
    dst[cellIndex] = tmp;
  }

  for (size_t i = 1; i < rowsize - 3; ++i) {
    uint_t cellIndex = DGEdge::indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE);
    ValueType tmp = scalars[0]*edge.getData(srcIds[0])->getPointer( Level )[cellIndex];
    for (uint_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( Level )[cellIndex];
    }
    dst[cellIndex] = tmp;
  }

}

SPECIALIZE_WITH_VALUETYPE(void, assignTmpl, assign)


}//namespace DGEdge
}//namespace hhg