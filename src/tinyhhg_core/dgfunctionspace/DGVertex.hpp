#pragma once

namespace hhg {
namespace DGVertex {

template< typename ValueType >
inline void enumerate(Vertex &vertex,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex> &dstId,
                      const uint_t level, uint_t& num)
{
  auto dst = vertex.getData(dstId)->getPointer( level );
  //for each adjacent edge there are two DoF where the first one is owned by the vertex
  for( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i){
    dst[i * 2] = static_cast< ValueType >(num++);
  }

}

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Vertex &vertex,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Vertex>& vertexMemoryId,
                            std::function<ValueType(const hhg::Point3D &)> &expr,
                            const std::shared_ptr< PrimitiveStorage >& storage ) {

  auto vertexMemory = vertex.getData(vertexMemoryId)->getPointer( Level );
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  Point3D dir1;
  Point3D dir2;
  Point3D x;

  for(auto faceIt : vertex.neighborFaces()){
    Face* face = storage->getFace(faceIt.getID());
    uint_t vertexIdOnFace = face->vertex_index(vertex.getID());
    dir1 = (face->getCoordinates()[(vertexIdOnFace + 1) % 3] - vertex.getCoordinates())/(walberla::real_c(rowsize - 1));
    dir2 = (face->getCoordinates()[(vertexIdOnFace + 2) % 3] - vertex.getCoordinates())/(walberla::real_c(rowsize - 1));
    x = vertex.getCoordinates() + 1.0 / 3.0 * (dir1 + dir2);
    vertexMemory[vertex.face_index(face->getID()) * 2] = expr(x);
  }
}

SPECIALIZE_WITH_VALUETYPE(void, interpolateTmpl, interpolate)

template< typename ValueType, uint_t Level >
inline void assignTmpl(Vertex &vertex,
                       const std::vector<ValueType>& scalars,
                       const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &srcIds,
                       const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId) {

  auto dst = vertex.getData(dstId)->getPointer(Level);

  for(uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i){
    uint_t index = i * 2;
    //tmp is necessary since dstId can also be in srcIds
    ValueType tmp = scalars[0] * *vertex.getData(srcIds[0])->getPointer( Level )[index];
    for(uint_t k = 1; k < srcIds.size(); ++k){
      tmp += scalars[k] * *vertex.getData(srcIds[k])->getPointer( Level )[index];
    }
    dst[index] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, assignTmpl, assign)




}//namespace DGVertex
}//namespace hhg