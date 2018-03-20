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

template< typename ValueType >
inline void interpolate(const uint_t & Level, Vertex &vertex,
                        const PrimitiveDataID<FunctionMemory< ValueType >, Vertex>& vertexMemoryId,
                        const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>>& srcMemoryIds,
                        std::function<ValueType(const hhg::Point3D &, const std::vector<ValueType>& f)> &expr,
                        const std::shared_ptr< PrimitiveStorage >& storage ) {

  auto vertexMemory = vertex.getData(vertexMemoryId)->getPointer( Level );
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  Point3D dir1;
  Point3D dir2;
  Point3D x;

  std::vector<ValueType*> srcPtr;
  for(auto src : srcMemoryIds){
    srcPtr.push_back(vertex.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVector(srcMemoryIds.size());



  for(auto faceIt : vertex.neighborFaces()){
    Face* face = storage->getFace(faceIt.getID());
    uint_t vertexIdOnFace = face->vertex_index(vertex.getID());
    dir1 = (face->getCoordinates()[(vertexIdOnFace + 1) % 3] - vertex.getCoordinates())/(walberla::real_c(rowsize - 1));
    dir2 = (face->getCoordinates()[(vertexIdOnFace + 2) % 3] - vertex.getCoordinates())/(walberla::real_c(rowsize - 1));
    x = vertex.getCoordinates() + 1.0 / 3.0 * (dir1 + dir2);
    for (size_t k = 0; k < srcPtr.size(); ++k) {
      srcVector[k] = srcPtr[k][vertex.face_index(face->getID()) * 2];
    }

    vertexMemory[vertex.face_index(face->getID()) * 2] = expr(x, srcVector);
  }
}


template< typename ValueType >
inline void assign(const uint_t & Level, Vertex &vertex,
                   const std::vector<ValueType>& scalars,
                   const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &srcIds,
                   const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId) {

  auto dst = vertex.getData(dstId)->getPointer(Level);

  for(uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i){
    uint_t index = i * 2;
    //tmp is necessary since dstId can also be in srcIds
    ValueType tmp = scalars[0] * vertex.getData(srcIds[0])->getPointer( Level )[index];
    for(uint_t k = 1; k < srcIds.size(); ++k){
      tmp += scalars[k] * vertex.getData(srcIds[k])->getPointer( Level )[index];
    }
    dst[index] = tmp;
  }
}



template< typename ValueType >
inline void upwind(const uint_t & Level, Vertex &vertex,
                   const std::shared_ptr< PrimitiveStorage >& storage,
                   const PrimitiveDataID < FunctionMemory< ValueType >, Vertex> &srcId,
                   const PrimitiveDataID < FunctionMemory< ValueType >, Vertex> &dstId,
                   const std::array<PrimitiveDataID< FunctionMemory< ValueType >, Vertex>, 2> &velocityIds,
                   UpdateType updateType) {

  auto src = vertex.getData(srcId)->getPointer( Level );
  auto dst = vertex.getData(dstId)->getPointer( Level );
  auto u = vertex.getData(velocityIds[0])->getPointer(Level);
  auto v = vertex.getData(velocityIds[1])->getPointer(Level);

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  ValueType tmp;
  Point2D u_0, u_1, u_2;
  real_t un_0, un_1, un_2;
  real_t c_up_0, c_up_1, c_up_2;


  for(auto faceIt : vertex.neighborFaces()) {
    Face* face = storage->getFace(faceIt.getID());

    real_t faceArea = std::pow(4.0, -walberla::real_c(Level))*face->area;
    real_t faceAreaInv = 1.0/faceArea;

    uint_t localFaceId = vertex.face_index(face->getID());

    uint_t faceMemoryIndex = 2 * localFaceId;
    uint_t blueMemoryIndex = faceMemoryIndex + 1;

    std::vector<PrimitiveID> adjEdgeIds = face->adjacent_edges(vertex.getID());
    std::vector<Edge*> adjEdges;
    adjEdges.push_back(storage->getEdge(adjEdgeIds[0]));
    adjEdges.push_back(storage->getEdge(adjEdgeIds[1]));

    uint_t v0 = face->vertex_index(vertex.getID());
    uint_t v1 = face->vertex_index(adjEdges[0]->get_opposite_vertex(vertex.getID()));
    uint_t v2 = face->vertex_index(adjEdges[1]->get_opposite_vertex(vertex.getID()));

    uint_t p1EdgeId0 = vertex.edge_index(adjEdgeIds[0]) + 1;
    uint_t p1EdgeId1 = vertex.edge_index(adjEdgeIds[1]) + 1;

    // compute edge directions
    auto d0 = (face->coords[v1] - face->coords[v0])/walberla::real_c(rowsize - 1);
    auto d1 = (face->coords[v0] - face->coords[v2])/walberla::real_c(rowsize - 1);
    auto d2 = (face->coords[v2] - face->coords[v1])/walberla::real_c(rowsize - 1);

    // compute edge lengths
    real_t d0Length = d0.norm();
    real_t d1Length = d1.norm();
    real_t d2Length = d2.norm();

    // compute normals
    auto n_0 = d0.normal2D()/d0Length;
    auto n_1 = d1.normal2D()/d1Length;
    auto n_2 = d2.normal2D()/d2Length;

    real_t faceOrientation = math::faceOrientation2D(face->coords[v0], face->coords[v1], face->coords[v2]);
    n_0 *= faceOrientation;
    n_1 *= faceOrientation;
    n_2 *= faceOrientation;


    u_0[0] = 0.5*(u[0] + u[p1EdgeId0]);
    u_0[1] = 0.5*(v[0] + v[p1EdgeId0]);

    u_1[0] = 0.5*(u[0] + u[p1EdgeId1]);
    u_1[1] = 0.5*(v[0] + v[p1EdgeId1]);

    u_2[0] = 0.5*(u[p1EdgeId0] + u[p1EdgeId1]);
    u_2[1] = 0.5*(v[p1EdgeId0] + v[p1EdgeId1]);

    un_0 = d0Length*u_0.dot(n_0);
    un_1 = d1Length*u_1.dot(n_1);
    un_2 = d2Length*u_2.dot(n_2);


    if (un_0 >= 0) {
      c_up_0 = src[faceMemoryIndex];
    } else {

      // check if neighbor exists and get its id
      if (adjEdges[0]->opposite_face_exists(face->getID())) {
        PrimitiveID oppositeFaceId = adjEdges[0]->get_opposite_face(face->getID());
        uint_t localOppositeFaceId = vertex.face_index(oppositeFaceId);
        uint_t oppositeFaceMemoryIndex = 2 * localOppositeFaceId;
        c_up_0 = src[oppositeFaceMemoryIndex];
      } else {
        // TODO: Handle boundary conditions in this case?
        c_up_0 = 0.0;
      }
    }

    if (un_1 >= 0) {
      c_up_1 = src[faceMemoryIndex];
    } else {

      // check if neighbor exists and get its id
      if (adjEdges[1]->opposite_face_exists(face->getID())) {
        PrimitiveID oppositeFaceId = adjEdges[1]->get_opposite_face(face->getID());
        uint_t localOppositeFaceId = vertex.face_index(oppositeFaceId);
        uint_t oppositeFaceMemoryIndex = 2 * localOppositeFaceId;
        c_up_1 = src[oppositeFaceMemoryIndex];
      } else {
        // TODO: Handle boundary conditions in this case?
        c_up_1 = 0.0;
      }

    }

    if (un_2 >= 0) {
      c_up_2 = src[faceMemoryIndex];
    } else {
      c_up_2 = src[blueMemoryIndex];
    }

    tmp = un_0*c_up_0 + un_1*c_up_1 + un_2*c_up_2;
    tmp *= faceAreaInv;

    if (updateType == Replace) {
      dst[faceMemoryIndex] = tmp;
    } else if (updateType == Add) {
      dst[faceMemoryIndex] += tmp;
    }
  }
}


template< typename ValueType >
inline void projectP1(const uint_t & Level, Vertex &vertex,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      const PrimitiveDataID < FunctionMemory< ValueType >, Vertex> &srcId,
                      const PrimitiveDataID < FunctionMemory< ValueType >, Vertex> &dstId,
                      UpdateType updateType) {

  auto src = vertex.getData(srcId)->getPointer( Level );
  auto dst = vertex.getData(dstId)->getPointer( Level );

  ValueType tmp;

  for(auto faceIt : vertex.neighborFaces()) {
    Face* face = storage->getFace(faceIt.getID());

    uint_t localFaceId = vertex.face_index(face->getID());

    uint_t faceMemoryIndex = 2 * localFaceId;

    std::vector<PrimitiveID> adjEdgeIds = face->adjacent_edges(vertex.getID());
    std::vector<Edge*> adjEdges;
    adjEdges.push_back(storage->getEdge(adjEdgeIds[0]));
    adjEdges.push_back(storage->getEdge(adjEdgeIds[1]));

    uint_t p1EdgeId0 = vertex.edge_index(adjEdgeIds[0]) + 1;
    uint_t p1EdgeId1 = vertex.edge_index(adjEdgeIds[1]) + 1;

    tmp = 1.0/3.0 * (src[0] + src[p1EdgeId0] + src[p1EdgeId1]);

    if (updateType == Replace) {
      dst[faceMemoryIndex] = tmp;
    } else if (updateType == Add) {
      dst[faceMemoryIndex] += tmp;
    }
  }
}


}//namespace DGVertex
}//namespace hhg