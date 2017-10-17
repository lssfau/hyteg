#pragma once
#include "DGEdgeIndex.hpp"
#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"

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
inline void upwindTmpl(Edge &edge,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       const PrimitiveDataID < FunctionMemory< ValueType >, Edge> &srcId,
                       const PrimitiveDataID < FunctionMemory< ValueType >, Edge> &dstId,
                       const std::array<PrimitiveDataID< FunctionMemory< ValueType >, Edge>, 2> &velocityIds,
                       UpdateType updateType) {

  using namespace P1Edge;

  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto u = edge.getData(velocityIds[0])->getPointer(Level);
  auto v = edge.getData(velocityIds[1])->getPointer(Level);

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  ValueType tmp;
  Point2D u_0, u_1, u_2;
  real_t un_0, un_1, un_2;
  real_t c_up_0, c_up_1, c_up_2;

  // first face (south)
  {
    Face *face = storage->getFace(edge.neighborFaces()[0]);
    real_t faceArea = std::pow(4.0, -walberla::real_c(Level))*face->area;
    real_t faceAreaInv = 1.0/faceArea;

    auto oppositeVertex = face->get_vertex_opposite_to_edge(edge.getID());

    uint_t v0 = face->vertex_index(edge.getVertexID0());
    uint_t v1 = face->vertex_index(edge.getVertexID1());
    uint_t v2 = face->vertex_index(oppositeVertex);

    // get edge directions
    auto d0 = (face->coords[v1] - face->coords[v0])/walberla::real_c(rowsize - 1);
    auto d1 = (face->coords[v2] - face->coords[v1])/walberla::real_c(rowsize - 1);
    auto d2 = (face->coords[v0] - face->coords[v2])/walberla::real_c(rowsize - 1);

    // compute edge lengths
    real_t d0Length = d0.norm();
    real_t d1Length = d1.norm();
    real_t d2Length = d2.norm();

    // compute normals
    auto n_0 = -1.0*d0.normal2D()/d0Length;
    auto n_1 = -1.0*d1.normal2D()/d1Length;
    auto n_2 = -1.0*d2.normal2D()/d2Length;

    for (uint_t i = 1; i < rowsize - 2; ++i) {
      u_0[0] = 0.5*(u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]);
      u_0[1] = 0.5*(v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]);

      u_1[0] = 0.5*(u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]
          + u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_SE)]);
      u_1[1] = 0.5*(v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]
          + v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_SE)]);

      u_2[0] = 0.5*(u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_SE)]);
      u_2[1] = 0.5*(v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_SE)]);

      // CENTER <-- CELL_GRAY_SE (i)
      // NORTH  <-- CELL_GRAY_NE (i)
      // WEST   <-- CELL_BLUE_SE (i)
      // EAST   <-- CELL_BLUE_SE (i+1)

      un_0 = d0Length*u_0.dot(n_0);
      un_1 = d1Length*u_1.dot(n_1);
      un_2 = d2Length*u_2.dot(n_2);

      if (un_0 >= 0) {
        c_up_0 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE)];
      } else {
        c_up_0 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE)];
      }

      if (un_1 >= 0) {
        c_up_1 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE)];
      } else {
        c_up_1 = src[indexDGFaceFromVertex<Level>(i + 1, stencilDirection::CELL_BLUE_SE)];
      }

      if (un_2 >= 0) {
        c_up_2 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE)];
      } else {
        c_up_2 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_BLUE_SE)];
      }

      tmp = un_0*c_up_0 + un_1*c_up_1 + un_2*c_up_2;
      tmp *= faceAreaInv;

      if (updateType==Replace) {
        dst[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE)] = tmp;
      } else if (updateType==Add) {
        dst[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE)] += tmp;
      }
    }
  }

  // second face (north)
  if (edge.getNumNeighborFaces() == 2)
  {
    Face *face = storage->getFace(edge.neighborFaces()[1]);
    real_t faceArea = std::pow(4.0, -walberla::real_c(Level))*face->area;
    real_t faceAreaInv = 1.0/faceArea;

    auto oppositeVertex = face->get_vertex_opposite_to_edge(edge.getID());

    uint_t v0 = face->vertex_index(edge.getVertexID0());
    uint_t v1 = face->vertex_index(edge.getVertexID1());
    uint_t v2 = face->vertex_index(oppositeVertex);

    // get edge directions
    auto d0 = (face->coords[v1] - face->coords[v0])/walberla::real_c(rowsize - 1);
    auto d1 = (face->coords[v2] - face->coords[v1])/walberla::real_c(rowsize - 1);
    auto d2 = (face->coords[v0] - face->coords[v2])/walberla::real_c(rowsize - 1);

    // compute edge lengths
    real_t d0Length = d0.norm();
    real_t d1Length = d1.norm();
    real_t d2Length = d2.norm();

    // compute normals
    auto n_0 = 1.0*d0.normal2D()/d0Length;
    auto n_1 = 1.0*d1.normal2D()/d1Length;
    auto n_2 = 1.0*d2.normal2D()/d2Length;

    for (uint_t i = 1; i < rowsize - 2; ++i) {
      u_0[0] = 0.5*(u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]);
      u_0[1] = 0.5*(v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]);

      u_1[0] = 0.5*(u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]
          + u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_N)]);
      u_1[1] = 0.5*(v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_E)]
          + v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_N)]);

      u_2[0] = 0.5*(u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + u[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_N)]);
      u_2[1] = 0.5*(v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_C)]
          + v[EdgeCoordsVertex::index<Level>(i, EdgeCoordsVertex::VERTEX_N)]);

      // CENTER <-- CELL_GRAY_NE (i)
      // SOUTH  <-- CELL_GRAY_SE (i)
      // WEST   <-- CELL_BLUE_NW (i+1)
      // EAST   <-- CELL_BLUE_NW (i)

      un_0 = d0Length*u_0.dot(n_0);
      un_1 = d1Length*u_1.dot(n_1);
      un_2 = d2Length*u_2.dot(n_2);

      if (un_0 >= 0) {
        c_up_0 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE)];
      } else {
        c_up_0 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE)];
      }

      if (un_1 >= 0) {
        c_up_1 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE)];
      } else {
        c_up_1 = src[indexDGFaceFromVertex<Level>(i+1, stencilDirection::CELL_BLUE_NW)];
      }

      if (un_2 >= 0) {
        c_up_2 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE)];
      } else {
        c_up_2 = src[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_BLUE_NW)];
      }

      tmp = un_0*c_up_0 + un_1*c_up_1 + un_2*c_up_2;
      tmp *= faceAreaInv;

      if (updateType==Replace) {
        dst[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE)] = tmp;
      } else if (updateType==Add) {
        dst[indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE)] += tmp;
      }
    }
  }

}

SPECIALIZE_WITH_VALUETYPE( void, upwindTmpl, upwind )

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


  Face* face0 = storage->getFace(edge.neighborFaces()[0].getID());
  uint_t outerVertexOnFace0 = face0->vertex_index(face0->get_vertex_opposite_to_edge(edge.getID()));


  Point3D face0d0 = dx;
  Point3D face0d2 = (face0->getCoordinates()[outerVertexOnFace0] - x0) /
                    (walberla::real_c(rowsize - 1));


  // gray south cells
  x = x0 + 1.0/3.0 * (face0d0 + face0d2) + dx;
  for (size_t i = 1; i < rowsize - 2; ++i) {
    edgeMemory[DGEdge::indexDGFaceFromVertex< Level >(i,stencilDirection::CELL_GRAY_SE)] = expr(x);
    x += dx;
  }

  if(edge.getNumNeighborFaces() == 2) {
    // gray north cells
    Face *face1 = storage->getFace(edge.neighborFaces()[1].getID());
    uint_t outerVertexOnFace1 = face1->vertex_index(face1->get_vertex_opposite_to_edge(edge.getID()));

    Point3D face1d0 = dx;
    Point3D face1d2 = (face1->getCoordinates()[outerVertexOnFace1] - x0) /
                      (walberla::real_c(rowsize - 1));

    x = x0 + 1.0 / 3.0 * (face1d0 + face1d2) + dx;
    for (size_t i = 1; i < rowsize - 2; ++i) {
      edgeMemory[DGEdge::indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_NE)] = expr(x);
      x += dx;
    }
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
  for (size_t i = 1; i < rowsize - 2; ++i) {
    uint_t cellIndex = DGEdge::indexDGFaceFromVertex<Level>(i, stencilDirection::CELL_GRAY_SE);
    ValueType tmp = scalars[0]*edge.getData(srcIds[0])->getPointer( Level )[cellIndex];
    for (uint_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( Level )[cellIndex];
    }
    dst[cellIndex] = tmp;
  }

  for (size_t i = 1; i < rowsize - 2; ++i) {
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