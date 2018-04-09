#pragma once

#include "core/debug/all.h"

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFaceIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/polynomial/PolynomialEvaluator.hpp"

namespace hhg {
namespace vertexdof {
namespace blendingnew {

template<class P1Form>
inline void assembleLocalStencil(const P1Form& form,
                                 const std::array<Point3D,3>& coords,
                                 const std::array<stencilDirection,3>& directions,
                                 std::vector<real_t>& opr_data)
{
  Point3D matrixRow;

  form.integrate(coords, matrixRow);

  opr_data[vertexdof::stencilIndexFromVertex(directions[0])] += matrixRow[0];
  opr_data[vertexdof::stencilIndexFromVertex(directions[1])] += matrixRow[1];
  opr_data[vertexdof::stencilIndexFromVertex(directions[2])] += matrixRow[2];
}

namespace macroface {

template<typename ValueType, class P1Form>
inline void applyBlending(uint_t Level, Face &face,
                          P1Form& form,
                          const PrimitiveDataID<FunctionMemory<ValueType>, Face> &srcId,
                          const PrimitiveDataID<FunctionMemory<ValueType>, Face> &dstId,
                          UpdateType update) {
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);

  Point3D x0(face.coords[0]), x;
  real_t h = 1.0 / (walberla::real_c(rowsize - 1));

  Point3D d0 = h * (face.coords[1] - face.coords[0]);
  Point3D d2 = h * (face.coords[2] - face.coords[0]);

  form.geometryMap = face.getGeometryMap();

  ValueType tmp;

  Point3D dirS = -1.0 * d2;
  Point3D dirSE = d0 - 1.0 * d2;
  Point3D dirE = d0;
  Point3D dirW = -1.0 * d0;
  Point3D dirNW = -1.0 * d0 + d2;
  Point3D dirN = d2;

  std::vector<real_t> opr_data(7);

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    x = x0;
    x += walberla::real_c(j)*d2 + d0;

    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      std::fill(opr_data.begin(), opr_data.end(), 0.0);

      assembleLocalStencil<P1Form>(form, {x, x + dirW, x + dirS}, P1Elements::FaceVertexDoF::elementSW, opr_data);
      assembleLocalStencil<P1Form>(form, {x, x + dirS, x + dirSE}, P1Elements::FaceVertexDoF::elementS, opr_data);
      assembleLocalStencil<P1Form>(form, {x, x + dirSE, x + dirE}, P1Elements::FaceVertexDoF::elementSE, opr_data);
      assembleLocalStencil<P1Form>(form, {x, x + dirE, x + dirN}, P1Elements::FaceVertexDoF::elementNE, opr_data);
      assembleLocalStencil<P1Form>(form, {x, x + dirN, x + dirNW}, P1Elements::FaceVertexDoF::elementN, opr_data);
      assembleLocalStencil<P1Form>(form, {x, x + dirNW, x + dirW}, P1Elements::FaceVertexDoF::elementNW, opr_data);

//      PointND<real_t, 7> test(opr_data.data());
//      WALBERLA_LOG_INFO("stencil = " << test);

      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex(Level, i, j, SD::VERTEX_C)];
      }

      tmp += opr_data[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex(Level, i, j, SD::VERTEX_C)];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex(Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

      dst[vertexdof::macroface::indexFromVertex(Level, i, j, SD::VERTEX_C)] = tmp;

      x += d0;
    }
    --inner_rowsize;
  }
}

} // macroface

namespace macroedge {

template< typename ValueType, class P1Form>
inline void applyBlending(uint_t Level, Edge &edge,
                          P1Form& form,
                          const std::shared_ptr< PrimitiveStorage >& storage,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                          UpdateType update) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  ValueType tmp;

  Face* faceS = storage->getFace(edge.neighborFaces()[0]);
  Face* faceN;
  uint_t s_south = faceS->vertex_index(edge.neighborVertices()[0]);
  uint_t e_south = faceS->vertex_index(edge.neighborVertices()[1]);
  uint_t o_south = faceS->vertex_index(faceS->get_vertex_opposite_to_edge(edge.getID()));

  real_t h = 1.0 / (walberla::real_c(rowsize - 1));

  Point3D dS_se = h * (faceS->coords[e_south] - faceS->coords[s_south]);
  Point3D dS_so = h * (faceS->coords[o_south] - faceS->coords[s_south]);
  Point3D dS_oe = h * (faceS->coords[e_south] - faceS->coords[o_south]);

  Point3D dir_S = -1.0 * dS_oe;
  Point3D dir_E = dS_se;
  Point3D dir_SE = dS_so;
  Point3D dir_W = -1.0 * dS_se;

  Point3D x = edge.getCoordinates()[0];
  Point3D dx = h * edge.getDirection();
  x += dx;

  uint_t s_north, e_north, o_north;
  Point3D dir_N;
  Point3D dir_NW;

  if (edge.getNumNeighborFaces() == 2) {
    faceN = storage->getFace(edge.neighborFaces()[1]);
    s_north = faceN->vertex_index(edge.neighborVertices()[0]);
    e_north = faceN->vertex_index(edge.neighborVertices()[1]);
    o_north = faceN->vertex_index(faceN->get_vertex_opposite_to_edge(edge.getID()));

    Point3D dN_so = h * (faceN->coords[o_north] - faceN->coords[s_north]);
    Point3D dN_oe = h * (faceN->coords[e_north] - faceN->coords[o_north]);

    dir_N = dN_so;
    dir_NW = -1.0 * dN_oe;
  }

  std::vector<real_t> opr_data(7);

  for (size_t i = 1; i < rowsize - 1; ++i) {

    std::fill(opr_data.begin(), opr_data.end(), 0.0);

    // assemble south
    form.geometryMap = faceS->getGeometryMap();
    assembleLocalStencil<P1Form>(form, {x, x + dir_W, x + dir_S}, P1Elements::FaceVertexDoF::elementSW, opr_data);
    assembleLocalStencil<P1Form>(form, {x, x + dir_S, x + dir_SE}, P1Elements::FaceVertexDoF::elementS, opr_data);
    assembleLocalStencil<P1Form>(form, {x, x + dir_SE, x + dir_E}, P1Elements::FaceVertexDoF::elementSE, opr_data);

    if (edge.getNumNeighborFaces() == 2) {
      form.geometryMap = faceN->getGeometryMap();
      assembleLocalStencil<P1Form>(form, {x, x + dir_E, x + dir_N}, P1Elements::FaceVertexDoF::elementNE, opr_data);
      assembleLocalStencil<P1Form>(form, {x, x + dir_N, x + dir_NW}, P1Elements::FaceVertexDoF::elementN, opr_data);
      assembleLocalStencil<P1Form>(form, {x, x + dir_NW, x + dir_W}, P1Elements::FaceVertexDoF::elementNW, opr_data);
    }

    tmp = opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] * src[ vertexdof::macroedge::indexFromVertex(Level, i, stencilDirection::VERTEX_C ) ];

    // neighbors on edge
    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      tmp += opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex(Level, i, neighbor ) ];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      tmp += opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex(Level, i, neighbor ) ];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        tmp += opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex(Level, i, neighbor ) ];
      }
    }

    if (update == Replace) {
      dst[ vertexdof::macroedge::indexFromVertex(Level, i, stencilDirection::VERTEX_C ) ] = tmp;
    } else if (update == Add) {
      dst[ vertexdof::macroedge::indexFromVertex(Level, i, stencilDirection::VERTEX_C ) ] += tmp;
    }

    x += dx;
  }
}
}

namespace macrovertex {

template< typename ValueType, class P1Form >
inline void applyBlending(uint_t level, Vertex &vertex,
                          P1Form& form,
                          const std::shared_ptr< PrimitiveStorage >& storage,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                          UpdateType update) {

  auto src = vertex.getData(srcId)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );

  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);

  std::vector<real_t> opr_data(1 + vertex.getNumNeighborEdges());
  std::fill(opr_data.begin(), opr_data.end(), 0.0);
  Point3D x;
  Point3D d0;
  Point3D d2;

  real_t h = 1.0 / (walberla::real_c(rowsize - 1));

  uint_t neighborId = 0;
  for (auto& faceId : vertex.neighborFaces()) {

    Face* face = storage->getFace(faceId);
    form.geometryMap = face->getGeometryMap();

    uint_t v_i = face->vertex_index(vertex.getID());
    std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

    x = face->coords[v_i];
    d0 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[0])->get_opposite_vertex(vertex.getID()))] - x) * h;
    d2 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[1])->get_opposite_vertex(vertex.getID()))] - x) * h;

    Point3D matrixRow;
    form.integrate({{x, x + d0, x + d2}}, matrixRow);

    uint_t i = 1;
    // iterate over adjacent edges
    for (auto &edgeId : adj_edges) {
      uint_t edge_idx = vertex.edge_index(edgeId) + 1;
      Edge *edge = storage->getEdge(edgeId);
      PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

      uint_t v_j = face->vertex_index(vertex_j);

      opr_data[edge_idx] += matrixRow[i];
      i += 1;
    }

    // add contribution of center vertex
    opr_data[0] += matrixRow[0];

    ++neighborId;
  }

  if (update==Replace) {
    dst[0] = opr_data[0]*src[0];
  } else if (update==Add) {
    dst[0] += opr_data[0]*src[0];
  }

  for (size_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
    dst[0] += opr_data[i + 1]*src[i + 1];
  }
}

} // macrovertex

} // blending
} // vertexdof
} // hhg
