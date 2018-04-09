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

  Point3D x0({0,0,0}), x;
  real_t h = 1.0 / (walberla::real_c(rowsize - 1));

  form.geometryMap = face.getGeometryMap();

  ValueType tmp;

  Point3D dirS({0, -h, 0});
  Point3D dirSE({h, -h, 0});
  Point3D dirE({h, 0, 0});
  Point3D dirW({-h, 0, 0});
  Point3D dirNW({-h, h, 0});
  Point3D dirN({0, h, 0});

  std::vector<real_t> opr_data(7);

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    x = x0;
    x[0] = h;
    x[1] = j * h;

    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      std::fill(opr_data.begin(), opr_data.end(), 0.0);

      assembleLocalStencil<P1Form>(form, {x0, x0 + dirW, x0 + dirS}, P1Elements::FaceVertexDoF::elementSW, opr_data);
      assembleLocalStencil<P1Form>(form, {x0, x0 + dirS, x0 + dirSE}, P1Elements::FaceVertexDoF::elementS, opr_data);
      assembleLocalStencil<P1Form>(form, {x0, x0 + dirSE, x0 + dirE}, P1Elements::FaceVertexDoF::elementSE, opr_data);
      assembleLocalStencil<P1Form>(form, {x0, x0 + dirE, x0 + dirN}, P1Elements::FaceVertexDoF::elementNE, opr_data);
      assembleLocalStencil<P1Form>(form, {x0, x0 + dirN, x0 + dirNW}, P1Elements::FaceVertexDoF::elementN, opr_data);
      assembleLocalStencil<P1Form>(form, {x0, x0 + dirNW, x0 + dirW}, P1Elements::FaceVertexDoF::elementNW, opr_data);

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

      x[0] += h;
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

  Point3D dS_se = (faceS->coords[e_south] - faceS->coords[s_south]) * h;
  Point3D dS_so = (faceS->coords[o_south] - faceS->coords[s_south]) * h;
  Point3D dS_oe = (faceS->coords[e_south] - faceS->coords[o_south]) * h;

  Point3D dirS_S = -1.0 * dS_oe;
  Point3D dirS_E = dS_se;
  Point3D dirS_SE = dS_so;
  Point3D dirS_W = -1.0 * dS_se;

  Point3D xS_0 = faceS->coords[s_south] + dS_se;

  uint_t s_north, e_north, o_north;
  Point3D dN_se;
  Point3D dN_so;
  Point3D dN_oe;
  Point3D dirN_E;
  Point3D dirN_W;
  Point3D dirN_NW;
  Point3D dirN_N;
  Point3D xN_0;

  if (edge.getNumNeighborFaces() == 2) {
    faceN = storage->getFace(edge.neighborFaces()[1]);
    s_north = faceN->vertex_index(edge.neighborVertices()[0]);
    e_north = faceN->vertex_index(edge.neighborVertices()[1]);
    o_north = faceN->vertex_index(faceN->get_vertex_opposite_to_edge(edge.getID()));

    dN_se = (faceN->coords[e_north] - faceN->coords[s_north]) * h;
    dN_so = (faceN->coords[o_north] - faceN->coords[s_north]) * h;
    dN_oe = (faceN->coords[e_north] - faceN->coords[o_north]) * h;

    dirN_E = dN_se;
    dirN_W = -1.0 * dirN_E;
    dirN_N = dN_so;
    dirN_NW = -1.0 * dN_oe;

    xN_0 = faceN->coords[s_north] + dN_se;
  }

  std::vector<real_t> opr_data(7);

  for (size_t i = 1; i < rowsize - 1; ++i) {

    std::fill(opr_data.begin(), opr_data.end(), 0.0);

    // assemble south
    form.geometryMap = faceS->getGeometryMap();
    assembleLocalStencil<P1Form>(form, {xS_0, xS_0 + dirS_W, xS_0 + dirS_S}, P1Elements::FaceVertexDoF::elementSW, opr_data);
    assembleLocalStencil<P1Form>(form, {xS_0, xS_0 + dirS_S, xS_0 + dirS_SE}, P1Elements::FaceVertexDoF::elementS, opr_data);
    assembleLocalStencil<P1Form>(form, {xS_0, xS_0 + dirS_SE, xS_0 + dirS_E}, P1Elements::FaceVertexDoF::elementSE, opr_data);

    if (edge.getNumNeighborFaces() == 2) {
      form.geometryMap = faceN->getGeometryMap();
      assembleLocalStencil<P1Form>(form, {xN_0, xN_0 + dirN_E, xN_0 + dirN_N}, P1Elements::FaceVertexDoF::elementNE, opr_data);
      assembleLocalStencil<P1Form>(form, {xN_0, xN_0 + dirN_N, xN_0 + dirN_NW}, P1Elements::FaceVertexDoF::elementN, opr_data);
      assembleLocalStencil<P1Form>(form, {xN_0, xN_0 + dirN_NW, xN_0 + dirN_W}, P1Elements::FaceVertexDoF::elementNW, opr_data);
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

    xS_0 += dS_se;
    xN_0 += dN_se;
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

    Point3D x0 = face->coords[v_i];
    d0 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[0])->get_opposite_vertex(vertex.getID()))] - x0) * h;
    d2 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[1])->get_opposite_vertex(vertex.getID()))] - x0) * h;

    Point3D matrixRow;
    form.integrate({{x0, x0 + d0, x0 + d2}}, matrixRow);

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
