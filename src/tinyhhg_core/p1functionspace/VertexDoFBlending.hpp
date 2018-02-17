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
namespace blending {
namespace macroface {

template<typename ValueType, uint_t Level>
inline void assembleLocalStencil(uint_t i, uint_t j, const Matrix3r& localMatrix,
                                 double* opr_data,
                                 real_t coeff,
                                 const std::array<stencilDirection,3>& vertices,
                                 const std::array<uint_t,3>& idx)
{
  opr_data[vertexdof::stencilIndexFromVertex(vertices[0])] += coeff * localMatrix(idx[0],idx[0]);
  opr_data[vertexdof::stencilIndexFromVertex(vertices[1])] += coeff * localMatrix(idx[0],idx[1]);
  opr_data[vertexdof::stencilIndexFromVertex(vertices[2])] += coeff * localMatrix(idx[0],idx[2]);
}

template<typename ValueType, uint_t Level>
inline void applyBlendingTmpl(Face &face,
                                 const std::vector<PrimitiveDataID<FaceP1LocalMatrixMemory, Face>> &operatorIds,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Face> &dstId,
                                 UpdateType update) {
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);

  Point3D x, x0;
  x0 = face.coords[0];
  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  Matrix2r mappingTensor;

  std::vector<FaceP1LocalMatrixMemory *> localMatricesVector;
  for (auto operatorId : operatorIds) {
    localMatricesVector.push_back(face.getData(operatorId));
  }

  ValueType tmp;

  std::array<SD, 3> triangleBlueSW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S};
  std::array<SD, 3> triangleGrayS = {SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE};
  std::array<SD, 3> triangleBlueSE = {SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E};
  std::array<SD, 3> triangleGrayNW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_NW};
  std::array<SD, 3> triangleBlueN = {SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N};
  std::array<SD, 3> triangleGrayNE = {SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_E};

  Point3D offsetSW = -1.0/3.0 * d0 - 1.0/3.0 * d2;
  Point3D offsetS = 1.0/3.0 * d0 - 2.0/3.0 * d2;
  Point3D offsetSE = 2.0/3.0 * d0 - 1.0/3.0 * d2;

  Point3D offsetNE = 1.0/3.0 * d0 + 1.0/3.0 * d2;
  Point3D offsetN = -1.0/3.0 * d0 + 2.0/3.0 * d2;
  Point3D offsetNW = -2.0/3.0 * d0 + 1.0/3.0 * d2;

  real_t* coeffs[] = { &mappingTensor(0,0), &mappingTensor(0,1), &mappingTensor(1,1) };
  std::vector<real_t> opr_data(vertexDoFMacroFaceStencilMemorySize(Level, 0));

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    x = x0;
    x += real_c(j)*d2 + d0;

    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      std::fill(opr_data.begin(), opr_data.end(), 0.0);

      face.blendingMap->evalTensorCoeff(x + offsetS, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayS,
                                               {2, 0, 1});
      }

      face.blendingMap->evalTensorCoeff(x + offsetSE, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleBlueSE,
                                               {1, 2, 0});
      }

      face.blendingMap->evalTensorCoeff(x + offsetSW, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleBlueSW,
                                               {0, 1, 2});
      }

      face.blendingMap->evalTensorCoeff(x + offsetNW, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayNW,
                                               {1, 0, 2});
      }

      face.blendingMap->evalTensorCoeff(x + offsetN, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleBlueN,
                                               {2, 1, 0});
      }

      face.blendingMap->evalTensorCoeff(x + offsetNE, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayNE,
                                               {0, 2, 1});
      }

      auto stencil = PointND<real_t,7>(opr_data.data());

      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      }

      tmp += opr_data[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[0])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[1])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[2])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[3])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[4])];
      tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex<Level>(i, j, vertexdof::macroface::neighborsWithoutCenter[5])];

      dst[vertexdof::macroface::indexFromVertex<Level>(i, j, SD::VERTEX_C)] = tmp;

      x += d0;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyBlendingTmpl, applyBlending)

} // macroface

namespace macroedge {

template<typename ValueType, uint_t Level>
inline void assembleLocalStencil(uint_t pos, const Matrix3r& localMatrix,
                                 real_t* opr_data,
                                 real_t coeff,
                                 const std::array< stencilDirection, 3 >& vertices,
                                 const std::array<uint_t,3>& idx)
{
  opr_data[vertexdof::stencilIndexFromVertex( vertices[0] )] += coeff * localMatrix(idx[0],idx[0]);
  opr_data[vertexdof::stencilIndexFromVertex( vertices[1] )] += coeff * localMatrix(idx[0],idx[1]);
  opr_data[vertexdof::stencilIndexFromVertex( vertices[2] )] += coeff * localMatrix(idx[0],idx[2]);
}

template< typename ValueType, uint_t Level >
inline void applyBlendingTmpl(Edge &edge,
                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                 const std::vector<PrimitiveDataID<EdgeP1LocalMatrixMemory, Edge>> &operatorIds,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                                 UpdateType update) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  Point3D x = edge.getCoordinates()[0];
  Point3D dx = edge.getDirection()/(real_t) (rowsize - 1);

  x += dx;

  std::vector<EdgeP1LocalMatrixMemory*> localMatricesVector;
  for(auto operatorId : operatorIds) {
    localMatricesVector.push_back(edge.getData(operatorId));
  }

  ValueType tmp;

  std::array< stencilDirection, 3 > triangleGraySW = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_W,  stencilDirection::VERTEX_S  };
  std::array< stencilDirection, 3 > triangleBlueS  = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_S,  stencilDirection::VERTEX_SE };
  std::array< stencilDirection, 3 > triangleGraySE = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_E  };
  std::array< stencilDirection, 3 > triangleGrayNW = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_W,  stencilDirection::VERTEX_NW };
  std::array< stencilDirection, 3 > triangleBlueN  = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_N  };
  std::array< stencilDirection, 3 > triangleGrayNE = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_N,  stencilDirection::VERTEX_E  };

  Face* faceS = storage->getFace(edge.neighborFaces()[0]);
  Face* faceN;
  uint_t s_south = faceS->vertex_index(edge.neighborVertices()[0]);
  uint_t e_south = faceS->vertex_index(edge.neighborVertices()[1]);
  uint_t o_south = faceS->vertex_index(faceS->get_vertex_opposite_to_edge(edge.getID()));

  Point3D d0S = (faceS->coords[e_south] - faceS->coords[s_south])/(walberla::real_c(rowsize - 1));
  Point3D d2S = (faceS->coords[o_south] - faceS->coords[e_south])/(walberla::real_c(rowsize - 1));

  Point3D offsetSW = -1.0/3.0 * d0S + 1.0/3.0 * d2S;
  Point3D offsetS = 1.0/3.0 * d0S + 2.0/3.0 * d2S;
  Point3D offsetSE = 2.0/3.0 * d0S + 1.0/3.0 * d2S;

  uint_t s_north, e_north, o_north;
  Point3D d0N;
  Point3D d2N;
  Point3D offsetNE;
  Point3D offsetN;
  Point3D offsetNW;

  if (edge.getNumNeighborFaces() == 2) {
    faceN = storage->getFace(edge.neighborFaces()[1]);
    s_north = faceN->vertex_index(edge.neighborVertices()[0]);
    e_north = faceN->vertex_index(edge.neighborVertices()[1]);
    o_north = faceN->vertex_index(faceN->get_vertex_opposite_to_edge(edge.getID()));

    d0N = (faceN->coords[e_north] - faceN->coords[s_north])/(walberla::real_c(rowsize - 1));
    d2N = (faceN->coords[o_north] - faceN->coords[e_north])/(walberla::real_c(rowsize - 1));

    offsetNE = 1.0/3.0 * d0N + 1.0/3.0 * d2N;
    offsetN = -1.0/3.0 * d0N + 2.0/3.0 * d2N;
    offsetNW = -2.0/3.0 * d0N + 1.0/3.0 * d2N;
  }

  Matrix2r mappingTensor;
  real_t* coeffs[] = { &mappingTensor(0,0), &mappingTensor(0,1), &mappingTensor(1,1) };
  std::vector<real_t> opr_data(7);

  for (size_t i = 1; i < rowsize - 1; ++i) {

    std::fill(opr_data.begin(), opr_data.end(), 0.0);

    faceS->blendingMap->evalTensorCoeff(x + offsetSW, mappingTensor);
    for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
      assembleLocalStencil<ValueType, Level>(i,
                                             localMatricesVector[coeffIdx]->getGrayMatrix(Level, 0),
                                             opr_data.data(),
                                             *coeffs[coeffIdx],
                                             triangleGraySW,
                                             {e_south, s_south, o_south});
    }

    faceS->blendingMap->evalTensorCoeff(x + offsetS, mappingTensor);
    for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
      assembleLocalStencil<ValueType, Level>(i,
                                             localMatricesVector[coeffIdx]->getBlueMatrix(Level, 0),
                                             opr_data.data(),
                                             *coeffs[coeffIdx],
                                             triangleBlueS,
                                             {o_south, e_south, s_south});
    }

    faceS->blendingMap->evalTensorCoeff(x + offsetSE, mappingTensor);
    for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
      assembleLocalStencil<ValueType, Level>(i,
                                             localMatricesVector[coeffIdx]->getGrayMatrix(Level, 0),
                                             opr_data.data(),
                                             *coeffs[coeffIdx],
                                             triangleGraySE,
                                             {s_south, o_south, e_south});
    }

    if (edge.getNumNeighborFaces() == 2) {
      faceN->blendingMap->evalTensorCoeff(x + offsetNW, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level, 1),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayNW,
                                               {e_north, s_north, o_north});
      }

      faceN->blendingMap->evalTensorCoeff(x + offsetN, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               localMatricesVector[coeffIdx]->getBlueMatrix(Level, 1),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleBlueN,
                                               {o_north, e_north, s_north});
      }

      faceN->blendingMap->evalTensorCoeff(x + offsetNE, mappingTensor);
      for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {
        assembleLocalStencil<ValueType, Level>(i,
                                               localMatricesVector[coeffIdx]->getGrayMatrix(Level, 1),
                                               opr_data.data(),
                                               *coeffs[coeffIdx],
                                               triangleGrayNE,
                                               {s_north, o_north, e_north});
      }
    }

    tmp = opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];

    // neighbors on edge
    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      tmp += opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>( i, neighbor ) ];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      tmp += opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>( i, neighbor ) ];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        tmp += opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>( i, neighbor ) ];
      }
    }

    if (update == Replace) {
      dst[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ] = tmp;
    } else if (update == Add) {
      dst[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ] += tmp;
    }

    x += dx;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, applyBlendingTmpl, applyBlending )

} // macroedge

namespace macrovertex {

template< typename ValueType >
inline void applyBlending(uint_t level, Vertex &vertex,
                              const std::shared_ptr< PrimitiveStorage >& storage,
                              const std::vector<PrimitiveDataID<VertexP1LocalMatrixMemory, Vertex>> &operatorIds,
                              const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId,
                              const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                              UpdateType update) {

  auto src = vertex.getData(srcId)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );

  std::vector<VertexP1LocalMatrixMemory*> localMatricesVector;
  for(auto operatorId : operatorIds) {
    localMatricesVector.push_back(vertex.getData(operatorId));
  }

  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);

  Matrix2r mappingTensor;
  real_t* coeffs[] = { &mappingTensor(0,0), &mappingTensor(0,1), &mappingTensor(1,1) };
  std::vector<real_t> opr_data(1 + vertex.getNumNeighborEdges());
  std::fill(opr_data.begin(), opr_data.end(), 0.0);
  Point3D x;
  Point3D d0;
  Point3D d2;

  uint_t neighborId = 0;
  for (auto& faceId : vertex.neighborFaces()) {

    Face* face = storage->getFace(faceId);

    uint_t v_i = face->vertex_index(vertex.getID());
    std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

    d0 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[0])->get_opposite_vertex(vertex.getID()))] - face->coords[v_i])/(walberla::real_c(rowsize - 1));
    d2 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[1])->get_opposite_vertex(vertex.getID()))] - face->coords[v_i])/(walberla::real_c(rowsize - 1));

    x = face->coords[v_i] + 1.0/3.0 * d0 + 1.0/3.0 * d2;
    face->blendingMap->evalTensorCoeff(x, mappingTensor);

    for (uint_t coeffIdx = 0; coeffIdx < 3; ++coeffIdx) {

      Matrix3r& local_stiffness = localMatricesVector[coeffIdx]->getGrayMatrix(level, neighborId);

      // iterate over adjacent edges
      for (auto &edgeId : adj_edges) {
        uint_t edge_idx = vertex.edge_index(edgeId) + 1;
        Edge *edge = storage->getEdge(edgeId);
        PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

        uint_t v_j = face->vertex_index(vertex_j);

        opr_data[edge_idx] += *coeffs[coeffIdx] * local_stiffness(v_i, v_j);
      }

      // add contribution of center vertex
      opr_data[0] += *coeffs[coeffIdx] * local_stiffness(v_i, v_i);
    }

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