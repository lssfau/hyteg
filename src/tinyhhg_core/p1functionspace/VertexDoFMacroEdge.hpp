#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/dgfunctionspace/DGEdgeIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

#include "core/DataTypes.h"

namespace hhg {
namespace vertexdof {
namespace macroedge {

template<typename ValueType, uint_t Level>
inline ValueType assembleLocal(uint_t pos, const Matrix3r& localMatrix,
                               double* src,
                               double* coeff,
                               const std::array< stencilDirection, 3 >& vertices,
                               const std::array<uint_t,3>& idx)
{

  ValueType meanCoeff = 1.0/3.0 * (coeff[ vertexdof::macroedge::indexFromVertex<Level>( pos, vertices[ 0 ] ) ]
                                 + coeff[ vertexdof::macroedge::indexFromVertex<Level>( pos, vertices[ 1 ] ) ]
                                 + coeff[ vertexdof::macroedge::indexFromVertex<Level>( pos, vertices[ 2 ] ) ]);

  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[ vertexdof::macroedge::indexFromVertex<Level>(pos, vertices[0]) ]
         + localMatrix(idx[0],idx[1]) * src[ vertexdof::macroedge::indexFromVertex<Level>(pos, vertices[1]) ]
         + localMatrix(idx[0],idx[2]) * src[ vertexdof::macroedge::indexFromVertex<Level>(pos, vertices[2]) ];
  return meanCoeff * tmp;
}

template<uint_t Level>
inline void fillLocalCoords(uint_t i, const std::array< stencilDirection, 3>& element, const std::array<real_t*, 2>& coords, real_t localCoords[6] )
{
  localCoords[0] = coords[0][ vertexdof::macroedge::indexFromVertex<Level>(i, element[0]) ];
  localCoords[1] = coords[1][ vertexdof::macroedge::indexFromVertex<Level>(i, element[0]) ];
  localCoords[2] = coords[0][ vertexdof::macroedge::indexFromVertex<Level>(i, element[1]) ];
  localCoords[3] = coords[1][ vertexdof::macroedge::indexFromVertex<Level>(i, element[1]) ];
  localCoords[4] = coords[0][ vertexdof::macroedge::indexFromVertex<Level>(i, element[2]) ];
  localCoords[5] = coords[1][ vertexdof::macroedge::indexFromVertex<Level>(i, element[2]) ];
}

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Edge &edge,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Edge> &edgeMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &srcIds,
                            std::function<ValueType(const hhg::Point3D &, const std::vector<ValueType>&)> &expr) {

  FunctionMemory< ValueType > *edgeMemory = edge.getData(edgeMemoryId);

  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(edge.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVector(srcIds.size());

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x = edge.getCoordinates()[0];
  Point3D dx = edge.getDirection()/(real_t) (rowsize - 1);

  x += dx;

  for (size_t i = 1; i < rowsize - 1; ++i) {

    for (uint_t k = 0; k < srcPtr.size(); ++k) {
      srcVector[k] = srcPtr[k][ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];
    }

    edgeMemory->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ] = expr(x, srcVector);
    x += dx;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate )

template< typename ValueType, uint_t Level >
inline void assignTmpl(Edge &edge,
                   const std::vector<ValueType> &scalars,
                   const std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge>> &srcIds,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge> &dstId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    ValueType tmp = scalars[0]*edge.getData(srcIds[0])->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];

    for (size_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];
    }

    edge.getData(dstId)->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, assignTmpl, assign )

template< typename ValueType, uint_t Level >
inline void addTmpl(Edge &edge,
                const std::vector<ValueType> &scalars,
                const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &srcIds,
                const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    ValueType tmp = 0.0;

    for (size_t k = 0; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];
    }

    edge.getData(dstId)->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ] += tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add )

template< typename ValueType, uint_t Level >
inline real_t dotTmpl(Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &lhsMemoryId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsMemoryId) {

  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    sp += edge.getData(lhsMemoryId)->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ]
        * edge.getData(rhsMemoryId)->getPointer( Level )[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];
  }

  return sp;
}

SPECIALIZE_WITH_VALUETYPE( real_t, dotTmpl, dot )

template< typename ValueType, uint_t Level >
inline void applyTmpl(Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId, UpdateType update) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto opr_data = edge.getData(operatorId)->getPointer( Level );
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  ValueType tmp;

  for (size_t i = 1; i < rowsize - 1; ++i) {

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
  }
}

SPECIALIZE_WITH_VALUETYPE( void, applyTmpl, apply )

template< typename ValueType, uint_t Level >
inline void applyCoefficientTmpl(Edge &edge,
                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                 const PrimitiveDataID<EdgeP1LocalMatrixMemory, Edge> &operatorId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &coeffId,
                                 UpdateType update) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto localMatrices = edge.getData(operatorId);
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto coeff = edge.getData(coeffId)->getPointer( Level );

  ValueType tmp;

  std::array< stencilDirection, 3 > triangleGraySW = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_W,  stencilDirection::VERTEX_S  };
  std::array< stencilDirection, 3 > triangleBlueS  = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_S,  stencilDirection::VERTEX_SE };
  std::array< stencilDirection, 3 > triangleGraySE = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_E  };
  std::array< stencilDirection, 3 > triangleGrayNW = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_W,  stencilDirection::VERTEX_NW };
  std::array< stencilDirection, 3 > triangleBlueN  = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_N  };
  std::array< stencilDirection, 3 > triangleGrayNE = { stencilDirection::VERTEX_C, stencilDirection::VERTEX_N,  stencilDirection::VERTEX_E  };

  Face* face = storage->getFace(edge.neighborFaces()[0]);
  uint_t s_south = face->vertex_index(edge.neighborVertices()[0]);
  uint_t e_south = face->vertex_index(edge.neighborVertices()[1]);
  uint_t o_south = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

  uint_t s_north, e_north, o_north;

  if (edge.getNumNeighborFaces() == 2) {
    face = storage->getFace(edge.neighborFaces()[1]);
    s_north = face->vertex_index(edge.neighborVertices()[0]);
    e_north = face->vertex_index(edge.neighborVertices()[1]);
    o_north = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));
  }

  for (size_t i = 1; i < rowsize - 1; ++i) {

    if (update == Replace) {
      tmp = ValueType(0);
    }
    else {
      tmp = dst[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];
    }

    tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 0), src, coeff, triangleGraySW, {e_south, s_south, o_south});
    tmp += assembleLocal<ValueType, Level>(i, localMatrices->getBlueMatrix(Level, 0), src, coeff, triangleBlueS, {o_south, e_south, s_south});
    tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 0), src, coeff, triangleGraySE, {s_south, o_south, e_south});

    if (edge.getNumNeighborFaces() == 2)
    {
      tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 1), src, coeff, triangleGrayNW, {e_north, s_north, o_north});
      tmp += assembleLocal<ValueType, Level>(i, localMatrices->getBlueMatrix(Level, 1), src, coeff, triangleBlueN, {o_north, e_north, s_north});
      tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 1), src, coeff, triangleGrayNE, {s_north, o_north, e_north});
    }

    dst[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, applyCoefficientTmpl, applyCoefficient )

template< typename ValueType, uint_t Level >
inline void applyElementwiseTmpl(Edge &edge,
                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                 std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                                 std::array<const PrimitiveDataID<FunctionMemory< ValueType >, Edge>, 2> &coordIds,
                                 UpdateType update) {
  using namespace P1Elements;
  typedef std::array< stencilDirection, 3 > Element;
  typedef std::array<uint_t, 3> DoFMap;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = edge.getData(srcId)->getPointer(Level);
  auto dst = edge.getData(dstId)->getPointer(Level);
  std::array<ValueType*, 2> globalCoords{{edge.getData(coordIds[0])->getPointer(Level), edge.getData(coordIds[1])->getPointer(Level)}};

  Face* face = storage->getFace(edge.neighborFaces()[0]);
  uint_t start = face->vertex_index(edge.neighborVertices()[0]);
  uint_t end = face->vertex_index(edge.neighborVertices()[1]);
  uint_t opposite = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

  Element elementSW = {{stencilDirection::VERTEX_C, stencilDirection::VERTEX_W, stencilDirection::VERTEX_S}};
  Element elementS = {{stencilDirection::VERTEX_C, stencilDirection::VERTEX_S, stencilDirection::VERTEX_SE}};
  Element elementSE = {{stencilDirection::VERTEX_C, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_E}};
  Element elementNE = {{stencilDirection::VERTEX_C, stencilDirection::VERTEX_E, stencilDirection::VERTEX_N}};
  Element elementN = {{stencilDirection::VERTEX_C, stencilDirection::VERTEX_N, stencilDirection::VERTEX_NW}};
  Element elementNW = {{stencilDirection::VERTEX_C, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_W}};

  DoFMap dofMapSW = {{end, start, opposite}};
  DoFMap dofMapS = {{opposite, end, start}};
  DoFMap dofMapSE = {{start, opposite, end}};

  DoFMap dofMapNE, dofMapN, dofMapNW;

  if (edge.getNumNeighborFaces() == 2) {
    face = storage->getFace(edge.neighborFaces()[1]);
    start = face->vertex_index(edge.neighborVertices()[0]);
    end = face->vertex_index(edge.neighborVertices()[1]);
    opposite = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

    dofMapNE = {{start, end, opposite}};
    dofMapN = {{opposite, start, end}};
    dofMapNW = {{end, opposite, start}};
  }

  ValueType tmp;
  real_t localCoords[6];
  Matrix3r localStiffness;
  std::vector<real_t> edgeStencil(7);

  for (size_t i = 1; i < rowsize - 1; ++i) {

    std::fill(edgeStencil.begin(), edgeStencil.end(), walberla::real_c(0.0));

    fillLocalCoords<Level>(i, elementSW, globalCoords, localCoords);
    computeElementMatrix(localStiffness, localCoords);
    assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSW ), {{0,1,2}}, localStiffness, edgeStencil);

    fillLocalCoords<Level>(i, elementS, globalCoords, localCoords);
    computeElementMatrix(localStiffness, localCoords);
    assembleP1LocalStencil( convertStencilDirectionsToIndices( elementS ), {{0,1,2}}, localStiffness, edgeStencil);

    fillLocalCoords<Level>(i, elementSE, globalCoords, localCoords);
    computeElementMatrix(localStiffness, localCoords);
    assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSE ), {{0,1,2}}, localStiffness, edgeStencil);

    if (edge.getNumNeighborFaces() == 2)
    {
      fillLocalCoords<Level>(i, elementNE, globalCoords, localCoords);
      computeElementMatrix(localStiffness, localCoords);
      assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNE ), {{0,1,2}}, localStiffness, edgeStencil);

      fillLocalCoords<Level>(i, elementN, globalCoords, localCoords);
      computeElementMatrix(localStiffness, localCoords);
      assembleP1LocalStencil( convertStencilDirectionsToIndices( elementN ), {{0,1,2}}, localStiffness, edgeStencil);

      fillLocalCoords<Level>(i, elementNW, globalCoords, localCoords);
      computeElementMatrix(localStiffness, localCoords);
      assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNW ), {{0,1,2}}, localStiffness, edgeStencil);
    }

//    WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("FACE.id = {}:edgeStencil = {}", edge.getID().getID(), PointND<real_t, 7>(&edgeStencil[0])));

    tmp = edgeStencil[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>( i, stencilDirection::VERTEX_C ) ];

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      tmp += edgeStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>(i, neighbor) ];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      tmp += edgeStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>(i, neighbor) ];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        tmp += edgeStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ] * src[ vertexdof::macroedge::indexFromVertex<Level>(i, neighbor) ];
      }
    }

    if (update == Replace) {
      dst[ vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] = tmp;
    } else if (update == Add) {
      dst[ vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] += tmp;
    }
  }
}

SPECIALIZE_WITH_VALUETYPE(void, applyElementwiseTmpl, applyElementwise)

template< typename ValueType, uint_t Level >
inline void smoothGSTmpl(Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto opr_data = edge.getData(operatorId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto rhs = edge.getData(rhsId)->getPointer( Level );

  for (size_t i = 1; i < rowsize - 1; ++i) {

    dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C) ] = rhs[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)];

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      dst[ vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C) ] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C) ] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        dst[ vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C) ] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
      }
    }

    dst[ vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C) ] /= opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ];
  }
}

SPECIALIZE_WITH_VALUETYPE( void, smoothGSTmpl, smooth_gs )

template< typename ValueType, uint_t Level >
inline void smoothSORTmpl(Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                          ValueType relax) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto opr_data = edge.getData(operatorId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto rhs = edge.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (size_t i = 1; i < rowsize - 1; ++i) {

    tmp = rhs[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)];

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      tmp -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      tmp -= opr_data[vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        tmp -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
      }
    }
    
    dst[ vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C ) ] = (1.0-relax) * dst[ vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C) ]
                                                                                                             + relax * tmp/opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ];
  }
}

SPECIALIZE_WITH_VALUETYPE( void, smoothSORTmpl, smooth_sor )

template< typename ValueType, uint_t Level >
inline void smoothJacTmpl(Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &tmpId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto opr_data = edge.getData(operatorId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto rhs = edge.getData(rhsId)->getPointer( Level );
  auto tmp = edge.getData(tmpId)->getPointer( Level );

  for (size_t i = 1; i < rowsize - 1; ++i) {

    dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] = rhs[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)];

    for (auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * tmp[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
    }

    for (auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * tmp[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for (auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * tmp[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
      }
    }

    dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] /= opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ];
  }
}

SPECIALIZE_WITH_VALUETYPE( void, smoothJacTmpl, smooth_jac )

template< typename ValueType, uint_t SourceLevel >
inline void prolongateTmpl(Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &memoryId) {

  size_t rowsize_c = levelinfo::num_microvertices_per_edge(SourceLevel);

  auto edge_data_f = edge.getData(memoryId)->getPointer( SourceLevel + 1 );
  auto edge_data_c = edge.getData(memoryId)->getPointer( SourceLevel );

  size_t i_c;
  for (i_c = 1; i_c < rowsize_c - 1; ++i_c) {

    edge_data_f[vertexdof::macroedge::indexFromVertex<SourceLevel+1>(2*i_c, stencilDirection::VERTEX_C)] = edge_data_c[vertexdof::macroedge::indexFromVertex<SourceLevel>(i_c, stencilDirection::VERTEX_C)];
    edge_data_f[vertexdof::macroedge::indexFromVertex<SourceLevel+1>(2*i_c-1, stencilDirection::VERTEX_C)] = 0.5*(edge_data_c[vertexdof::macroedge::indexFromVertex<SourceLevel>(i_c-1, stencilDirection::VERTEX_C)]
                                                              + edge_data_c[vertexdof::macroedge::indexFromVertex<SourceLevel>(i_c, stencilDirection::VERTEX_C)]);
  }

  edge_data_f[vertexdof::macroedge::indexFromVertex<SourceLevel+1>(2*i_c-1, stencilDirection::VERTEX_C)] = 0.5*(edge_data_c[vertexdof::macroedge::indexFromVertex<SourceLevel>(i_c-1, stencilDirection::VERTEX_C)]
                                                            + edge_data_c[vertexdof::macroedge::indexFromVertex<SourceLevel>(i_c, stencilDirection::VERTEX_C)]);
}

SPECIALIZE_WITH_VALUETYPE( void, prolongateTmpl, prolongate )

template< typename ValueType, uint_t Level >
inline void prolongateQuadraticTmpl(Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &memoryId) {

  //TODO: rewrite using index function possible? maybe more generalized notion of Operator between different levels
  // is required.

  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(Level);
  size_t i_fine = 1;
  ValueType invtemp = 1/8.;
  const ValueType s1[3] = {3*invtemp, 6*invtemp, -invtemp};
  const ValueType s2[3] = {-invtemp, 6*invtemp, 3*invtemp};

  auto edge_data_f = edge.getData(memoryId)->getPointer( Level + 1 );
  auto edge_data_c = edge.getData(memoryId)->getPointer( Level );
  size_t i_coarse;
  for (i_coarse = 0; i_coarse < rowsize_coarse - 2; ++i_coarse) {
    edge_data_f[i_fine] =
        (s1[0]*edge_data_c[i_coarse] + s1[1]*edge_data_c[i_coarse + 1] + s1[2]*edge_data_c[i_coarse + 2]);
    edge_data_f[i_fine + 1] = edge_data_c[i_coarse + 1];
    i_fine += 2;
  }
  i_coarse--;
  edge_data_f[i_fine] =
      (s2[0]*edge_data_c[i_coarse] + s2[1]*edge_data_c[i_coarse + 1] + s2[2]*edge_data_c[i_coarse + 2]);

}

SPECIALIZE_WITH_VALUETYPE( void, prolongateQuadraticTmpl, prolongateQuadratic )

template< typename ValueType, uint_t Level >
inline void restrictTmpl(Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &memoryId) {

  size_t rowsize_c = levelinfo::num_microvertices_per_edge(Level - 1);

  auto edge_data_f = edge.getData(memoryId)->getPointer( Level );
  auto edge_data_c = edge.getData(memoryId)->getPointer( Level - 1 );

  uint_t i_c;
  for ( i_c = 1; i_c < rowsize_c - 1; ++i_c)
  {
    edge_data_c[vertexdof::macroedge::indexFromVertex<Level-1>(i_c, stencilDirection::VERTEX_C)] = 1.0 * edge_data_f[vertexdof::macroedge::indexFromVertex<Level>(2*i_c, stencilDirection::VERTEX_C)];

    for (auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      edge_data_c[vertexdof::macroedge::indexFromVertex<Level-1>(i_c, stencilDirection::VERTEX_C)] += 0.5 * edge_data_f[vertexdof::macroedge::indexFromVertex<Level>(2*i_c, neighbor)];
    }

    for (auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      edge_data_c[vertexdof::macroedge::indexFromVertex<Level-1>(i_c, stencilDirection::VERTEX_C)] += 0.5 * edge_data_f[vertexdof::macroedge::indexFromVertex<Level>(2*i_c, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for (auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        edge_data_c[vertexdof::macroedge::indexFromVertex<Level-1>(i_c, stencilDirection::VERTEX_C)] += 0.5 * edge_data_f[vertexdof::macroedge::indexFromVertex<Level>(2*i_c, neighbor)];
      }
    }
  }
}

SPECIALIZE_WITH_VALUETYPE( void, restrictTmpl, restrict )

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId, uint_t& num) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    edge.getData(dstId)->getPointer( Level )[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] = walberla::real_c(num++);
  }
}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )

template< typename ValueType, uint_t Level >
inline void integrateDGTmpl(Edge &edge,
                            const std::shared_ptr< PrimitiveStorage >& storage,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsP1Id,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId) {

  typedef stencilDirection sD;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto rhs = edge.getData(rhsId)->getPointer( Level );
  auto rhsP1 = edge.getData(rhsP1Id)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  ValueType tmp;

  Face* face = storage->getFace(edge.neighborFaces()[0]);
  real_t weightedFaceArea0, weightedFaceArea1;

  weightedFaceArea0 = std::pow(4.0, -walberla::real_c(Level)) * face->area / 3.0;

  uint_t s_north, e_north, o_north;

  if (edge.getNumNeighborFaces() == 2) {
    face = storage->getFace(edge.neighborFaces()[1]);
    weightedFaceArea1 = std::pow(4.0, -walberla::real_c(Level)) * face->area / 3.0;
  }

  for (size_t i = 1; i < rowsize - 1; ++i) {

    tmp =  weightedFaceArea0 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_SW)] * (0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_W)]) + 0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_S)]));
    tmp += weightedFaceArea0 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_BLUE_SE)] * (0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_S)]) + 0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_SE)]));
    tmp += weightedFaceArea0 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_SE)] * (0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_SE)]) + 0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_E)]));

    if (edge.getNumNeighborFaces() == 2) {

      tmp += weightedFaceArea1 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_NW)] * (0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_W)]) + 0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_NW)]));
      tmp += weightedFaceArea1 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_BLUE_NW)] * (0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_NW)]) + 0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_N)]));
      tmp += weightedFaceArea1 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_NE)] * (0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_N)]) + 0.5 * 0.5 * (rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] + rhsP1[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_E)]));
    }

    dst[vertexdof::macroedge::indexFromVertex<Level>(i, sD::VERTEX_C)] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, integrateDGTmpl, integrateDG )

#ifdef HHG_BUILD_WITH_PETSC
template<uint_t Level>
inline void saveOperatorTmpl(Edge &edge, const PrimitiveDataID< StencilMemory< real_t >, Edge> &operatorId,
                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &srcId,
                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &dstId, Mat& mat) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto opr_data = edge.getData(operatorId)->getPointer( Level );
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 1; ++i) {
    PetscInt dstint = dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)];
    PetscInt srcint = src[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)];
    //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, VERTEX_C)], opr_data[VERTEX_C]);
    MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] ,INSERT_VALUES);         //TODO: Make this more efficient by grouping all of them in an array

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF ) {
      srcint = src[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] ,INSERT_VALUES);
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF ) {
      srcint = src[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] ,INSERT_VALUES);
    }

    if (edge.getNumNeighborFaces() == 2) {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF ) {
        srcint = src[vertexdof::macroedge::indexFromVertex<Level>(i, neighbor)];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] ,INSERT_VALUES);
      }
    }
  }
}

SPECIALIZE(void, saveOperatorTmpl, saveOperator)

template< typename ValueType, uint_t Level >
inline void createVectorFromFunctionTmpl(Edge &edge,
                                     const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                     const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId,
                                     Vec& vec) {
  PetscInt rowsize = (PetscInt) levelinfo::num_microvertices_per_edge(Level);

  auto src = edge.getData(srcId)->getPointer( Level );
  auto numerator = edge.getData(numeratorId)->getPointer( Level );

  VecSetValues(vec,rowsize-2,&numerator[1],&src[1],INSERT_VALUES);
}

SPECIALIZE_WITH_VALUETYPE(void, createVectorFromFunctionTmpl, createVectorFromFunction)

template< typename ValueType, uint_t Level >
inline void createFunctionFromVectorTmpl(Edge &edge,
                                         const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId,
                                         Vec& vec) {
  PetscInt rowsize = (PetscInt) levelinfo::num_microvertices_per_edge(Level);

  auto numerator = edge.getData(numeratorId)->getPointer( Level );

  VecGetValues(vec,rowsize-2,&numerator[1],&edge.getData(srcId)->getPointer( Level )[1]);
}

SPECIALIZE_WITH_VALUETYPE(void, createFunctionFromVectorTmpl, createFunctionFromVector)

template< uint_t Level >
inline void applyDirichletBCTmpl(Edge &edge,std::vector<PetscInt> &mat,
                                 const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId){

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for(uint_t i = 1;i<rowsize-1; i++)
  {
    mat.push_back(edge.getData(numeratorId)->getPointer( Level )[i]);
  }

}
SPECIALIZE(void, applyDirichletBCTmpl, applyDirichletBC)
#endif

} // namespace macroedge
} // namespace vertexdof
} // namespace hhg
