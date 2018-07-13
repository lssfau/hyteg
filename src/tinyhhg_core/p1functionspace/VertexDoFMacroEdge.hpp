#pragma once

#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/facedofspace/FaceDoFIndexing.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/indexing/Common.hpp"

#include "core/math/KahanSummation.h"
#include "core/DataTypes.h"

#ifdef DEBUG_ELEMENTWISE
#include "tinyhhg_core/format.hpp"
#endif

namespace hhg {
namespace vertexdof {
namespace macroedge {

using walberla::real_c;
using indexing::Index;


inline Point3D coordinateFromIndex( const uint_t & level, const Edge & edge, const Index & index )
{
  const real_t  stepFrequency = 1.0 / levelinfo::num_microedges_per_edge( level );
  const Point3D step         = ( edge.getCoordinates()[1] - edge.getCoordinates()[0] ) * stepFrequency;
  return edge.getCoordinates()[0] + step * real_c( index.x() );
}

template<typename ValueType >
inline ValueType assembleLocal(const uint_t & level, uint_t pos, const Matrix3r& localMatrix,
                               double* src,
                               double* coeff,
                               const std::array< stencilDirection, 3 >& vertices,
                               const std::array<uint_t,3>& idx)
{

  ValueType meanCoeff = 1.0/3.0 * (coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[0] )]
                                 + coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[1] )]
                                 + coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[2] )]);

  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[vertexdof::macroedge::indexFromVertex( level, pos, vertices[0] )]
         + localMatrix(idx[0],idx[1]) * src[vertexdof::macroedge::indexFromVertex( level, pos, vertices[1] )]
         + localMatrix(idx[0],idx[2]) * src[vertexdof::macroedge::indexFromVertex( level, pos, vertices[2] )];
  return meanCoeff * tmp;
}

template<typename ValueType>
inline void assembleLocalStencil(uint_t level, uint_t pos, const Matrix3r& localMatrix,
                                 real_t* opr_data,
                                 real_t* coeff,
                                 const std::array< stencilDirection, 3 >& vertices,
                                 const std::array<uint_t,3>& idx)
{

  ValueType meanCoeff = 1.0/3.0 * (coeff[ vertexdof::macroedge::indexFromVertex( level, pos, vertices[ 0 ] ) ]
                                   + coeff[ vertexdof::macroedge::indexFromVertex( level, pos, vertices[ 1 ] ) ]
                                   + coeff[ vertexdof::macroedge::indexFromVertex( level, pos, vertices[ 2 ] ) ]);

  opr_data[vertexdof::stencilIndexFromVertex( vertices[0] )] += meanCoeff * localMatrix(idx[0],idx[0]);
  opr_data[vertexdof::stencilIndexFromVertex( vertices[1] )] += meanCoeff * localMatrix(idx[0],idx[1]);
  opr_data[vertexdof::stencilIndexFromVertex( vertices[2] )] += meanCoeff * localMatrix(idx[0],idx[2]);
}

template< typename ValueType >
inline void interpolate(const uint_t & level, Edge &edge,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge> &edgeMemoryId,
                        const ValueType & scalar )
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  auto edgeData =  edge.getData(edgeMemoryId)->getPointer( level );
  for (size_t i = 1; i < rowsize - 1; ++i)
  {
    edgeData[i] = scalar;
  }
}

template< typename ValueType >
inline void interpolate(const uint_t & level, Edge &edge,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Edge> &edgeMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &srcIds,
                            const std::function<ValueType(const hhg::Point3D &, const std::vector<ValueType>&)> &expr)
{
  ValueType * edgeData = edge.getData( edgeMemoryId )->getPointer( level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( edge.getData(src)->getPointer( level ) );
  }

  std::vector<ValueType> srcVector( srcIds.size() );

  Point3D xBlend;

  for ( const auto & it : vertexdof::macroedge::Iterator( level, 1 ) )
  {
    const Point3D coordinate = coordinateFromIndex( level, edge, it );
    const uint_t  idx        = vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C );

    for ( uint_t k = 0; k < srcPtr.size(); ++k )
    {
      srcVector[ k ] = srcPtr[ k ][ idx ];
    }
    edge.getGeometryMap()->evalF(coordinate, xBlend);
    edgeData[ idx ] = expr( xBlend, srcVector );
  }
}


template< typename ValueType >
inline void assign( const uint_t & level, Edge &edge,
                   const std::vector<ValueType> &scalars,
                   const std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge>> &srcIds,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge> &dstId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    ValueType tmp = scalars[0]*edge.getData(srcIds[0])->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];

    for (size_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
    }

    edge.getData(dstId)->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = tmp;
  }
}


template< typename ValueType >
inline void add(const uint_t & level,
                const Edge & edge,
                const ValueType & scalar,
                const PrimitiveDataID<FunctionMemory< ValueType >, Edge> & dstId )
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i)
  {
    edge.getData(dstId)->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] += scalar;
  }
}

template< typename ValueType >
inline void add( const uint_t & level, Edge &edge,
                const std::vector<ValueType> &scalars,
                const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &srcIds,
                const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    ValueType tmp = 0.0;

    for (size_t k = 0; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
    }

    edge.getData(dstId)->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] += tmp;
  }
}


template< typename ValueType >
inline real_t dot( const uint_t & level, Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &lhsMemoryId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsMemoryId) {

  walberla::math::KahanAccumulator< ValueType > scalarProduct;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    scalarProduct += edge.getData(lhsMemoryId)->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )]
        * edge.getData(rhsMemoryId)->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
  }

  return scalarProduct.get();
}


template< typename ValueType >
inline void apply( const uint_t & level, Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId, UpdateType update) {

  typedef stencilDirection sD;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto opr_data = edge.getData(operatorId)->getPointer( level );
  auto src = edge.getData(srcId)->getPointer( level );
  auto dst = edge.getData(dstId)->getPointer( level );

  ValueType tmp;

  for (size_t i = 1; i < rowsize - 1; ++i)
  {
    const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
    const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
    const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

    const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
    const auto dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
    const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

    tmp = opr_data[ stencilIdxW ] * src[ dofIdxW ]
      + opr_data[ stencilIdxC ] * src[ dofIdxC ]
      + opr_data[ stencilIdxE ] * src[ dofIdxE ];

    for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
    {
      const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
      const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
      const auto stencilWeightW = opr_data[ stencilIdxW ];
      const auto stencilWeightE = opr_data[ stencilIdxE ];
      const auto dofIdxW = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
      const auto dofIdxE = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
      tmp += stencilWeightW * src[dofIdxW] + stencilWeightE * src[dofIdxE];
    }

    for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
    {
      const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
      const auto dofIdx = vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
      tmp += opr_data[ stencilIdx ] * src[ dofIdx ];
    }

    if (update == Replace) {
      dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = tmp;
    } else if (update == Add) {
      dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] += tmp;
    }
  }
}


template< typename ValueType >
inline void applyCoefficient( const uint_t & level, Edge &edge,
                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                 const std::vector<PrimitiveDataID<EdgeP1LocalMatrixMemory, Edge>> &operatorIds,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                                 const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &coeffIds,
                                 UpdateType update) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto src = edge.getData(srcId)->getPointer( level );
  auto dst = edge.getData(dstId)->getPointer( level );

  std::vector<EdgeP1LocalMatrixMemory*> localMatricesVector;
  for(auto operatorId : operatorIds) {
    localMatricesVector.push_back(edge.getData(operatorId));
  }

  std::vector<real_t*> coeffs;
  for(auto coeffId : coeffIds) {
    coeffs.push_back(edge.getData(coeffId)->getPointer( level ));
  }

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
      tmp = dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
    }

    for (uint_t coeffIdx = 0; coeffIdx < coeffIds.size(); ++coeffIdx) {

      tmp += assembleLocal<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 0), src, coeffs[coeffIdx], triangleGraySW,
                                             {e_south, s_south, o_south});
      tmp += assembleLocal<ValueType>(level, i, localMatricesVector[coeffIdx]->getBlueMatrix(level, 0), src, coeffs[coeffIdx], triangleBlueS,
                                             {o_south, e_south, s_south});
      tmp += assembleLocal<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 0), src, coeffs[coeffIdx], triangleGraySE,
                                             {s_south, o_south, e_south});

      if (edge.getNumNeighborFaces() == 2) {
        tmp += assembleLocal<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 1), src, coeffs[coeffIdx], triangleGrayNW,
                                               {e_north, s_north, o_north});
        tmp += assembleLocal<ValueType>(level, i, localMatricesVector[coeffIdx]->getBlueMatrix(level, 1), src, coeffs[coeffIdx], triangleBlueN,
                                               {o_north, e_north, s_north});
        tmp += assembleLocal<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 1), src, coeffs[coeffIdx], triangleGrayNE,
                                               {s_north, o_north, e_north});
      }
    }

    dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = tmp;
  }
}


template< typename ValueType >
inline void smooth_gs( const uint_t & level, Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId) {

  typedef stencilDirection sD;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto opr_data = edge.getData(operatorId)->getPointer( level );
  auto rhs = edge.getData(rhsId)->getPointer( level );
  auto dst = edge.getData(dstId)->getPointer( level );

  const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
  const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
  const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

  const auto invCenterWeight = 1.0 / opr_data[ stencilIdxC ];

  ValueType tmp;

  for (size_t i = 1; i < rowsize - 1; ++i)
  {
    const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
    const auto dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
    const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

    tmp = rhs[ dofIdxC ];

    tmp -= opr_data[ stencilIdxW ] * dst[ dofIdxW ] + opr_data[ stencilIdxE ] * dst[ dofIdxE ];

    for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
    {
      const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
      const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
      const auto stencilWeightW = opr_data[ stencilIdxW ];
      const auto stencilWeightE = opr_data[ stencilIdxE ];
      const auto dofIdxW = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
      const auto dofIdxE = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
      tmp -= stencilWeightW * dst[dofIdxW] + stencilWeightE * dst[dofIdxE];
    }

    for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
    {
      const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
      const auto dofIdx = vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
      tmp -= opr_data[ stencilIdx ] * dst[ dofIdx ];
    }

    dst[ dofIdxC ] = invCenterWeight * tmp;
  }
}


template< typename ValueType >
inline void smooth_gs_coefficient(uint_t level, Edge &edge,
                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                  const std::vector<PrimitiveDataID<EdgeP1LocalMatrixMemory, Edge>> &operatorIds,
                                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                                  const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                                  const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &coeffIds) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto dst = edge.getData(dstId)->getPointer( level );
  auto rhs = edge.getData(rhsId)->getPointer( level );

  std::vector<EdgeP1LocalMatrixMemory*> localMatricesVector;
  for(auto operatorId : operatorIds) {
    localMatricesVector.push_back(edge.getData(operatorId));
  }

  std::vector<real_t*> coeffs;
  for(auto coeffId : coeffIds) {
    coeffs.push_back(edge.getData(coeffId)->getPointer( level ));
  }

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

  std::vector<real_t> opr_data(7);

  for (size_t i = 1; i < rowsize - 1; ++i) {

    std::fill(opr_data.begin(), opr_data.end(), 0.0);

    for (uint_t coeffIdx = 0; coeffIdx < coeffIds.size(); ++coeffIdx) {

      assembleLocalStencil<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 0), opr_data.data(), coeffs[coeffIdx],
                                             triangleGraySW, {e_south, s_south, o_south});
      assembleLocalStencil<ValueType>(level, i, localMatricesVector[coeffIdx]->getBlueMatrix(level, 0), opr_data.data(), coeffs[coeffIdx],
                                             triangleBlueS, {o_south, e_south, s_south});
      assembleLocalStencil<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 0), opr_data.data(), coeffs[coeffIdx],
                                             triangleGraySE, {s_south, o_south, e_south});

      if (edge.getNumNeighborFaces() == 2) {
        assembleLocalStencil<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 1), opr_data.data(), coeffs[coeffIdx],
                                               triangleGrayNW, {e_north, s_north, o_north});
        assembleLocalStencil<ValueType>(level, i, localMatricesVector[coeffIdx]->getBlueMatrix(level, 1), opr_data.data(), coeffs[coeffIdx],
                                               triangleBlueN, {o_north, e_north, s_north});
        assembleLocalStencil<ValueType>(level, i, localMatricesVector[coeffIdx]->getGrayMatrix(level, 1), opr_data.data(), coeffs[coeffIdx],
                                               triangleGrayNE, {s_north, o_north, e_north});
      }
    }

    dst[vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C) ] = rhs[vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C)];

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      dst[ vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C) ] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex(level, i, neighbor)];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C) ] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex(level, i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        dst[ vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C) ] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex(level, i, neighbor)];
      }
    }

    dst[ vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C) ] /= opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
inline void smooth_sor(const uint_t & level, Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                          ValueType relax) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto opr_data = edge.getData(operatorId)->getPointer( level );
  auto dst = edge.getData(dstId)->getPointer( level );
  auto rhs = edge.getData(rhsId)->getPointer( level );

  ValueType tmp;

  for (size_t i = 1; i < rowsize - 1; ++i) {

    tmp = rhs[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      tmp -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      tmp -= opr_data[vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        tmp -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * dst[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
      }
    }

    dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = ( 1.0 - relax) * dst[vertexdof::macroedge::indexFromVertex( level, i,
                                                                                                                                                                                 stencilDirection::VERTEX_C )]
                                                                                                            + relax * tmp/opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ];
  }
}


template< typename ValueType >
inline void smooth_jac(const uint_t & level, Edge &edge, const PrimitiveDataID< StencilMemory< ValueType >, Edge> &operatorId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                          const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &tmpId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto opr_data = edge.getData(operatorId)->getPointer( level );
  auto dst = edge.getData(dstId)->getPointer( level );
  auto rhs = edge.getData(rhsId)->getPointer( level );
  auto tmp = edge.getData(tmpId)->getPointer( level );

  for (size_t i = 1; i < rowsize - 1; ++i) {

    dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = rhs[vertexdof::macroedge::indexFromVertex( level, i,
                                                                                                                                                                stencilDirection::VERTEX_C )];

    for (auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * tmp[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
    }

    for (auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * tmp[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
    }

    if (edge.getNumNeighborFaces() == 2)
    {
      for (auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] -= opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] * tmp[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
      }
    }

    dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] /= opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ];
  }
}


template< typename ValueType >
inline void enumerate( const uint_t & level, Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId, uint_t& num) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    edge.getData(dstId)->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = walberla::real_c( num++);
  }
}


template< typename ValueType >
inline void integrateDG(const uint_t & level, Edge &edge,
                            const std::shared_ptr< PrimitiveStorage >& storage,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsP1Id,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId) {

  typedef stencilDirection sD;

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto rhs = edge.getData(rhsId)->getPointer( level );
  auto rhsP1 = edge.getData(rhsP1Id)->getPointer( level );
  auto dst = edge.getData(dstId)->getPointer( level );

  ValueType tmp;

  Face* face = storage->getFace(edge.neighborFaces()[0]);
  real_t weightedFaceArea0, weightedFaceArea1;

  weightedFaceArea0 = std::pow(4.0, -walberla::real_c(level)) * face->area / 3.0;

  uint_t s_north, e_north, o_north;

  if (edge.getNumNeighborFaces() == 2) {
    face = storage->getFace(edge.neighborFaces()[1]);
    weightedFaceArea1 = std::pow(4.0, -walberla::real_c(level)) * face->area / 3.0;
  }

  for (size_t i = 1; i < rowsize - 1; ++i) {

    tmp = weightedFaceArea0 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_SW )] * ( 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i,
                                                                                                                                                                        sD::VERTEX_C )]
                                                                                                                           + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W )] )
                                                                                                             + 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                                                                                                                             rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_S )] ));
    tmp += weightedFaceArea0 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_BLUE_SE )] * ( 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                                                                                                                            rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_S )] ) +
                                                                                                              0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                                                                                                                            rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_SE )] ));
    tmp += weightedFaceArea0 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_SE )] *
      ( 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_SE )] ) +
        0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E )] ));

    if (edge.getNumNeighborFaces() == 2) {

      tmp += weightedFaceArea1 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_NW )] * ( 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )]
                                                                                                                              + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W )])
                                                                                                                + 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )]
                                                                                                                                + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_NW )]));
      tmp += weightedFaceArea1 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_BLUE_NW )] * ( 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )]
                                                                                                                              + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_NW )])
                                                                                                                + 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )]
                                                                                                                                + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_N )]));
      tmp += weightedFaceArea1 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_NE )] * ( 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )]
                                                                                                                              + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_N )])
                                                                                                                + 0.5 * 0.5 * ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )]
                                                                                                                                + rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E )]));
    }

    dst[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] = tmp;
  }
}

template< typename ValueType >
inline void printFunctionMemory( const uint_t & level, const Edge& edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId){
  ValueType* edgeMemory = edge.getData(dstId)->getPointer( level );
  using namespace std;
  cout << setfill('=') << setw(100) << "" << endl;
  cout << edge << std::left << setprecision(1) << fixed << setfill(' ') << endl;
  uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

  if(edge.getNumNeighborFaces() == 2) {
    for (uint_t i = 0; i < (rowsize - 1); ++i) {
      cout << setw(5) << edgeMemory[hhg::vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_N )] << "|";
    }
    cout << endl;
  }
  for(uint_t i = 0; i < rowsize; ++i){
    cout << setw(5) << edgeMemory[hhg::vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] << "|";
  }
  cout << endl << "     |";
  for(uint_t i = 1; i < rowsize; ++i){
    cout << setw(5) << edgeMemory[hhg::vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_S )] << "|";
  }
  cout << endl << setfill('=') << setw(100) << "" << endl << setfill(' ');

}


template< typename ValueType >
inline real_t getMaxValue( const uint_t & level, Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId ) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

  auto src = edge.getData( srcId )->getPointer( level );
  real_t localMax = -std::numeric_limits<real_t>::max();

  for( size_t i = 1; i < rowsize - 1; ++i ) {
    localMax = std::max( localMax, src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C ) ] );
  }

  return localMax;
}


template< typename ValueType >
inline real_t getMaxMagnitude( const uint_t & level, Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId ) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

  auto src = edge.getData( srcId )->getPointer( level );
  real_t localMax = real_t(0.0);

  for( size_t i = 1; i < rowsize - 1; ++i ) {
    localMax = std::max( localMax, std::abs( src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C ) ] ));
  }

  return localMax;
}


template< typename ValueType >
inline real_t getMinValue( const uint_t & level, Edge &edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId ) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

  auto src = edge.getData( srcId )->getPointer( level );
  real_t localMin = std::numeric_limits<real_t>::max();

  for( size_t i = 1; i < rowsize - 1; ++i ) {
    localMin = std::min( localMin, src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C ) ] );
  }

  return localMin;
}



#ifdef HHG_BUILD_WITH_PETSC
inline void saveOperator( const uint_t & level, Edge &edge, const PrimitiveDataID< StencilMemory< real_t >, Edge> &operatorId,
                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &srcId,
                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &dstId, Mat& mat) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  auto opr_data = edge.getData(operatorId)->getPointer( level );
  auto src = edge.getData(srcId)->getPointer( level );
  auto dst = edge.getData(dstId)->getPointer( level );


  for (uint_t i = 1; i < rowsize - 1; ++i) {
    PetscInt dstint = dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
    PetscInt srcint = src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
    //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, VERTEX_C)], opr_data[VERTEX_C]);
    MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C ) ] ,INSERT_VALUES);         //TODO: Make this more efficient by grouping all of them in an array

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF ) {
      srcint = src[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] ,INSERT_VALUES);
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF ) {
      srcint = src[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] ,INSERT_VALUES);
    }

    if (edge.getNumNeighborFaces() == 2) {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF ) {
        srcint = src[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[ vertexdof::stencilIndexFromVertex( neighbor ) ] ,INSERT_VALUES);
      }
    }
  }
}


template< typename ValueType >
inline void createVectorFromFunction(const uint_t & level, Edge &edge,
                                     const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                     const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId,
                                     Vec& vec) {
  PetscInt rowsize = (PetscInt) levelinfo::num_microvertices_per_edge(level);

  auto src = edge.getData(srcId)->getPointer( level );
  auto numerator = edge.getData(numeratorId)->getPointer( level );

  VecSetValues(vec,rowsize-2,&numerator[1],&src[1],INSERT_VALUES);
}


template< typename ValueType >
inline void createFunctionFromVector(const uint_t & level, Edge &edge,
                                         const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId,
                                         Vec& vec) {
  PetscInt rowsize = (PetscInt) levelinfo::num_microvertices_per_edge(level);

  auto numerator = edge.getData(numeratorId)->getPointer( level );

  VecGetValues(vec,rowsize-2,&numerator[1],&edge.getData(srcId)->getPointer( level )[1]);
}



inline void applyDirichletBC(const uint_t & level, Edge &edge,std::vector<PetscInt> &mat,
                                 const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId){

  size_t rowsize = levelinfo::num_microvertices_per_edge(level);

  for(uint_t i = 1;i<rowsize-1; i++)
  {
    mat.push_back(edge.getData(numeratorId)->getPointer( level )[i]);
  }

}
#endif

} // namespace macroedge
} // namespace vertexdof
} // namespace hhg
