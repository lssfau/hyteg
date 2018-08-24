#pragma once

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/primitives/Cell.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements3D.hpp"

namespace hhg{
namespace EdgeDoFToVertexDoF {

inline void applyVertex(uint_t level,
                  Vertex &vertex,
                  const PrimitiveDataID<StencilMemory < real_t >, Vertex> &operatorId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &srcId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstId,
                  UpdateType update)
{

  real_t * opr_data = vertex.getData(operatorId)->getPointer( level );
  real_t * src      = vertex.getData(srcId)->getPointer( level );
  real_t * dst      = vertex.getData(dstId)->getPointer( level );

  real_t tmp = 0;
  for(uint_t i = 0; i < vertex.getData(operatorId)->getSize( level ); ++i){
    tmp += src[i] * opr_data[i];
  }

  if (update==Replace) {
    dst[0] = tmp;
  } else if (update==Add) {
    dst[0] += tmp;
  }
}


inline void applyEdge(const uint_t & Level, Edge &edge,
                  const PrimitiveDataID<StencilMemory < real_t >, Edge> &operatorId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                  UpdateType update)
{
  using namespace hhg::edgedof;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  real_t * src      = edge.getData(srcId)->getPointer( Level );
  real_t * dst      = edge.getData(dstId)->getPointer( Level );

  real_t tmp;

  for(uint_t i = 1; i < rowsize - 1; ++i){
    tmp = 0.0;
    for(uint_t k = 0; k < macroedge::neighborsOnEdgeFromVertex.size(); ++k){
      tmp += opr_data[stencilIndexFromVertex(macroedge::neighborsOnEdgeFromVertex[k])] *
             src[macroedge::indexFromVertex( Level, i, macroedge::neighborsOnEdgeFromVertex[k] )];
    }
    for(uint_t k = 0; k < macroedge::neighborsOnSouthFaceFromVertex.size(); ++k){
      tmp += opr_data[stencilIndexFromVertex(macroedge::neighborsOnSouthFaceFromVertex[k])] *
             src[macroedge::indexFromVertex( Level, i, macroedge::neighborsOnSouthFaceFromVertex[k] )];
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < macroedge::neighborsOnNorthFaceFromVertex.size(); ++k){
        tmp += opr_data[stencilIndexFromVertex(macroedge::neighborsOnNorthFaceFromVertex[k])] *
               src[macroedge::indexFromVertex( Level, i, macroedge::neighborsOnNorthFaceFromVertex[k] )];
      }
    }
    if (update==Replace) {
      dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] = tmp;
    } else if (update==Add) {
      dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] += tmp;
    }
  }
}



inline void applyFace(const uint_t & Level, Face &face,
                       const PrimitiveDataID<StencilMemory < real_t >, Face> &operatorId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstId,
                       UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  real_t * src      = face.getData(srcId)->getPointer( Level );
  real_t * dst      = face.getData(dstId)->getPointer( Level );

  real_t tmp;

  using namespace edgedof::macroface;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = 0.0;

      for(uint_t k = 0; k < neighborsFromVertex.size(); ++k){
        tmp += opr_data[edgedof::stencilIndexFromVertex(neighborsFromVertex[k])] *
               src[indexFromVertex( Level, i, j, neighborsFromVertex[k] )];
      }

      if (update==Replace) {
        dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
      } else if (update==Add) {
        dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
      }
    }
    --inner_rowsize;
  }
}


inline void applyCell(const uint_t & Level, Cell & cell,
                      const PrimitiveDataID<StencilMemory < real_t >, Cell> &operatorId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Cell> &srcId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Cell> &dstId,
                      UpdateType update)
{

  real_t * opr_data = cell.getData(operatorId)->getPointer( Level );
  real_t * src      = cell.getData(srcId)->getPointer( Level );
  real_t * dst      = cell.getData(dstId)->getPointer( Level );

  real_t tmp;

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
    tmp = 0.0;
    for ( const auto & orientation : edgedof::EdgeDoFOrientation )
    {
      const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromVertexDoFInMacroCell( orientation );
      for ( const auto & neighbor : edgeDoFNeighbors )
      {
        const uint_t stencilIdx  = edgedof::stencilIndexFromVertex3D( neighbor, orientation );
        const auto   srcIdx      = it + neighbor;
        const auto   srcArrayIdx = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), orientation );
        tmp += opr_data[stencilIdx] * src[srcArrayIdx];
      }
    }

    const auto   dstArrayIdx = vertexdof::macrocell::index( Level, it.x(), it.y(), it.z() );

    if ( update == Replace )
    {
      dst[dstArrayIdx] = tmp;
    }
    else if ( update == Add )
    {
      dst[dstArrayIdx] += tmp;
    }
  }
}


} /// EdgeDoFToVertexDoF
} /// namespace hhg
