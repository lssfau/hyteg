#pragma once

#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"

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


template<uint_t Level>
inline void applyEdgeTmpl(Edge &edge,
                  const PrimitiveDataID<StencilMemory < real_t >, Edge> &operatorId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                  const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                  UpdateType update)
{
  using namespace hhg::indexing::edgedof;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  real_t * src      = edge.getData(srcId)->getPointer( Level );
  real_t * dst      = edge.getData(dstId)->getPointer( Level );

  real_t tmp;

  for(uint_t i = 1; i < rowsize - 1; ++i){
    tmp = 0.0;
    for(uint_t k = 0; k < macroedge::neighborsOnEdgeFromVertex.size(); ++k){
      tmp += opr_data[stencilIndexFromVertex(macroedge::neighborsOnEdgeFromVertex[k])] *
             src[macroedge::indexFromVertex< Level >(i, macroedge::neighborsOnEdgeFromVertex[k])];
    }
    for(uint_t k = 0; k < macroedge::neighborsOnSouthFaceFromVertex.size(); ++k){
      tmp += opr_data[stencilIndexFromVertex(macroedge::neighborsOnSouthFaceFromVertex[k])] *
             src[macroedge::indexFromVertex< Level >(i, macroedge::neighborsOnSouthFaceFromVertex[k])];
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < macroedge::neighborsOnNorthFaceFromVertex.size(); ++k){
        tmp += opr_data[stencilIndexFromVertex(macroedge::neighborsOnNorthFaceFromVertex[k])] *
               src[macroedge::indexFromVertex< Level >(i, macroedge::neighborsOnNorthFaceFromVertex[k])];
      }
    }
    if (update==Replace) {
      dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] = tmp;
    } else if (update==Add) {
      dst[vertexdof::macroedge::indexFromVertex<Level>(i, stencilDirection::VERTEX_C)] += tmp;
    }
  }
}

SPECIALIZE(void, applyEdgeTmpl, applyEdge)


template<uint_t Level>
inline void applyFaceTmpl(Face &face,
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

  using namespace indexing::edgedof::macroface;

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = 0.0;

      for(uint_t k = 0; k < neighborsFromVertex.size(); ++k){
        tmp += opr_data[indexing::edgedof::stencilIndexFromVertex(neighborsFromVertex[k])] *
               src[indexFromVertex< Level >(i, j, neighborsFromVertex[k])];
      }

      if (update==Replace) {
        dst[ vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp;
      } else if (update==Add) {
        dst[ vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, applyFaceTmpl, applyFace)

} /// EdgeDoFToVertexDoF
} /// namespace hhg
