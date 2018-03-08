#pragma once

#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/levelinfo.hpp"

namespace hhg{
namespace VertexDoFToEdgeDoF{


template<uint_t Level>
inline void applyEdgeTmpl(Edge &edge,
                          const PrimitiveDataID<StencilMemory < real_t >, Edge> &operatorId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                          UpdateType update) {

  size_t rowsize = levelinfo::num_microedges_per_edge(Level);

  real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  real_t * src      = edge.getData(srcId)->getPointer( Level );
  real_t * dst      = edge.getData(dstId)->getPointer( Level );

  real_t tmp;

  for(uint_t i = 0; i < rowsize; ++i){
    tmp = 0.0;
    for(uint_t k = 0; k < hhg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF.size(); ++k){
      tmp += opr_data[hhg::vertexdof::stencilIndexFromHorizontalEdge(hhg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF[k])] *
             src[hhg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hhg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF[k] )];
    }
    for(uint_t k = 0; k < hhg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF.size(); ++k){
      tmp += opr_data[hhg::vertexdof::stencilIndexFromHorizontalEdge(hhg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF[k])] *
             src[hhg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hhg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF[k] )];
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < hhg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF.size(); ++k){
        tmp += opr_data[hhg::vertexdof::stencilIndexFromHorizontalEdge(hhg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF[k])] *
               src[hhg::vertexdof::macroedge::indexFromHorizontalEdge( Level, i, hhg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF[k] )];
      }
    }

    if (update==Replace) {
      dst[hhg::edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] = tmp;
    } else if (update==Add) {
      dst[hhg::edgedof::macroedge::indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] += tmp;
    }
  }

}


SPECIALIZE(void, applyEdgeTmpl, applyEdge)

template<uint_t Level>
inline void applyFaceTmpl(Face &face,
                          const PrimitiveDataID<StencilMemory < real_t >, Face> &operatorId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Face> &srcId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstId,
                          UpdateType update){

  real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  real_t * src      = face.getData(srcId)->getPointer( Level );
  real_t * dst      = face.getData(dstId)->getPointer( Level );

  real_t tmp;

  using namespace vertexdof::macroface;

  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromHorizontalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromHorizontalEdge(neighborsFromHorizontalEdge[k])] *
               src[indexFromHorizontalEdge( Level, it.col(), it.row(), neighborsFromHorizontalEdge[k] )];
      }
      if (update==Replace) {
        dst[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] += tmp;
      }
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromDiagonalEdge(neighborsFromDiagonalEdge[k])] *
               src[indexFromDiagonalEdge( Level, it.col(), it.row(), neighborsFromDiagonalEdge[k] )];
      }
      if (update==Replace) {
        dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] += tmp;
      }
    }
    if( it.col() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k){
        tmp += opr_data[vertexdof::stencilIndexFromVerticalEdge(neighborsFromVerticalEdge[k])] *
               src[indexFromVerticalEdge( Level, it.col(), it.row(), neighborsFromVerticalEdge[k] )];
      }

      if (update==Replace) {
        dst[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] = tmp;
      } else if ( update==Add ) {
        dst[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] += tmp;
      }
    }
  }

}

SPECIALIZE(void, applyFaceTmpl, applyFace)


} // VertexDoFToEdgeDoFOperator
} // namespace hhg
