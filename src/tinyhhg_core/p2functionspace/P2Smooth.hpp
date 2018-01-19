#pragma once

#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/macros.hpp"

namespace hhg{
namespace P2{

namespace vertex {

void smoothGSvertexDoF(Vertex vertex,
                       const PrimitiveDataID<StencilMemory< real_t >, Vertex> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory< real_t >, Vertex> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &getVertexDoFID,
                       uint_t level);

} /// namespace vertex

namespace edge {

void smoothGSvertexDoF(Edge edge,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &rhsVertexDoFID,
                       uint_t level);

void smoothGSedgeDoF(Edge edge,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &rhsEdgeDoFID,
                       uint_t level);

} /// namespace edge

namespace face {


template< uint_t Level>
void smoothGSvertexDoFTmpl(Face &face, const PrimitiveDataID<StencilMemory<double>, Face> &vertexDoFStencilID,
                            const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstVertexDoFID,
                            const PrimitiveDataID<StencilMemory<double>, Face> &edgeDoFStencilID,
                            const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstEdgeDoFID,
                            const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsVertexDoFID) {

  size_t rowsize = levelinfo::num_microvertices_per_edge( Level );
  size_t inner_rowsize = rowsize;

  real_t * vertexDoFStencil = face.getData(vertexDoFStencilID)->getPointer( Level );
  real_t * dstVertexDoF = face.getData(dstVertexDoFID)->getPointer( Level );
  real_t * edgeDoFStencil = face.getData(edgeDoFStencilID)->getPointer( Level );
  real_t * dstEdgeDoF = face.getData(dstEdgeDoFID)->getPointer( Level );
  real_t * rhs = face.getData(rhsVertexDoFID)->getPointer( Level );


  real_t tmp;

  for (size_t row = 1; row < rowsize - 2; ++row) {
    for (size_t col = 1; col < inner_rowsize - 2; ++col) {
      tmp = rhs[vertexdof::macroface::indexFromVertex< Level >(col, row, stencilDirection::VERTEX_C)];

      /// update from vertex dofs
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        size_t tmp2 = vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k]);
        uint_t tmp3 = vertexdof::macroface::indexFromVertex< Level >(col, row, vertexdof::macroface::neighborsWithoutCenter[k]);
        tmp += vertexDoFStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])] *
               dstVertexDoF[vertexdof::macroface::indexFromVertex< Level >(col, row, vertexdof::macroface::neighborsWithoutCenter[k])];
      }
      /// update from edge dofs
      for(uint_t k = 0; k < edgedof::macroface::neighborsFromVertex.size(); ++k){
        tmp += edgeDoFStencil[edgedof::stencilIndexFromVertex(edgedof::macroface::neighborsFromVertex[k])] *
               dstEdgeDoF[edgedof::macroface::indexFromVertex< Level >(col, row, edgedof::macroface::neighborsFromVertex[k])];
      }
      uint_t tmp5 = vertexdof::macroface::indexFromVertex<Level>(col, row, stencilDirection::VERTEX_C);
      dstVertexDoF[vertexdof::macroface::indexFromVertex<Level>(col, row, stencilDirection::VERTEX_C)] = tmp / vertexDoFStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];

    }
    --inner_rowsize;
  }

}

SPECIALIZE(void, smoothGSvertexDoFTmpl, smoothGSvertexDoF)

void smoothGSedgeDoF(Face face,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &vertexDoFStencil,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstVertexDoFID,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &edgeDoFStencil,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstEdgeDoFID,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &rhsEdgeDoFID,
                     uint_t level);

} /// namespace face

} /// namespace P2
} /// namespace hhg