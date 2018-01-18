#include "P2SmoothTest.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

namespace hhg {
namespace P2 {

namespace vertex {

void smoothGSvertexDof(Vertex vertex, const PrimitiveDataID<StencilMemory<double>, Vertex> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory<real_t>, Vertex> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory<double>, Vertex> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory<real_t>, Vertex> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Vertex> &getVertexDoFID, uint_t level) {

}
} /// namespace vertex

namespace edge {

void smoothGSvertexDof(Edge edge, const PrimitiveDataID<StencilMemory<double>, Edge> &vertexDoFStencil,
                       const PrimitiveDataID<FunctionMemory<real_t>, Edge> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory<double>, Edge> &edgeDoFStencil,
                       const PrimitiveDataID<FunctionMemory<real_t>, Edge> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Edge> &rhsVertexDoFID, uint_t level) {

}

void smoothGSedgeDof(Edge edge, const PrimitiveDataID<StencilMemory<double>, Edge> &vertexDoFStencil,
                     const PrimitiveDataID<FunctionMemory<real_t>, Edge> &dstVertexDoFID,
                     const PrimitiveDataID<StencilMemory<double>, Edge> &edgeDoFStencil,
                     const PrimitiveDataID<FunctionMemory<real_t>, Edge> &dstEdgeDoFID,
                     const PrimitiveDataID<FunctionMemory<real_t>, Edge> &rhsEdgeDoFID, uint_t level) {

}
} /// namespace edge

namespace face {

template<uint_t Level>
void smoothGSvertexDof_tmpl(Face face, const PrimitiveDataID<StencilMemory<double>, Face> &vertexDoFStencilID,
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

  for (size_t i = 1; i < rowsize - 2; ++i) {
    for (size_t j = 1; j < inner_rowsize - 2; ++j) {
      tmp = rhs[vertexdof::macroface::indexFromVertex< Level >(i, j, stencilDirection::VERTEX_C)];

      /// update from vertex dofs
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp += vertexDoFStencil[edgedof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])] *
               dstVertexDoF[vertexdof::macroface::indexFromVertex< Level >(i, j, vertexdof::macroface::neighborsWithoutCenter[k])];
      }
      /// update from vertex dofs
      for(uint_t k = 0; k < edgedof::macroface::neighborsFromVertex.size(); ++k){
        tmp += vertexDoFStencil[edgedof::stencilIndexFromVertex(edgedof::macroface::neighborsFromVertex[k])] *
               dstVertexDoF[edgedof::macroface::indexFromVertex< Level >(i, j, edgedof::macroface::neighborsFromVertex[k])];
      }
      dstEdgeDoF[vertexdof::macroface::indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] = tmp / vertexDoFStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];

    }
    --inner_rowsize;
  }

}

void smoothGSedgeDof(Face face, const PrimitiveDataID<StencilMemory<double>, Face> &vertexDoFStencil,
                     const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstVertexDoFID,
                     const PrimitiveDataID<StencilMemory<double>, Face> &edgeDoFStencil,
                     const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstEdgeDoFID,
                     const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsEdgeDoFID, uint_t level) {

}
} /// namespace face

} /// namespace P2
} /// namespace hhg