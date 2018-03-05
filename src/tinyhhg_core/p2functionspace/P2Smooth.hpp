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

void smoothGSvertexDoF(uint_t level, Vertex &vertex,
                           const PrimitiveDataID<StencilMemory< real_t >, Vertex> &vertexDoFStencilID,
                           const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstVertexDoFID,
                           const PrimitiveDataID<StencilMemory< real_t >, Vertex> &edgeDoFStencilID,
                           const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &dstEdgeDoFID,
                           const PrimitiveDataID<FunctionMemory< real_t >, Vertex> &rhsVertexDoFID){
  real_t * vertexDoFStencil = vertex.getData(vertexDoFStencilID)->getPointer( level );
  real_t * dstVertexDoF = vertex.getData(dstVertexDoFID)->getPointer( level );
  real_t * edgeDoFStencil = vertex.getData(edgeDoFStencilID)->getPointer( level );
  real_t * dstEdgeDoF = vertex.getData(dstEdgeDoFID)->getPointer( level );
  real_t * rhs = vertex.getData(rhsVertexDoFID)->getPointer( level );

  real_t tmp = 0;
  tmp = rhs[0];
  for(uint_t i = 0; i < vertex.getData(edgeDoFStencilID)->getSize( level ); ++i){
    tmp -= dstEdgeDoF[i] * edgeDoFStencil[i];
  }
  for(uint_t i = 1; i < vertex.getData(vertexDoFStencilID)->getSize( level); ++i){
    tmp -= dstVertexDoF[i] * vertexDoFStencil[i];
  }

  dstVertexDoF[0] = tmp / vertexDoFStencil[0];

}

} /// namespace vertex

namespace edge {

template< uint_t Level >
void smoothGSvertexDoFTmpl(Edge &edge,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &vertexDoFStencilID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory< real_t >, Edge> &edgeDoFStencilID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &rhsVertexDoFID){
  size_t rowsize = levelinfo::num_microvertices_per_edge( Level );

  real_t * vertexDoFStencil = edge.getData(vertexDoFStencilID)->getPointer( Level );
  real_t * dstVertexDoF = edge.getData(dstVertexDoFID)->getPointer( Level );
  real_t * edgeDoFStencil = edge.getData(edgeDoFStencilID)->getPointer( Level );
  real_t * dstEdgeDoF = edge.getData(dstEdgeDoFID)->getPointer( Level );
  real_t * rhs = edge.getData(rhsVertexDoFID)->getPointer( Level );

  real_t tmp;

  for(uint_t i = 1; i < rowsize - 1; ++i){
    tmp = rhs[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )];

    // neighbors on edge
    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF ) {
      tmp -= vertexDoFStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ] *
        dstVertexDoF[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
    }

    for(uint_t k = 0; k < edgedof::macroedge::neighborsOnEdgeFromVertex.size(); ++k){
      tmp -= edgeDoFStencil[edgedof::stencilIndexFromVertex(edgedof::macroedge::neighborsOnEdgeFromVertex[k])] *
             dstEdgeDoF[edgedof::macroedge::indexFromVertex< Level >(i, edgedof::macroedge::neighborsOnEdgeFromVertex[k])];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF ) {
      tmp -= vertexDoFStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ] *
        dstVertexDoF[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
    }
    for(uint_t k = 0; k < edgedof::macroedge::neighborsOnSouthFaceFromVertex.size(); ++k){
      tmp -= edgeDoFStencil[edgedof::stencilIndexFromVertex(edgedof::macroedge::neighborsOnSouthFaceFromVertex[k])] *
             dstEdgeDoF[edgedof::macroedge::indexFromVertex< Level >(i, edgedof::macroedge::neighborsOnSouthFaceFromVertex[k])];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF ) {
        tmp -= vertexDoFStencil[ vertexdof::stencilIndexFromVertex( neighbor ) ] *
          dstVertexDoF[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
      }
      for(uint_t k = 0; k < edgedof::macroedge::neighborsOnNorthFaceFromVertex.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromVertex(edgedof::macroedge::neighborsOnNorthFaceFromVertex[k])] *
               dstEdgeDoF[edgedof::macroedge::indexFromVertex< Level >(i, edgedof::macroedge::neighborsOnNorthFaceFromVertex[k])];
      }
    }

    dstVertexDoF[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] =
      tmp / vertexDoFStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];

    }
}

SPECIALIZE(void, smoothGSvertexDoFTmpl, smoothGSvertexDoF)

template< uint_t Level >
void smoothGSedgeDoFTmpl(Edge &edge,
                         const PrimitiveDataID<StencilMemory< real_t >, Edge> &vertexDoFStencilID,
                         const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstVertexDoFID,
                         const PrimitiveDataID<StencilMemory< real_t >, Edge> &edgeDoFStencilID,
                         const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstEdgeDoFID,
                         const PrimitiveDataID<FunctionMemory< real_t >, Edge> &rhsEdgeDoFID){
  size_t rowsize = levelinfo::num_microedges_per_edge( Level );

  real_t * vertexDoFStencil = edge.getData(vertexDoFStencilID)->getPointer( Level );
  real_t * dstVertexDoF = edge.getData(dstVertexDoFID)->getPointer( Level );
  real_t * edgeDoFStencil = edge.getData(edgeDoFStencilID)->getPointer( Level );
  real_t * dstEdgeDoF = edge.getData(dstEdgeDoFID)->getPointer( Level );
  real_t * rhs = edge.getData(rhsEdgeDoFID)->getPointer( Level );

  real_t tmp;

  for(uint_t i = 0; i < rowsize ; ++i){
    tmp = rhs[edgedof::macroedge::indexFromHorizontalEdge< Level >(i, stencilDirection::EDGE_HO_C)];


    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF) {
      tmp -= vertexDoFStencil[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ] *
             dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( Level, i, neighbor )];
    }

    for ( const auto & neighbor : edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge) {
      tmp -= edgeDoFStencil[ edgedof::stencilIndexFromHorizontalEdge( neighbor ) ] *
             dstEdgeDoF[ edgedof::macroedge::indexFromHorizontalEdge<Level>( i, neighbor ) ];
    }

    for ( const auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF) {
      tmp -= vertexDoFStencil[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ] *
             dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( Level, i, neighbor )];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (const auto &neighbor : edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge) {
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge(neighbor)] *
               dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge<Level>(i, neighbor)];
      }

      for ( const auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF) {
        tmp -= vertexDoFStencil[ vertexdof::stencilIndexFromHorizontalEdge( neighbor ) ] *
               dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( Level, i, neighbor )];
      }
    }

    dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge<Level>(i, stencilDirection::EDGE_HO_C)] =
      tmp / edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)];

  }
}

SPECIALIZE(void, smoothGSedgeDoFTmpl, smoothGSedgeDoF)

} /// namespace edge

namespace face {


template< uint_t Level>
void smoothGSvertexDoFTmpl(Face &face, const PrimitiveDataID<StencilMemory<real_t>, Face> &vertexDoFStencilID,
                            const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstVertexDoFID,
                            const PrimitiveDataID<StencilMemory<real_t>, Face> &edgeDoFStencilID,
                            const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstEdgeDoFID,
                            const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsVertexDoFID)
{
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
      tmp = rhs[vertexdof::macroface::indexFromVertex( Level, col, row, stencilDirection::VERTEX_C )];

      /// update from vertex dofs
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])] *
               dstVertexDoF[vertexdof::macroface::indexFromVertex( Level, col, row,
                                                                            vertexdof::macroface::neighborsWithoutCenter[k] )];
      }
      /// update from edge dofs
      for(uint_t k = 0; k < edgedof::macroface::neighborsFromVertex.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromVertex(edgedof::macroface::neighborsFromVertex[k])] *
               dstEdgeDoF[edgedof::macroface::indexFromVertex< Level >(col, row, edgedof::macroface::neighborsFromVertex[k])];
      }
      dstVertexDoF[vertexdof::macroface::indexFromVertex( Level, col, row,
                                                                   stencilDirection::VERTEX_C )] = tmp / vertexDoFStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C)];

    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, smoothGSvertexDoFTmpl, smoothGSvertexDoF)

template< uint_t Level>
void smoothGSedgeDoFTmpl(Face &face,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &vertexDoFStencilID,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstVertexDoFID,
                     const PrimitiveDataID<StencilMemory< real_t >, Face> &edgeDoFStencilID,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstEdgeDoFID,
                     const PrimitiveDataID<FunctionMemory< real_t >, Face> &rhsEdgeDoFID)
{
  real_t * vertexDoFStencil = face.getData(vertexDoFStencilID)->getPointer( Level );
  real_t * dstVertexDoF = face.getData(dstVertexDoFID)->getPointer( Level );
  real_t * edgeDoFStencil = face.getData(edgeDoFStencilID)->getPointer( Level );
  real_t * dstEdgeDoF = face.getData(dstEdgeDoFID)->getPointer( Level );
  real_t * rhs = face.getData(rhsEdgeDoFID)->getPointer( Level );

  real_t tmp;

  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0) {
      tmp = rhs[edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)];
      for(uint_t k = 1; k < edgedof::macroface::neighborsFromHorizontalEdge.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge(edgedof::macroface::neighborsFromHorizontalEdge[k])] *
               dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge< Level >(it.col(), it.row(), edgedof::macroface::neighborsFromHorizontalEdge[k])];
      }

      for(uint_t k = 0; k < vertexdof::macroface::neighborsFromHorizontalEdge.size(); ++k){
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromHorizontalEdge(vertexdof::macroface::neighborsFromHorizontalEdge[k])] *
               dstVertexDoF[vertexdof::macroface::indexFromHorizontalEdge< Level >(it.col(), it.row(), vertexdof::macroface::neighborsFromHorizontalEdge[k])];
      }

      dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge<Level>(it.col(), it.row(), stencilDirection::EDGE_HO_C)] =
        tmp / edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)];

    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      tmp = rhs[edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)];
      for(uint_t k = 1; k < edgedof::macroface::neighborsFromDiagonalEdge.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromDiagonalEdge(edgedof::macroface::neighborsFromDiagonalEdge[k])] *
               dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge< Level >(it.col(), it.row(), edgedof::macroface::neighborsFromDiagonalEdge[k])];
      }

      for(uint_t k = 0; k < vertexdof::macroface::neighborsFromDiagonalEdge.size(); ++k){
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromDiagonalEdge(vertexdof::macroface::neighborsFromDiagonalEdge[k])] *
          dstVertexDoF[vertexdof::macroface::indexFromDiagonalEdge< Level >(it.col(), it.row(), vertexdof::macroface::neighborsFromDiagonalEdge[k])];
      }

      dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge<Level>(it.col(), it.row(), stencilDirection::EDGE_DI_C)] =
        tmp / edgeDoFStencil[edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)];
    }
    if( it.col() != 0) {
      tmp = rhs[edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)];
      for(uint_t k = 1; k < edgedof::macroface::neighborsFromVerticalEdge.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromVerticalEdge(edgedof::macroface::neighborsFromVerticalEdge[k])] *
               dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge< Level >(it.col(), it.row(), edgedof::macroface::neighborsFromVerticalEdge[k])];
      }

      for(uint_t k = 0; k < vertexdof::macroface::neighborsFromVerticalEdge.size(); ++k){
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromVerticalEdge(vertexdof::macroface::neighborsFromVerticalEdge[k])] *
               dstVertexDoF[vertexdof::macroface::indexFromVerticalEdge< Level >(it.col(), it.row(), vertexdof::macroface::neighborsFromVerticalEdge[k])];
      }

      dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge<Level>(it.col(), it.row(), stencilDirection::EDGE_VE_C)] =
        tmp / edgeDoFStencil[edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)];
    }
  }
}

SPECIALIZE(void, smoothGSedgeDoFTmpl, smoothGSedgeDoF)

} /// namespace face

} /// namespace P2
} /// namespace hhg