#include "P2MacroFace.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

namespace hhg{
namespace P2{
namespace macroface{

void smoothJacobiVertexDoF(const uint_t & level, const Face &face,
                           const PrimitiveDataID<StencilMemory<real_t>, Face> &vertexDoFStencilID,
                           const PrimitiveDataID<FunctionMemory<real_t>, Face> &srcVertexDoFID,
                           const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstVertexDoFID,
                           const PrimitiveDataID<StencilMemory<real_t>, Face> &edgeDoFStencilID,
                           const PrimitiveDataID<FunctionMemory<real_t>, Face> &srcEdgeDoFID,
                           const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsVertexDoFID){
///TODO: remove hardcoded damping factor
  real_t * vertexDoFStencil = face.getData(vertexDoFStencilID)->getPointer( level );
  real_t * dstVertexDoF = face.getData(dstVertexDoFID)->getPointer( level );
  real_t * srcVertexDoF = face.getData(srcVertexDoFID)->getPointer( level );
  real_t * edgeDoFStencil = face.getData(edgeDoFStencilID)->getPointer( level );
  real_t * srcEdgeDoF = face.getData(srcEdgeDoFID)->getPointer( level );
  real_t * rhs = face.getData(rhsVertexDoFID)->getPointer( level );

//  vertexDoFStencil[0] = 1.0/3.0;
//  vertexDoFStencil[1] = 0.0;
//  vertexDoFStencil[2] = 1.0/3.0;
//  vertexDoFStencil[3] = 4.0;
//  vertexDoFStencil[4] = 1.0/3.0;
//  vertexDoFStencil[5] = 0.0;
//  vertexDoFStencil[6] = 1.0/3.0;
//
//  edgeDoFStencil[0] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[1] = 0.0;
//  edgeDoFStencil[2] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[3] = 0.0;
//  edgeDoFStencil[4] = 0.0;
//  edgeDoFStencil[5] = 0.0;
//  edgeDoFStencil[6] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[7] = 0.0;
//  edgeDoFStencil[8] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[9] = 0.0;
//  edgeDoFStencil[10] = 0.0;
//  edgeDoFStencil[11] = 0.0;

  real_t tmp;


  for ( const auto & it : hhg::vertexdof::macroface::Iterator( level, 1 ) ){
      tmp = rhs[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), stencilDirection::VERTEX_C )];

      /// update from vertex dofs
      for (auto k : vertexdof::macroface::neighborsWithoutCenter)
      {
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromVertex(k)] *
               srcVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(),
                                                                   k)];
      }
      /// update from edge dofs
      for (auto k : edgedof::macroface::neighborsFromVertex)
      {
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromVertex(k)] *
               srcEdgeDoF[edgedof::macroface::indexFromVertex( level, it.col(), it.row(), k)];
      }
      dstVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(),
                                                          stencilDirection::VERTEX_C )] = 0.66 * ( tmp / vertexDoFStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C)]);

  }
}


void smoothJacobiEdgeDoF(const uint_t & Level, const Face &face,
                         const PrimitiveDataID<StencilMemory<real_t>, Face> &vertexDoFStencilID,
                         const PrimitiveDataID<FunctionMemory<real_t>, Face> &srcVertexDoFID,
                         const PrimitiveDataID<StencilMemory<real_t>, Face> &edgeDoFStencilID,
                         const PrimitiveDataID<FunctionMemory<real_t>, Face> &srcEdgeDoFID,
                         const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstEdgeDoFID,
                         const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsEdgeDoFID){
///TODO: remove hardcoded damping factor
  real_t * vertexDoFStencil = face.getData(vertexDoFStencilID)->getPointer( Level );
  real_t * srcVertexDoF = face.getData(srcVertexDoFID)->getPointer( Level );
  real_t * edgeDoFStencil = face.getData(edgeDoFStencilID)->getPointer( Level );
  real_t * srcEdgeDoF = face.getData(srcEdgeDoFID)->getPointer( Level );
  real_t * dstEdgeDoF = face.getData(dstEdgeDoFID)->getPointer( Level );
  real_t * rhs = face.getData(rhsEdgeDoFID)->getPointer( Level );
//
//  vertexDoFStencil[0] = - (1.0 + 1.0/3.0);
//  vertexDoFStencil[1] = - (1.0 + 1.0/3.0);
//  vertexDoFStencil[2] = 0.0;
//  vertexDoFStencil[3] = 0.0;
//  vertexDoFStencil[4] = 0.0;
//  vertexDoFStencil[5] = 0.0;
//  vertexDoFStencil[6] = 0.0;
//  vertexDoFStencil[7] = 0.0;
//  vertexDoFStencil[8] = - (1.0 + 1.0/3.0);
//  vertexDoFStencil[9] = 0.0;
//  vertexDoFStencil[10] = - (1.0 + 1.0/3.0);;
//  vertexDoFStencil[11] = 0.0;
//
//  edgeDoFStencil[0] = (5.0 + 1.0/3.0);
//  edgeDoFStencil[1] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[2] = 0.0;
//  edgeDoFStencil[3] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[4] = 0.0;
//  edgeDoFStencil[5] = (5.0 + 1.0/3.0);
//  edgeDoFStencil[6] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[7] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[8] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[9] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[10] = (5.0 + 1.0/3.0);
//  edgeDoFStencil[11] = 0.0;
//  edgeDoFStencil[12] = - (1.0 + 1.0/3.0);
//  edgeDoFStencil[13] = 0.0;
//  edgeDoFStencil[14] = - (1.0 + 1.0/3.0);


  real_t tmp;
  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0) {
      tmp = rhs[edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)];
      for(uint_t k = 1; k < edgedof::macroface::neighborsFromHorizontalEdge.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge(edgedof::macroface::neighborsFromHorizontalEdge[k])] *
               srcEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), edgedof::macroface::neighborsFromHorizontalEdge[k] )];
      }

      for(uint_t k = 0; k < vertexdof::macroface::neighborsFromHorizontalEdge.size(); ++k){
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromHorizontalEdge(vertexdof::macroface::neighborsFromHorizontalEdge[k])] *
               srcVertexDoF[vertexdof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), vertexdof::macroface::neighborsFromHorizontalEdge[k] )];
      }

      dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] =
          0.66 * (tmp / edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)]);

    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      tmp = rhs[edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)];
      for(uint_t k = 1; k < edgedof::macroface::neighborsFromDiagonalEdge.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromDiagonalEdge(edgedof::macroface::neighborsFromDiagonalEdge[k])] *
               srcEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), edgedof::macroface::neighborsFromDiagonalEdge[k] )];
      }

      for(uint_t k = 0; k < vertexdof::macroface::neighborsFromDiagonalEdge.size(); ++k){
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromDiagonalEdge(vertexdof::macroface::neighborsFromDiagonalEdge[k])] *
               srcVertexDoF[vertexdof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), vertexdof::macroface::neighborsFromDiagonalEdge[k] )];
      }

      dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] =
          0.66 * (tmp / edgeDoFStencil[edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)]);
    }
    if( it.col() != 0) {
      tmp = rhs[edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)];
      for(uint_t k = 1; k < edgedof::macroface::neighborsFromVerticalEdge.size(); ++k){
        tmp -= edgeDoFStencil[edgedof::stencilIndexFromVerticalEdge(edgedof::macroface::neighborsFromVerticalEdge[k])] *
               srcEdgeDoF[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), edgedof::macroface::neighborsFromVerticalEdge[k] )];
      }

      for(uint_t k = 0; k < vertexdof::macroface::neighborsFromVerticalEdge.size(); ++k){
        tmp -= vertexDoFStencil[vertexdof::stencilIndexFromVerticalEdge(vertexdof::macroface::neighborsFromVerticalEdge[k])] *
               srcVertexDoF[vertexdof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), vertexdof::macroface::neighborsFromVerticalEdge[k] )];
      }

      dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] =
          0.66 * (tmp / edgeDoFStencil[edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)]);
    }
  }
}

}
}
}

