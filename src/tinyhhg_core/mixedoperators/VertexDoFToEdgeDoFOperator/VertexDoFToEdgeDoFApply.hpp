#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements3D.hpp"

namespace hhg{
namespace VertexDoFToEdgeDoF{


inline void applyEdge(const uint_t & Level, Edge &edge,
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


inline void applyFace(const uint_t & Level, Face &face,
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


inline void applyCell(const uint_t & Level, Cell & cell,
                      const PrimitiveDataID<StencilMemory < real_t >, Cell> &operatorId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Cell> &srcId,
                      const PrimitiveDataID<FunctionMemory< real_t >, Cell> &dstId,
                      UpdateType update){

  real_t * opr_data = cell.getData(operatorId)->getPointer( Level );
  real_t * src      = cell.getData(srcId)->getPointer( Level );
  real_t * dst      = cell.getData(dstId)->getPointer( Level );

  for ( const auto & it : hhg::edgedof::macrocell::Iterator( Level, 0 ) )
  {
    std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

    if ( edgedof::macrocell::isInnerXEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
    if ( edgedof::macrocell::isInnerYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
    if ( edgedof::macrocell::isInnerZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
    if ( edgedof::macrocell::isInnerXYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
    if ( edgedof::macrocell::isInnerXZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
    if ( edgedof::macrocell::isInnerYZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

    for ( const auto & centerOrientation : innerOrientations )
    {
      real_t tmp = 0.0;

      const auto vertexDoFNeighbors = P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );
      for ( const auto & neighbor : vertexDoFNeighbors )
      {
        const uint_t stencilIdx  = vertexdof::stencilIndexFromEdge3D( vertexdof::stencilDirectionFromLogicalOffset( neighbor ), centerOrientation );
        const auto   srcIdx      = it + neighbor;
        const auto   srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
        tmp += opr_data[stencilIdx] * src[srcArrayIdx];
      }

      const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );

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

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    real_t tmp = 0.0;
    const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;

    const auto vertexDoFNeighbors = P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );
    for ( const auto & neighbor : vertexDoFNeighbors )
    {
      const uint_t stencilIdx  = vertexdof::stencilIndexFromEdge3D( vertexdof::stencilDirectionFromLogicalOffset( neighbor ), centerOrientation );
      const auto   srcIdx      = it + neighbor;
      const auto   srcArrayIdx = vertexdof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
      tmp += opr_data[stencilIdx] * src[srcArrayIdx];
    }

    const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );

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


} // VertexDoFToEdgeDoFOperator
} // namespace hhg
