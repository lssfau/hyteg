#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>
#include <tinyhhg_core/edgedofspace/EdgeDoFMacroEdge.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp>

#include "P2MacroFace.hpp"

namespace hhg {
namespace P2 {
namespace macroedge {

void smoothGaussSeidl( const uint_t&                                            level,
                       const Edge&                                              edge,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToVertexStencilID,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToVertexStencilID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToEdgeStencilID,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToEdgeStencilID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsVertexDoFID,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsEdgeDoFID )
{
   real_t* vertexToVertexStencil = edge.getData( vertexToVertexStencilID )->getPointer( level );
   real_t* edgeToVertexStencil   = edge.getData( edgeToVertexStencilID )->getPointer( level );
   real_t* dstVertexDoF          = edge.getData( dstVertexDoFID )->getPointer( level );
   real_t* vertexToEdgeStencil   = edge.getData( vertexToEdgeStencilID )->getPointer( level );
   real_t* edgeToEdgeStencil     = edge.getData( edgeToEdgeStencilID )->getPointer( level );
   real_t* dstEdgeDoF            = edge.getData( dstEdgeDoFID )->getPointer( level );
   real_t* rhsVertexDoF          = edge.getData( rhsVertexDoFID )->getPointer( level );
   real_t* rhsEdgeDoF            = edge.getData( rhsEdgeDoFID )->getPointer( level );

   real_t tmpVertex = 0, tmpEdgeHO = 0;
   for( const auto& it : hhg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      ////////// VERTEX //////////
      if( it.col() != 0 )
      {
         tmpVertex = rhsVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.col(), stencilDirection::VERTEX_C )];
         /// on edge vertex dof
         for( const auto& dir : hhg::vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
         {
            tmpVertex -= dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.col(), dir )] *
                         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
         }
         /// on edge edge dof
         for( const auto& dir : hhg::edgedof::macroedge::neighborsOnEdgeFromVertex )
         {
            tmpVertex -= dstEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.col(), dir )] *
                         edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
         }
         /// south face vertex dof
         for( const auto& dir : hhg::vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
         {
            tmpVertex -= dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.col(), dir )] *
                         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
         }
         /// south face edge
         for( const auto& dir : hhg::edgedof::macroedge::neighborsOnSouthFaceFromVertex )
         {
            tmpVertex -= dstEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.col(), dir )] *
                         edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
         }
         if( edge.getNumNeighborFaces() == 2 )
         {
            /// north face vertex dof
            for( const auto& dir : hhg::vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
            {
               tmpVertex -= dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.col(), dir )] *
                            vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
            }
            /// north face edge
            for( const auto& dir : hhg::edgedof::macroedge::neighborsOnNorthFaceFromVertex )
            {
               tmpVertex -= dstEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.col(), dir )] *
                            edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
            }
         }
      }
      ////////// HORIZONTAL EDGE //////////
      tmpEdgeHO = rhsEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.col(), stencilDirection::EDGE_HO_C )];
      /// on edge
      for( const auto& dir : hhg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF )
      {
         tmpEdgeHO -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.col(), dir )] *
                      vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
      }
      /// on south face
      for( const auto& dir : hhg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF )
      {
         tmpEdgeHO -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.col(), dir )] *
                      vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
      }

      for( const auto& dir : hhg::edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge )
      {
         tmpEdgeHO -= dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.col(), dir )] *
                      edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( dir )];
      }
      /// on north face
      if( edge.getNumNeighborFaces() == 2 )
      {
         for( const auto& dir : hhg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF )
         {
            tmpEdgeHO -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.col(), dir )] *
                         vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
         }
         for( const auto& dir : hhg::edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge )
         {
            tmpEdgeHO -= dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.col(), dir )] *
                         edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( dir )];
         }
      }
      ////////// VERTEX //////////
      if( it.col() != 0 )
      {
         dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.col(), stencilDirection::VERTEX_C )] =
             tmpVertex / vertexToVertexStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
      }
      dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.col(), stencilDirection::EDGE_HO_C )] =
          tmpEdgeHO / edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C )];
   }
}

} // namespace macroedge
} // namespace P2
} // namespace hhg
