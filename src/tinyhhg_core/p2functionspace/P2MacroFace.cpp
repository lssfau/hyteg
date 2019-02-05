#include "P2MacroFace.hpp"

#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
void smoothJacobiVertexDoF( const uint_t&                                            level,
                            const Face&                                              face,
                            const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexDoFStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcVertexDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstVertexDoFID,
                            const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeDoFStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcEdgeDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsVertexDoFID,
                            const real_t                                             dampingFactor )
{
   typedef stencilDirection sD;

   real_t* vertexDoFStencil = face.getData( vertexDoFStencilID )->getPointer( level );
   real_t* dstVertexDoF     = face.getData( dstVertexDoFID )->getPointer( level );
   real_t* srcVertexDoF     = face.getData( srcVertexDoFID )->getPointer( level );
   real_t* edgeDoFStencil   = face.getData( edgeDoFStencilID )->getPointer( level );
   real_t* srcEdgeDoF       = face.getData( srcEdgeDoFID )->getPointer( level );
   real_t* rhs              = face.getData( rhsVertexDoFID )->getPointer( level );

   real_t tmp;

   for( const auto& it : hhg::vertexdof::macroface::Iterator( level, 1 ) )
   {
      tmp = rhs[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )];

      /// update from vertex dofs
      for( auto k : vertexdof::macroface::neighborsWithoutCenter )
      {
         tmp -= vertexDoFStencil[vertexdof::stencilIndexFromVertex( k )] *
                srcVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), k )];
      }
      /// update from edge dofs
      for( auto k : edgedof::macroface::neighborsFromVertex )
      {
         tmp -= edgeDoFStencil[edgedof::stencilIndexFromVertex( k )] *
                srcEdgeDoF[edgedof::macroface::indexFromVertex( level, it.col(), it.row(), k )];
      }
      dstVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )] =
          dampingFactor * ( tmp / vertexDoFStencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] ) +
          ( 1.0 - dampingFactor ) *
              srcVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )];
   }
}

void smoothJacobiEdgeDoF( const uint_t&                                            Level,
                          const Face&                                              face,
                          const PrimitiveDataID< StencilMemory< real_t >, Face >&  vertexDoFStencilID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcVertexDoFID,
                          const PrimitiveDataID< StencilMemory< real_t >, Face >&  edgeDoFStencilID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcEdgeDoFID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstEdgeDoFID,
                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsEdgeDoFID,
                          const real_t                                             dampingFactor )
{
   typedef stencilDirection sD;

   real_t* vertexDoFStencil = face.getData( vertexDoFStencilID )->getPointer( Level );
   real_t* srcVertexDoF     = face.getData( srcVertexDoFID )->getPointer( Level );
   real_t* edgeDoFStencil   = face.getData( edgeDoFStencilID )->getPointer( Level );
   real_t* srcEdgeDoF       = face.getData( srcEdgeDoFID )->getPointer( Level );
   real_t* dstEdgeDoF       = face.getData( dstEdgeDoFID )->getPointer( Level );
   real_t* rhs              = face.getData( rhsEdgeDoFID )->getPointer( Level );

   real_t tmp;
   for( const auto& it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
   {
      if( it.row() != 0 )
      {
         tmp = rhs[edgedof::stencilIndexFromHorizontalEdge( sD::EDGE_HO_C )];
         for( uint_t k = 1; k < edgedof::macroface::neighborsFromHorizontalEdge.size(); ++k )
         {
            tmp -= edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge( edgedof::macroface::neighborsFromHorizontalEdge[k] )] *
                   srcEdgeDoF[edgedof::macroface::indexFromHorizontalEdge(
                       Level, it.col(), it.row(), edgedof::macroface::neighborsFromHorizontalEdge[k] )];
         }

         for( uint_t k = 0; k < vertexdof::macroface::neighborsFromHorizontalEdge.size(); ++k )
         {
            tmp -= vertexDoFStencil[vertexdof::stencilIndexFromHorizontalEdge(
                       vertexdof::macroface::neighborsFromHorizontalEdge[k] )] *
                   srcVertexDoF[vertexdof::macroface::indexFromHorizontalEdge(
                       Level, it.col(), it.row(), vertexdof::macroface::neighborsFromHorizontalEdge[k] )];
         }

         dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), sD::EDGE_HO_C )] =
             dampingFactor * ( tmp / edgeDoFStencil[edgedof::stencilIndexFromHorizontalEdge( sD::EDGE_HO_C )] ) +
             ( 1 - dampingFactor ) *
                 srcEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), sD::EDGE_HO_C )];
      }
      if( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         tmp = rhs[edgedof::stencilIndexFromDiagonalEdge( sD::EDGE_DI_C )];
         for( uint_t k = 1; k < edgedof::macroface::neighborsFromDiagonalEdge.size(); ++k )
         {
            tmp -= edgeDoFStencil[edgedof::stencilIndexFromDiagonalEdge( edgedof::macroface::neighborsFromDiagonalEdge[k] )] *
                   srcEdgeDoF[edgedof::macroface::indexFromDiagonalEdge(
                       Level, it.col(), it.row(), edgedof::macroface::neighborsFromDiagonalEdge[k] )];
         }

         for( uint_t k = 0; k < vertexdof::macroface::neighborsFromDiagonalEdge.size(); ++k )
         {
            tmp -=
                vertexDoFStencil[vertexdof::stencilIndexFromDiagonalEdge( vertexdof::macroface::neighborsFromDiagonalEdge[k] )] *
                srcVertexDoF[vertexdof::macroface::indexFromDiagonalEdge(
                    Level, it.col(), it.row(), vertexdof::macroface::neighborsFromDiagonalEdge[k] )];
         }

         dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), sD::EDGE_DI_C )] =
             dampingFactor * ( tmp / edgeDoFStencil[edgedof::stencilIndexFromDiagonalEdge( sD::EDGE_DI_C )] ) +
             ( 1.0 - dampingFactor ) *
                 srcEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), sD::EDGE_DI_C )];
      }
      if( it.col() != 0 )
      {
         tmp = rhs[edgedof::stencilIndexFromVerticalEdge( sD::EDGE_VE_C )];
         for( uint_t k = 1; k < edgedof::macroface::neighborsFromVerticalEdge.size(); ++k )
         {
            tmp -= edgeDoFStencil[edgedof::stencilIndexFromVerticalEdge( edgedof::macroface::neighborsFromVerticalEdge[k] )] *
                   srcEdgeDoF[edgedof::macroface::indexFromVerticalEdge(
                       Level, it.col(), it.row(), edgedof::macroface::neighborsFromVerticalEdge[k] )];
         }

         for( uint_t k = 0; k < vertexdof::macroface::neighborsFromVerticalEdge.size(); ++k )
         {
            tmp -=
                vertexDoFStencil[vertexdof::stencilIndexFromVerticalEdge( vertexdof::macroface::neighborsFromVerticalEdge[k] )] *
                srcVertexDoF[vertexdof::macroface::indexFromVerticalEdge(
                    Level, it.col(), it.row(), vertexdof::macroface::neighborsFromVerticalEdge[k] )];
         }

         dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), sD::EDGE_VE_C )] =
             dampingFactor * ( tmp / edgeDoFStencil[edgedof::stencilIndexFromVerticalEdge( sD::EDGE_VE_C )] ) +
             ( 1 - dampingFactor ) *
                 srcEdgeDoF[edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), sD::EDGE_VE_C )];
      }
   }
}


void smoothSOR(const uint_t &level,
                       const Face &face,
                       const real_t & relax,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &vertexToVertexStencilID,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &edgeToVertexStencilID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &vertexToEdgeStencilID,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &edgeToEdgeStencilID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsVertexDoFID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsEdgeDoFID)
{
   typedef stencilDirection sD;

   real_t* vertexToVertexStencil = face.getData( vertexToVertexStencilID )->getPointer( level );
   real_t* edgeToVertexStencil   = face.getData( edgeToVertexStencilID )->getPointer( level );
   real_t* dstVertexDoF          = face.getData( dstVertexDoFID )->getPointer( level );
   real_t* vertexToEdgeStencil   = face.getData( vertexToEdgeStencilID )->getPointer( level );
   real_t* edgeToEdgeStencil     = face.getData( edgeToEdgeStencilID )->getPointer( level );
   real_t* dstEdgeDoF            = face.getData( dstEdgeDoFID )->getPointer( level );
   real_t* rhsVertexDoF          = face.getData( rhsVertexDoFID )->getPointer( level );
   real_t* rhsEdgeDoF            = face.getData( rhsEdgeDoFID )->getPointer( level );

   real_t tmpVertex = 0, tmpEdgeHO = 0, tmpEdgeDI = 0, tmpEdgeVE = 0;

   // invert center weights
   const real_t invVertexCenter = 1.0 / vertexToVertexStencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )];
   const real_t invEdgeXCenter  = 1.0 / edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( sD::EDGE_HO_C )];
   const real_t invEdgeXYCenter = 1.0 / edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge( sD::EDGE_DI_C )];
   const real_t invEdgeYCenter  = 1.0 / edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge( sD::EDGE_VE_C )];

   /// sum up weighted values first for vertex and edges and write to corresponding dof
   for( const auto& it : hhg::edgedof::macroface::Iterator( level, 0 ) )
   {
      ////////// VERTEX //////////
      if( !vertexdof::macroface::isVertexOnBoundary( level, it ) )
      {
         tmpVertex = rhsVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )];
         /// vertex to vertex
         for( const auto& dir : vertexdof::macroface::neighborsWithoutCenter )
         {
            tmpVertex -= dstVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), dir )] *
                         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
         }
         /// edge to vertex
         for( const auto& dir : edgedof::macroface::neighborsFromVertex )
         {
            tmpVertex -= dstEdgeDoF[edgedof::macroface::indexFromVertex( level, it.col(), it.row(), dir )] *
                         edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
         }
         dstVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )] =
             (1.0 - relax) * dstVertexDoF[vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )] +
             relax * invVertexCenter * tmpVertex;
      }
      ////////// HORIZONTAL EDGE //////////
      if( !edgedof::isHorizontalEdgeOnBoundary( level, it ) )
      {
         tmpEdgeHO = rhsEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), sD::EDGE_HO_C )];
         /// vertex to edge
         for( const auto& dir : vertexdof::macroface::neighborsFromHorizontalEdge )
         {
            tmpEdgeHO -= dstVertexDoF[vertexdof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), dir )] *
                         vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
         }
         /// edge to edge
         for( const auto& dir : edgedof::macroface::neighborsFromHorizontalEdgeWithoutCenter )
         {
            tmpEdgeHO -= dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), dir )] *
                         edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( dir )];
         }
         dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), sD::EDGE_HO_C )] =
             (1.0 - relax) * dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), sD::EDGE_HO_C )] +
             relax * invEdgeXCenter * tmpEdgeHO;      
      }
      ////////// VERTICAL EDGE //////////
      if( !edgedof::isVerticalEdgeOnBoundary( level, it ) )
      {
         tmpEdgeVE = rhsEdgeDoF[edgedof::macroface::indexFromVerticalEdge( level, it.col(), it.row(), sD::EDGE_VE_C )];
         /// vertex to edge
         for( const auto& dir : vertexdof::macroface::neighborsFromVerticalEdge )
         {
            tmpEdgeVE -= dstVertexDoF[vertexdof::macroface::indexFromVerticalEdge( level, it.col(), it.row(), dir )] *
                         vertexToEdgeStencil[vertexdof::stencilIndexFromVerticalEdge( dir )];
         }
         /// edge to edge
         for( const auto& dir : edgedof::macroface::neighborsFromVerticalEdgeWithoutCenter )
         {
            tmpEdgeVE -= dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge( level, it.col(), it.row(), dir )] *
                         edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge( dir )];
         }
         dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge( level, it.col(), it.row(), sD::EDGE_VE_C )] =
             (1.0 - relax) * dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge( level, it.col(), it.row(), sD::EDGE_VE_C )] +
             relax * invEdgeYCenter * tmpEdgeVE;      
      }
      ////////// DIAGONAL EDGE //////////
      if( !edgedof::isDiagonalEdgeOnBoundary( level, it ) )
      {
         tmpEdgeDI = rhsEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( level, it.col(), it.row(), sD::EDGE_DI_C )];
         /// vertex to edge
         for( const auto& dir : vertexdof::macroface::neighborsFromDiagonalEdge )
         {
            tmpEdgeDI -= dstVertexDoF[vertexdof::macroface::indexFromDiagonalEdge( level, it.col(), it.row(), dir )] *
                         vertexToEdgeStencil[vertexdof::stencilIndexFromDiagonalEdge( dir )];
         }
         /// edge to edge
         for( const auto& dir : edgedof::macroface::neighborsFromDiagonalEdgeWithoutCenter )
         {
            tmpEdgeDI -= dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( level, it.col(), it.row(), dir )] *
                         edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge( dir )];
         }
         dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( level, it.col(), it.row(), sD::EDGE_DI_C )] =
             (1.0 - relax) * dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge( level, it.col(), it.row(), sD::EDGE_DI_C )] +
             relax * invEdgeXYCenter * tmpEdgeDI;
      }
   }
}

void smoothGaussSeidel(const uint_t &level,
                       const Face &face,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &vertexToVertexStencilID,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &edgeToVertexStencilID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstVertexDoFID,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &vertexToEdgeStencilID,
                       const PrimitiveDataID<StencilMemory<real_t>, Face> &edgeToEdgeStencilID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &dstEdgeDoFID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsVertexDoFID,
                       const PrimitiveDataID<FunctionMemory<real_t>, Face> &rhsEdgeDoFID)
{
  smoothSOR(level,
            face,
            1.0,
           vertexToVertexStencilID,
           edgeToVertexStencilID,
           dstVertexDoFID,
           vertexToEdgeStencilID,
           edgeToEdgeStencilID,
           dstEdgeDoFID,
           rhsVertexDoFID,
           rhsEdgeDoFID);
}

} // namespace macroface
} // namespace P2
} // namespace hhg
