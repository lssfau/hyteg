/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "P2MacroVertex.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitives/all.hpp"

using walberla::real_t;

namespace hyteg {
namespace P2 {

namespace macrovertex {

void smoothSORVertexDoF( uint_t                                                     level,
                         Vertex&                                                    vertex,
                         const real_t&                                              relax,
                         const PrimitiveDataID< StencilMemory< real_t >, Vertex >&  vertexDoFStencilID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstVertexDoFID,
                         const PrimitiveDataID< StencilMemory< real_t >, Vertex >&  edgeDoFStencilID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstEdgeDoFID,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& rhsVertexDoFID )
{
   real_t* vertexDoFStencil = vertex.getData( vertexDoFStencilID )->getPointer( level );
   real_t* dstVertexDoF     = vertex.getData( dstVertexDoFID )->getPointer( level );
   real_t* edgeDoFStencil   = vertex.getData( edgeDoFStencilID )->getPointer( level );
   real_t* dstEdgeDoF       = vertex.getData( dstEdgeDoFID )->getPointer( level );
   real_t* rhs              = vertex.getData( rhsVertexDoFID )->getPointer( level );

   real_t tmp = 0;
   tmp        = rhs[0];
   for( uint_t i = 0; i < vertex.getData( edgeDoFStencilID )->getSize( level ); ++i )
   {
      tmp -= dstEdgeDoF[i] * edgeDoFStencil[i];
   }
   for( uint_t i = 1; i < vertex.getData( vertexDoFStencilID )->getSize( level ); ++i )
   {
      tmp -= dstVertexDoF[i] * vertexDoFStencil[i];
   }

   dstVertexDoF[0] = (real_c(1.0) - relax) * dstVertexDoF[0] + (relax * tmp) / vertexDoFStencil[0];
}

void smoothSOR3D(
    const uint_t&                                                                                    level,
    const PrimitiveStorage&                                                                          storage,
    Vertex&                                                                                          vertex,
    const real_t&                                                                                    relax,
    const PrimitiveDataID< StencilMemory< real_t >, Vertex >&                                        vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex >& edgeToVertexOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                                       vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                                       vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                                       edgeDoFDstId )
{
   using edgedof::EdgeDoFOrientation;
   using indexing::IndexIncrement;

   auto v2v_operator = vertex.getData( vertexToVertexOperatorId )->getPointer( level );
   auto e2v_operator = vertex.getData( edgeToVertexOperatorId )->getData( level );

   real_t* vertexDoFDst = vertex.getData( vertexDoFDstId )->getPointer( level );
   real_t* vertexDoFRhs = vertex.getData( vertexDoFRhsId )->getPointer( level );
   real_t* edgeDoFDst   = vertex.getData( edgeDoFDstId )->getPointer( level );

   const real_t vertexDoFRelaxOverCenter = relax / v2v_operator[0];
   const real_t oneMinusRelax            = real_c( 1 ) - relax;

   real_t tmp = vertexDoFRhs[0];

   // vertex leaves
   for ( size_t i = 1; i < vertex.getNumNeighborEdges() + 1; i++ )
   {
      tmp -= v2v_operator[i] * vertexDoFDst[i];
   }

   // edge leaves
   const auto centerIndexOnVertex = indexing::Index( 0, 0, 0 );

   for ( uint_t neighborCellID = 0; neighborCellID < vertex.getNumNeighborCells(); neighborCellID++ )
   {
      const Cell& neighborCell      = *( storage.getCell( vertex.neighborCells().at( neighborCellID ) ) );
      auto        cellLocalVertexID = neighborCell.getLocalVertexID( vertex.getID() );

      const auto basisInCell = algorithms::getMissingIntegersAscending< 1, 4 >( {cellLocalVertexID} );

      const auto centerIndexInCell = indexing::basisConversion(
          centerIndexOnVertex, basisInCell, {0, 1, 2, 3}, levelinfo::num_microvertices_per_edge( level ) );

      for ( const auto& leafOrientationInCell : edgedof::allEdgeDoFOrientationsWithoutXYZ )
      {
         for ( const auto& stencilIt : e2v_operator[neighborCellID][leafOrientationInCell] )
         {
            const auto stencilOffset = stencilIt.first;
            const auto stencilWeight = stencilIt.second;

            const auto leafIndexInCell = centerIndexInCell + stencilOffset.cast< idx_t >();

            const auto onCellFacesSet = edgedof::macrocell::isOnCellFaces( level, leafIndexInCell, leafOrientationInCell );
            const auto onCellEdgesSet = edgedof::macrocell::isOnCellEdges( level, leafIndexInCell, leafOrientationInCell );

            uint_t leafArrayIndexOnVertex = std::numeric_limits< uint_t >::max();

            WALBERLA_ASSERT_GREATER( onCellFacesSet.size(), 0 );

            if ( onCellFacesSet.size() == 1 )
            {
               // on macro-face
               WALBERLA_ASSERT_EQUAL( onCellEdgesSet.size(), 0 );
               const auto faceID            = neighborCell.neighborFaces().at( *onCellFacesSet.begin() );
               const auto vertexLocalFaceID = vertex.face_index( faceID );
               leafArrayIndexOnVertex       = vertex.getNumNeighborEdges() + vertexLocalFaceID;
            }
            else
            {
               // on macro-edge
               WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), 2 );
               WALBERLA_ASSERT_EQUAL( onCellEdgesSet.size(), 1 );
               const auto edgeID            = neighborCell.neighborEdges().at( *onCellEdgesSet.begin() );
               const auto vertexLocalEdgeID = vertex.edge_index( edgeID );
               leafArrayIndexOnVertex       = vertexLocalEdgeID;
            }

            tmp -= edgeDoFDst[leafArrayIndexOnVertex] * stencilWeight;
         }
      }
   }

   vertexDoFDst[0] = oneMinusRelax * vertexDoFDst[0] + vertexDoFRelaxOverCenter * tmp;
}

void smoothJacobiVertexDoF( uint_t                                                     level,
                            Vertex&                                                    vertex,
                            const real_t&                                              relax,
                            const PrimitiveDataID< StencilMemory < real_t >, Vertex >& vertexToVertexStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& srcVertexDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstVertexDoFID,
                            const PrimitiveDataID< StencilMemory < real_t >, Vertex >& edgeToVertexStencilID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& srcEdgeDoFID,
                            const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& rhsVertexDoFID )
{
   real_t* vertexDoFStencil = vertex.getData( vertexToVertexStencilID )->getPointer( level );
   real_t* edgeDoFStencil   = vertex.getData( edgeToVertexStencilID   )->getPointer( level );

   real_t* srcVertexDoF     = vertex.getData( srcVertexDoFID )->getPointer( level );
   real_t* dstVertexDoF     = vertex.getData( dstVertexDoFID )->getPointer( level );
   real_t* rhs              = vertex.getData( rhsVertexDoFID )->getPointer( level );

   real_t* srcEdgeDoF       = vertex.getData( srcEdgeDoFID )->getPointer( level );

   real_t tmp = 0;
   tmp        = rhs[0];
   for( uint_t i = 0; i < vertex.getData( edgeToVertexStencilID )->getSize( level ); ++i )
   {
      tmp -= srcEdgeDoF[i] * edgeDoFStencil[i];
   }
   for( uint_t i = 1; i < vertex.getData( vertexToVertexStencilID )->getSize( level ); ++i )
   {
      tmp -= srcVertexDoF[i] * vertexDoFStencil[i];
   }

   dstVertexDoF[0] = (real_c(1.0) - relax) * srcVertexDoF[0] + (relax * tmp) / vertexDoFStencil[0];
}

} // namespace vertex

} // namespace P2
} // namespace hyteg
