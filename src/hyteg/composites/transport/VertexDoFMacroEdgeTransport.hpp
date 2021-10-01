/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#pragma once

#include "core/DataTypes.h"
#include "core/debug/all.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {
namespace vertexdof {
namespace transport {
namespace macroedge {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using indexing::Index;

template < typename ValueType, bool AlgebraicUpwind >
inline void apply( const uint_t&                                               level,
                   Edge&                                                       edge,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& uxId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& uyId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& uzId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Edge >&  xOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Edge >&  yOprId,
                   const PrimitiveDataID< StencilMemory< ValueType >, Edge >&  zOprId )
{
   typedef stencilDirection sD;
   size_t                   rowsize = levelinfo::num_microvertices_per_edge( level );

   const ValueType* src = edge.getData( srcId )->getPointer( level );
   ValueType*       dst = edge.getData( dstId )->getPointer( level );

   const ValueType* ux = edge.getData( uxId )->getPointer( level );
   const ValueType* uy = edge.getData( uyId )->getPointer( level );
   const ValueType* uz = edge.getData( uzId )->getPointer( level );

   const ValueType* xOperatorData = edge.getData( xOprId )->getPointer( level );
   const ValueType* yOperatorData = edge.getData( yOprId )->getPointer( level );
   const ValueType* zOperatorData = edge.getData( zOprId )->getPointer( level );

   std::vector< ValueType > stencil( edge.getData( xOprId )->getSize( level ) );
   real_t                   dTmp;

   ValueType tmp;

   for( size_t i = 1; i < rowsize - 1; ++i )
   {
      const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
      const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
      const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

      const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
      const auto dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
      const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

      // fill stencil
      stencil[stencilIdxC] = 0.0;

      stencil[stencilIdxW] = 0.5 * ( ux[dofIdxC] + ux[dofIdxW] ) * xOperatorData[stencilIdxW];
      stencil[stencilIdxW] += 0.5 * ( uy[dofIdxC] + uy[dofIdxW] ) * yOperatorData[stencilIdxW];
      stencil[stencilIdxW] += 0.5 * ( uz[dofIdxC] + uz[dofIdxW] ) * zOperatorData[stencilIdxW];
      stencil[stencilIdxC] -= stencil[stencilIdxW];

      // algebraic upwind
      if( AlgebraicUpwind )
      {
         dTmp = std::abs( stencil[stencilIdxW] );
         stencil[stencilIdxC] += dTmp;
         stencil[stencilIdxW] -= dTmp;
      }

      stencil[stencilIdxE] = 0.5 * ( ux[dofIdxC] + ux[dofIdxE] ) * xOperatorData[stencilIdxE];
      stencil[stencilIdxE] += 0.5 * ( uy[dofIdxC] + uy[dofIdxE] ) * yOperatorData[stencilIdxE];
      stencil[stencilIdxE] += 0.5 * ( uz[dofIdxC] + uz[dofIdxE] ) * zOperatorData[stencilIdxE];
      stencil[stencilIdxC] -= stencil[stencilIdxE];

      // algebraic upwind
      if( AlgebraicUpwind )
      {
         dTmp = std::abs( stencil[stencilIdxE] );
         stencil[stencilIdxC] += dTmp;
         stencil[stencilIdxE] -= dTmp;
      }

      for( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
      {
         const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
         const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
         const auto dofIdxWNeighborFace     = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
         const auto dofIdxENeighborFace     = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
         stencil[stencilIdxWNeighborFace]   = 0.5 * ( ux[dofIdxC] + ux[dofIdxWNeighborFace] ) * xOperatorData[stencilIdxWNeighborFace];
         stencil[stencilIdxWNeighborFace] += 0.5 * ( uy[dofIdxC] + uy[dofIdxWNeighborFace] ) * yOperatorData[stencilIdxWNeighborFace];
         stencil[stencilIdxWNeighborFace] += 0.5 * ( uz[dofIdxC] + uz[dofIdxWNeighborFace] ) * zOperatorData[stencilIdxWNeighborFace];
         stencil[stencilIdxC] -= stencil[stencilIdxWNeighborFace];

         // algebraic upwind
         if( AlgebraicUpwind )
         {
            dTmp = std::abs( stencil[stencilIdxWNeighborFace] );
            stencil[stencilIdxC] += dTmp;
            stencil[stencilIdxWNeighborFace] -= dTmp;
         }

         stencil[stencilIdxENeighborFace] = 0.5 * ( ux[dofIdxC] + ux[dofIdxENeighborFace] ) * xOperatorData[stencilIdxENeighborFace];
         stencil[stencilIdxENeighborFace] += 0.5 * ( uy[dofIdxC] + uy[dofIdxENeighborFace] ) * yOperatorData[stencilIdxENeighborFace];
         stencil[stencilIdxENeighborFace] += 0.5 * ( uz[dofIdxC] + uz[dofIdxENeighborFace] ) * zOperatorData[stencilIdxENeighborFace];
         stencil[stencilIdxC] -= stencil[stencilIdxENeighborFace];

         // algebraic upwind
         if( AlgebraicUpwind )
         {
            dTmp = std::abs( stencil[stencilIdxENeighborFace] );
            stencil[stencilIdxC] += dTmp;
            stencil[stencilIdxENeighborFace] -= dTmp;
         }
      }

      for( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
      {
         const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
         const auto dofIdx =
             vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
         stencil[stencilIdx] = 0.5 * ( ux[dofIdxC] + ux[dofIdx] ) * xOperatorData[stencilIdx];
         stencil[stencilIdx] += 0.5 * ( uy[dofIdxC] + uy[dofIdx] ) * yOperatorData[stencilIdx];
         stencil[stencilIdx] += 0.5 * ( uz[dofIdxC] + uz[dofIdx] ) * zOperatorData[stencilIdx];
         stencil[stencilIdxC] -= stencil[stencilIdx];

         // algebraic upwind
         if( AlgebraicUpwind )
         {
            dTmp = std::abs( stencil[stencilIdx] );
            stencil[stencilIdxC] += dTmp;
            stencil[stencilIdx] -= dTmp;
         }
      }

      // apply stencil
      tmp = stencil[stencilIdxW] * src[dofIdxW] + stencil[stencilIdxC] * src[dofIdxC] + stencil[stencilIdxE] * src[dofIdxE];

      for( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
      {
         const auto stencilIdxWNeighborFace    = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
         const auto stencilIdxENeighborFace    = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
         const auto stencilWeightW = stencil[stencilIdxWNeighborFace];
         const auto stencilWeightE = stencil[stencilIdxENeighborFace];
         const auto dofIdxWNeighborFace        = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
         const auto dofIdxENeighborFace        = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
         tmp += stencilWeightW * src[dofIdxWNeighborFace] + stencilWeightE * src[dofIdxENeighborFace];
      }

      for( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
      {
         const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
         const auto dofIdx =
             vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
         tmp += stencil[stencilIdx] * src[dofIdx];
      }

      dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = tmp;
   }
}

} // namespace macroedge
} // namespace transport
} // namespace vertexdof
} // namespace hyteg
