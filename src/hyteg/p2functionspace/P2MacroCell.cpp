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

#include "hyteg/p2functionspace/P2MacroCell.hpp"

#include "hyteg/indexing/Common.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "hyteg/p2functionspace/P2Elements3D.hpp"
#include "hyteg/primitives/Cell.hpp"

namespace hyteg {
namespace P2 {
namespace macrocell {

using edgedof::EdgeDoFOrientation;
using indexing::Index;

void smoothSOR(
    const uint_t&                                                                                level,
    Cell&                                                                                        cell,
    const real_t&                                                                                relax,
    const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >&        vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     edgeDoFRhsId )
{
   auto v2v_operator = cell.getData( vertexToVertexOperatorId )->getData( level );
   auto e2v_operator = cell.getData( edgeToVertexOperatorId )->getData( level );
   auto v2e_operator = cell.getData( vertexToEdgeOperatorId )->getData( level );
   auto e2e_operator = cell.getData( edgeToEdgeOperatorId )->getData( level );

   real_t* vertexDoFDst = cell.getData( vertexDoFDstId )->getPointer( level );
   real_t* vertexDoFRhs = cell.getData( vertexDoFRhsId )->getPointer( level );
   real_t* edgeDoFDst   = cell.getData( edgeDoFDstId )->getPointer( level );
   real_t* edgeDoFRhs   = cell.getData( edgeDoFRhsId )->getPointer( level );

   const real_t vertexDoFRelaxOverCenter = relax / v2v_operator[ { 0, 0, 0 } ];
   const real_t oneMinusRelax            = real_c( 1 ) - relax;

   std::map< EdgeDoFOrientation, real_t > edgeDoFRelaxOverCenter;
   edgeDoFRelaxOverCenter[EdgeDoFOrientation::X] =
       relax / e2e_operator[EdgeDoFOrientation::X][EdgeDoFOrientation::X][Index( 0, 0, 0 )];
   edgeDoFRelaxOverCenter[EdgeDoFOrientation::Y] =
       relax / e2e_operator[EdgeDoFOrientation::Y][EdgeDoFOrientation::Y][Index( 0, 0, 0 )];
   edgeDoFRelaxOverCenter[EdgeDoFOrientation::Z] =
       relax / e2e_operator[EdgeDoFOrientation::Z][EdgeDoFOrientation::Z][Index( 0, 0, 0 )];
   edgeDoFRelaxOverCenter[EdgeDoFOrientation::XY] =
       relax / e2e_operator[EdgeDoFOrientation::XY][EdgeDoFOrientation::XY][Index( 0, 0, 0 )];
   edgeDoFRelaxOverCenter[EdgeDoFOrientation::XZ] =
       relax / e2e_operator[EdgeDoFOrientation::XZ][EdgeDoFOrientation::XZ][Index( 0, 0, 0 )];
   edgeDoFRelaxOverCenter[EdgeDoFOrientation::YZ] =
       relax / e2e_operator[EdgeDoFOrientation::YZ][EdgeDoFOrientation::YZ][Index( 0, 0, 0 )];
   edgeDoFRelaxOverCenter[EdgeDoFOrientation::XYZ] =
       relax / e2e_operator[EdgeDoFOrientation::XYZ][EdgeDoFOrientation::XYZ][Index( 0, 0, 0 )];

   real_t tmp;

   // update vertex unknowns
   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const auto dstArrayIdx = vertexdof::macrocell::index( level, it.x(), it.y(), it.z() );
      tmp                    = vertexDoFRhs[dstArrayIdx];

      // vertex leaves
      for ( const auto& vertexLeaf : vertexdof::macrocell::neighborsWithoutCenter )
      {
         const auto idx         = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), vertexLeaf );
         tmp -= v2v_operator[ vertexdof::logicalIndexOffsetFromVertex( vertexLeaf ) ] * vertexDoFDst[idx];
      }

      // edge leaves
      for ( const auto& orientation : edgedof::allEdgeDoFOrientations )
      {
         const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromVertexDoFInMacroCell( orientation );
         for ( const auto& neighbor : edgeDoFNeighbors )
         {
            const auto srcIdx      = it + neighbor;
            const auto srcArrayIdx = edgedof::macrocell::index( level, srcIdx.x(), srcIdx.y(), srcIdx.z(), orientation );
            tmp -= e2v_operator[orientation][neighbor] * edgeDoFDst[srcArrayIdx];
         }
      }

      vertexDoFDst[dstArrayIdx] = oneMinusRelax * vertexDoFDst[dstArrayIdx] + vertexDoFRelaxOverCenter * tmp;
   }

   // update edge unknowns part I (all but XYZ)
   for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) )
   {
      std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

      for ( const auto& orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
      {
         if ( edgedof::macrocell::isInnerEdgeDoF( level, it, orientation ) )
            innerOrientations.push_back( orientation );
      }

      for ( const auto& centerOrientation : innerOrientations )
      {
         const auto dstArrayIdx = edgedof::macrocell::index( level, it.x(), it.y(), it.z(), centerOrientation );
         tmp                    = edgeDoFRhs[dstArrayIdx];

         // vertex leaves
         const auto vertexDoFNeighbors =
             P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );
         for ( const auto& neighbor : vertexDoFNeighbors )
         {
            const auto srcIdx      = it + neighbor;
            const auto srcArrayIdx = vertexdof::macrocell::index( level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
            tmp -= v2e_operator[centerOrientation][neighbor] * vertexDoFDst[srcArrayIdx];
         }

         // edge leaves
         for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
         {
            const auto edgeDoFNeighbors =
                P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
            for ( const auto& neighbor : edgeDoFNeighbors )
            {
               // skip center
               if ( centerOrientation == leafOrientation && neighbor == Index( 0, 0, 0 ) )
                  continue;

               const auto   srcIdx      = it + neighbor;
               const auto   srcArrayIdx = edgedof::macrocell::index( level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
               const real_t stencilWeight = e2e_operator[centerOrientation][leafOrientation][neighbor];
               tmp -= stencilWeight * edgeDoFDst[srcArrayIdx];
            }
         }

         edgeDoFDst[dstArrayIdx] = oneMinusRelax * edgeDoFDst[dstArrayIdx] + edgeDoFRelaxOverCenter[centerOrientation] * tmp;
      }
   }

   // update edge unknowns part II (only XYZ)
   for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) )
   {
      const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;
      const auto dstArrayIdx       = edgedof::macrocell::index( level, it.x(), it.y(), it.z(), centerOrientation );
      tmp                          = edgeDoFRhs[dstArrayIdx];

      // vertex leaves
      const auto vertexDoFNeighbors =
          P2Elements::P2Elements3D::getAllVertexDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation );
      for ( const auto& neighbor : vertexDoFNeighbors )
      {
         const auto srcIdx      = it + neighbor;
         const auto srcArrayIdx = vertexdof::macrocell::index( level, srcIdx.x(), srcIdx.y(), srcIdx.z() );
         tmp -= v2e_operator[centerOrientation][neighbor] * vertexDoFDst[srcArrayIdx];
      }

      // edge leaves
      for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
      {
         const auto edgeDoFNeighbors =
             P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
         for ( const auto& neighbor : edgeDoFNeighbors )
         {
            // skip center
            if ( centerOrientation == leafOrientation && neighbor == Index( 0, 0, 0 ) )
               continue;

            const auto   srcIdx        = it + neighbor;
            const auto   srcArrayIdx   = edgedof::macrocell::index( level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
            const real_t stencilWeight = e2e_operator[centerOrientation][leafOrientation][neighbor];
            tmp -= stencilWeight * edgeDoFDst[srcArrayIdx];
         }
      }

      edgeDoFDst[dstArrayIdx] = oneMinusRelax * edgeDoFDst[dstArrayIdx] + edgeDoFRelaxOverCenter[centerOrientation] * tmp;
   }
}

} // namespace macrocell
} // namespace P2
} // namespace hyteg
