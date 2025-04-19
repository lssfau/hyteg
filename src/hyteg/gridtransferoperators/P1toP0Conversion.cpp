
/*
 * Copyright (c) 2025 Ponsuganth Ilangovan.
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

#include "hyteg/gridtransferoperators/P1toP0Conversion.hpp"

#include "hyteg/averaging/CellAveraging.hpp"
#include "hyteg/averaging/FaceAveraging.hpp"

namespace hyteg {

void P1toP0Conversion( const P1Function< real_t >& src,
                       const P0Function< real_t >& dst,
                       const uint_t                level,
                       const hyteg::AveragingType  averagingType )

{
   if ( dst.getStorage()->hasGlobalCells() )
   {
      for ( auto it : dst.getStorage()->getCells() )
      {
         PrimitiveID cellId = it.first;
         Cell&       cell   = *( it.second );

         const auto p0DofMemory = dst.getDGFunction()->volumeDoFFunction()->dofMemory( cellId, level );

         auto       p1FuncId   = src.getCellDataID();
         const auto p1FuncData = cell.getData( p1FuncId )->getPointer( level );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               uint_t p0DofIdx = volumedofspace::indexing::index(
                   idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA );

               const std::array< indexing::Index, 4 > vertexIndices =
                   celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );

               auto microTet0 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[0] );
               auto microTet1 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[1] );
               auto microTet2 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[2] );
               auto microTet3 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[3] );

               auto valueTet0 = p1FuncData[vertexdof::macrocell::index(
                   level, vertexIndices[0].x(), vertexIndices[0].y(), vertexIndices[0].z() )];
               auto valueTet1 = p1FuncData[vertexdof::macrocell::index(
                   level, vertexIndices[1].x(), vertexIndices[1].y(), vertexIndices[1].z() )];
               auto valueTet2 = p1FuncData[vertexdof::macrocell::index(
                   level, vertexIndices[2].x(), vertexIndices[2].y(), vertexIndices[2].z() )];
               auto valueTet3 = p1FuncData[vertexdof::macrocell::index(
                   level, vertexIndices[3].x(), vertexIndices[3].y(), vertexIndices[3].z() )];

               real_t sampledAverage = evaluateSampledAverage( { microTet0, microTet1, microTet2, microTet3 },
                                                               { valueTet0, valueTet1, valueTet2, valueTet3 },
                                                               averagingType );

               p0DofMemory[p0DofIdx] = sampledAverage;
            }
         }
      }
   }
   else
   {
      for ( auto& it : dst.getStorage()->getFaces() )
      {
         PrimitiveID faceId = it.first;
         Face&       face   = *( it.second );

         const auto p0DofMemory = dst.getDGFunction()->volumeDoFFunction()->dofMemory( faceId, level );

         auto       p1FuncId   = src.getFaceDataID();
         const auto p1FuncData = face.getData( p1FuncId )->getPointer( level );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               uint_t p0DofIdx = volumedofspace::indexing::index(
                   idxIt.x(), idxIt.y(), faceType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA );

               const std::array< indexing::Index, 3 > vertexIndices =
                   facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
               std::array< Point3D, 3 > elementVertices;
               for ( uint_t i = 0; i < 3; i++ )
               {
                  const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                  elementVertices[i]( 0 )  = elementVertex[0];
                  elementVertices[i]( 1 )  = elementVertex[1];
                  elementVertices[i]( 2 )  = 0.0;
               }

               auto valueTri0 = p1FuncData[vertexdof::macroface::indexFromVertex(
                   level, vertexIndices[0].x(), vertexIndices[0].y(), stencilDirection::VERTEX_C )];
               auto valueTri1 = p1FuncData[vertexdof::macroface::indexFromVertex(
                   level, vertexIndices[1].x(), vertexIndices[1].y(), stencilDirection::VERTEX_C )];
               auto valueTri2 = p1FuncData[vertexdof::macroface::indexFromVertex(
                   level, vertexIndices[2].x(), vertexIndices[2].y(), stencilDirection::VERTEX_C )];

               real_t sampledAverage =
                   evaluateSampledAverage( elementVertices, { valueTri0, valueTri1, valueTri2 }, averagingType );

               p0DofMemory[p0DofIdx] = sampledAverage;
            }
         }
      }
   }
}

} // namespace hyteg