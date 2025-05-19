/*
* Copyright (c) 2017-2024 Ponsuganth Ilangovan, Nils Kohl, Marcus Mohr.
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

#include "hyteg/p0functionspace/P0Function.hpp"

namespace hyteg {

template < typename ValueType >
void P0Function< ValueType >::writeElementVolumesToDoFs( uint_t level )
{
   if ( this->storage_->hasGlobalCells() )
   {
      for ( auto it : this->storage_->getCells() )
      {
         PrimitiveID cellId = it.first;
         Cell&       cell   = *( it.second );

         const auto p0DofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( cellId, level );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               uint_t p0DofIdx = volumedofspace::indexing::index(
                   idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA );

               const indexing::Index& coarseElementIdx = idxIt;

               const std::array< indexing::Index, 4 > coarseVertexIndices =
                   celldof::macrocell::getMicroVerticesFromMicroCell( coarseElementIdx, cellType );

               std::array< Point3D, 4 > coarseVertexCoordinates;
               for ( uint_t i = 0; i < 4; i++ )
               {
                  const auto elementVertex = vertexdof::macrocell::coordinateFromIndex( level, cell, coarseVertexIndices[i] );

                  Point3D elementVertexMapped;
                  cell.getGeometryMap()->evalF(elementVertex, elementVertexMapped);

                  coarseVertexCoordinates[i]( 0 ) = elementVertexMapped[0];
                  coarseVertexCoordinates[i]( 1 ) = elementVertexMapped[1];
                  coarseVertexCoordinates[i]( 2 ) = elementVertexMapped[2];
               }

               for ( uint_t i = 1; i < 4; i++ )
               {
                  coarseVertexCoordinates[i] -= coarseVertexCoordinates[0];
               }

               real_t coarseTetVolume =
                   ( 1.0 / 6.0 ) *
                   std::abs( coarseVertexCoordinates[1].dot( coarseVertexCoordinates[2].cross( coarseVertexCoordinates[3] ) ) );

               p0DofMemory[p0DofIdx] = coarseTetVolume;
            }
         }
      }
   }
   else
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         PrimitiveID faceId = it.first;
         Face&       face   = *( it.second );

         const auto p0DofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( faceId, level );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               uint_t p0DofIdx = volumedofspace::indexing::index(
                   idxIt.x(), idxIt.y(), faceType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA );

               const indexing::Index& coarseElementIdx = idxIt;

               const std::array< indexing::Index, 3u > coarseVertexIndices =
                   facedof::macroface::getMicroVerticesFromMicroFace( coarseElementIdx, faceType );

               std::array< Point3D, 3u > coarseVertexCoordinates;
               for ( uint_t i = 0; i < 3u; i++ )
               {
                  const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, coarseVertexIndices[i] );

                  Point3D elementVertexMapped;
                  face.getGeometryMap()->evalF(elementVertex, elementVertexMapped);

                  coarseVertexCoordinates[i]( 0 ) = elementVertexMapped[0];
                  coarseVertexCoordinates[i]( 1 ) = elementVertexMapped[1];
                  coarseVertexCoordinates[i]( 2 ) = elementVertexMapped[2];
               }

               auto ab = coarseVertexCoordinates[1] - coarseVertexCoordinates[0];
               auto ac = coarseVertexCoordinates[2] - coarseVertexCoordinates[0];

               real_t faceVolume     = ( 1.0 / 2.0 ) * ( ab.cross( ac ) ).norm();
               p0DofMemory[p0DofIdx] = faceVolume;
            }
         }
      }
   }
}

template class P0Function< real_t >;

} // namespace hyteg
