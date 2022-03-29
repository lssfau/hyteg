/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"

namespace hyteg {
namespace volumedofspace {
namespace indexing {

ElementNeighborInfo::ElementNeighborInfo( Index                                      elementIdx,
                                          FaceType                                   faceType,
                                          uint_t                                     level,
                                          BoundaryCondition                          boundaryCondition,
                                          PrimitiveID                                faceID,
                                          const std::shared_ptr< PrimitiveStorage >& storage )
: dim_( 2 )
, elementIdx_( elementIdx )
, volumeID_( faceID )
, storage_( storage )
, faceType_( faceType )
, level_( level )
{
   WALBERLA_ASSERT( storage->faceExistsLocally( faceID ) );
   const auto face = storage->getFace( faceID );

   vertexCoordsVolume_.resize( 3 );
   neighborElementIndices_.resize( 3 );

   neighborElementVertexCoords_.resize( 3 );
   for ( uint_t i = 0; i < 3; i++ )
   {
      neighborElementVertexCoords_[i].resize( 3 );
   }

   neighborFaceElementTypes_.resize( 3 );

   interfaceVertexIndices_.resize( 3 );
   for ( uint_t i = 0; i < 3; i++ )
   {
      interfaceVertexIndices_[i].resize( 2 );
   }

   interfaceVertexCoords_.resize( 3 );
   for ( uint_t i = 0; i < 3; i++ )
   {
      interfaceVertexCoords_[i].resize( 2 );
   }

   oppositeVertexIndex_.resize( 3 );
   oppositeVertexCoords_.resize( 3 );
   neighborOppositeVertexIndex_.resize( 3 );
   neighborOppositeVertexCoords_.resize( 3 );

   onMacroBoundary_.resize( 3 );
   neighborBoundaryType_.resize( 3 );

   outwardNormal_.resize( 3 );

   const auto vertexIndicesVolume = facedof::macroface::getMicroVerticesFromMicroFace( elementIdx, faceType );

   for ( uint_t i = 0; i < 3; i++ )
   {
      const auto coord            = vertexdof::macroface::coordinateFromIndex( level, *face, vertexIndicesVolume[i] );
      vertexCoordsVolume_[i]( 0 ) = coord[0];
      vertexCoordsVolume_[i]( 1 ) = coord[1];
      vertexCoordsVolume_[i]( 2 ) = 0;
   }

   if ( faceType == FaceType::GRAY )
   {
      neighborElementIndices_[0] = Index( elementIdx.x(), elementIdx.y() - 1, 0 );
      neighborElementIndices_[1] = Index( elementIdx.x() - 1, elementIdx.y(), 0 );
      neighborElementIndices_[2] = Index( elementIdx.x(), elementIdx.y(), 0 );

      std::fill( neighborFaceElementTypes_.begin(), neighborFaceElementTypes_.end(), FaceType::BLUE );

      // bottom neighbor
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x(), elementIdx.y(), 0 );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      oppositeVertexIndex_[0]         = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + 1, elementIdx.y() - 1, 0 );

      // left neighbor
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x(), elementIdx.y(), 0 );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() - 1, elementIdx.y() + 1, 0 );

      // diagonal neighbor
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      oppositeVertexIndex_[2]         = Index( elementIdx.x(), elementIdx.y(), 0 );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
   }
   else
   {
      neighborElementIndices_[0] = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      neighborElementIndices_[1] = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      neighborElementIndices_[2] = Index( elementIdx.x(), elementIdx.y(), 0 );

      std::fill( neighborFaceElementTypes_.begin(), neighborFaceElementTypes_.end(), FaceType::GRAY );

      // right neighbor
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
      oppositeVertexIndex_[0]         = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + 2, elementIdx.y(), 0 );

      // top neighbor
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x(), elementIdx.y() + 2, 0 );

      // diagonal neighbor
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      oppositeVertexIndex_[2]         = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x(), elementIdx.y(), 0 );
   }

   onMacroBoundary_[0] = faceType == facedof::FaceType::GRAY && elementIdx.y() == 0;
   onMacroBoundary_[1] = faceType == facedof::FaceType::GRAY && elementIdx.x() == 0;
   onMacroBoundary_[2] = faceType == facedof::FaceType::GRAY &&
                         elementIdx.x() + elementIdx.y() == idx_t( levelinfo::num_microedges_per_edge( level ) - 1 );

   for ( uint_t n = 0; n < 3; n++ )
   {
      neighborBoundaryType_[n] =
          boundaryCondition.getBoundaryType( storage->getEdge( face->neighborEdges()[n] )->getMeshBoundaryFlag() );
   }

   // Looping over neighbor elements.
   for ( uint_t n = 0; n < 3; n++ )
   {
      const auto vertexIndices =
          facedof::macroface::getMicroVerticesFromMicroFace( neighborElementIndices_[n], neighborFaceElementTypes_[n] );

      for ( uint_t i = 0; i < 3; i++ )
      {
         const auto coord                        = vertexdof::macroface::coordinateFromIndex( level, *face, vertexIndices[i] );
         neighborElementVertexCoords_[n][i]( 0 ) = coord[0];
         neighborElementVertexCoords_[n][i]( 1 ) = coord[1];
         neighborElementVertexCoords_[n][i]( 2 ) = 0;
      }

      for ( uint_t i = 0; i < 2; i++ )
      {
         const auto coord = vertexdof::macroface::coordinateFromIndex( level, *face, interfaceVertexIndices_[n][i] );
         interfaceVertexCoords_[n][i]( 0 ) = coord[0];
         interfaceVertexCoords_[n][i]( 1 ) = coord[1];
         interfaceVertexCoords_[n][i]( 2 ) = 0;
      }

      const auto oppCoord           = vertexdof::macroface::coordinateFromIndex( level, *face, oppositeVertexIndex_[n] );
      oppositeVertexCoords_[n]( 0 ) = oppCoord[0];
      oppositeVertexCoords_[n]( 1 ) = oppCoord[1];
      oppositeVertexCoords_[n]( 2 ) = 0;

      const auto nOppCoord = vertexdof::macroface::coordinateFromIndex( level, *face, neighborOppositeVertexIndex_[n] );
      neighborOppositeVertexCoords_[n]( 0 ) = nOppCoord[0];
      neighborOppositeVertexCoords_[n]( 1 ) = nOppCoord[1];
      neighborOppositeVertexCoords_[n]( 2 ) = 0;

      // TODO: improve normal computation!
      Point      outerPoint = ( 1 / 3. ) * ( neighborElementVertexCoords_[n][0] + neighborElementVertexCoords_[n][1] +
                                        neighborElementVertexCoords_[n][2] );
      const auto s          = ( outerPoint - interfaceVertexCoords_[n][0] )
                         .template dot( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] ) /
                     ( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] )
                         .template dot( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] );
      const Point proj  = interfaceVertexCoords_[n][0] + s * ( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] );
      outwardNormal_[n] = ( outerPoint - proj );
      outwardNormal_[n].normalize();
   }
}

ElementNeighborInfo::ElementNeighborInfo( Index                                      elementIdx,
                                          CellType                                   cellType,
                                          uint_t                                     level,
                                          BoundaryCondition                          boundaryCondition,
                                          PrimitiveID                                cellID,
                                          const std::shared_ptr< PrimitiveStorage >& storage )
: dim_( 3 )
, elementIdx_( elementIdx )
, volumeID_( cellID )
, storage_( storage )
, cellType_( cellType )
, level_( level )
{
   WALBERLA_ASSERT( storage->cellExistsLocally( cellID ) );
   const auto cell = storage->getCell( cellID );

   vertexCoordsVolume_.resize( 4 );
   neighborElementIndices_.resize( 4 );

   neighborElementVertexCoords_.resize( 4 );
   for ( uint_t i = 0; i < 4; i++ )
   {
      neighborElementVertexCoords_[i].resize( 4 );
   }

   neighborCellElementTypes_.resize( 4 );

   interfaceVertexIndices_.resize( 4 );
   for ( uint_t i = 0; i < 4; i++ )
   {
      interfaceVertexIndices_[i].resize( 3 );
   }

   interfaceVertexCoords_.resize( 4 );
   for ( uint_t i = 0; i < 4; i++ )
   {
      interfaceVertexCoords_[i].resize( 3 );
   }

   oppositeVertexIndex_.resize( 4 );
   oppositeVertexCoords_.resize( 4 );
   neighborOppositeVertexIndex_.resize( 4 );
   neighborOppositeVertexCoords_.resize( 4 );

   onMacroBoundary_.resize( 4 );
   neighborBoundaryType_.resize( 4 );

   outwardNormal_.resize( 4 );

   const auto vertexIndicesVolume = celldof::macrocell::getMicroVerticesFromMicroCell( elementIdx, cellType );

   for ( uint_t i = 0; i < 4; i++ )
   {
      const auto coord            = vertexdof::macrocell::coordinateFromIndex( level, *cell, vertexIndicesVolume[i] );
      vertexCoordsVolume_[i]( 0 ) = coord[0];
      vertexCoordsVolume_[i]( 1 ) = coord[1];
      vertexCoordsVolume_[i]( 2 ) = coord[2];
   }

   if ( cellType == CellType::WHITE_UP )
   {
      neighborCellElementTypes_[0]    = CellType::BLUE_UP;
      neighborElementIndices_[0]      = Index( elementIdx.x() + ( -1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[0][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[0]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + ( -1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[1]    = CellType::BLUE_DOWN;
      neighborElementIndices_[1]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( -1 ) );
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( -1 ) );

      neighborCellElementTypes_[2]    = CellType::GREEN_DOWN;
      neighborElementIndices_[2]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( -1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[2]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( -1 ), elementIdx.z() + ( 1 ) );

      neighborCellElementTypes_[3]    = CellType::GREEN_UP;
      neighborElementIndices_[3]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[3]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[3] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
   }

   if ( cellType == CellType::WHITE_DOWN )
   {
      neighborCellElementTypes_[0]    = CellType::BLUE_UP;
      neighborElementIndices_[0]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[0][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[0]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 2 ) );

      neighborCellElementTypes_[1]    = CellType::GREEN_UP;
      neighborElementIndices_[1]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[1][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 2 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[2]    = CellType::GREEN_DOWN;
      neighborElementIndices_[2]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[2]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[3]    = CellType::BLUE_DOWN;
      neighborElementIndices_[3]      = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][2]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[3]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[3] = Index( elementIdx.x() + ( 2 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
   }

   if ( cellType == CellType::BLUE_UP )
   {
      neighborCellElementTypes_[0]    = CellType::WHITE_UP;
      neighborElementIndices_[0]      = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[0][2]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[0]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + ( 2 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[1]    = CellType::WHITE_DOWN;
      neighborElementIndices_[1]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( -1 ) );
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( -1 ) );

      neighborCellElementTypes_[2]    = CellType::GREEN_UP;
      neighborElementIndices_[2]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[2][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[2]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );

      neighborCellElementTypes_[3]    = CellType::GREEN_DOWN;
      neighborElementIndices_[3]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[3]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[3] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
   }

   if ( cellType == CellType::BLUE_DOWN )
   {
      neighborCellElementTypes_[0]    = CellType::WHITE_DOWN;
      neighborElementIndices_[0]      = Index( elementIdx.x() + ( -1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[0]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + ( -1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );

      neighborCellElementTypes_[1]    = CellType::GREEN_UP;
      neighborElementIndices_[1]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[1][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[2]    = CellType::GREEN_DOWN;
      neighborElementIndices_[2]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[2]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[3]    = CellType::WHITE_UP;
      neighborElementIndices_[3]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[3]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[3] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 2 ) );
   }

   if ( cellType == CellType::GREEN_UP )
   {
      neighborCellElementTypes_[0]    = CellType::WHITE_UP;
      neighborElementIndices_[0]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[0][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[0]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[1]    = CellType::WHITE_DOWN;
      neighborElementIndices_[1]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( -1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[1][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( -1 ), elementIdx.z() + ( 1 ) );

      neighborCellElementTypes_[2]    = CellType::BLUE_UP;
      neighborElementIndices_[2]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[2][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[2]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[3]    = CellType::BLUE_DOWN;
      neighborElementIndices_[3]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[3]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[3] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
   }

   if ( cellType == CellType::GREEN_DOWN )
   {
      neighborCellElementTypes_[0]    = CellType::WHITE_UP;
      neighborElementIndices_[0]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[0][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[0]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[0] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 2 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[1]    = CellType::BLUE_UP;
      neighborElementIndices_[1]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[1][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );

      neighborCellElementTypes_[2]    = CellType::BLUE_DOWN;
      neighborElementIndices_[2]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[2][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[2]         = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );

      neighborCellElementTypes_[3]    = CellType::WHITE_DOWN;
      neighborElementIndices_[3]      = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][0]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 0 ), elementIdx.z() + ( 1 ) );
      interfaceVertexIndices_[3][1]   = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      interfaceVertexIndices_[3][2]   = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
      oppositeVertexIndex_[3]         = Index( elementIdx.x() + ( 0 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 0 ) );
      neighborOppositeVertexIndex_[3] = Index( elementIdx.x() + ( 1 ), elementIdx.y() + ( 1 ), elementIdx.z() + ( 1 ) );
   }

   for ( uint_t n = 0; n < 4; n++ )
   {
      auto interface = interfaceVertexIndices_[n];
      WALBERLA_ASSERT_GREATER_EQUAL( interface[0][0], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[0][1], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[0][2], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[1][0], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[1][1], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[1][2], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[2][0], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[2][1], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[2][2], 0 );
      bool onFace0 = interface[0][2] + interface[1][2] + interface[2][2] == 0;
      bool onFace1 = interface[0][1] + interface[1][1] + interface[2][1] == 0;
      bool onFace2 = interface[0][0] + interface[1][0] + interface[2][0] == 0;
      bool onFace3 =
          ( interface[0][0] + interface[0][1] + interface[0][2] == idx_t( levelinfo::num_microvertices_per_edge( level ) ) - 1 &&
            interface[1][0] + interface[1][1] + interface[1][2] == idx_t( levelinfo::num_microvertices_per_edge( level ) ) - 1 &&
            interface[2][0] + interface[2][1] + interface[2][2] == idx_t( levelinfo::num_microvertices_per_edge( level ) ) - 1 );
      onMacroBoundary_[n] = onFace0 || onFace1 || onFace2 || onFace3;
   }

   for ( uint_t n = 0; n < 4; n++ )
   {
      neighborBoundaryType_[n] =
          boundaryCondition.getBoundaryType( storage->getFace( cell->neighborFaces()[n] )->getMeshBoundaryFlag() );
   }

   // Looping over neighbor elements.
   for ( uint_t n = 0; n < 4; n++ )
   {
      const auto vertexIndices =
          celldof::macrocell::getMicroVerticesFromMicroCell( neighborElementIndices_[n], neighborCellElementTypes_[n] );

      for ( uint_t i = 0; i < 4; i++ )
      {
         const auto coord                        = vertexdof::macrocell::coordinateFromIndex( level, *cell, vertexIndices[i] );
         neighborElementVertexCoords_[n][i]( 0 ) = coord[0];
         neighborElementVertexCoords_[n][i]( 1 ) = coord[1];
         neighborElementVertexCoords_[n][i]( 2 ) = coord[2];
      }

      for ( uint_t i = 0; i < 3; i++ )
      {
         const auto coord = vertexdof::macrocell::coordinateFromIndex( level, *cell, interfaceVertexIndices_[n][i] );
         interfaceVertexCoords_[n][i]( 0 ) = coord[0];
         interfaceVertexCoords_[n][i]( 1 ) = coord[1];
         interfaceVertexCoords_[n][i]( 2 ) = coord[2];
      }

      const auto oppCoord           = vertexdof::macrocell::coordinateFromIndex( level, *cell, oppositeVertexIndex_[n] );
      oppositeVertexCoords_[n]( 0 ) = oppCoord[0];
      oppositeVertexCoords_[n]( 1 ) = oppCoord[1];
      oppositeVertexCoords_[n]( 2 ) = oppCoord[2];

      const auto nOppCoord = vertexdof::macrocell::coordinateFromIndex( level, *cell, neighborOppositeVertexIndex_[n] );
      neighborOppositeVertexCoords_[n]( 0 ) = nOppCoord[0];
      neighborOppositeVertexCoords_[n]( 1 ) = nOppCoord[1];
      neighborOppositeVertexCoords_[n]( 2 ) = nOppCoord[2];

      // TODO: improve normal computation!
      Point outerPoint = ( 1 / 4. ) * ( neighborElementVertexCoords_[n][0] + neighborElementVertexCoords_[n][1] +
                                        neighborElementVertexCoords_[n][2] + neighborElementVertexCoords_[n][3] );

      auto x = outerPoint( 0 );
      auto y = outerPoint( 1 );
      auto z = outerPoint( 2 );

      auto normal = ( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] )
                        .cross( interfaceVertexCoords_[n][2] - interfaceVertexCoords_[n][0] )
                        .normalized();

      auto a = normal( 0 );
      auto b = normal( 1 );
      auto c = normal( 2 );

      auto d = interfaceVertexCoords_[n][0]( 0 );
      auto e = interfaceVertexCoords_[n][0]( 1 );
      auto f = interfaceVertexCoords_[n][0]( 2 );

      auto t = ( a * d - a * x + b * e - b * y + c * f - c * z ) / ( a * a + b * b + c * c );

      auto proj = outerPoint + t * normal;

      outwardNormal_[n] = ( outerPoint - proj );
      outwardNormal_[n].normalize();

      //      WALBERLA_LOG_INFO_ON_ROOT( "eltype: " << celldof::CellTypeToStr.at( cellType )
      //                                            << ", neighbor: " << celldof::CellTypeToStr.at( neighborCellType( n ) )
      //                                            << ", normal: " << outwardNormal_[n] );
   }
}

void ElementNeighborInfo::macroBoundaryNeighborElementVertexCoords( uint_t                neighbor,
                                                                    std::vector< Point >& neighborElementVertexCoords,
                                                                    Point&                neighborOppositeVertexCoords ) const
{
   // Overwriting some stuff before computing the element matrix.

   // TODO: some of what follows can be precomputed

   // The only complicated thing we have to do here is to compute the micro-vertex coordinates
   // of the micro-vertex that is inside the neighboring macro-volume. Then we must insert the
   // micro-element coordinates in the _same_ order as this is done locally.

   // The idea is to perform a "trafo" of the index space to obtain the local logical micro-element
   // index in the neighboring macro-volume.

   if ( storage_->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented." );
   }
   else
   {
      neighborElementVertexCoords.clear();
      neighborElementVertexCoords.resize( 3 );

      WALBERLA_ASSERT( storage_->faceExistsLocally( volumeID_ ) );
      const auto face = storage_->getFace( volumeID_ );

      WALBERLA_ASSERT_GREATER(
          face->getIndirectNeighborFaceIDsOverEdges().count( neighbor ), 0, "Neighbor face should exist but doesn't ..." );
      const auto neighborFace = storage_->getFace( face->getIndirectNeighborFaceIDsOverEdges().at( neighbor ) );

      // PID of the opposite macro-vertex.
      const auto oppositeMacroVertexID = neighborFace->get_vertex_opposite_to_edge( face->neighborEdges()[neighbor] );

      Index                   pseudoLocalIndex;
      std::array< uint_t, 4 > srcBasis;

      switch ( neighbor )
      {
      case 0:
         pseudoLocalIndex = Index( elementIdx_.x(), 0, 0 );
         srcBasis[0]      = neighborFace->vertex_index( face->neighborVertices()[0] );
         srcBasis[1]      = neighborFace->vertex_index( face->neighborVertices()[1] );
         srcBasis[2]      = neighborFace->vertex_index( oppositeMacroVertexID );
         srcBasis[3]      = 3;
         break;
      case 1:
         pseudoLocalIndex = Index( elementIdx_.y(), 0, 0 );
         srcBasis[0]      = neighborFace->vertex_index( face->neighborVertices()[0] );
         srcBasis[1]      = neighborFace->vertex_index( face->neighborVertices()[2] );
         srcBasis[2]      = neighborFace->vertex_index( oppositeMacroVertexID );
         srcBasis[3]      = 3;
         break;
      case 2:
         pseudoLocalIndex = Index( elementIdx_.y(), 0, 0 );
         srcBasis[0]      = neighborFace->vertex_index( face->neighborVertices()[1] );
         srcBasis[1]      = neighborFace->vertex_index( face->neighborVertices()[2] );
         srcBasis[2]      = neighborFace->vertex_index( oppositeMacroVertexID );
         srcBasis[3]      = 3;
         break;
      default:
         WALBERLA_ABORT( "Invalid neighbor ID." );
      }

      const auto elementIdxNeighborMicro = hyteg::indexing::basisConversion(
          pseudoLocalIndex, srcBasis, { 0, 1, 2, 3 }, levelinfo::num_microedges_per_edge( level_ ) );

      // In 2D the neighboring micro-face is always of type GRAY since only those are on the boundary.
      const auto nFaceTypeLocalToNeighbor = facedof::FaceType::GRAY;

      const auto neighborElementVertexIndices =
          facedof::macroface::getMicroVerticesFromMicroFace( elementIdxNeighborMicro, nFaceTypeLocalToNeighbor );

      // Eventually we need to know the coordinates of the micro-vertex that is opposite to the interface.
      // First find out what the local interface ID of the other macro-volume is.
      uint_t neighborInterfaceID;
      for ( const auto& [nifID, pid] : neighborFace->getIndirectNeighborFaceIDsOverEdges() )
      {
         if ( pid == volumeID_ )
         {
            neighborInterfaceID = nifID;
            break;
         }
      }

      // Then find the idx of the opposing micro-vertex.
      Index opposingMicroVertexIdx;
      for ( const auto& nev : neighborElementVertexIndices )
      {
         const auto localInterfaceIDs = isOnCellEdge( nev, level_ );
         if ( !algorithms::contains( localInterfaceIDs, neighborInterfaceID ) )
         {
            opposingMicroVertexIdx = nev;
            break;
         }
      }

      for ( uint_t i = 0; i < 3; i++ )
      {
         const auto coord = vertexdof::macroface::coordinateFromIndex( level_, *neighborFace, neighborElementVertexIndices[i] );
         neighborElementVertexCoords[i]( 0 ) = coord[0];
         neighborElementVertexCoords[i]( 1 ) = coord[1];
         neighborElementVertexCoords[i]( 2 ) = 0;
      }

      const auto coord = vertexdof::macroface::coordinateFromIndex( level_, *neighborFace, opposingMicroVertexIdx );

      neighborOppositeVertexCoords( 0 ) = coord[0];
      neighborOppositeVertexCoords( 1 ) = coord[1];
      neighborOppositeVertexCoords( 2 ) = 0;
   }
}

} // namespace indexing
} // namespace volumedofspace
} // namespace hyteg
