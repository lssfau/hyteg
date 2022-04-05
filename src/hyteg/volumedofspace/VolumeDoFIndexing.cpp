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

const uint_t ElementNeighborInfo::NOT_AT_BOUNDARY = std::numeric_limits< uint_t >::max();

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

   macroBoundary_.resize( 3 );
   std::fill( macroBoundary_.begin(), macroBoundary_.end(), NOT_AT_BOUNDARY );

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
      if ( elementIdx.y() == 0 )
      {
         macroBoundary_[0] = 0;
      }

      // left neighbor
      interfaceVertexIndices_[1][0]   = Index( elementIdx.x(), elementIdx.y(), 0 );
      interfaceVertexIndices_[1][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      oppositeVertexIndex_[1]         = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      neighborOppositeVertexIndex_[1] = Index( elementIdx.x() - 1, elementIdx.y() + 1, 0 );
      if ( elementIdx.x() == 0 )
      {
         macroBoundary_[1] = 1;
      }

      // diagonal neighbor
      interfaceVertexIndices_[2][0]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
      interfaceVertexIndices_[2][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
      oppositeVertexIndex_[2]         = Index( elementIdx.x(), elementIdx.y(), 0 );
      neighborOppositeVertexIndex_[2] = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
      if ( elementIdx.x() + elementIdx.y() == idx_t( levelinfo::num_microedges_per_edge( level ) ) - 1 )
      {
         macroBoundary_[2] = 2;
      }
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

   for ( uint_t n = 0; n < 3; n++ )
   {
      if ( atMacroBoundary( n ) )
      {
         neighborBoundaryType_[n] = boundaryCondition.getBoundaryType(
             storage->getEdge( face->neighborEdges()[macroBoundary_[n]] )->getMeshBoundaryFlag() );
      }
   }

   // Looping over neighbor elements.
   for ( uint_t n = 0; n < 3; n++ )
   {
      if ( !atMacroBoundary( n ) )
      {
         // Certain properties are only meaningful if we are not at a macro-boundary.
         // Neighbor in local macro.
         const auto vertexIndices =
             facedof::macroface::getMicroVerticesFromMicroFace( neighborElementIndices_[n], neighborFaceElementTypes_[n] );

         for ( uint_t i = 0; i < 3; i++ )
         {
            const auto coord                        = vertexdof::macroface::coordinateFromIndex( level, *face, vertexIndices[i] );
            neighborElementVertexCoords_[n][i]( 0 ) = coord[0];
            neighborElementVertexCoords_[n][i]( 1 ) = coord[1];
            neighborElementVertexCoords_[n][i]( 2 ) = 0;
         }

         const auto nOppCoord = vertexdof::macroface::coordinateFromIndex( level, *face, neighborOppositeVertexIndex_[n] );
         neighborOppositeVertexCoords_[n]( 0 ) = nOppCoord[0];
         neighborOppositeVertexCoords_[n]( 1 ) = nOppCoord[1];
         neighborOppositeVertexCoords_[n]( 2 ) = 0;
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

      // TODO: improve normal computation!
      Point      innerPoint = ( 1 / 3. ) * ( elementVertexCoords()[0] + elementVertexCoords()[1] + elementVertexCoords()[2] );
      const auto s =
          ( innerPoint - interfaceVertexCoords_[n][0] ).dot( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] ) /
          ( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] )
              .dot( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] );
      const Point proj  = interfaceVertexCoords_[n][0] + s * ( interfaceVertexCoords_[n][1] - interfaceVertexCoords_[n][0] );
      outwardNormal_[n] = ( innerPoint - proj );
      outwardNormal_[n].normalize();
      outwardNormal_[n] = -outwardNormal_[n];
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

   macroBoundary_.resize( 4 );
   std::fill( macroBoundary_.begin(), macroBoundary_.end(), NOT_AT_BOUNDARY );

   neighborBoundaryType_.resize( 4 );
   std::fill( neighborBoundaryType_.begin(), neighborBoundaryType_.end(), Inner );

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

   // Checking which macro-boundary interfaces the micro-element touches.
   for ( uint_t n = 0; n < 4; n++ )
   {
      auto interface = interfaceVertexIndices_[n];

      // These vertex indices better be contained in the macro-volume ...
      WALBERLA_ASSERT_GREATER_EQUAL( interface[0][0], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[0][1], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[0][2], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[1][0], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[1][1], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[1][2], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[2][0], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[2][1], 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( interface[2][2], 0 );

      std::vector< bool > onFace( 4 );

      onFace[0] = interface[0][2] + interface[1][2] + interface[2][2] == 0;
      onFace[1] = interface[0][1] + interface[1][1] + interface[2][1] == 0;
      onFace[2] = interface[0][0] + interface[1][0] + interface[2][0] == 0;
      onFace[3] =
          ( interface[0][0] + interface[0][1] + interface[0][2] == idx_t( levelinfo::num_microvertices_per_edge( level ) ) - 1 &&
            interface[1][0] + interface[1][1] + interface[1][2] == idx_t( levelinfo::num_microvertices_per_edge( level ) ) - 1 &&
            interface[2][0] + interface[2][1] + interface[2][2] == idx_t( levelinfo::num_microvertices_per_edge( level ) ) - 1 );

      WALBERLA_ASSERT_LESS_EQUAL( (int) onFace[0] + (int) onFace[1] + (int) onFace[2] + (int) onFace[3],
                                  1,
                                  "Each micro-interface should only be located at a single (or no) macro-boundary." );

      uint_t mBound = NOT_AT_BOUNDARY;
      for ( uint_t i = 0; i < 4; i++ )
      {
         if ( onFace[i] )
         {
            mBound = i;
            break;
         }
      }

      // Most cell types can only be inner or at a certain boundary.
      WALBERLA_ASSERT( !( cellType == celldof::CellType::WHITE_DOWN && mBound != NOT_AT_BOUNDARY ) );
      WALBERLA_ASSERT( !( cellType == celldof::CellType::BLUE_UP && !( mBound == NOT_AT_BOUNDARY || mBound == 0 ) ) );
      WALBERLA_ASSERT( !( cellType == celldof::CellType::GREEN_UP && !( mBound == NOT_AT_BOUNDARY || mBound == 1 ) ) );
      WALBERLA_ASSERT( !( cellType == celldof::CellType::BLUE_DOWN && !( mBound == NOT_AT_BOUNDARY || mBound == 2 ) ) );
      WALBERLA_ASSERT( !( cellType == celldof::CellType::GREEN_DOWN && !( mBound == NOT_AT_BOUNDARY || mBound == 3 ) ) );

      macroBoundary_[n] = mBound;
   }

   for ( uint_t n = 0; n < 4; n++ )
   {
      if ( atMacroBoundary( n ) )
      {
         neighborBoundaryType_[n] = boundaryCondition.getBoundaryType(
             storage->getFace( cell->neighborFaces()[macroBoundary_[n]] )->getMeshBoundaryFlag() );
      }
   }

   // Looping over neighbor elements.
   for ( uint_t n = 0; n < 4; n++ )
   {
      if ( !atMacroBoundary( n ) )
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

         const auto nOppCoord = vertexdof::macrocell::coordinateFromIndex( level, *cell, neighborOppositeVertexIndex_[n] );
         neighborOppositeVertexCoords_[n]( 0 ) = nOppCoord[0];
         neighborOppositeVertexCoords_[n]( 1 ) = nOppCoord[1];
         neighborOppositeVertexCoords_[n]( 2 ) = nOppCoord[2];
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

      // TODO: improve normal computation!
      Point innerPoint = ( 1 / 4. ) * ( elementVertexCoords()[0] + elementVertexCoords()[1] + elementVertexCoords()[2] +
                                        elementVertexCoords()[3] );

      auto x = innerPoint( 0 );
      auto y = innerPoint( 1 );
      auto z = innerPoint( 2 );

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

      auto proj = innerPoint + t * normal;

      outwardNormal_[n] = ( innerPoint - proj );
      outwardNormal_[n].normalize();
      outwardNormal_[n] = -outwardNormal_[n];
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

   WALBERLA_CHECK( atMacroBoundary( neighbor ), "Element is not located at boundary." );

   if ( storage_->hasGlobalCells() )
   {
      neighborElementVertexCoords.clear();
      neighborElementVertexCoords.resize( 4 );

      WALBERLA_ASSERT( storage_->cellExistsLocally( volumeID_ ) );
      const auto cell = storage_->getCell( volumeID_ );

      WALBERLA_ASSERT_GREATER( cell->getIndirectNeighborCellIDsOverFaces().count( macroBoundaryID( neighbor ) ),
                               0,
                               "Neighbor cell should exist but doesn't ..." );
      const auto neighborCell =
          storage_->getCell( cell->getIndirectNeighborCellIDsOverFaces().at( macroBoundaryID( neighbor ) ) );

      // PID of the opposite macro-vertex.
      const auto oppositeMacroVertexID = neighborCell->getOppositeVertexID( cell->neighborFaces()[macroBoundaryID( neighbor )] );

      Index                   pseudoLocalIndex;
      std::array< uint_t, 4 > srcBasis;

      switch ( macroBoundaryID( neighbor ) )
      {
      case 0:
         pseudoLocalIndex = Index( elementIdx_.x(), elementIdx_.y(), 0 );
         srcBasis[0]      = neighborCell->getLocalVertexID( cell->neighborVertices()[0] );
         srcBasis[1]      = neighborCell->getLocalVertexID( cell->neighborVertices()[1] );
         srcBasis[2]      = neighborCell->getLocalVertexID( cell->neighborVertices()[2] );
         srcBasis[3]      = neighborCell->getLocalVertexID( oppositeMacroVertexID );
         break;
      case 1:
         pseudoLocalIndex = Index( elementIdx_.x(), elementIdx_.z(), 0 );
         srcBasis[0]      = neighborCell->getLocalVertexID( cell->neighborVertices()[0] );
         srcBasis[1]      = neighborCell->getLocalVertexID( cell->neighborVertices()[1] );
         srcBasis[2]      = neighborCell->getLocalVertexID( cell->neighborVertices()[3] );
         srcBasis[3]      = neighborCell->getLocalVertexID( oppositeMacroVertexID );
         break;
      case 2:
         pseudoLocalIndex = Index( elementIdx_.y(), elementIdx_.z(), 0 );
         srcBasis[0]      = neighborCell->getLocalVertexID( cell->neighborVertices()[0] );
         srcBasis[1]      = neighborCell->getLocalVertexID( cell->neighborVertices()[2] );
         srcBasis[2]      = neighborCell->getLocalVertexID( cell->neighborVertices()[3] );
         srcBasis[3]      = neighborCell->getLocalVertexID( oppositeMacroVertexID );
         break;
      case 3:
         pseudoLocalIndex = Index( elementIdx_.y(), elementIdx_.z(), 0 );
         srcBasis[0]      = neighborCell->getLocalVertexID( cell->neighborVertices()[1] );
         srcBasis[1]      = neighborCell->getLocalVertexID( cell->neighborVertices()[2] );
         srcBasis[2]      = neighborCell->getLocalVertexID( cell->neighborVertices()[3] );
         srcBasis[3]      = neighborCell->getLocalVertexID( oppositeMacroVertexID );
         break;
      default:
         WALBERLA_ABORT( "Invalid neighbor ID." );
      }

      auto cellWidth = levelinfo::num_microedges_per_edge( level_ );
      if ( cellType_ != celldof::CellType::WHITE_UP )
      {
         cellWidth -= 1;
      }
      const auto elementIdxNeighborMicro =
          hyteg::indexing::basisConversion( pseudoLocalIndex, srcBasis, { 0, 1, 2, 3 }, cellWidth );

      // We need to find out the neighboring micro-cell type local to the neighboring macro-cell.
      auto nCellTypeLocalToNeighbor = celldof::CellType::WHITE_UP;
      // It's WHITE_UP if the local micro-cell is WHITE_UP.
      WALBERLA_ASSERT_UNEQUAL( cellType_, celldof::CellType::WHITE_DOWN, "Invalid cell type." );
      // For the remaining four types it remains to find out the local macro-face ID of the interface from the view of the
      // neighboring macro-cell.
      const auto neighborMacroLocalFaceID =
          neighborCell->getLocalFaceID( cell->neighborFaces().at( macroBoundaryID( neighbor ) ) );
      if ( cellType_ != celldof::CellType::WHITE_UP )
      {
         switch ( neighborMacroLocalFaceID )
         {
         case 0:
            nCellTypeLocalToNeighbor = celldof::CellType::BLUE_UP;
            break;
         case 1:
            nCellTypeLocalToNeighbor = celldof::CellType::GREEN_UP;
            break;
         case 2:
            nCellTypeLocalToNeighbor = celldof::CellType::BLUE_DOWN;
            break;
         case 3:
            nCellTypeLocalToNeighbor = celldof::CellType::GREEN_DOWN;
            break;
         }
      }

      const auto neighborElementVertexIndices =
          celldof::macrocell::getMicroVerticesFromMicroCell( elementIdxNeighborMicro, nCellTypeLocalToNeighbor );

      // Find the idx of the opposing micro-vertex. It should not be located on the interface macro.
      Index opposingMicroVertexIdx;
      for ( const auto& nev : neighborElementVertexIndices )
      {
         const auto localInterfaceIDs = isOnCellFace( nev, level_ );
         if ( !algorithms::contains( localInterfaceIDs, neighborMacroLocalFaceID ) )
         {
            opposingMicroVertexIdx = nev;
            break;
         }
      }

      for ( uint_t i = 0; i < 4; i++ )
      {
         const auto coord = vertexdof::macrocell::coordinateFromIndex( level_, *neighborCell, neighborElementVertexIndices[i] );
         neighborElementVertexCoords[i]( 0 ) = coord[0];
         neighborElementVertexCoords[i]( 1 ) = coord[1];
         neighborElementVertexCoords[i]( 2 ) = coord[2];
      }

      const auto coord = vertexdof::macrocell::coordinateFromIndex( level_, *neighborCell, opposingMicroVertexIdx );

      neighborOppositeVertexCoords( 0 ) = coord[0];
      neighborOppositeVertexCoords( 1 ) = coord[1];
      neighborOppositeVertexCoords( 2 ) = coord[2];
   }
   else
   {
      neighborElementVertexCoords.clear();
      neighborElementVertexCoords.resize( 3 );

      WALBERLA_ASSERT( storage_->faceExistsLocally( volumeID_ ) );
      const auto face = storage_->getFace( volumeID_ );

      WALBERLA_ASSERT_GREATER( face->getIndirectNeighborFaceIDsOverEdges().count( macroBoundaryID( neighbor ) ),
                               0,
                               "Neighbor face should exist but doesn't ..." );
      const auto neighborFace =
          storage_->getFace( face->getIndirectNeighborFaceIDsOverEdges().at( macroBoundaryID( neighbor ) ) );

      // PID of the opposite macro-vertex.
      const auto oppositeMacroVertexID =
          neighborFace->get_vertex_opposite_to_edge( face->neighborEdges()[macroBoundaryID( neighbor )] );

      Index                   pseudoLocalIndex;
      std::array< uint_t, 4 > srcBasis;

      switch ( macroBoundaryID( neighbor ) )
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
      const auto neighborMacroLocalEdgeID = neighborFace->edge_index( face->neighborEdges().at( neighbor ) );

      // Then find the idx of the opposing micro-vertex.
      Index opposingMicroVertexIdx;
      for ( const auto& nev : neighborElementVertexIndices )
      {
         const auto localInterfaceIDs = isOnCellEdge( nev, level_ );
         if ( !algorithms::contains( localInterfaceIDs, neighborMacroLocalEdgeID ) )
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
