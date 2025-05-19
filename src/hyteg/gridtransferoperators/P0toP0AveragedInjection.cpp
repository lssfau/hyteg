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

#include "hyteg/gridtransferoperators/P0toP0AveragedInjection.hpp"

#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

void P0toP0AveragedInjection::restrict( const P0Function< real_t >& function,
                                        const uint_t&               level,
                                        const DoFType&              flag ) const
{
  WALBERLA_UNUSED( flag );

   if ( function.getStorage()->hasGlobalCells() )
   {
      for ( auto& it : function.getStorage()->getCells() )
      {
         PrimitiveID cellID = it.first;
         Cell&       cell   = *( it.second );

         uint_t coarseLevel = level - 1;
         uint_t fineLevel   = level;

         const auto coarseDofMemory = function.getDGFunction()->volumeDoFFunction()->dofMemory( cellID, coarseLevel );
         const auto fineDofMemory   = function.getDGFunction()->volumeDoFFunction()->dofMemory( cellID, fineLevel );

         for ( auto coarseCellType : celldof::allCellTypes )
         {
            for ( const auto& itSrc : celldof::macrocell::Iterator( coarseLevel, coarseCellType ) )
            {
               const indexing::Index& coarseElementIdx = itSrc;

               const std::array< indexing::Index, 4 > coarseVertexIndices =
                   celldof::macrocell::getMicroVerticesFromMicroCell( coarseElementIdx, coarseCellType );

               std::array< Point3D, 4 > coarseVertexCoordinates;
               for ( uint_t i = 0; i < 4; i++ )
               {
                  const auto elementVertex = vertexdof::macrocell::coordinateFromIndex( level, cell, coarseVertexIndices[i] );
                  coarseVertexCoordinates[i]( 0 ) = elementVertex[0];
                  coarseVertexCoordinates[i]( 1 ) = elementVertex[1];
                  coarseVertexCoordinates[i]( 2 ) = elementVertex[2];
               }

               for ( uint_t i = 1; i < 4; i++ )
               {
                  coarseVertexCoordinates[i] -= coarseVertexCoordinates[0];
               }

               real_t coarseTetVolume =
                   ( 1.0 / 6.0 ) *
                   std::abs( coarseVertexCoordinates[1].dot( coarseVertexCoordinates[2].cross( coarseVertexCoordinates[3] ) ) );

               std::vector< hyteg::indexing::Index > fineElementIndices;
               std::vector< celldof::CellType >      fineCellTypes;

               real_t fineAverage    = 0.0;
               real_t fineInvAverage = 0.0;
               real_t fineGeoAverage = 1.0;

               volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                   coarseElementIdx, coarseCellType, fineElementIndices, fineCellTypes );

               WALBERLA_CHECK_EQUAL( fineElementIndices.size(), fineCellTypes.size() );

               uint_t nFineElements = fineElementIndices.size();

               for ( uint_t fineIdx = 0; fineIdx < fineElementIndices.size(); fineIdx += 1 )
               {
                  auto fineElementIdx = fineElementIndices[fineIdx];
                  auto fineCellType   = fineCellTypes[fineIdx];

                  const std::array< indexing::Index, 4 > fineVertexIndices =
                      celldof::macrocell::getMicroVerticesFromMicroCell( fineElementIdx, fineCellType );

                  std::array< Point3D, 4 > fineVertexCoordinates;
                  for ( uint_t i = 0; i < 4; i++ )
                  {
                     const auto elementVertex = vertexdof::macrocell::coordinateFromIndex( level, cell, fineVertexIndices[i] );
                     fineVertexCoordinates[i]( 0 ) = elementVertex[0];
                     fineVertexCoordinates[i]( 1 ) = elementVertex[1];
                     fineVertexCoordinates[i]( 2 ) = elementVertex[2];
                  }

                  for ( uint_t i = 1; i < 4; i++ )
                  {
                     fineVertexCoordinates[i] -= fineVertexCoordinates[0];
                  }

                  real_t fineTetVolume =
                      ( 1.0 / 6.0 ) *
                      std::abs( fineVertexCoordinates[1].dot( fineVertexCoordinates[2].cross( fineVertexCoordinates[3] ) ) );

                  real_t fineVal =
                      fineDofMemory[volumedofspace::indexing::index( fineElementIdx.x(),
                                                                     fineElementIdx.y(),
                                                                     fineElementIdx.z(),
                                                                     fineCellType,
                                                                     0u,
                                                                     1u,
                                                                     fineLevel,
                                                                     volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )];

                  real_t wi = ( fineTetVolume / coarseTetVolume );

                  if ( volumeWeighted_ )
                  {
                     fineAverage += wi * fineVal;
                     fineInvAverage += ( wi / fineVal );
                     fineGeoAverage *= std::pow( fineVal, wi );
                  }
                  else
                  {
                     fineAverage += fineVal;
                     fineInvAverage += ( 1.0 / fineVal );
                     fineGeoAverage *= std::pow( fineVal, 1.0 / nFineElements );
                  }
               }

               real_t finalAvgVal = 0.0;

               if ( averagingType_ == hyteg::AveragingType::ARITHMETIC )
               {
                  finalAvgVal = fineAverage / nFineElements;
               }
               else if ( averagingType_ == hyteg::AveragingType::HARMONIC )
               {
                  finalAvgVal = nFineElements / fineInvAverage;
               }
               else if ( averagingType_ == hyteg::AveragingType::GEOMETRIC )
               {
                  finalAvgVal = fineGeoAverage;
               }
               else
               {
                  WALBERLA_ABORT( "Not implemented" );
               }

               coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
                                                                coarseElementIdx.y(),
                                                                coarseElementIdx.z(),
                                                                coarseCellType,
                                                                0u,
                                                                1u,
                                                                coarseLevel,
                                                                volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )] =
                   finalAvgVal;
            }
         }
      }
   }
   else
   {
      for ( auto& it : function.getStorage()->getFaces() )
      {
         PrimitiveID faceID = it.first;
         Face&       face   = *( it.second );

         uint_t coarseLevel = level - 1;
         uint_t fineLevel   = level;

         const auto coarseDofMemory = function.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, coarseLevel );
         const auto fineDofMemory   = function.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, fineLevel );

         for ( auto coarseFaceType : facedof::allFaceTypes )
         {
            for ( const auto& itSrc : facedof::macroface::Iterator( coarseLevel, coarseFaceType ) )
            {
               const indexing::Index& coarseElementIdx = itSrc;

               std::vector< hyteg::indexing::Index > fineElementIndices;
               std::vector< facedof::FaceType >      fineFaceTypes;

               real_t fineAverage    = 0.0;
               real_t fineInvAverage = 0.0;
               real_t fineGeoAverage = 1.0;

               volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                   coarseElementIdx, coarseFaceType, fineElementIndices, fineFaceTypes );

               WALBERLA_CHECK_EQUAL( fineElementIndices.size(), fineFaceTypes.size() );

               uint_t nFineElements = fineElementIndices.size();

               for ( uint_t fineIdx = 0; fineIdx < fineElementIndices.size(); fineIdx += 1 )
               {
                  auto fineElementIdx = fineElementIndices[fineIdx];
                  auto fineFaceType   = fineFaceTypes[fineIdx];

                  real_t fineVal =
                      fineDofMemory[volumedofspace::indexing::index( fineElementIdx.x(),
                                                                     fineElementIdx.y(),
                                                                     fineFaceType,
                                                                     0u,
                                                                     1u,
                                                                     fineLevel,
                                                                     volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )];

                  if ( volumeWeighted_ )
                  {
                     WALBERLA_ABORT( "Not implemented" );
                  }
                  else
                  {
                     fineAverage += fineVal;
                     fineInvAverage += ( 1.0 / fineVal );
                     fineGeoAverage *= std::pow( fineVal, 1.0 / nFineElements );
                  }
               }

               real_t finalAvgVal = 0.0;

               if ( averagingType_ == hyteg::AveragingType::ARITHMETIC )
               {
                  finalAvgVal = fineAverage / nFineElements;
               }
               else if ( averagingType_ == hyteg::AveragingType::HARMONIC )
               {
                  finalAvgVal = nFineElements / fineInvAverage;
               }
               else
               {
                  WALBERLA_ABORT( "Not implemented" );
               }

               coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
                                                                coarseElementIdx.y(),
                                                                coarseFaceType,
                                                                0u,
                                                                1u,
                                                                coarseLevel,
                                                                volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )] =
                   finalAvgVal;
            }
         }
      }
   }
}

void P0toP0AveragedInjection::restrictToAllLowerLevels( const P0Function< real_t >& function, const uint_t& sourceLevel )
{
  for ( uint_t level = sourceLevel; level > function.getMinLevel(); level-- )
   {
      this->restrict( function, level );
   }
}

} // namespace hyteg
