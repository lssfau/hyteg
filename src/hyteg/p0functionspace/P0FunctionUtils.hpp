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
#pragma once

#include "hyteg/p0functionspace/P0Function.hpp"

namespace hyteg {

namespace p0averaging {

inline real_t evaluateSampledAverage( std::array< Point3D, 3 >             microTriangles,
                                      std::array< real_t, 3 >              valueTriangles,
                                      hyteg::p0averaging::AVERAGING_METHOD averagingMethod )
{
   Point3D microTet0 = microTriangles[0];
   Point3D microTet1 = microTriangles[1];
   Point3D microTet2 = microTriangles[2];

   real_t valueTet0 = valueTriangles[0];
   real_t valueTet1 = valueTriangles[1];
   real_t valueTet2 = valueTriangles[2];

   Point3D coordinates = ( microTet0 + microTet1 + microTet2 ) / 3.0;

   auto xLocal = vertexdof::macroface::transformToLocalTri( microTet0, microTet1, microTet2, coordinates );

   auto value = valueTet0 * ( real_c( 1.0 ) - xLocal[0] - xLocal[1] ) + valueTet1 * xLocal[0] + valueTet2 * xLocal[1];

   if ( averagingMethod == p0averaging::AVERAGING_METHOD::ARITHMETIC )
   {
      return ( valueTet0 + valueTet1 + valueTet2 ) / 3.0;
   }
   else
   {
      WALBERLA_ABORT( "Not implemented" );
   }
}

inline real_t evaluateSampledAverage( std::array< Point3D, 4 >             microTets,
                                      std::array< real_t, 4 >              valueTets,
                                      hyteg::p0averaging::AVERAGING_METHOD averagingMethod )
{
   Point3D microTet0 = microTets[0];
   Point3D microTet1 = microTets[1];
   Point3D microTet2 = microTets[2];
   Point3D microTet3 = microTets[3];

   real_t valueTet0 = valueTets[0];
   real_t valueTet1 = valueTets[1];
   real_t valueTet2 = valueTets[2];
   real_t valueTet3 = valueTets[3];

   Point3D coordinates = ( microTet0 + microTet1 + microTet2 + microTet3 ) / 4.0;

   // auto xLocal = vertexdof::macrocell::detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );

   std::function< real_t( const Point3D& ) > locallyEvaluate = [&]( const Point3D& x ) {
      return valueTet0 * ( real_c( 1.0 ) - x[0] - x[1] - x[2] ) + valueTet1 * x[0] + valueTet2 * x[1] + valueTet3 * x[2];
   };

   // This 5 point quadrature is hard coded for now
   Point3D qp1 = Point3D( 0.25, 0.25, 0.25 );
   Point3D qp2 = Point3D( 0.16666667, 0.16666667, 0.5 );
   Point3D qp3 = Point3D( 0.16666667, 0.5, 0.16666667 );
   Point3D qp4 = Point3D( 0.5, 0.16666667, 0.16666667 );
   Point3D qp5 = Point3D( 0.16666667, 0.16666667, 0.16666667 );

   real_t val1 = locallyEvaluate( qp1 );
   real_t val2 = locallyEvaluate( qp2 );
   real_t val3 = locallyEvaluate( qp3 );
   real_t val4 = locallyEvaluate( qp4 );
   real_t val5 = locallyEvaluate( qp5 );

   auto xLocal = vertexdof::macrocell::detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );
   auto value  = valueTet0 * ( real_c( 1.0 ) - xLocal[0] - xLocal[1] - xLocal[2] ) + valueTet1 * xLocal[0] +
                valueTet2 * xLocal[1] + valueTet3 * xLocal[2];

   if ( averagingMethod == p0averaging::AVERAGING_METHOD::ARITHMETIC )
   {
      return ( value + valueTet0 + valueTet1 + valueTet2 + valueTet3 ) / 5.0;
   }
   else if ( averagingMethod == p0averaging::AVERAGING_METHOD::ARITHMETIC_QP )
   {
      return ( val1 + val2 + val3 + val4 + val5 ) / 5.0;
   }
   else if ( averagingMethod == p0averaging::AVERAGING_METHOD::HARMONIC )
   {
      return 5.0 / ( ( 1.0 / value ) + ( 1.0 / valueTet0 ) + ( 1.0 / valueTet1 ) + ( 1.0 / valueTet2 ) + ( 1.0 / valueTet3 ) );
   }
   else if ( averagingMethod == p0averaging::AVERAGING_METHOD::HARMONIC_QP )
   {
      return 5.0 / ( ( 1.0 / val1 ) + ( 1.0 / val2 ) + ( 1.0 / val3 ) + ( 1.0 / val4 ) + ( 1.0 / val5 ) );
   }
   else if ( averagingMethod == p0averaging::AVERAGING_METHOD::GEOMETRIC )
   {
      return std::pow( value, 1.0 / 5.0 ) * std::pow( valueTet0, 1.0 / 5.0 ) * std::pow( valueTet1, 1.0 / 5.0 ) *
             std::pow( valueTet2, 1.0 / 5.0 ) * std::pow( valueTet3, 1.0 / 5.0 );
   }
   else if ( averagingMethod == p0averaging::AVERAGING_METHOD::GEOMETRIC_QP )
   {
      return std::pow( val1, 1.0 / 5.0 ) * std::pow( val2, 1.0 / 5.0 ) * std::pow( val3, 1.0 / 5.0 ) *
             std::pow( val4, 1.0 / 5.0 ) * std::pow( val5, 1.0 / 5.0 );
   }
   else
   {
      WALBERLA_ABORT( "Not implemented" );
   }
}
} // namespace p0averaging

template < typename ValueType >
void P0Function< ValueType >::transferToLowerLevel( uint_t                               level,
                                                    hyteg::p0averaging::AVERAGING_METHOD averagingMethod,
                                                    bool                                 volumeWeighted = false )
{
   if ( this->storage_->hasGlobalCells() )
   {
      for ( auto& it : this->storage_->getCells() )
      {
         PrimitiveID cellID = it.first;
         Cell&       cell   = *( it.second );

         uint_t coarseLevel = level - 1;
         uint_t fineLevel   = level;

         const auto coarseDofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( cellID, coarseLevel );
         const auto fineDofMemory   = this->getDGFunction()->volumeDoFFunction()->dofMemory( cellID, fineLevel );

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

                  if ( volumeWeighted )
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

               if ( averagingMethod == p0averaging::AVERAGING_METHOD::ARITHMETIC )
               {
                  finalAvgVal = fineAverage / nFineElements;
               }
               else if ( averagingMethod == p0averaging::AVERAGING_METHOD::HARMONIC )
               {
                  finalAvgVal = nFineElements / fineInvAverage;
               }
               else if ( averagingMethod == p0averaging::AVERAGING_METHOD::GEOMETRIC )
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
      for ( auto& it : this->storage_->getFaces() )
      {
         PrimitiveID faceID = it.first;
         Face&       face   = *( it.second );

         uint_t coarseLevel = level - 1;
         uint_t fineLevel   = level;

         const auto coarseDofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( faceID, coarseLevel );
         const auto fineDofMemory   = this->getDGFunction()->volumeDoFFunction()->dofMemory( faceID, fineLevel );

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

                  if ( volumeWeighted )
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

               if ( averagingMethod == p0averaging::AVERAGING_METHOD::ARITHMETIC )
               {
                  finalAvgVal = fineAverage / nFineElements;
               }
               else if ( averagingMethod == p0averaging::AVERAGING_METHOD::HARMONIC )
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

template < typename ValueType >
void P0Function< ValueType >::transferToAllLowerLevels( uint_t                               sourceLevel,
                                                        hyteg::p0averaging::AVERAGING_METHOD averagingMethod,
                                                        bool                                 volumeWeighted )
{
   for ( uint_t level = sourceLevel; level > this->getMinLevel(); level-- )
   {
      transferToLowerLevel( level, averagingMethod, volumeWeighted );
   }
}

template < typename ValueType >
void P0Function< ValueType >::averageFromP1( P1Function< real_t >                 src,
                                             uint_t                               level,
                                             hyteg::p0averaging::AVERAGING_METHOD averagingMethod )
{
   if ( this->storage_->hasGlobalCells() )
   {
      for ( auto it : this->storage_->getCells() )
      {
         PrimitiveID cellId = it.first;
         Cell&       cell   = *( it.second );

         const auto p0DofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( cellId, level );

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
                                                               averagingMethod );

               p0DofMemory[p0DofIdx] = sampledAverage;
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
                   evaluateSampledAverage( elementVertices, { valueTri0, valueTri1, valueTri2 }, averagingMethod );

               p0DofMemory[p0DofIdx] = sampledAverage;
            }
         }
      }
   }
}
} // namespace hyteg
