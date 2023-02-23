/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {
namespace edgedof {
namespace macroface {

using indexing::Index;
using walberla::real_c;
using walberla::uint_t;

inline Point3D
    coordinateFromIndex( const uint_t& Level, const Face& face, const Index& index, const EdgeDoFOrientation& orientation )
{
   const real_t  stepFrequency = 1.0 / levelinfo::num_microedges_per_edge( Level );
   const Point3D xStep         = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) * stepFrequency;
   const Point3D yStep         = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) * stepFrequency;
   Point3D       increment;
   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
      increment = 0.5 * xStep;
      break;
   case EdgeDoFOrientation::Y:
      increment = 0.5 * yStep;
      break;
   case EdgeDoFOrientation::XY:
      increment = 0.5 * xStep + 0.5 * yStep;
      break;
   default:
      WALBERLA_ABORT( "Invalid edgedof orientation." )
   }
   return face.getCoordinates()[0] + xStep * real_c( index.x() ) + yStep * real_c( index.y() ) + increment;
}

inline indexing::Index getIndexInNeighboringMacroCell( const indexing::Index&  edgeDoFIndexInMacroFace,
                                                       const Face&             face,
                                                       const uint_t&           neighborCellID,
                                                       const PrimitiveStorage& storage,
                                                       const uint_t&           level )
{
   const Cell&  neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
   const uint_t localFaceID  = neighborCell.getLocalFaceID( face.getID() );

   const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 3, 4 >(
       { neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
         neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
         neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ) } );

   const auto indexInMacroCell = indexing::basisConversion(
       edgeDoFIndexInMacroFace, localVertexIDsAtCell, { 0, 1, 2, 3 }, levelinfo::num_microedges_per_edge( level ) );
   return indexInMacroCell;
}

inline edgedof::EdgeDoFOrientation getOrientattionInNeighboringMacroCell( const EdgeDoFOrientation& orientationInMacroFace,
                                                                          const Face&               face,
                                                                          const uint_t&             neighborCellID,
                                                                          const PrimitiveStorage&   storage )
{
   const Cell&  neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
   const uint_t localFaceID  = neighborCell.getLocalFaceID( face.getID() );

   const auto orientationInMacroCell = edgedof::convertEdgeDoFOrientationFaceToCell(
       orientationInMacroFace,
       neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
       neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
       neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ) );
   return orientationInMacroCell;
}

template < typename ValueType >
inline void getLocalElementDoFIndicesFromCoordinates( const uint_t&                                               level,
                                                      Face&                                                       face,
                                                      const Point3D&                                              coordinates,
                                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcID,
                                                      Point2D&  localCoordinates,
                                                      Matrix2r& transform,
                                                      Point3D&  dofs )
{
   // Get the element local coordinates and DoFs from physical coordinates
   // The local DoFs are sorted in following order
   // x
   // | \
  // 1  0
   // |    \
  // x--2--x

   // Transform absolute coordinates to macro element relative coordinates
   Matrix2r A;
   A( 0, 0 ) = ( face.getCoordinates()[1] - face.getCoordinates()[0] )[0];
   A( 0, 1 ) = ( face.getCoordinates()[2] - face.getCoordinates()[0] )[0];
   A( 1, 0 ) = ( face.getCoordinates()[1] - face.getCoordinates()[0] )[1];
   A( 1, 1 ) = ( face.getCoordinates()[2] - face.getCoordinates()[0] )[1];
   transform = A.inverse();

   Point2D x( coordinates[0] - face.getCoordinates()[0][0], coordinates[1] - face.getCoordinates()[0][1] );

   Point2D xRelMacro = transform * x;

   // Determine lower-left corner index of the quad where the evaluation point lies in
   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );
   real_t hInv    = walberla::real_c( rowsize - 1 );
   real_t h       = walberla::real_c( 1.0 / hInv );

   int binX = static_cast< int >( std::floor( xRelMacro[0] * ( rowsize - 1 ) ) );
   int binY = static_cast< int >( std::floor( xRelMacro[1] * ( rowsize - 1 ) ) );

   if ( binX < 0 )
   {
      binX = 0;
   }

   if ( binY < 0 )
   {
      binY = 0;
   }

   if ( binX >= static_cast< int >( rowsize - 1 ) )
   {
      binX = static_cast< int >( rowsize - 2 );
   }

   if ( binY >= static_cast< int >( rowsize - 1 ) )
   {
      binY = static_cast< int >( rowsize - 2 );
   }

   if ( binX + binY >= static_cast< int >( rowsize - 1 ) )
   {
      int binXDec = ( binX + binY - static_cast< int >( rowsize - 2 ) ) / 2;
      binXDec += ( binX + binY - static_cast< int >( rowsize - 2 ) ) % 2;
      int binYDec = ( binX + binY - static_cast< int >( rowsize - 2 ) ) / 2;
      binX -= binXDec;
      binY -= binYDec;
   }

   uint_t col = uint_c( binX );
   uint_t row = uint_c( binY );

   WALBERLA_ASSERT_LESS( col, rowsize - 1 );
   WALBERLA_ASSERT_LESS( row, rowsize - 1 );
   WALBERLA_ASSERT_LESS( col + row, rowsize - 1, "index.x(): " << col << ", index.y()" << row );

   localCoordinates[0] = xRelMacro[0] - col * h;
   localCoordinates[1] = xRelMacro[1] - row * h;
   localCoordinates *= hInv;

   auto srcData = face.getData( srcID )->getPointer( level );

   transform *= hInv;
   transform.transposeInPlace();

   // decide if up or down triangle
   // clamp to macro-face if the corresponding down-triangle would be out of the macro-face
   // otherwise check floating point distance
   bool upTriangle = ( col + row == rowsize - 2 ) || ( localCoordinates[0] + localCoordinates[1] <= 1.0 );

   if ( upTriangle )
   {
      // Up triangle
      dofs[0] = srcData[edgedof::macroface::diagonalIndex( level, col, row )];
      dofs[1] = srcData[edgedof::macroface::verticalIndex( level, col, row )];
      dofs[2] = srcData[edgedof::macroface::horizontalIndex( level, col, row )];
   }
   else
   {
      // Down triangle
      localCoordinates[0] = 1.0 - localCoordinates[0];
      localCoordinates[1] = 1.0 - localCoordinates[1];
      transform *= -1.0;
      dofs[0] = srcData[edgedof::macroface::diagonalIndex( level, col, row )];
      dofs[1] = srcData[edgedof::macroface::verticalIndex( level, col + 1, row )];
      dofs[2] = srcData[edgedof::macroface::horizontalIndex( level, col, row + 1 )];
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                               Level,
                         Face&                                                       face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceMemoryId,
                         const ValueType&                                            constant )
{
   auto faceData = face.getData( faceMemoryId )->getPointer( Level );

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         faceData[edgedof::macroface::horizontalIndex( Level, it.col(), it.row() )] = constant;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         faceData[edgedof::macroface::verticalIndex( Level, it.col(), it.row() )] = constant;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         faceData[edgedof::macroface::diagonalIndex( Level, it.col(), it.row() )] = constant;
      }
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                                                               Level,
                         Face&                                                                                       face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                                 faceMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >&                  srcIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& ) >& expr )
{
   auto faceData = face.getData( faceMemoryId )->getPointer( Level );

   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   std::vector< ValueType > srcVectorHorizontal( srcIds.size() );
   std::vector< ValueType > srcVectorVertical( srcIds.size() );
   std::vector< ValueType > srcVectorDiagonal( srcIds.size() );

   const Point3D faceBottomLeftCoords  = face.getCoordinates()[0];
   const Point3D faceBottomRightCoords = face.getCoordinates()[1];
   const Point3D faceTopLeftCoords     = face.getCoordinates()[2];

   const Point3D horizontalMicroEdgeOffset =
       ( ( faceBottomRightCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( Level ) ) ) * 0.5;
   const Point3D verticalMicroEdgeOffset =
       ( ( faceTopLeftCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( Level ) ) ) * 0.5;

   Point3D xBlend;

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      const Point3D horizontalMicroEdgePosition =
          faceBottomLeftCoords +
          ( ( real_c( it.col() ) * 2 + 1 ) * horizontalMicroEdgeOffset + ( real_c( it.row() ) * 2 ) * verticalMicroEdgeOffset );
      const Point3D verticalMicroEdgePosition =
          faceBottomLeftCoords +
          ( ( real_c( it.col() ) * 2 ) * horizontalMicroEdgeOffset + ( real_c( it.row() ) * 2 + 1 ) * verticalMicroEdgeOffset );
      const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         for ( uint_t k = 0; k < srcPtr.size(); ++k )
         {
            srcVectorHorizontal[k] = srcPtr[k][edgedof::macroface::horizontalIndex( Level, it.col(), it.row() )];
         }

         face.getGeometryMap()->evalF( horizontalMicroEdgePosition, xBlend );
         faceData[edgedof::macroface::horizontalIndex( Level, it.col(), it.row() )] = expr( xBlend, srcVectorHorizontal );
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         for ( uint_t k = 0; k < srcPtr.size(); ++k )
         {
            srcVectorVertical[k] = srcPtr[k][edgedof::macroface::verticalIndex( Level, it.col(), it.row() )];
         }

         face.getGeometryMap()->evalF( verticalMicroEdgePosition, xBlend );
         faceData[edgedof::macroface::verticalIndex( Level, it.col(), it.row() )] = expr( xBlend, srcVectorVertical );
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         for ( uint_t k = 0; k < srcPtr.size(); ++k )
         {
            srcVectorDiagonal[k] = srcPtr[k][edgedof::macroface::diagonalIndex( Level, it.col(), it.row() )];
         }

         face.getGeometryMap()->evalF( diagonalMicroEdgePosition, xBlend );
         faceData[edgedof::macroface::diagonalIndex( Level, it.col(), it.row() )] = expr( xBlend, srcVectorDiagonal );
      }
   }
}

template < typename ValueType >
inline void swap( const uint_t&                                               level,
                  Face&                                                       face,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstID )
{
   auto srcData = face.getData( srcID );
   auto dstData = face.getData( dstID );
   srcData->swap( *dstData, level );
}

template < typename ValueType >
inline void add( const uint_t&                                               Level,
                 Face&                                                       face,
                 const ValueType&                                            scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   auto dstData = face.getData( dstId )->getPointer( Level );

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      const uint_t idxHorizontal = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
      const uint_t idxVertical   = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
      const uint_t idxDiagonal   = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         dstData[idxHorizontal] += scalar;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         dstData[idxVertical] += scalar;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         dstData[idxDiagonal] += scalar;
      }
   }
}

template < typename ValueType >
inline void add( const uint_t&                                                              Level,
                 Face&                                                                      face,
                 const std::vector< ValueType >&                                            scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" )
   WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" )

   auto dstData = face.getData( dstId )->getPointer( Level );

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      auto tmpHorizontal = static_cast< ValueType >( 0.0 );
      auto tmpVertical   = static_cast< ValueType >( 0.0 );
      auto tmpDiagonal   = static_cast< ValueType >( 0.0 );

      const uint_t idxHorizontal = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
      const uint_t idxVertical   = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
      const uint_t idxDiagonal   = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );

      for ( uint_t i = 0; i < scalars.size(); i++ )
      {
         const ValueType scalar  = scalars[i];
         const auto      srcData = face.getData( srcIds[i] )->getPointer( Level );

         // Do not update horizontal DoFs at bottom
         if ( it.row() != 0 )
         {
            tmpHorizontal += scalar * srcData[idxHorizontal];
         }

         // Do not update vertical DoFs at left border
         if ( it.col() != 0 )
         {
            tmpVertical += scalar * srcData[idxVertical];
         }

         // Do not update diagonal DoFs at diagonal border
         if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
         {
            tmpDiagonal += scalar * srcData[idxDiagonal];
         }
      }

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         dstData[idxHorizontal] += tmpHorizontal;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         dstData[idxVertical] += tmpVertical;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         dstData[idxDiagonal] += tmpDiagonal;
      }
   }
}

template < typename ValueType >
inline void assign( const uint_t&                                                              Level,
                    Face&                                                                      face,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" )
   WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" )

   auto dstData = face.getData( dstId )->getPointer( Level );

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      auto tmpHorizontal = static_cast< ValueType >( 0.0 );
      auto tmpVertical   = static_cast< ValueType >( 0.0 );
      auto tmpDiagonal   = static_cast< ValueType >( 0.0 );

      const uint_t idxHorizontal = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
      const uint_t idxVertical   = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
      const uint_t idxDiagonal   = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );

      for ( uint_t i = 0; i < scalars.size(); i++ )
      {
         const ValueType scalar  = scalars[i];
         const auto      srcData = face.getData( srcIds[i] )->getPointer( Level );

         // Do not update horizontal DoFs at bottom
         if ( it.row() != 0 )
         {
            tmpHorizontal += scalar * srcData[idxHorizontal];
         }

         // Do not update vertical DoFs at left border
         if ( it.col() != 0 )
         {
            tmpVertical += scalar * srcData[idxVertical];
         }

         // Do not update diagonal DoFs at diagonal border
         if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
         {
            tmpDiagonal += scalar * srcData[idxDiagonal];
         }
      }

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         dstData[idxHorizontal] = tmpHorizontal;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         dstData[idxVertical] = tmpVertical;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         dstData[idxDiagonal] = tmpDiagonal;
      }
   }
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                              level,
                             Face&                                                                      face,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   auto dstData = face.getData( dstId )->getPointer( level );

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
   {
      const uint_t idxHorizontal = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
      const uint_t idxVertical   = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
      const uint_t idxDiagonal   = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );

      ValueType tmpHorizontal = face.getData( srcIds[0] )->getPointer( level )[idxHorizontal];
      ValueType tmpVertical   = face.getData( srcIds[0] )->getPointer( level )[idxVertical];
      ValueType tmpDiagonal   = face.getData( srcIds[0] )->getPointer( level )[idxDiagonal];

      for ( uint_t i = 1; i < srcIds.size(); ++i )
      {
         // Do not update horizontal DoFs at bottom
         if ( it.row() != 0 )
         {
            tmpHorizontal *= face.getData( srcIds[i] )->getPointer( level )[idxHorizontal];
         }

         // Do not update vertical DoFs at left border
         if ( it.col() != 0 )
         {
            tmpVertical *= face.getData( srcIds[i] )->getPointer( level )[idxVertical];
         }

         // Do not update diagonal DoFs at diagonal border
         if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
         {
            tmpDiagonal *= face.getData( srcIds[i] )->getPointer( level )[idxDiagonal];
         }
      }

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         dstData[idxHorizontal] = tmpHorizontal;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         dstData[idxVertical] = tmpVertical;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         dstData[idxDiagonal] = tmpDiagonal;
      }
   }
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                               Level,
                      Face&                                                       face,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& lhsId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId )
{
   auto lhsData = face.getData( lhsId )->getPointer( Level );
   auto rhsData = face.getData( rhsId )->getPointer( Level );

   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
         scalarProduct += lhsData[idx] * rhsData[idx];
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
         scalarProduct += lhsData[idx] * rhsData[idx];
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );
         scalarProduct += lhsData[idx] * rhsData[idx];
      }
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline ValueType sum( const uint_t&                                               Level,
                      Face&                                                       face,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dataId,
                      const bool&                                                 absolute )
{
   auto data = face.getData( dataId )->getPointer( Level );

   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
         if ( absolute )
            scalarProduct += std::abs( data[idx] );
         else
            scalarProduct += data[idx];
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
         if ( absolute )
            scalarProduct += std::abs( data[idx] );
         else
            scalarProduct += data[idx];
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );
         if ( absolute )
            scalarProduct += std::abs( data[idx] );
         else
            scalarProduct += data[idx];
      }
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline void enumerate( const uint_t&                                               Level,
                       Face&                                                       face,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                       ValueType&                                                  num )
{
   ValueType* dst            = face.getData( dstId )->getPointer( Level );
   auto       horizontal_num = static_cast< size_t >( num );
   auto       diagonal_num   = static_cast< size_t >( num ) + hyteg::edgedof::levelToFaceSizeAnyEdgeDoF( Level ) -
                       hyteg::levelinfo::num_microedges_per_edge( Level );
   auto vertical_num =
       static_cast< size_t >( num ) +
       ( hyteg::edgedof::levelToFaceSizeAnyEdgeDoF( Level ) - hyteg::levelinfo::num_microedges_per_edge( Level ) ) * 2;
   for ( const auto& it : hyteg::edgedof::macroface::Iterator( Level, 0 ) )
   {
      /// the border edge DoFs belong to the corresponding edges
      if ( it.row() != 0 )
      {
         dst[hyteg::edgedof::macroface::horizontalIndex( Level, it.col(), it.row() )] = ValueType( horizontal_num );
         ++horizontal_num;
         ++num;
      }
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         dst[hyteg::edgedof::macroface::diagonalIndex( Level, it.col(), it.row() )] = ValueType( diagonal_num );
         ++diagonal_num;
         ++num;
      }
      if ( it.col() != 0 )
      {
         dst[hyteg::edgedof::macroface::verticalIndex( Level, it.col(), it.row() )] = ValueType( vertical_num );
         ++vertical_num;
         ++num;
      }
   }
}

inline void apply( const uint_t&                                            Level,
                   Face&                                                    face,
                   const PrimitiveDataID< StencilMemory< real_t >, Face >&  operatorId,
                   const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                   const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                   UpdateType                                               update )
{
   real_t* opr_data = face.getData( operatorId )->getPointer( Level );
   real_t* src      = face.getData( srcId )->getPointer( Level );
   real_t* dst      = face.getData( dstId )->getPointer( Level );

   real_t tmp;

   using namespace edgedof::macroface;

   for ( const auto& it : hyteg::edgedof::macroface::Iterator( Level, 0 ) )
   {
      if ( it.row() != 0 )
      {
         tmp = 0.0;
         for ( auto k : neighborsFromHorizontalEdge )
         {
            tmp += opr_data[edgedof::stencilIndexFromHorizontalEdge( k )] *
                   src[indexFromHorizontalEdge( Level, it.col(), it.row(), k )];
         }
         if ( update == Replace )
         {
            dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] = tmp;
         }
         else if ( update == Add )
         {
            dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] += tmp;
         }
      }
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         tmp = 0.0;
         for ( auto k : neighborsFromDiagonalEdge )
         {
            tmp +=
                opr_data[edgedof::stencilIndexFromDiagonalEdge( k )] * src[indexFromDiagonalEdge( Level, it.col(), it.row(), k )];
         }
         if ( update == Replace )
         {
            dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] = tmp;
         }
         else if ( update == Add )
         {
            dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] += tmp;
         }
      }
      if ( it.col() != 0 )
      {
         tmp = 0.0;
         for ( auto k : neighborsFromVerticalEdge )
         {
            tmp +=
                opr_data[edgedof::stencilIndexFromVerticalEdge( k )] * src[indexFromVerticalEdge( Level, it.col(), it.row(), k )];
         }

         if ( update == Replace )
         {
            dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] = tmp;
         }
         else if ( update == Add )
         {
            dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] += tmp;
         }
      }
   }
}

inline void apply3D( const uint_t&                                                   level,
                     Face&                                                           face,
                     const PrimitiveStorage&                                         storage,
                     const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Face >& operatorId,
                     const PrimitiveDataID< FunctionMemory< real_t >, Face >&        srcId,
                     const PrimitiveDataID< FunctionMemory< real_t >, Face >&        dstId,
                     UpdateType                                                      update )
{
   auto    opr_data = face.getData( operatorId )->getData( level );
   real_t* src      = face.getData( srcId )->getPointer( level );
   real_t* dst      = face.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexInFace : hyteg::edgedof::macroface::Iterator( level, 0 ) )
   {
      std::map< edgedof::EdgeDoFOrientation, real_t > tmpResults = {
          { edgedof::EdgeDoFOrientation::X, real_c( 0 ) },
          { edgedof::EdgeDoFOrientation::Y, real_c( 0 ) },
          { edgedof::EdgeDoFOrientation::XY, real_c( 0 ) },
      };

      for ( const auto& faceCenterOrientation : edgedof::faceLocalEdgeDoFOrientations )
      {
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::X &&
              edgedof::macroface::isHorizontalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::Y &&
              edgedof::macroface::isVerticalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::XY &&
              edgedof::macroface::isDiagonalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;

         for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
         {
            const Cell&  neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
            const uint_t localFaceID  = neighborCell.getLocalFaceID( face.getID() );

            const auto centerIndexInCell =
                getIndexInNeighboringMacroCell( centerIndexInFace, face, neighborCellID, storage, level );
            const auto cellCenterOrientation =
                getOrientattionInNeighboringMacroCell( faceCenterOrientation, face, neighborCellID, storage );

            for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
               for ( const auto& stencilIt : opr_data[neighborCellID][cellCenterOrientation][leafOrientation] )
               {
                  const auto stencilOffset = stencilIt.first;
                  const auto stencilWeight = stencilIt.second;

                  const auto leafOrientationInFace =
                      macrocell::getOrientattionInNeighboringMacroFace( leafOrientation, neighborCell, localFaceID, storage );

                  const auto leafIndexInCell = centerIndexInCell + stencilOffset;
                  const auto leafIndexInFace =
                      leafOrientation == edgedof::EdgeDoFOrientation::XYZ ?
                          macrocell::getIndexInNeighboringMacroFaceXYZ(
                              leafIndexInCell, neighborCell, localFaceID, storage, level ) :
                          macrocell::getIndexInNeighboringMacroFace( leafIndexInCell, neighborCell, localFaceID, storage, level );

                  WALBERLA_ASSERT_LESS_EQUAL( leafIndexInFace.z(), 1 )

                  uint_t leafArrayIndexInFace;
                  if ( algorithms::contains( edgedof::faceLocalEdgeDoFOrientations, leafOrientationInFace ) &&
                       leafIndexInFace.z() == 0 )
                  {
                     leafArrayIndexInFace =
                         edgedof::macroface::index( level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace );
                  }
                  else
                  {
                     leafArrayIndexInFace = edgedof::macroface::index(
                         level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace, neighborCellID );
                  }

                  tmpResults[faceCenterOrientation] += stencilWeight * src[leafArrayIndexInFace];
               }
            }
         }

         const auto dstIdx =
             edgedof::macroface::index( level, centerIndexInFace.x(), centerIndexInFace.y(), faceCenterOrientation );
         if ( update == Replace )
         {
            dst[dstIdx] = tmpResults[faceCenterOrientation];
         }
         else
         {
            dst[dstIdx] += tmpResults[faceCenterOrientation];
         }
      }
   }
}

template < typename ValueType >
inline void
    printFunctionMemory( const uint_t& Level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   ValueType* faceMemory = face.getData( dstId )->getPointer( Level );
   using namespace std;
   cout << setfill( '=' ) << setw( 100 ) << "" << endl;
   cout << face << std::left << setprecision( 1 ) << fixed << setfill( ' ' ) << endl;
   cout << "Horizontal Edge";
   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      if ( it.col() == 0 )
         std::cout << std::endl;
      cout << setw( 5 )
           << faceMemory[hyteg::edgedof::macroface::indexFromHorizontalEdge(
                  Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )]
           << "|";
   }
   cout << endl << "Diagonal Edge";
   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      if ( it.col() == 0 )
         std::cout << std::endl;
      cout << setw( 5 )
           << faceMemory[hyteg::edgedof::macroface::indexFromDiagonalEdge(
                  Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )]
           << "|";
   }
   cout << endl << "Vertical Edge";
   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      if ( it.col() == 0 )
         std::cout << std::endl;
      cout << setw( 5 )
           << faceMemory[hyteg::edgedof::macroface::indexFromVerticalEdge(
                  Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )]
           << "|";
   }
   cout << endl << setfill( '=' ) << setw( 100 ) << "" << endl << setfill( ' ' );
}

template < typename ValueType >
inline ValueType getMaxValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   auto src      = face.getData( srcId )->getPointer( level );
   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
         localMax         = std::max( localMax, src[idx] );
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
         localMax         = std::max( localMax, src[idx] );
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );
         localMax         = std::max( localMax, src[idx] );
      }
   }

   return localMax;
}

template < typename ValueType >
inline ValueType getMinValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   auto src      = face.getData( srcId )->getPointer( level );
   auto localMin = std::numeric_limits< ValueType >::max();

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
         localMin         = std::min( localMin, src[idx] );
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
         localMin         = std::min( localMin, src[idx] );
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );
         localMin         = std::min( localMin, src[idx] );
      }
   }

   return localMin;
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   auto src      = face.getData( srcId )->getPointer( level );
   auto localMax = ValueType( 0.0 );

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
         localMax         = std::max( localMax, std::abs( src[idx] ) );
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
         localMax         = std::max( localMax, std::abs( src[idx] ) );
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );
         localMax         = std::max( localMax, std::abs( src[idx] ) );
      }
   }

   return localMax;
}

inline void
    invertElementwise( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< real_t >, Face >& faceDataID )
{
   real_t* data = face.getData( faceDataID )->getPointer( level );

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
   {
      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         uint_t idx = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
         data[idx]  = real_c( 1.0 ) / data[idx];
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         uint_t idx = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
         data[idx]  = real_c( 1.0 ) / data[idx];
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         uint_t idx = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );
         data[idx]  = real_c( 1.0 ) / data[idx];
      }
   }
}

template < typename ValueType >
inline void createVectorFromFunction( const uint_t&                                               Level,
                                      Face&                                                       face,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Face >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   auto src       = face.getData( srcId )->getPointer( Level );
   auto numerator = face.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
         vec->setValue( uint_c( numerator[idx] ), src[idx] );
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
         vec->setValue( uint_c( numerator[idx] ), src[idx] );
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );
         vec->setValue( uint_c( numerator[idx] ), src[idx] );
      }
   }
}

template < typename ValueType >
inline void createFunctionFromVector( const uint_t&                                               Level,
                                      Face&                                                       face,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Face >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   auto src       = face.getData( srcId )->getPointer( Level );
   auto numerator = face.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
         src[idx]         = vec->getValue( uint_c( numerator[idx] ) );
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
         src[idx]         = vec->getValue( uint_c( numerator[idx] ) );
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );
         src[idx]         = vec->getValue( uint_c( numerator[idx] ) );
      }
   }
}

inline void applyDirichletBC( const uint_t&                                           Level,
                              Face&                                                   face,
                              std::vector< idx_t >&                                   mat,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Face >& numeratorId )
{
   auto numerator = face.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : edgedof::macroface::Iterator( Level, 0 ) )
   {
      // Do not read horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
         mat.push_back( numerator[idx] );
      }

      // Do not read vertical DoFs at left border
      if ( it.col() != 0 )
      {
         const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
         mat.push_back( numerator[idx] );
      }

      // Do not read diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );
         mat.push_back( numerator[idx] );
      }
   }
}

template < typename ValueType >
inline void setBoundaryToZero( const uint_t&                                               level,
                               const Face&                                                 face,
                               const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceDataID )
{
   real_t* data = face.getData( faceDataID )->getPointer( level );

   for ( const auto& idx : edgedof::macroface::BoundaryIterator( level, indexing::FaceBoundaryDirection::BOTTOM_LEFT_TO_RIGHT ) )
   {
      data[edgedof::macroface::horizontalIndex( level, idx.col(), idx.row() )] = 0;
   }

   for ( const auto& idx : edgedof::macroface::BoundaryIterator( level, indexing::FaceBoundaryDirection::LEFT_BOTTOM_TO_TOP ) )
   {
      data[edgedof::macroface::verticalIndex( level, idx.col(), idx.row() )] = 0;
   }

   for ( const auto& idx :
         edgedof::macroface::BoundaryIterator( level, indexing::FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP ) )
   {
      data[edgedof::macroface::diagonalIndex( level, idx.col(), idx.row() )] = 0;
   }
}

} // namespace macroface
} // namespace edgedof
} // namespace hyteg
