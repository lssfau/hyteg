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

// #include "DGFaceIndex.hpp"
// #include "hyteg/bubblefunctionspace/BubbleFaceIndex.hpp"

#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"

namespace hyteg {
namespace facedof {
namespace macroface {

template < typename ValueType >
inline void
    enumerate( const uint_t& Level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId, uint_t& num )
{
   ValueType* dst = face.getData( dstId )->getPointer( Level );
   //the outermost gray cells belong to the edge
   const uint_t grayRowsize   = levelinfo::num_microvertices_per_edge( Level ) - 1;
   uint_t       inner_rowsize = grayRowsize - 2;
   for ( uint_t col = 1; col < grayRowsize - 2; ++col )
   {
      for ( uint_t row = 1; row < inner_rowsize; ++row )
      {
         dst[facedof::macroface::indexFaceFromVertex( Level, col, row, stencilDirection::CELL_GRAY_NE )] = num;
         ++num;
      }
      --inner_rowsize;
   }

   const uint_t blueRowsize = levelinfo::num_microvertices_per_edge( Level ) - 2;
   inner_rowsize            = blueRowsize;
   for ( uint_t col = 0; col < blueRowsize; ++col )
   {
      for ( uint_t row = 0; row < inner_rowsize; ++row )
      {
         dst[facedof::macroface::indexFaceFromGrayFace( Level, col, row, stencilDirection::CELL_BLUE_E )] = num;
         ++num;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                                              Level,
                         Face&                                                                      face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                faceMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcMemoryIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& f ) >& expr )
{
   auto faceMemory = face.getData( faceMemoryId )->getPointer( Level );

   uint_t  rowsize = levelinfo::num_microvertices_per_edge( Level );
   Point3D x, x0, xBlend;

   x0 = face.coords[0];

   Point3D d0 = ( face.coords[1] - face.coords[0] ) / ( walberla::real_c( rowsize - 1 ) );
   Point3D d2 = ( face.coords[2] - face.coords[0] ) / ( walberla::real_c( rowsize - 1 ) );

   uint_t inner_rowsize = rowsize;

   std::vector< ValueType* > srcPtr;
   for ( auto src : srcMemoryIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   std::vector< ValueType > srcVector( srcMemoryIds.size() );

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      x = x0 + 1.0 / 3.0 * ( d0 + d2 ) + walberla::real_c( j ) * d2 + d0;

      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         for ( size_t k = 0; k < srcPtr.size(); ++k )
         {
            srcVector[k] = srcPtr[k][facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )];
         }

         face.getGeometryMap()->evalF( x, xBlend );
         faceMemory[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )] =
             expr( xBlend, srcVector );
         x += d0;
      }
      --inner_rowsize;
   }

   inner_rowsize = rowsize;

   // blue cells
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      x = x0 + 2.0 / 3.0 * ( d0 + d2 ) + walberla::real_c( j ) * d2;

      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         for ( size_t k = 0; k < srcPtr.size(); ++k )
         {
            srcVector[k] = srcPtr[k][facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )];
         }

         face.getGeometryMap()->evalF( x, xBlend );
         faceMemory[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )] =
             expr( xBlend, srcVector );
         x += d0;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                                              Level,
                 Face&                                                                      face,
                 const std::vector< ValueType >&                                            scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   size_t inner_rowsize = rowsize;

   auto dst = face.getData( dstId )->getPointer( Level );

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C );

         ValueType tmp = 0.0;

         for ( uint_t k = 0; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * face.getData( srcIds[k] )->getPointer( Level )[cellIndex];
         }

         dst[cellIndex] += tmp;
      }
      --inner_rowsize;
   }

   inner_rowsize = rowsize;

   // blue cells
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C );

         ValueType tmp = 0.0;

         for ( uint_t k = 0; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * face.getData( srcIds[k] )->getPointer( Level )[cellIndex];
         }

         dst[cellIndex] += tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void assign( const uint_t&                                                              Level,
                    Face&                                                                      face,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   size_t inner_rowsize = rowsize;

   auto                      dst = face.getData( dstId )->getPointer( Level );
   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }
   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C );

         ValueType tmp = scalars[0] * srcPtr[0][cellIndex];

         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * srcPtr[k][cellIndex];
         }

         dst[cellIndex] = tmp;
      }
      --inner_rowsize;
   }

   inner_rowsize = rowsize;

   // blue cells
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C );

         ValueType tmp = scalars[0] * srcPtr[0][cellIndex];

         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * srcPtr[k][cellIndex];
         }

         dst[cellIndex] = tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline ValueType getMaxValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   size_t inner_rowsize = rowsize;

   auto   src      = face.getData( srcId )->getPointer( level );
   ValueType localMax = -std::numeric_limits< ValueType >::max();

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C );
         localMax       = std::max( localMax, src[cellIndex] );
      }
      --inner_rowsize;
   }

   // blue cells
   inner_rowsize = rowsize;
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C );
         localMax       = std::max( localMax, src[cellIndex] );
      }
      --inner_rowsize;
   }

   return localMax;
}

template < typename ValueType >
inline ValueType getMinValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   size_t inner_rowsize = rowsize;

   auto   src      = face.getData( srcId )->getPointer( level );
   ValueType localMin = std::numeric_limits< ValueType >::max();

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C );
         localMin       = std::min( localMin, src[cellIndex] );
      }
      --inner_rowsize;
   }

   // blue cells
   inner_rowsize = rowsize;
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C );
         localMin       = std::min( localMin, src[cellIndex] );
      }
      --inner_rowsize;
   }

   return localMin;
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   size_t inner_rowsize = rowsize;

   auto   src      = face.getData( srcId )->getPointer( level );
   ValueType localMax = ValueType( 0.0 );

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C );
         localMax       = std::max( localMax, std::abs( src[cellIndex] ) );
      }
      --inner_rowsize;
   }

   // blue cells
   inner_rowsize = rowsize;
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C );
         localMax       = std::max( localMax, std::abs( src[cellIndex] ) );
      }
      --inner_rowsize;
   }

   return localMax;
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                              level,
                             Face&                                                                      face,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   size_t inner_rowsize = rowsize;

   ValueType*                dstPtr = face.getData( dstId )->getPointer( level );
   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( level ) );
   }

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C );

         ValueType tmp = srcPtr[0][cellIndex];

         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp *= srcPtr[k][cellIndex];
         }

         dstPtr[cellIndex] = tmp;
      }
      --inner_rowsize;
   }

   inner_rowsize = rowsize;

   // blue cells
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C );

         ValueType tmp = srcPtr[0][cellIndex];

         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp *= srcPtr[k][cellIndex];
         }

         dstPtr[cellIndex] = tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                               level,
                 Face&                                                       face,
                 const ValueType                                             scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   ValueType* dstPtr = face.getData( dstId )->getPointer( level );

   size_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   size_t inner_rowsize = rowsize;

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C );
         dstPtr[cellIndex] += scalar;
      }
      --inner_rowsize;
   }

   // blue cells
   inner_rowsize = rowsize;
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C );
         dstPtr[cellIndex] += scalar;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                               level,
                      Face&                                                       face,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& lhsMemoryId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsMemoryId )
{
   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   ValueType* lhsPtr = face.getData( lhsMemoryId )->getPointer( level );
   ValueType* rhsPtr = face.getData( rhsMemoryId )->getPointer( level );

   size_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   size_t inner_rowsize = rowsize;

   // gray cells
   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C );
         scalarProduct += lhsPtr[cellIndex] * rhsPtr[cellIndex];
      }
      --inner_rowsize;
   }

   // blue cells
   inner_rowsize = rowsize;
   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         auto cellIndex = facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C );
         scalarProduct += lhsPtr[cellIndex] * rhsPtr[cellIndex];
      }
      --inner_rowsize;
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline void swap( const uint_t&                                               level,
                  Face&                                                       face,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstID )
{
   face.getData( srcID )->swap( *face.getData( dstID ), level );
}

} // namespace macroface
} // namespace facedof
} //namespace hyteg
