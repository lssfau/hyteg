/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
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

#include <algorithm>
#include <cstdlib>

#include "core/debug/all.h"
#include "core/math/KahanSummation.h"

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {

using indexing::Index;
using walberla::real_c;
using walberla::uint_t;

inline indexing::Index getIndexInNeighboringMacroCell( const indexing::Index&  vertexDoFIndexInMacroFace,
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

   Index indexInMacroCell = indexing::basisConversion(
       vertexDoFIndexInMacroFace, localVertexIDsAtCell, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ) );
   return indexInMacroCell;
}

inline Point3D coordinateFromIndex( const uint_t& Level, const Face& face, const Index& index )
{
   const real_t  stepFrequency = 1.0 / real_c(levelinfo::num_microedges_per_edge( Level ));
   const Point3D xStep         = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) * stepFrequency;
   const Point3D yStep         = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) * stepFrequency;
   return face.getCoordinates()[0] + xStep * real_c( index.x() ) + yStep * real_c( index.y() );
}

template < typename ValueType >
inline ValueType assembleLocal( const uint_t&                            Level,
                                idx_t                                    i,
                                idx_t                                    j,
                                const Matrix3r&                          localMatrix,
                                double*                                  src,
                                double*                                  coeff,
                                const std::array< stencilDirection, 3 >& vertices,
                                const std::array< uint_t, 3 >&           idx )
{
   ValueType meanCoeff = 1.0 / 3.0 *
                         ( coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )] );

   ValueType tmp;
   tmp = localMatrix( idx[0], idx[0] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
         localMatrix( idx[0], idx[1] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
         localMatrix( idx[0], idx[2] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )];
   return meanCoeff * tmp;
}

template < typename ValueType >
inline void assembleLocalStencil( uint_t                                   Level,
                                  idx_t                                    i,
                                  idx_t                                    j,
                                  const Matrix3r&                          localMatrix,
                                  double*                                  opr_data,
                                  double*                                  coeff,
                                  const std::array< stencilDirection, 3 >& vertices,
                                  const std::array< uint_t, 3 >&           idx )
{
   ValueType meanCoeff = 1.0 / 3.0 *
                         ( coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )] );

   opr_data[vertexdof::stencilIndexFromVertex( vertices[0] )] += meanCoeff * localMatrix( idx[0], idx[0] );
   opr_data[vertexdof::stencilIndexFromVertex( vertices[1] )] += meanCoeff * localMatrix( idx[0], idx[1] );
   opr_data[vertexdof::stencilIndexFromVertex( vertices[2] )] += meanCoeff * localMatrix( idx[0], idx[2] );
}

template < typename ValueType >
inline ValueType assembleLocalDG( const uint_t&                            Level,
                                  uint_t                                   i,
                                  uint_t                                   j,
                                  const Matrix3r&                          localMatrix,
                                  double*                                  src,
                                  const std::array< stencilDirection, 3 >& vertices,
                                  const std::array< uint_t, 3 >&           idx )
{
   ValueType tmp;
   tmp = localMatrix( idx[0], idx[0] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
         localMatrix( idx[0], idx[1] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
         localMatrix( idx[0], idx[2] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )];
   return tmp;
}

template < typename ValueType >
inline void getLocalElementDoFIndicesFromCoordinates( const uint_t&                                               level,
                                                      const Face&                                                 face,
                                                      const Point3D                                               coordinates,
                                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcID,
                                                      Point2D&  localCoordinates,
                                                      Matrix2r& transform,
                                                      Point3D&  dofs )
{
   /** Get the element local coordinates and DoFs from physical coordinates
    * The local DoFs are sorted in following order
    * 2
    * | \
    * |  \
    * 0 - 1
   **/

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
   int    binX    = static_cast< int >( std::floor( xRelMacro[0] * ( rowsize - 1 ) ) );
   int    binY    = static_cast< int >( std::floor( xRelMacro[1] * ( rowsize - 1 ) ) );

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

   Index index;
   index.x() = binX;
   index.y() = binY;

   WALBERLA_ASSERT_LESS( index.x(), rowsize - 1 );
   WALBERLA_ASSERT_LESS( index.y(), rowsize - 1 );
   WALBERLA_ASSERT_LESS( index.x() + index.y(), rowsize - 1, "index.x(): " << index.x() << ", index.y()" << index.y() );

   localCoordinates[0] = xRelMacro[0] - real_c( index.x() ) * h;
   localCoordinates[1] = xRelMacro[1] - real_c( index.y() ) * h;
   localCoordinates *= hInv;

   auto srcData = face.getData( srcID )->getPointer( level );

   transform *= hInv;
   transform.transposeInPlace();

   // decide if up or down triangle
   // clamp to macro-face if the corresponding down-triangle would be out of the macro-face
   // otherwise check floating point distance
   bool upTriangle = ( index.x() + index.y() == rowsize - 2 ) || ( localCoordinates[0] + localCoordinates[1] <= 1.0 );

   // Combine solution from linear basis functions depending on triangle orientation
   if ( upTriangle )
   {
      dofs[0] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_C )];
      dofs[1] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_E )];
      dofs[2] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_N )];
   }
   else
   {
      WALBERLA_ASSERT_LESS( index.x(), rowsize - 2 );
      WALBERLA_ASSERT_LESS( index.y(), rowsize - 2 );
      WALBERLA_ASSERT_LESS( index.x() + index.y(), rowsize - 2 );
      index.y() += 1;
      localCoordinates[0] = 1.0 - localCoordinates[0];
      localCoordinates[1] = 1.0 - localCoordinates[1];
      transform *= -1.0;
      dofs[0] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_E )];
      dofs[1] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_C )];
      dofs[2] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_SE )];
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                               Level,
                         Face&                                                       face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceMemoryId,
                         const ValueType&                                            scalar,
                         const uint_t&                                               offset = 1 )
{
   ValueType* faceData = face.getData( faceMemoryId )->getPointer( Level );

   for ( const auto& it : vertexdof::macroface::Iterator( Level, offset ) )
   {
      const uint_t idx = vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(), stencilDirection::VERTEX_C );
      faceData[idx]    = scalar;
   }
}

template < typename ValueType >
inline void interpolate( const std::shared_ptr< PrimitiveStorage >&                                                  storage,
                         const uint_t&                                                                               Level,
                         Face&                                                                                       face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                                 faceMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >&                  srcIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& ) >& expr,
                         const uint_t&                                                                               offset = 1 )
{
   ValueType* faceData = face.getData( faceMemoryId )->getPointer( Level );

   std::vector< ValueType* > srcPtr;
   for ( const auto& src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   std::vector< ValueType > srcVector( srcIds.size() );

   for ( const auto& it : vertexdof::macroface::Iterator( Level, offset ) )
   {
      const auto   coordinate = micromesh::microVertexPosition( storage, face.getID(), Level, it );
      const uint_t idx        = vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(), stencilDirection::VERTEX_C );

      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
         srcVector[k] = srcPtr[k][idx];
      }

      faceData[idx] = expr( coordinate, srcVector );
   }
}

template < typename ValueType >
inline real_t evaluate( const uint_t&                                               level,
                        const Face&                                                 face,
                        const Point3D                                               coordinates,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcID )
{
   Point2D  localCoordinates;
   Matrix2r transform;
   Point3D  localDoFs;
   getLocalElementDoFIndicesFromCoordinates( level, face, coordinates, srcID, localCoordinates, transform, localDoFs );
   real_t value = ( 1.0 - localCoordinates[0] - localCoordinates[1] ) * localDoFs[0];
   value += localCoordinates[0] * localDoFs[1];
   value += localCoordinates[1] * localDoFs[2];

   return value;
}

template < typename ValueType >
inline void evaluateGradient( const uint_t&                                               level,
                              Face&                                                       face,
                              const Point3D                                               coordinates,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcID,
                              Point3D&                                                    gradient )
{
   Point2D  localCoordinates;
   Matrix2r transform;
   Point3D  localDoFs;
   getLocalElementDoFIndicesFromCoordinates( level, face, coordinates, srcID, localCoordinates, transform, localDoFs );
   Point2D gradient_;
   gradient_[0] = localDoFs[1] - localDoFs[0];
   gradient_[1] = localDoFs[2] - localDoFs[0];
   gradient_    = transform * gradient_;
   gradient[0]  = gradient_[0];
   gradient[1]  = gradient_[1];
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
inline void assign( const uint_t&                                                              Level,
                    Face&                                                                      face,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType*                dst = face.getData( dstId )->getPointer( Level );
   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }
   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         ValueType tmp = scalars[0] * srcPtr[0][vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * srcPtr[k][vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         }
         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                               level,
                 const Face&                                                 face,
                 const ValueType&                                            scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   ValueType* dstPtr = face.getData( dstId )->getPointer( level );

   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         dstPtr[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] += scalar;
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
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType*                dstPtr = face.getData( dstId )->getPointer( Level );
   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         auto tmp = ValueType( 0 );

         for ( uint_t k = 0; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * srcPtr[k][vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         }
         dstPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                              Level,
                             Face&                                                                      face,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType*                dstPtr = face.getData( dstId )->getPointer( Level );
   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   for ( idx_t j = 1; j < rowsize - 2; ++j )
   {
      for ( idx_t i = 1; i < inner_rowsize - 2; ++i )
      {
         const uint_t idx = vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C );
         ValueType    tmp = srcPtr[0][idx];

         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp *= srcPtr[k][idx];
         }
         dstPtr[idx] = tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                               Level,
                      Face&                                                       face,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& lhsId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                      const uint_t&                                               offset = 1 )
{
   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   uint_t rowsizeY = levelinfo::num_microvertices_per_edge( Level );
   uint_t rowsizeX;

   ValueType* lhsPtr = face.getData( lhsId )->getPointer( Level );
   ValueType* rhsPtr = face.getData( rhsId )->getPointer( Level );

   for ( uint_t j = offset; j < rowsizeY - offset; ++j )
   {
      rowsizeX = rowsizeY - j;
      for ( uint_t i = offset; i < rowsizeX - offset; ++i )
      {
         scalarProduct += lhsPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] *
                          rhsPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
      }
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline ValueType sum( const uint_t&                                               level,
                      const Face&                                                 face,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dataID,
                      const bool&                                                 absolute )
{
   ValueType* faceData = face.getData( dataID )->getPointer( level );
   auto       sum      = ValueType( 0 );
   for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
   {
      const uint_t idx = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C );
      if ( absolute )
      {
         sum += std::abs( faceData[idx] );
      }
      else
      {
         sum += faceData[idx];
      }
   }
   return sum;
}

template < typename ValueType >
inline void apply( const uint_t&                                               Level,
                   Face&                                                       face,
                   const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                   UpdateType                                                  update )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType* opr_data = face.getData( operatorId )->getPointer( Level );
   ValueType* src      = face.getData( srcId )->getPointer( Level );
   ValueType* dst      = face.getData( dstId )->getPointer( Level );

   ValueType tmp = real_c( 0 );

   for ( int j = 1; j < rowsize - 2; ++j )
   {
      for ( int i = 1; i < inner_rowsize - 2; ++i )
      {
         if ( face.getNumNeighborCells() == 0 )
         {
            static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
            tmp = real_c( 0 );
            for ( const auto direction : vertexdof::macroface::neighborsWithCenter )
            {
               tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }
         else if ( face.getNumNeighborCells() == 1 )
         {
            tmp = real_c( 0 );
            for ( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter )
            {
               tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }
         else if ( face.getNumNeighborCells() == 2 )
         {
            tmp = real_c( 0 );
            for ( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter )
            {
               tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }

         WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

         if ( update == Replace )
         {
            dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
         }
         else
         {
            dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
         }
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void apply3D( const uint_t&                                                   Level,
                     Face&                                                           face,
                     const PrimitiveStorage&                                         storage,
                     const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Face >& operatorId,
                     const PrimitiveDataID< FunctionMemory< ValueType >, Face >&     srcId,
                     const PrimitiveDataID< FunctionMemory< ValueType >, Face >&     dstId,
                     UpdateType                                                      update )
{
   auto       opr_data = face.getData( operatorId )->getData( Level );
   ValueType* src      = face.getData( srcId )->getPointer( Level );
   ValueType* dst      = face.getData( dstId )->getPointer( Level );

   for ( const auto& idxIt : Iterator( Level, 1 ) )
   {
      ValueType tmp = real_c( 0 );

      for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
      {
         auto neighborCell = storage.getCell( face.neighborCells().at( neighborCellIdx ) );
         auto centerIndexInCell =
             vertexdof::macroface::getIndexInNeighboringMacroCell( idxIt, face, neighborCellIdx, storage, Level );
         for ( const auto& stencilIt : opr_data[neighborCellIdx] )
         {
            auto  weight = stencilIt.second;
            Index leafIndexInMacroCell( centerIndexInCell );
            leafIndexInMacroCell += stencilIt.first;
            auto leafIndexInMacroFace = macrocell::getIndexInNeighboringMacroFace(
                leafIndexInMacroCell, *neighborCell, neighborCell->getLocalFaceID( face.getID() ), storage, Level );

            uint_t leafArrayIndexInMacroFace;
            if ( leafIndexInMacroFace.z() == 0 )
            {
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y() );
            }
            else
            {
               WALBERLA_ASSERT_EQUAL( leafIndexInMacroFace.z(), 1 );
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y(), neighborCellIdx );
            }

            tmp += weight * src[leafArrayIndexInMacroFace];
         }
      }

      if ( update == Replace )
      {
         dst[vertexdof::macroface::index( Level, idxIt.x(), idxIt.y() )] = tmp;
      }
      else if ( update == Add )
      {
         dst[vertexdof::macroface::index( Level, idxIt.x(), idxIt.y() )] += tmp;
      }
   }
}

template < typename ValueType >
inline void smooth_gs( const uint_t&                                               Level,
                       Face&                                                       face,
                       const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto opr_data = face.getData( operatorId )->getPointer( Level );
   auto dst      = face.getData( dstId )->getPointer( Level );
   auto rhs      = face.getData( rhsId )->getPointer( Level );

   const auto invCenterWeight = 1.0 / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];

   ValueType tmp;

   for ( int j = 1; j < rowsize - 2; ++j )
   {
      for ( int i = 1; i < inner_rowsize - 2; ++i )
      {
         tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         if ( face.getNumNeighborCells() == 0 )
         {
            for ( const auto direction : vertexdof::macroface::neighborsWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }
         else if ( face.getNumNeighborCells() == 1 )
         {
            for ( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }
         else if ( face.getNumNeighborCells() == 2 )
         {
            for ( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }

         WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp * invCenterWeight;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void smooth_sor( const uint_t&                                               Level,
                        Face&                                                       face,
                        const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                        ValueType                                                   relax )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto opr_data = face.getData( operatorId )->getPointer( Level );
   auto dst      = face.getData( dstId )->getPointer( Level );
   auto rhs      = face.getData( rhsId )->getPointer( Level );

   const auto invCenterWeight = 1.0 / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];

   ValueType tmp;

   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         if ( face.getNumNeighborCells() == 0 )
         {
            for ( const auto direction : vertexdof::macroface::neighborsWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }
         else if ( face.getNumNeighborCells() == 1 )
         {
            for ( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }
         else if ( face.getNumNeighborCells() == 2 )
         {
            for ( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }

         WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] =
             ( 1.0 - relax ) * dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] +
             relax * tmp * invCenterWeight;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void smoothSOR3D( const uint_t&                                                   Level,
                         Face&                                                           face,
                         const PrimitiveStorage&                                         storage,
                         const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Face >& operatorId,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >&     dstId,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >&     rhsId,
                         ValueType                                                       relax )
{
   auto       opr_data = face.getData( operatorId )->getData( Level );
   ValueType* rhs      = face.getData( rhsId )->getPointer( Level );
   ValueType* dst      = face.getData( dstId )->getPointer( Level );

   real_t centerWeight = real_c( 0 );
   for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
   {
      centerWeight += opr_data[neighborCellIdx][{ 0, 0, 0 }];
   }
   const auto invCenterWeight = 1.0 / centerWeight;

   for ( const auto& idxIt : Iterator( Level, 1 ) )
   {
      ValueType tmp = rhs[vertexdof::macroface::index( Level, idxIt.x(), idxIt.y() )];

      for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
      {
         auto  neighborCell = storage.getCell( face.neighborCells().at( neighborCellIdx ) );
         Index centerIndexInCell =
             vertexdof::macroface::getIndexInNeighboringMacroCell( idxIt, face, neighborCellIdx, storage, Level );
         for ( const auto& stencilIt : opr_data[neighborCellIdx] )
         {
            if ( stencilIt.first == indexing::Index( { 0, 0, 0 } ) )
               continue;

            auto  weight = stencilIt.second;
            Index leafIndexInMacroCell( centerIndexInCell );
            leafIndexInMacroCell += stencilIt.first;
            Index leafIndexInMacroFace = macrocell::getIndexInNeighboringMacroFace(
                leafIndexInMacroCell, *neighborCell, neighborCell->getLocalFaceID( face.getID() ), storage, Level );

            uint_t leafArrayIndexInMacroFace;
            if ( leafIndexInMacroFace.z() == 0 )
            {
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y() );
            }
            else
            {
               WALBERLA_ASSERT_EQUAL( leafIndexInMacroFace.z(), 1 );
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y(), neighborCellIdx );
            }

            tmp -= weight * dst[leafArrayIndexInMacroFace];
         }
      }

      dst[vertexdof::macroface::index( Level, idxIt.x(), idxIt.y() )] =
          ( 1.0 - relax ) * dst[vertexdof::macroface::index( Level, idxIt.x(), idxIt.y() )] + relax * tmp * invCenterWeight;
   }
}

template < typename ValueType >
inline void smooth_jac( const uint_t&                                               Level,
                        Face&                                                       face,
                        const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& tmpId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto opr_data = face.getData( operatorId )->getPointer( Level );
   auto dst      = face.getData( dstId )->getPointer( Level );
   auto rhs      = face.getData( rhsId )->getPointer( Level );
   auto tmpVar   = face.getData( tmpId )->getPointer( Level );

   ValueType tmp;

   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         for ( auto neighbor : vertexdof::macroface::neighborsWithoutCenter )
         {
            tmp -= opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
                   tmpVar[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
         }

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] =
             tmp / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
      }
      --inner_rowsize;
   }
}

/// Checks if a given index is a the boundary of the face
/// \param index The index which should be checked
/// \param length Size of the triangle in the first dimension
inline bool is_boundary( uint_t index, uint_t length )
{
   if ( index < length )
      return true;
   while ( index >= length )
   {
      index -= length;
      length--;
   }
   return ( index == 0 || index == ( length - 1 ) );
}

template < typename ValueType >
inline void enumerate( const uint_t&                                               Level,
                       Face&                                                       face,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                       ValueType&                                                  num )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   uint_t mr = 1 + rowsize;

   ValueType* dstPtr = face.getData( dstId )->getPointer( Level );

   for ( uint_t i = 0; i < rowsize - 3; ++i )
   {
      for ( uint_t j = 0; j < inner_rowsize - 3; ++j )
      {
         dstPtr[mr] = num;
         num++;

         mr += 1;
      }

      mr += 2;
      --inner_rowsize;
   }
}

template < typename ValueType >
inline ValueType
    getMaxDoFValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   auto src      = face.getData( srcId )->getPointer( level );
   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         localMax = std::max( localMax, src[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] );
      }
      --inner_rowsize;
   }

   return localMax;
}

template < typename ValueType >
inline ValueType
    getMaxDoFMagnitude( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   auto src      = face.getData( srcId )->getPointer( level );
   auto localMax = ValueType( 0.0 );

   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         localMax = std::max( localMax,
                              std::abs( src[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] ) );
      }
      --inner_rowsize;
   }

   return localMax;
}

template < typename ValueType >
inline ValueType
    getMinDoFValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   auto src      = face.getData( srcId )->getPointer( level );
   auto localMin = std::numeric_limits< ValueType >::max();

   for ( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         localMin = std::min( localMin, src[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] );
      }
      --inner_rowsize;
   }

   return localMin;
}

template < typename ValueType >
inline ValueType reduce( uint_t                                                      level,
                         std::function< ValueType( ValueType, ValueType ) >&         reduceOperation,
                         ValueType                                                   initialValue,
                         Face&                                                       face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   auto src = face.getData( srcId )->getPointer( level );
   for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
   {
      initialValue = reduceOperation(
          initialValue, src[vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C )] );
   }
   return initialValue;
}

template < typename ValueType >
inline void saveOperator( const uint_t&                                              Level,
                          Face&                                                      face,
                          const PrimitiveDataID< StencilMemory< ValueType >, Face >& operatorId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Face >&    srcId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Face >&    dstId,
                          const std::shared_ptr< SparseMatrixProxy >&                mat )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto opr_data = face.getData( operatorId )->getPointer( Level );
   auto src      = face.getData( srcId )->getPointer( Level );
   auto dst      = face.getData( dstId )->getPointer( Level );

   for ( idx_t i = 1; i < rowsize - 2; ++i )
   {
      for ( idx_t j = 1; j < inner_rowsize - 2; ++j )
      {
         idx_t srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         idx_t dstInt = dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         mat->addValue(
             uint_c( dstInt ), uint_c( srcInt ), opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] );

         if ( face.getNumNeighborCells() == 0 )
         {
            for ( const auto& neighbor : vertexdof::macroface::neighborsWithoutCenter )
            {
               srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
               mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[vertexdof::stencilIndexFromVertex( neighbor )] );
            }
         }
         else if ( face.getNumNeighborCells() == 1 )
         {
            for ( const auto neighbor : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
            {
               srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
               mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[vertexdof::stencilIndexFromVertex( neighbor )] );
            }
         }
         else if ( face.getNumNeighborCells() == 2 )
         {
            for ( const auto neighbor : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
            {
               srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
               mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[vertexdof::stencilIndexFromVertex( neighbor )] );
            }
         }
      }
      --inner_rowsize;
   }
}

inline void saveIdentityOperator( const uint_t&                                           Level,
                                  Face&                                                   face,
                                  const PrimitiveDataID< FunctionMemory< idx_t >, Face >& dstId,
                                  const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto dst = face.getData( dstId )->getPointer( Level );

   for ( idx_t i = 1; i < rowsize - 2; ++i )
   {
      for ( idx_t j = 1; j < inner_rowsize - 2; ++j )
      {
         idx_t dstInt = dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
      }
      --inner_rowsize;
   }
}

inline void saveOperator3D( const uint_t&                                                   Level,
                            Face&                                                           face,
                            const PrimitiveStorage&                                         storage,
                            const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Face >& operatorId,
                            const PrimitiveDataID< FunctionMemory< idx_t >, Face >&         srcId,
                            const PrimitiveDataID< FunctionMemory< idx_t >, Face >&         dstId,
                            const std::shared_ptr< SparseMatrixProxy >&                     mat )
{
   auto opr_data = face.getData( operatorId )->getData( Level );
   auto src      = face.getData( srcId )->getPointer( Level );
   auto dst      = face.getData( dstId )->getPointer( Level );

   for ( const auto& idxIt : Iterator( Level, 1 ) )
   {
      const idx_t dstInt = dst[vertexdof::macroface::index( Level, idxIt.x(), idxIt.y() )];

      for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
      {
         auto  neighborCell = storage.getCell( face.neighborCells().at( neighborCellIdx ) );
         Index centerIndexInCell =
             vertexdof::macroface::getIndexInNeighboringMacroCell( idxIt, face, neighborCellIdx, storage, Level );
         for ( const auto& stencilIt : opr_data[neighborCellIdx] )
         {
            auto  weight = stencilIt.second;
            Index leafIndexInMacroCell( centerIndexInCell + stencilIt.first );
            auto  leafIndexInMacroFace = macrocell::getIndexInNeighboringMacroFace(
                leafIndexInMacroCell, *neighborCell, neighborCell->getLocalFaceID( face.getID() ), storage, Level );

            uint_t leafArrayIndexInMacroFace;
            if ( leafIndexInMacroFace.z() == 0 )
            {
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y() );
            }
            else
            {
               WALBERLA_ASSERT_EQUAL( leafIndexInMacroFace.z(), 1 );
               leafArrayIndexInMacroFace =
                   vertexdof::macroface::index( Level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y(), neighborCellIdx );
            }

            const idx_t srcInt = src[leafArrayIndexInMacroFace];
            mat->addValue( uint_c( dstInt ), uint_c( srcInt ), weight );
         }
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
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto src       = face.getData( srcId )->getPointer( Level );
   auto numerator = face.getData( numeratorId )->getPointer( Level );

   for ( idx_t i = 1; i < rowsize - 2; ++i )
   {
      for ( idx_t j = 1; j < inner_rowsize - 2; ++j )
      {
         idx_t numeratorInt = numerator[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         vec->setValue( uint_c( numeratorInt ),
                        src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] );
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void createFunctionFromVector( const uint_t&                                               Level,
                                      Face&                                                       face,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Face >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto src       = face.getData( srcId )->getPointer( Level );
   auto numerator = face.getData( numeratorId )->getPointer( Level );

   for ( idx_t i = 1; i < rowsize - 2; ++i )
   {
      for ( idx_t j = 1; j < inner_rowsize - 2; ++j )
      {
         idx_t numeratorInt = numerator[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] =
             vec->getValue( uint_c( numeratorInt ) );
      }
      --inner_rowsize;
   }
}

inline void applyDirichletBC( const uint_t&                                           level,
                              Face&                                                   face,
                              std::vector< idx_t >&                                   mat,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Face >& numeratorId )
{
   for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
   {
      mat.push_back( face.getData( numeratorId )->getPointer( level )[vertexdof::macroface::index( level, it.x(), it.y() )] );
   }
}

template < typename ValueType >
inline void
    printFunctionMemory( const uint_t Level, const Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   ValueType* faceMemory = face.getData( dstId )->getPointer( Level );
   using namespace std;
   cout << setfill( '=' ) << setw( 100 ) << "" << endl;
   cout << face << std::left << setprecision( 1 ) << fixed << setfill( ' ' ) << endl << "Vertex DoFs: ";
   for ( const auto& it : vertexdof::macroface::Iterator( Level, 0 ) )
   {
      if ( it.x() == 0 )
         std::cout << std::endl;
      cout << setw( 5 )
           << faceMemory[hyteg::vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(), stencilDirection::VERTEX_C )]
           << "|";
   }
   cout << endl << setfill( '=' ) << setw( 100 ) << "" << endl << setfill( ' ' );
}

} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
