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

#include "core/DataTypes.h"
#include "core/Format.hpp"
#include "core/debug/all.h"
#include "core/math/Matrix3.h"

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/types/types.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg::vertexdof::macrocell {

using walberla::int64_c;
using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using indexing::Index;

inline indexing::Index getIndexInNeighboringMacroFace( const indexing::Index&  vertexDoFIndexInMacroCell,
                                                       const Cell&             cell,
                                                       const uint_t&           neighborFaceID,
                                                       const PrimitiveStorage& storage,
                                                       const uint_t&           level )
{
   const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 3, 4 >(
       { cell.getFaceLocalVertexToCellLocalVertexMaps().at( neighborFaceID ).at( 0 ),
         cell.getFaceLocalVertexToCellLocalVertexMaps().at( neighborFaceID ).at( 1 ),
         cell.getFaceLocalVertexToCellLocalVertexMaps().at( neighborFaceID ).at( 2 ) } );

   const auto indexInMacroFace = indexing::basisConversion(
       vertexDoFIndexInMacroCell, { 0, 1, 2, 3 }, localVertexIDsAtCell, levelinfo::num_microvertices_per_edge( level ) );
   return indexInMacroFace;
}

inline Point3D coordinateFromIndex( const uint_t& level, const Cell& cell, const Index& index )
{
   const real_t  stepFrequency = real_c( 1.0 ) / real_c( levelinfo::num_microedges_per_edge( level ) );
   const Point3D xStep         = ( cell.getCoordinates()[1] - cell.getCoordinates()[0] ) * stepFrequency;
   const Point3D yStep         = ( cell.getCoordinates()[2] - cell.getCoordinates()[0] ) * stepFrequency;
   const Point3D zStep         = ( cell.getCoordinates()[3] - cell.getCoordinates()[0] ) * stepFrequency;
   return cell.getCoordinates()[0] + xStep * real_c( index.x() ) + yStep * real_c( index.y() ) + zStep * real_c( index.z() );
}

template < typename ValueType >
inline void interpolate( const uint_t&                                               level,
                         const Cell&                                                 cell,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& cellMemoryId,
                         const ValueType&                                            scalar,
                         const uint_t&                                               offset = 1 )
{
   ValueType* cellData = cell.getData( cellMemoryId )->getPointer( level );

   for ( const auto& it : vertexdof::macrocell::Iterator( level, offset ) )
   {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
      cellData[idx]    = scalar;
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                                                               level,
                         const Cell&                                                                                 cell,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&                                 cellMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >&                  srcIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& ) >& expr,
                         const uint_t&                                                                               offset = 1 )
{
   ValueType* cellData = cell.getData( cellMemoryId )->getPointer( level );

   std::vector< ValueType* > srcPtr;
   for ( const auto& src : srcIds )
   {
      srcPtr.push_back( cell.getData( src )->getPointer( level ) );
   }

   std::vector< ValueType > srcVector( srcIds.size() );

   Point3D xBlend;

   for ( const auto& it : vertexdof::macrocell::Iterator( level, offset ) )
   {
      const Point3D coordinate = coordinateFromIndex( level, cell, it );
      const uint_t  idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
         srcVector[k] = srcPtr[k][idx];
      }
      cell.getGeometryMap()->evalF( coordinate, xBlend );
      cellData[idx] = expr( xBlend, srcVector );
   }
}

template < typename ValueType >
inline void swap( const uint_t&                                               level,
                  Cell&                                                       cell,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dstID )
{
   auto srcData = cell.getData( srcID );
   auto dstData = cell.getData( dstID );
   srcData->swap( *dstData, level );
}

template < typename ValueType >
inline real_t evaluate( const uint_t&, const Cell&, const Point3D&, const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& )
{
   WALBERLA_ABORT( "VertexDoF macro-cell evaluate not implemented for this data type." )
}

namespace detail {

// we should find better locations for alle functions here in detail namespace

inline constexpr int clamp( int x, int min, int max )
{
   if ( x < min )
      return min;
   if ( x > max )
      return max;
   return x;
}

inline Point3D transformToLocalTet( const Point3D& tet0,
                                    const Point3D& tet1,
                                    const Point3D& tet2,
                                    const Point3D& tet3,
                                    const Point3D& globalPoint )
{
   walberla::math::Matrix3< real_t > A;

   A( 0, 0 ) = ( tet1 - tet0 )[0];
   A( 0, 1 ) = ( tet2 - tet0 )[0];
   A( 0, 2 ) = ( tet3 - tet0 )[0];

   A( 1, 0 ) = ( tet1 - tet0 )[1];
   A( 1, 1 ) = ( tet2 - tet0 )[1];
   A( 1, 2 ) = ( tet3 - tet0 )[1];

   A( 2, 0 ) = ( tet1 - tet0 )[2];
   A( 2, 1 ) = ( tet2 - tet0 )[2];
   A( 2, 2 ) = ( tet3 - tet0 )[2];

   A.invert();

   walberla::Vector3< real_t > x( globalPoint[0] - tet0[0], globalPoint[1] - tet0[1], globalPoint[2] - tet0[2] );

   auto result = A * x;
   return Point3D( result[0], result[1], result[2] );
}

inline std::array< Index, 4 > findLocalMicroCell( const uint_t& level, const Cell& cell, const Point3D& coordinates )
{
   // Assuming now that the passed coordinates are in the cell.
   // Otherwise they are clamped.

   // 1. Affine transformation to local macro-tet.

   auto xRelMacro = detail::transformToLocalTet(
       cell.getCoordinates()[0], cell.getCoordinates()[1], cell.getCoordinates()[2], cell.getCoordinates()[3], coordinates );

   // 2. Find micro-cube in macro-cell. Each micro-cube is composed of 6 cells of all 6 cell-types.

   const int    numMicroEdges = (int) levelinfo::num_microedges_per_edge( level );
   const real_t microEdgeSize = real_c( 1.0 ) / real_c( numMicroEdges );

   int planeX = (int) ( xRelMacro[0] / microEdgeSize );
   int planeY = (int) ( xRelMacro[1] / microEdgeSize );
   int planeZ = (int) ( xRelMacro[2] / microEdgeSize );

   // clip to prevent element outside macro-cell. std::clamp is C++17 ...
   planeX = detail::clamp( planeX, 0, numMicroEdges - 1 );
   planeY = detail::clamp( planeY, 0, numMicroEdges - 1 - planeX );
   planeZ = detail::clamp( planeZ, 0, numMicroEdges - 1 - planeX - planeY );

   // 3. In the interior of the macro-cell, the micro-cubes are contained entirely.
   //    On the boundary, some micro-tets in the micro-cube are outside of the macro-cell.
   //    Let's check that.

   // check if the point is located in a sub-cube, i.e. all cell types are possible
   const bool inFullCube = planeX + planeY + planeZ < numMicroEdges - 2;
   // check if cube is at boundary and lacks the WHITE_DOWN cell
   const bool inCutCube = planeX + planeY + planeZ == numMicroEdges - 2;
   // if both are false, the point is located in a white cell near the macro-edges / macro-vertices

   // 4. Now we got through all <= 6 micro-tets in the micro-cube and check if the point is contained.
   //    We can perform some shortcuts, in general, however, we perform the coordinate transformation
   //    for all cells.
   //    If no micro-tet can be chosen directly (due to floating-point math errors) we choose the "closest"
   //    micro-tet by minimal summed distance to all faces.
   //
   //    Possible optimization: the transformation to the 6 local tets can be optimized by transforming to the
   //                           local micro-cube space first, and then apply the precomputed transforms to the point.
   //                           That could save up to 5 3x3 matrix solves(!)

   auto cellType = celldof::CellType::WHITE_UP;

   if ( inFullCube || inCutCube )
   {
      std::vector< celldof::CellType > possibleCellTypes = { celldof::CellType::WHITE_UP,
                                                             celldof::CellType::BLUE_UP,
                                                             celldof::CellType::GREEN_UP,
                                                             celldof::CellType::BLUE_DOWN,
                                                             celldof::CellType::GREEN_DOWN };
      if ( inFullCube )
      {
         possibleCellTypes.push_back( celldof::CellType::WHITE_DOWN );
      }

      real_t maxDistSum = std::numeric_limits< real_t >::max();

      for ( auto ct : possibleCellTypes )
      {
         auto mci = celldof::macrocell::getMicroVerticesFromMicroCell(
             Index( int64_c( planeX ), int64_c( planeY ), int64_c( planeZ ) ), ct );
         auto mt0 = coordinateFromIndex( level, cell, mci[0] );
         auto mt1 = coordinateFromIndex( level, cell, mci[1] );
         auto mt2 = coordinateFromIndex( level, cell, mci[2] );
         auto mt3 = coordinateFromIndex( level, cell, mci[3] );

         auto xl = detail::transformToLocalTet( mt0, mt1, mt2, mt3, coordinates );
         auto s  = xl[0] + xl[1] + xl[2];

         Point4D rel( xl[0], xl[1], xl[2], s );

         real_t distSum  = 0;
         bool   contains = true;
         for ( uint_t i = 0; i < 4; i++ )
         {
            if ( rel[i] < 0 )
            {
               distSum += std::abs( rel[i] );
               contains = false;
            }
            else if ( rel[i] > 1 )
            {
               distSum += std::abs( rel[i] - 1 );
               contains = false;
            }
         }

         if ( contains )
         {
            cellType = ct;
            break;
         }

         if ( distSum < maxDistSum )
         {
            cellType   = ct;
            maxDistSum = distSum;
         }
      }
   }

   auto microCellIndices = celldof::macrocell::getMicroVerticesFromMicroCell(
       Index( int64_c( planeX ), int64_c( planeY ), int64_c( planeZ ) ), cellType );

   const uint_t numMicroVertices = levelinfo::num_microvertices_per_edge( level );
   WALBERLA_ASSERT_LESS( microCellIndices[0].x() + microCellIndices[0].y() + microCellIndices[0].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[1].x() + microCellIndices[1].y() + microCellIndices[1].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[2].x() + microCellIndices[2].y() + microCellIndices[2].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[3].x() + microCellIndices[3].y() + microCellIndices[3].z(), numMicroVertices );

   WALBERLA_DEBUG_SECTION()
   {
      auto microTet0 = coordinateFromIndex( level, cell, microCellIndices[0] );
      auto microTet1 = coordinateFromIndex( level, cell, microCellIndices[1] );
      auto microTet2 = coordinateFromIndex( level, cell, microCellIndices[2] );
      auto microTet3 = coordinateFromIndex( level, cell, microCellIndices[3] );

      auto xLocal = detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );

      auto sum = xLocal[0] + xLocal[1] + xLocal[2];

      if ( xLocal[0] < -1e-8 || xLocal[0] > 1.0 + 1e-8 || xLocal[1] < -1e-8 || xLocal[1] > 1.0 + 1e-8 || xLocal[2] < -1e-8 ||
           xLocal[2] > 1.0 + 1e-8 || sum < -1e-8 || sum > 1.0 + 1e-8 )
      {
         WALBERLA_LOG_DEVEL( "Bad cell choice:" )
         WALBERLA_LOG_DEVEL( "local " << xLocal << ", global " << coordinates
                                      << ", local sum: " << walberla::format( "%10.5e", sum ) << ", micro-cell type"
                                      << celldof::CellTypeToStr.at( cellType ) );
         for ( auto ct : celldof::allCellTypes )
         {
            auto mci = celldof::macrocell::getMicroVerticesFromMicroCell(
                Index( int64_c( planeX ), int64_c( planeY ), int64_c( planeZ ) ), ct );
            auto mt0 = coordinateFromIndex( level, cell, mci[0] );
            auto mt1 = coordinateFromIndex( level, cell, mci[1] );
            auto mt2 = coordinateFromIndex( level, cell, mci[2] );
            auto mt3 = coordinateFromIndex( level, cell, mci[3] );

            auto xl = detail::transformToLocalTet( mt0, mt1, mt2, mt3, coordinates );
            auto s  = xl[0] + xl[1] + xl[2];
            WALBERLA_LOG_DEVEL( "Cell type: " << celldof::CellTypeToStr.at( ct ) );
            WALBERLA_LOG_DEVEL( "x local:   " << xl << ", local sum: " << walberla::format( "%10.5e", s ) );
         }

         WALBERLA_LOG_DEVEL( "Breakpoint here ..." );
      }
   }

   return microCellIndices;
}

} // namespace detail

template <>
inline real_t evaluate( const uint_t&                                            level,
                        const Cell&                                              cell,
                        const Point3D&                                           coordinates,
                        const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dataID )
{
   auto microCellIndices = detail::findLocalMicroCell( level, cell, coordinates );

   // 3. perform interpolation

   auto microTet0 = coordinateFromIndex( level, cell, microCellIndices[0] );
   auto microTet1 = coordinateFromIndex( level, cell, microCellIndices[1] );
   auto microTet2 = coordinateFromIndex( level, cell, microCellIndices[2] );
   auto microTet3 = coordinateFromIndex( level, cell, microCellIndices[3] );

   auto cellData = cell.getData( dataID )->getPointer( level );

   auto valueTet0 =
       cellData[vertexdof::macrocell::index( level, microCellIndices[0].x(), microCellIndices[0].y(), microCellIndices[0].z() )];
   auto valueTet1 =
       cellData[vertexdof::macrocell::index( level, microCellIndices[1].x(), microCellIndices[1].y(), microCellIndices[1].z() )];
   auto valueTet2 =
       cellData[vertexdof::macrocell::index( level, microCellIndices[2].x(), microCellIndices[2].y(), microCellIndices[2].z() )];
   auto valueTet3 =
       cellData[vertexdof::macrocell::index( level, microCellIndices[3].x(), microCellIndices[3].y(), microCellIndices[3].z() )];

   auto xLocal = detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );

   auto value = valueTet0 * ( real_c( 1.0 ) - xLocal[0] - xLocal[1] - xLocal[2] ) + valueTet1 * xLocal[0] +
                valueTet2 * xLocal[1] + valueTet3 * xLocal[2];

   return value;
}

template < typename ValueType >
inline void assign( const uint_t&                                                              level,
                    const Cell&                                                                cell,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&                dstId )
{
   ValueType* dst = cell.getData( dstId )->getPointer( level );

   std::vector< ValueType* > srcPtr;
   for ( const auto& src : srcIds )
   {
      srcPtr.push_back( cell.getData( src )->getPointer( level ) );
   }

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

      ValueType tmp = scalars[0] * srcPtr[0][idx];

      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp += scalars[k] * srcPtr[k][idx];
      }
      dst[idx] = tmp;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                               level,
                 const Cell&                                                 cell,
                 const ValueType&                                            scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dstId )
{
   ValueType* dst = cell.getData( dstId )->getPointer( level );

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
      dst[idx] += scalar;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                                              level,
                 const Cell&                                                                cell,
                 const std::vector< ValueType >&                                            scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&                dstId )
{
   ValueType* dst = cell.getData( dstId )->getPointer( level );

   std::vector< ValueType* > srcPtr;
   for ( const auto& src : srcIds )
   {
      srcPtr.push_back( cell.getData( src )->getPointer( level ) );
   }

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

      ValueType tmp = scalars[0] * srcPtr[0][idx];

      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp += scalars[k] * srcPtr[k][idx];
      }
      dst[idx] += tmp;
   }
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                              level,
                             const Cell&                                                                cell,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&                dstId )
{
   ValueType* dst = cell.getData( dstId )->getPointer( level );

   std::vector< ValueType* > srcPtr;
   for ( const auto& src : srcIds )
   {
      srcPtr.push_back( cell.getData( src )->getPointer( level ) );
   }

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

      ValueType tmp = srcPtr[0][idx];

      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp *= srcPtr[k][idx];
      }
      dst[idx] = tmp;
   }
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                               level,
                      const Cell&                                                 cell,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& lhsId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& rhsId,
                      const uint_t&                                               offset = 1 )
{
   auto sp = ValueType( 0 );

   const ValueType* lhsPtr = cell.getData( lhsId )->getPointer( level );
   const ValueType* rhsPtr = cell.getData( rhsId )->getPointer( level );

   for ( const auto& it : vertexdof::macrocell::Iterator( level, offset ) )
   {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
      sp += lhsPtr[idx] * rhsPtr[idx];
   }

   return sp;
}

template < typename ValueType >
inline ValueType sum( const uint_t&                                               level,
                      const Cell&                                                 cell,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dataID,
                      const bool&                                                 absolute )
{
   auto sum = ValueType( 0 );

   const ValueType* data = cell.getData( dataID )->getPointer( level );

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
      if ( absolute )
      {
         sum += std::abs( data[idx] );
      }
      else
      {
         sum += data[idx];
      }
   }

   return sum;
}

template < typename ValueType >
inline void apply( const uint_t&                                                   level,
                   Cell&                                                           cell,
                   const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Cell >& operatorId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&     srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&     dstId,
                   const UpdateType                                                update )
{
   typedef stencilDirection sd;

   auto             operatorData = cell.getData( operatorId )->getData( level );
   const ValueType* src          = cell.getData( srcId )->getPointer( level );
   ValueType*       dst          = cell.getData( dstId )->getPointer( level );

   ValueType tmp;

   if ( update == Replace )
   {
      for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
      {
         const idx_t x = it.x();
         const idx_t y = it.y();
         const idx_t z = it.z();

         const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

         tmp = operatorData.at( { 0, 0, 0 } ) * src[centerIdx];

         for ( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
         {
            const uint_t idx = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
            WALBERLA_ASSERT_GREATER( operatorData.count( vertexdof::logicalIndexOffsetFromVertex( neighbor ) ), 0 );
            tmp += operatorData.at( vertexdof::logicalIndexOffsetFromVertex( neighbor ) ) * src[idx];
         }

         dst[centerIdx] = tmp;
      }
   }
   else
   {
      for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
      {
         const idx_t x = it.x();
         const idx_t y = it.y();
         const idx_t z = it.z();

         const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

         tmp = operatorData.at( { 0, 0, 0 } ) * src[centerIdx];

         for ( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
         {
            const uint_t idx = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
            tmp += operatorData.at( logicalIndexOffsetFromVertex( neighbor ) ) * src[idx];
         }

         dst[centerIdx] += tmp;
      }
   }
}

template < typename ValueType >
inline void smooth_gs( const uint_t&                                                   level,
                       Cell&                                                           cell,
                       const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Cell >& operatorId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&     dstId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&     rhsId )
{
   typedef stencilDirection sd;

   auto             operatorData = cell.getData( operatorId )->getData( level );
   const ValueType* rhs          = cell.getData( rhsId )->getPointer( level );
   ValueType*       dst          = cell.getData( dstId )->getPointer( level );

   ValueType tmp;

   const auto inverseCenterWeight = 1.0 / operatorData[{ 0, 0, 0 }];

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

      tmp = rhs[centerIdx];

      for ( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         const uint_t idx = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
         tmp -= operatorData[logicalIndexOffsetFromVertex( neighbor )] * dst[idx];
      }

      dst[centerIdx] = tmp * inverseCenterWeight;
   }
}

template < typename ValueType >
inline void smooth_sor( const uint_t&                                                   level,
                        Cell&                                                           cell,
                        const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Cell >& operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&     dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&     rhsId,
                        ValueType                                                       relax )
{
   typedef stencilDirection sd;

   auto             operatorData = cell.getData( operatorId )->getData( level );
   const ValueType* rhs          = cell.getData( rhsId )->getPointer( level );
   ValueType*       dst          = cell.getData( dstId )->getPointer( level );

   ValueType tmp;

   const auto inverseCenterWeight = 1.0 / operatorData[{ 0, 0, 0 }];

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

      tmp = rhs[centerIdx];

      for ( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         const uint_t idx = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
         tmp -= operatorData[logicalIndexOffsetFromVertex( neighbor )] * dst[idx];
      }

      dst[centerIdx] = ( 1.0 - relax ) * dst[centerIdx] + tmp * relax * inverseCenterWeight;
   }
}

template < typename ValueType >
inline void enumerate( const uint_t&                                               Level,
                       Cell&                                                       cell,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dstId,
                       ValueType&                                                  num )
{
   ValueType* dstPtr = cell.getData( dstId )->getPointer( Level );

   for ( const auto& it : vertexdof::macrocell::Iterator( Level, 1 ) )
   {
      const uint_t idx = vertexdof::macrocell::index( Level, it.x(), it.y(), it.z() );
      dstPtr[idx]      = num;
      num++;
   }
}

template < typename ValueType >
inline ValueType getMaxValue( const uint_t& level, Cell& cell, const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId )
{
   auto src      = cell.getData( srcId )->getPointer( level );
   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      localMax = std::max(
          localMax, src[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )] );
   }

   return localMax;
}

template < typename ValueType >
inline ValueType getMinValue( const uint_t& level, Cell& cell, const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId )
{
   auto src      = cell.getData( srcId )->getPointer( level );
   auto localMin = std::numeric_limits< ValueType >::max();

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      localMin = std::min(
          localMin, src[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )] );
   }

   return localMin;
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Cell& cell, const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId )
{
   auto src      = cell.getData( srcId )->getPointer( level );
   auto localMax = ValueType( 0.0 );

   for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
   {
      localMax = std::max(
          localMax,
          std::abs( src[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )] ) );
   }

   return localMax;
}

inline void saveOperator( const uint_t&                                                                         Level,
                          Cell&                                                                                 cell,
                          const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >& operatorId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Cell >&                               srcId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Cell >&                               dstId,
                          const std::shared_ptr< SparseMatrixProxy >&                                           mat )
{
   auto opr_data = cell.getData( operatorId )->getData( Level );
   auto src      = cell.getData( srcId )->getPointer( Level );
   auto dst      = cell.getData( dstId )->getPointer( Level );

   for ( const auto& it : vertexdof::macrocell::Iterator( Level, 1 ) )
   {
      idx_t srcInt = src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];
      idx_t dstInt = dst[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];

      mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[{ 0, 0, 0 }] );

      for ( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
         srcInt = src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), neighbor )];

         mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[vertexdof::logicalIndexOffsetFromVertex( neighbor )] );
      }
   }
}

inline void saveIdentityOperator( const uint_t&                                           Level,
                                  Cell&                                                   cell,
                                  const PrimitiveDataID< FunctionMemory< idx_t >, Cell >& dstId,
                                  const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   auto dst = cell.getData( dstId )->getPointer( Level );

   for ( const auto& it : vertexdof::macrocell::Iterator( Level, 1 ) )
   {
      idx_t dstInt = dst[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];

      mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
   }
}

template < typename ValueType >
inline void createVectorFromFunction( const uint_t&                                               Level,
                                      Cell&                                                       cell,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Cell >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   auto src       = cell.getData( srcId )->getPointer( Level );
   auto numerator = cell.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : vertexdof::macrocell::Iterator( Level, 1 ) )
   {
      idx_t numeratorInt =
          numerator[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];
      vec->setValue( uint_c( numeratorInt ),
                     src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )] );
   }
}

template < typename ValueType >
inline void createFunctionFromVector( const uint_t&                                               Level,
                                      Cell&                                                       cell,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Cell >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   auto src       = cell.getData( srcId )->getPointer( Level );
   auto numerator = cell.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : vertexdof::macrocell::Iterator( Level, 1 ) )
   {
      idx_t numeratorInt =
          numerator[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];
      src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )] =
          vec->getValue( uint_c( numeratorInt ) );
   }
}

} // namespace hyteg::vertexdof::macrocell
