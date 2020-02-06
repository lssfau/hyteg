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

#include "core/debug/all.h"
#include "core/DataTypes.h"
#include "core/math/Matrix3.h"

#include "hyteg/primitives/Cell.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/FunctionMemory.hpp"
#include "hyteg/StencilMemory.hpp"
#include "hyteg/types/flags.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/Algorithms.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/LevelWiseMemory.hpp"
#include "hyteg/celldofspace/CellDoFIndexing.hpp"

namespace hyteg {
namespace vertexdof {
namespace macrocell {

using walberla::uint_t;
using walberla::real_t;
using walberla::real_c;

using indexing::Index;

inline indexing::Index getIndexInNeighboringMacroFace( const indexing::Index  & vertexDoFIndexInMacroCell,
                                                       const Cell             & cell,
                                                       const uint_t           & neighborFaceID,
                                                       const PrimitiveStorage & storage,
                                                       const uint_t           & level )
{
  const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 3, 4 >(
  { cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(0),
    cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(1),
    cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(2) } );

  const auto indexInMacroFace = indexing::basisConversion( vertexDoFIndexInMacroCell, {0, 1, 2, 3},
                                                           localVertexIDsAtCell, levelinfo::num_microvertices_per_edge( level ) );
  return indexInMacroFace;
}

inline Point3D coordinateFromIndex( const uint_t & level, const Cell & cell, const Index & index )
{
  const real_t  stepFrequency = 1.0 / levelinfo::num_microedges_per_edge( level );
  const Point3D xStep         = ( cell.getCoordinates()[1] - cell.getCoordinates()[0] ) * stepFrequency;
  const Point3D yStep         = ( cell.getCoordinates()[2] - cell.getCoordinates()[0] ) * stepFrequency;
  const Point3D zStep         = ( cell.getCoordinates()[3] - cell.getCoordinates()[0] ) * stepFrequency;
  return cell.getCoordinates()[0] + xStep * real_c( index.x() ) + yStep * real_c( index.y() ) + zStep * real_c( index.z() );
}

template< typename ValueType >
inline void interpolate( const uint_t & level,
                         const Cell & cell,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& cellMemoryId,
                         const ValueType & scalar )
{
  ValueType * cellData = cell.getData( cellMemoryId )->getPointer( level );

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
    cellData[ idx ] = scalar;
  }
}


template< typename ValueType >
inline void interpolate( const uint_t & level,
                         const Cell & cell,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& cellMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                         const std::function< ValueType( const hyteg::Point3D &, const std::vector< ValueType > & )> & expr)
{
  ValueType * cellData = cell.getData( cellMemoryId )->getPointer( level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( cell.getData(src)->getPointer( level ) );
  }

  std::vector<ValueType> srcVector( srcIds.size() );

  Point3D xBlend;

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const Point3D coordinate = coordinateFromIndex( level, cell, it );
    const uint_t  idx        = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

    for ( uint_t k = 0; k < srcPtr.size(); ++k )
    {
      srcVector[ k ] = srcPtr[ k ][ idx ];
    }
    cell.getGeometryMap()->evalF( coordinate, xBlend );
    cellData[ idx ] = expr( xBlend, srcVector );
  }
}

template< typename ValueType >
inline void swap( const uint_t & level, Cell & cell,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstID )
{
  auto srcData = cell.getData( srcID );
  auto dstData = cell.getData( dstID );
  srcData->swap( *dstData, level );
}

template< typename ValueType >
inline real_t evaluate( const uint_t & level, const Cell & cell, const Point3D & coordinates,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dataID )
{
   WALBERLA_ABORT("VertexDoF macro-cell evaluate not implemented for this data type.")
}

namespace detail {

// we should find better locations for alle functions here in detail namespace

inline constexpr int clamp( int x, int min, int max )
{
   if ( x < min ) return min;
   if ( x > max ) return max;
   return x;
}

inline Point3D transformToLocalTet( const Point3D & tet0, const Point3D & tet1, const Point3D & tet2, const Point3D & tet3, const Point3D & globalPoint )
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

   auto transform = A.invert();

   walberla::Vector3< real_t > x( globalPoint[0] - tet0[0], globalPoint[1] - tet0[1], globalPoint[2] - tet0[2] );

   auto result = transform * x;
   return Point3D( { result[0], result[1], result[2] } );
}

}

template<>
inline real_t evaluate( const uint_t & level, const Cell & cell, const Point3D & coordinates,
                        const PrimitiveDataID< FunctionMemory< real_t >, Cell > & dataID )
{
   // Assuming now that the passed coordinates are in the cell.
   // Otherwise they are clamped.

   // 1. affine transformation to local macro-tet

   auto xRelMacro = detail::transformToLocalTet(
       cell.getCoordinates()[0], cell.getCoordinates()[1], cell.getCoordinates()[2], cell.getCoordinates()[3], coordinates );

   // 2. Clip along the six relevant half planes that intersect the macro-tetrahedron.
   //    Those planes have the following normals:
   //       (1, 0, 0) -> x-plane
   //       (0, 1, 0) -> y-plane
   //       (0, 0, 1) -> z-plane
   //       (1, 1, 1) -> xyz-plane
   //       (1, 1, 0) -> xy-plane
   //       (0, 1, 1) -> yz-plane
   //    All coordinates are clipped to the next lower plane (relative to origin at vertex 0, or (0, 0, 0)).
   //    Also they are clipped to the "interior" of the tetrahedron.
   //
   //    Depending on the index of the plane, the exact micro-tetrahedron can be defined.
   //    First we determine the "micro-cube" depending on the x-, y- and z-planes.
   //    Then the other plane-indices are evaluated.
   //    E.g.:
   //       x-plane = 0, y-plane = 0, z-plane = 0 (=> micro-cube at (0, 0, 0) - (1, 1, 1))
   //       xyz-plane == 0 -> white-up cell
   //       xyz-plane == 2 -> white-down cell
   //       xyz-plane == 1 ->
   //          xy-plane | yz-plane |     cell
   //          ---------+----------+---------
   //              0    |     0    | green-up
   //              0    |     1    | blue-down
   //              1    |     0    | blue-up
   //              1    |     1    | green-down
   //
   //    In general we can calculate the final cell directly using the following rules:
   //       xyz-plane - (x-plane + y-plane + z-plane) == 0 -> white-up cell
   //       xyz-plane - (x-plane + y-plane + z-plane) == 2 -> white-down cell
   //       xyz-plane - (x-plane + y-plane + z-plane) == 1 ->
   //
   //          xy-plane-tilde = xy-plane - (x-plane + y-plane)
   //          yz-plane-tilde = yz-plane - (y-plane + z-plane)
   //
   //          xy-plane-tilde | yz-plane-tilde |     cell
   //          ---------------+----------------+---------
   //                    0    |     0          | green-up
   //                    0    |     1          | blue-down
   //                    1    |     0          | blue-up
   //                    1    |     1          | green-down
   //    The final cell has the local cell-coordinates of the cube (x-plane, y-plane, z-plane).

   const int    numMicroEdges = (int) levelinfo::num_microedges_per_edge( level );
   const real_t microEdgeSize = 1.0 / real_c( numMicroEdges );

   const real_t planeWidthXYZ = microEdgeSize * std::sqrt( 3.0 ) / 3.0;
   const real_t planeWidthXY  = microEdgeSize * std::sqrt( 2.0 ) / 2.0;
   const real_t planeWidthYZ  = microEdgeSize * std::sqrt( 2.0 ) / 2.0;

   int planeX = (int) ( xRelMacro[0] / microEdgeSize );
   int planeY = (int) ( xRelMacro[1] / microEdgeSize );
   int planeZ = (int) ( xRelMacro[2] / microEdgeSize );

   const real_t lengthXYZ = xRelMacro.norm();
   const real_t lengthXY  = std::sqrt( xRelMacro[0] * xRelMacro[0] + xRelMacro[1] * xRelMacro[1] );
   const real_t lengthYZ  = std::sqrt( xRelMacro[1] * xRelMacro[1] + xRelMacro[2] * xRelMacro[2] );

   int planeXYZ = (int) ( lengthXYZ / planeWidthXYZ );
   int planeXY  = (int) ( lengthXY / planeWidthXY );
   int planeYZ  = (int) ( lengthYZ / planeWidthYZ );

   // clip to prevent element outside macro-cell. std::clamp is C++17 ...

   planeX = detail::clamp( planeX, 0, numMicroEdges - 1 );
   planeY = detail::clamp( planeY, 0, numMicroEdges - 1 );
   planeZ = detail::clamp( planeZ, 0, numMicroEdges - 1 );

   planeXYZ = detail::clamp( planeXYZ, planeX + planeY + planeZ, planeX + planeY + planeZ + 2 );
   planeXY  = detail::clamp( planeXY, planeX + planeY, planeX + planeY + 1 );
   planeYZ  = detail::clamp( planeYZ, planeY + planeZ, planeY + planeZ + 1 );

   celldof::CellType cellType;

   // check if the point is located in a sub-cube, i.e. all cell types are possible
   const bool inFullCube = planeX + planeY + planeZ < numMicroEdges - 2;
   // check if cube is at boundary and lacks the WHITE_DOWN cell
   const bool inCutCube  = planeX + planeY + planeZ == numMicroEdges - 2;
   // if both are false, the point is located in a white cell near the macro-edges / macro-vertices

   const int whiteCellDecider = planeXYZ - ( planeX + planeY + planeZ );

   if ( (!inFullCube && !inCutCube) || whiteCellDecider == 0 )
   {
      cellType = celldof::CellType::WHITE_UP;
   }
   else if ( inFullCube && whiteCellDecider == 2 )
   {
      cellType = celldof::CellType::WHITE_DOWN;
   }
   else
   {
      if ( planeXY - ( planeX + planeY ) == 0 )
      {
         if ( planeYZ - ( planeY + planeZ ) == 0 )
         {
            cellType = celldof::CellType::GREEN_UP;
         }
         else
         {
            cellType = celldof::CellType::BLUE_DOWN;
         }
      }
      else
      {
         if ( planeYZ - ( planeY + planeZ ) == 0 )
         {
            cellType = celldof::CellType::BLUE_UP;
         }
         else
         {
            cellType = celldof::CellType::GREEN_DOWN;
         }
      }
   }

   const auto microCellIndices = celldof::macrocell::getMicroVerticesFromMicroCell( Index( uint_c(planeX), uint_c(planeY), uint_c(planeZ) ), cellType );

   const uint_t numMicroVertices = levelinfo::num_microvertices_per_edge( level );
   WALBERLA_ASSERT_LESS( microCellIndices[0].x() + microCellIndices[0].y() + microCellIndices[0].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[1].x() + microCellIndices[1].y() + microCellIndices[1].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[2].x() + microCellIndices[2].y() + microCellIndices[2].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[3].x() + microCellIndices[3].y() + microCellIndices[3].z(), numMicroVertices );

   // ^^^ in extra function - code can be re-used for any element order ^^^

   // 3. perform interpolation

   auto microTet0 = coordinateFromIndex( level, cell, microCellIndices[0] );
   auto microTet1 = coordinateFromIndex( level, cell, microCellIndices[1] );
   auto microTet2 = coordinateFromIndex( level, cell, microCellIndices[2] );
   auto microTet3 = coordinateFromIndex( level, cell, microCellIndices[3] );

   auto cellData = cell.getData( dataID )->getPointer( level );

   auto valueTet0 = cellData[ vertexdof::macrocell::index( level, microCellIndices[0].x(), microCellIndices[0].y(), microCellIndices[0].z()) ];
   auto valueTet1 = cellData[ vertexdof::macrocell::index( level, microCellIndices[1].x(), microCellIndices[1].y(), microCellIndices[1].z()) ];
   auto valueTet2 = cellData[ vertexdof::macrocell::index( level, microCellIndices[2].x(), microCellIndices[2].y(), microCellIndices[2].z()) ];
   auto valueTet3 = cellData[ vertexdof::macrocell::index( level, microCellIndices[3].x(), microCellIndices[3].y(), microCellIndices[3].z()) ];

   auto xLocal = detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );

   auto value = valueTet0 * (1.0 - xLocal[0] - xLocal[1] - xLocal[2]) + valueTet1 * xLocal[0] + valueTet2 * xLocal[1] + valueTet3 * xLocal[2];

   return value;

}

template< typename ValueType >
inline void assign( const uint_t & level,
                    const Cell & cell,
                    const std::vector< ValueType > & scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  ValueType * dst = cell.getData( dstId )->getPointer( level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( cell.getData( src )->getPointer( level ));
  }

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

    ValueType tmp = scalars[0] * srcPtr[0][ idx ];

    for ( uint_t k = 1; k < srcIds.size(); ++k )
    {
      tmp += scalars[k] * srcPtr[k][ idx ];
    }
    dst[ idx ] = tmp;
  }
}


template< typename ValueType >
inline void add( const uint_t & level,
                 const Cell & cell,
                 const ValueType & scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  ValueType * dst = cell.getData( dstId )->getPointer( level );

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
    dst[ idx ] += scalar;
  }
}


template< typename ValueType >
inline void add( const uint_t & level,
                 const Cell & cell,
                 const std::vector< ValueType > & scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  ValueType * dst = cell.getData( dstId )->getPointer( level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( cell.getData( src )->getPointer( level ));
  }

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

    ValueType tmp = scalars[0] * srcPtr[0][ idx ];

    for ( uint_t k = 1; k < srcIds.size(); ++k )
    {
      tmp += scalars[k] * srcPtr[k][ idx ];
    }
    dst[ idx ] += tmp;
  }
}


template< typename ValueType >
inline void multElementwise( const uint_t & level,
                               const Cell & cell,
                               const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                               const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  ValueType * dst = cell.getData( dstId )->getPointer( level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( cell.getData( src )->getPointer( level ));
  }

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );

    ValueType tmp = srcPtr[0][ idx ];

    for ( uint_t k = 1; k < srcIds.size(); ++k )
    {
      tmp *= srcPtr[k][ idx ];
    }
    dst[ idx ] = tmp;
  }

}

template< typename ValueType >
inline ValueType dot( const uint_t & level,
                   const Cell & cell,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & lhsId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & rhsId)
{
   auto sp = ValueType( 0 );

  const ValueType * lhsPtr = cell.getData( lhsId )->getPointer( level );
  const ValueType * rhsPtr = cell.getData( rhsId )->getPointer( level );

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
    sp += lhsPtr[ idx ] * rhsPtr[ idx ];
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

template< typename ValueType >
inline void apply( const uint_t & level,
                   Cell & cell,
                   const PrimitiveDataID< LevelWiseMemory< StencilMap_T >,  Cell > & operatorId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId,
                   const UpdateType update )
{
  typedef stencilDirection sd;

  auto operatorData     = cell.getData( operatorId )->getData( level );
  const ValueType * src = cell.getData( srcId )->getPointer( level );
        ValueType * dst = cell.getData( dstId )->getPointer( level );

  ValueType tmp;

  if( update == Replace )
  {
    for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
    {
      const uint_t x = it.x();
      const uint_t y = it.y();
      const uint_t z = it.z();

      const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

      tmp = operatorData.at({0, 0, 0}) * src[ centerIdx ];

      for ( const auto & neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
        const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
        WALBERLA_ASSERT_GREATER( operatorData.count( vertexdof::logicalIndexOffsetFromVertex( neighbor ) ), 0 );
        tmp += operatorData.at( vertexdof::logicalIndexOffsetFromVertex( neighbor ) ) * src[ idx ];
      }

      dst[ centerIdx ] = tmp;
    }
  }
  else
  {
    for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
    {
      const uint_t x = it.x();
      const uint_t y = it.y();
      const uint_t z = it.z();

      const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

      tmp = operatorData.at( {0, 0, 0} ) * src[ centerIdx ];

      for ( const auto & neighbor : vertexdof::macrocell::neighborsWithoutCenter )
      {
        const uint_t stencilIdx = vertexdof::stencilIndexFromVertex( neighbor );
        const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
        tmp += operatorData.at( logicalIndexOffsetFromVertex( neighbor ) ) * src[ idx ];
      }

      dst[ centerIdx ] += tmp;
    }
  }
}


template< typename ValueType >
inline void smooth_gs( const uint_t & level,
                       Cell & cell,
                       const PrimitiveDataID<LevelWiseMemory< StencilMap_T >,  Cell > & operatorId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & rhsId )
{
  typedef stencilDirection sd;

  auto operatorData = cell.getData( operatorId )->getData( level );
  const ValueType * rhs          = cell.getData( rhsId )->getPointer( level );
        ValueType * dst          = cell.getData( dstId )->getPointer( level );

  ValueType tmp;

  const auto inverseCenterWeight = 1.0 / operatorData[ { 0, 0, 0 } ];

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t x = it.x();
    const uint_t y = it.y();
    const uint_t z = it.z();

    const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

    tmp = rhs[ centerIdx ];

    for ( const auto & neighbor : vertexdof::macrocell::neighborsWithoutCenter )
    {
      const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
      tmp -= operatorData[ logicalIndexOffsetFromVertex( neighbor ) ] * dst[ idx ];
    }

    dst[ centerIdx ] = tmp * inverseCenterWeight;
  }
}

template< typename ValueType >
inline void smooth_sor( const uint_t & level,
                       Cell & cell,
                       const PrimitiveDataID< LevelWiseMemory< StencilMap_T >,  Cell > & operatorId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & rhsId,
                       ValueType                                                    relax )
{
  typedef stencilDirection sd;

  auto operatorData = cell.getData( operatorId )->getData( level );
  const ValueType * rhs          = cell.getData( rhsId )->getPointer( level );
  ValueType * dst          = cell.getData( dstId )->getPointer( level );

  ValueType tmp;

  const auto inverseCenterWeight = 1.0 / operatorData[ { 0, 0, 0 } ];

  for ( const auto & it : vertexdof::macrocell::Iterator( level, 1 ) )
  {
    const uint_t x = it.x();
    const uint_t y = it.y();
    const uint_t z = it.z();

    const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, x, y, z, sd::VERTEX_C );

    tmp = rhs[ centerIdx ];

    for ( const auto & neighbor : vertexdof::macrocell::neighborsWithoutCenter )
    {
      const uint_t idx        = vertexdof::macrocell::indexFromVertex( level, x, y, z, neighbor );
      tmp -= operatorData[ logicalIndexOffsetFromVertex( neighbor ) ] * dst[ idx ];
    }

    dst[ centerIdx ] = ( 1.0 - relax ) * dst[ centerIdx ] + tmp * relax * inverseCenterWeight;
  }
}



template< typename ValueType >
inline void enumerate(const uint_t & Level, Cell & cell, const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &dstId, ValueType& num) {

  ValueType* dstPtr = cell.getData(dstId)->getPointer( Level );

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
    const uint_t idx = vertexdof::macrocell::index( Level, it.x(), it.y(), it.z() );
    dstPtr[idx] = num;
    num++;
  }
}


template< typename ValueType >
inline ValueType getMaxValue( const uint_t & level, Cell &cell, const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId ) {

  auto src = cell.getData( srcId )->getPointer( level );
  auto localMax = -std::numeric_limits< ValueType >::max();

  for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) ) {
    localMax = std::max( localMax, src[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ] );
  }

  return localMax;
}


template< typename ValueType >
inline ValueType getMinValue( const uint_t & level, Cell &cell, const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId ) {

  auto src = cell.getData( srcId )->getPointer( level );
  auto localMin = std::numeric_limits< ValueType >::max();

  for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) ) {
    localMin = std::min( localMin, src[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ] );
  }

  return localMin;
}


template< typename ValueType >
inline ValueType getMaxMagnitude( const uint_t & level, Cell &cell, const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId ) {

  auto src = cell.getData( srcId )->getPointer( level );
  auto localMax = ValueType(0.0);

  for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) ) {
    localMax = std::max( localMax, std::abs( src[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ] ));
  }

  return localMax;
}


#ifdef HYTEG_BUILD_WITH_PETSC

inline void saveOperator(const uint_t & Level, Cell & cell, const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T > , Cell>& operatorId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Cell> &srcId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Cell> &dstId, Mat& mat)
{

  auto opr_data = cell.getData(operatorId)->getData( Level );
  auto src = cell.getData(srcId)->getPointer( Level );
  auto dst = cell.getData(dstId)->getPointer( Level );

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
      PetscInt srcInt = src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C)];
      PetscInt dstInt = dst[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C)];

      MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[ { 0, 0, 0 } ] ,ADD_VALUES);

      for ( const auto & neighbor : vertexdof::macrocell::neighborsWithoutCenter ) {
        srcInt = src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), neighbor)];

        MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[vertexdof::logicalIndexOffsetFromVertex(neighbor)] ,ADD_VALUES);
      }
  }
}



template< typename ValueType >
inline void createVectorFromFunction(const uint_t & Level, Cell & cell,
                              const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Cell> &numeratorId,
                              Vec& vec)
{

  auto src = cell.getData(srcId)->getPointer( Level );
  auto numerator = cell.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
      PetscInt numeratorInt = numerator[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C)];
      VecSetValues(vec,1,&numeratorInt,&src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C)],INSERT_VALUES);
  }
}


template< typename ValueType >
inline void createFunctionFromVector(const uint_t & Level, Cell & cell,
                                         const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId,
                                         const PrimitiveDataID<FunctionMemory< PetscInt >, Cell> &numeratorId,
                                         Vec& vec)
{
  auto src = cell.getData(srcId)->getPointer( Level );
  auto numerator = cell.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : vertexdof::macrocell::Iterator( Level, 1 ) )
  {
      PetscInt numeratorInt = numerator[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C)];
      VecGetValues(vec,1,&numeratorInt,&src[vertexdof::macrocell::indexFromVertex( Level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C)]);
  }
}

#endif


} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg
