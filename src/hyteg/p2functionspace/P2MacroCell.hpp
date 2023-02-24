/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/primitives/Cell.hpp"

namespace hyteg {
namespace P2 {
namespace macrocell {

using walberla::real_t;
using walberla::uint_t;
using walberla::int_c;
using indexing::Index;
using indexing::Index;

template < typename ValueType >
inline void evaluate( const uint_t&,
                      const Cell&,
                      const Point3D&,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >&,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >&,
                      std::vector< ValueType >& )
{
   WALBERLA_ABORT( "P2 3D evaluate not implemented for this data type-" )
}



template < typename ValueType >
inline ValueType evaluate( const uint_t&                                               level,
                           const Cell&                                                 cell,
                           const Point3D&                                              coordinates,
                           const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& vertexDoFDataID,
                           const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& edgeDoFDataID )
{
   std::vector< ValueType > results( 1 );
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > vertexDoFDataIDs( { vertexDoFDataID } );
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > edgeDoFDataIDs( { edgeDoFDataID } );
   evaluate< ValueType >( level, cell, coordinates, vertexDoFDataIDs, edgeDoFDataIDs, results );
   return results[0];
}


template <>
inline void evaluate( const uint_t&                                            level,
                      const Cell&                                              cell,
                      const Point3D&                                           coordinates,
                      const std::vector< PrimitiveDataID< FunctionMemory< real_t >, Cell > >& vertexDoFDataIDs,
                      const std::vector< PrimitiveDataID< FunctionMemory< real_t >, Cell > >& edgeDoFDataIDs,
                      std::vector< real_t > & results )
{
   auto microCellIndices = vertexdof::macrocell::detail::findLocalMicroCell( level, cell, coordinates );

   // 3. perform interpolation

   auto microTet0 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[0] );
   auto microTet1 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[1] );
   auto microTet2 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[2] );
   auto microTet3 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[3] );

   Index vertexIndex0 =
       Index( int_c( microCellIndices[0].x() ), int_c( microCellIndices[0].y() ), int_c( microCellIndices[0].z() ) );
   Index vertexIndex1 =
       Index( int_c( microCellIndices[1].x() ), int_c( microCellIndices[1].y() ), int_c( microCellIndices[1].z() ) );
   Index vertexIndex2 =
       Index( int_c( microCellIndices[2].x() ), int_c( microCellIndices[2].y() ), int_c( microCellIndices[2].z() ) );
   Index vertexIndex3 =
       Index( int_c( microCellIndices[3].x() ), int_c( microCellIndices[3].y() ), int_c( microCellIndices[3].z() ) );

   auto edgeIndex0 = edgedof::calcEdgeDoFIndex( vertexIndex0, vertexIndex1 );
   auto edgeIndex1 = edgedof::calcEdgeDoFIndex( vertexIndex0, vertexIndex2 );
   auto edgeIndex2 = edgedof::calcEdgeDoFIndex( vertexIndex1, vertexIndex2 );
   auto edgeIndex3 = edgedof::calcEdgeDoFIndex( vertexIndex0, vertexIndex3 );
   auto edgeIndex4 = edgedof::calcEdgeDoFIndex( vertexIndex1, vertexIndex3 );
   auto edgeIndex5 = edgedof::calcEdgeDoFIndex( vertexIndex2, vertexIndex3 );

   auto edgeOrientation0 = edgedof::calcEdgeDoFOrientation( vertexIndex0, vertexIndex1 );
   auto edgeOrientation1 = edgedof::calcEdgeDoFOrientation( vertexIndex0, vertexIndex2 );
   auto edgeOrientation2 = edgedof::calcEdgeDoFOrientation( vertexIndex1, vertexIndex2 );
   auto edgeOrientation3 = edgedof::calcEdgeDoFOrientation( vertexIndex0, vertexIndex3 );
   auto edgeOrientation4 = edgedof::calcEdgeDoFOrientation( vertexIndex1, vertexIndex3 );
   auto edgeOrientation5 = edgedof::calcEdgeDoFOrientation( vertexIndex2, vertexIndex3 );

   auto xLocal = vertexdof::macrocell::detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );

   auto xi_1 = xLocal[0];
   auto xi_2 = xLocal[1];
   auto xi_3 = xLocal[2];

   const auto tetV0ArrayIdx =
       vertexdof::macrocell::index( level, microCellIndices[0].x(), microCellIndices[0].y(), microCellIndices[0].z() );
   const auto tetV1ArrayIdx =
       vertexdof::macrocell::index( level, microCellIndices[1].x(), microCellIndices[1].y(), microCellIndices[1].z() );
   const auto tetV2ArrayIdx =
       vertexdof::macrocell::index( level, microCellIndices[2].x(), microCellIndices[2].y(), microCellIndices[2].z() );
   const auto tetV3ArrayIdx =
       vertexdof::macrocell::index( level, microCellIndices[3].x(), microCellIndices[3].y(), microCellIndices[3].z() );

   const auto tetE0ArrayIdx = edgedof::macrocell::index( level, uint_c(edgeIndex0.x()), uint_c(edgeIndex0.y()), uint_c(edgeIndex0.z()), edgeOrientation0 );
   const auto tetE1ArrayIdx = edgedof::macrocell::index( level, uint_c(edgeIndex1.x()), uint_c(edgeIndex1.y()), uint_c(edgeIndex1.z()), edgeOrientation1 );
   const auto tetE2ArrayIdx = edgedof::macrocell::index( level, uint_c(edgeIndex2.x()), uint_c(edgeIndex2.y()), uint_c(edgeIndex2.z()), edgeOrientation2 );
   const auto tetE3ArrayIdx = edgedof::macrocell::index( level, uint_c(edgeIndex3.x()), uint_c(edgeIndex3.y()), uint_c(edgeIndex3.z()), edgeOrientation3 );
   const auto tetE4ArrayIdx = edgedof::macrocell::index( level, uint_c(edgeIndex4.x()), uint_c(edgeIndex4.y()), uint_c(edgeIndex4.z()), edgeOrientation4 );
   const auto tetE5ArrayIdx = edgedof::macrocell::index( level, uint_c(edgeIndex5.x()), uint_c(edgeIndex5.y()), uint_c(edgeIndex5.z()), edgeOrientation5 );


   for ( uint_t i = 0; i < results.size(); i++ )
   {
      auto vertexDoFDataID = vertexDoFDataIDs[i];
      auto edgeDoFDataID = edgeDoFDataIDs[i];

      auto vertexdofData = cell.getData( vertexDoFDataID )->getPointer( level );
      auto edgedofData   = cell.getData( edgeDoFDataID )->getPointer( level );

      auto valueTetV0 = vertexdofData[tetV0ArrayIdx];
      auto valueTetV1 = vertexdofData[tetV1ArrayIdx];
      auto valueTetV2 = vertexdofData[tetV2ArrayIdx];
      auto valueTetV3 = vertexdofData[tetV3ArrayIdx];

      auto valueTetE0 = edgedofData[tetE0ArrayIdx];
      auto valueTetE1 = edgedofData[tetE1ArrayIdx];
      auto valueTetE2 = edgedofData[tetE2ArrayIdx];
      auto valueTetE3 = edgedofData[tetE3ArrayIdx];
      auto valueTetE4 = edgedofData[tetE4ArrayIdx];
      auto valueTetE5 = edgedofData[tetE5ArrayIdx];

      // basis functions P2 N(xi_1, xi_2, xi_3):
      //   at [0. 0. 0.]: 2.0*xi_1**2 + 4.0*xi_1*xi_2 + 4.0*xi_1*xi_3 - 3.0*xi_1 + 2.0*xi_2**2 + 4.0*xi_2*xi_3 - 3.0*xi_2 + 2.0*xi_3**2 - 3.0*xi_3 + 1.0
      //   at [1. 0. 0.]: 2.0*xi_1**2 - 1.0*xi_1
      //   at [0. 1. 0.]: 2.0*xi_2**2 - 1.0*xi_2
      //   at [0. 0. 1.]: 2.0*xi_3**2 - 1.0*xi_3
      //   at [0.5 0.  0. ]: -4.0*xi_1**2 - 4.0*xi_1*xi_2 - 4.0*xi_1*xi_3 + 4.0*xi_1
      //   at [0.5 0.5 0. ]: 4.0*xi_1*xi_2
      //   at [0.  0.5 0. ]: -4.0*xi_1*xi_2 - 4.0*xi_2**2 - 4.0*xi_2*xi_3 + 4.0*xi_2
      //   at [0.  0.  0.5]: -4.0*xi_1*xi_3 - 4.0*xi_2*xi_3 - 4.0*xi_3**2 + 4.0*xi_3
      //   at [0.5 0.  0.5]: 4.0*xi_1*xi_3
      //   at [0.  0.5 0.5]: 4.0*xi_2*xi_3

      auto scaleV0 = valueTetV0 * ( 2.0 * xi_1 * xi_1 + 4.0 * xi_1 * xi_2 + 4.0 * xi_1 * xi_3 - 3.0 * xi_1 + 2.0 * xi_2 * xi_2 +
                                    4.0 * xi_2 * xi_3 - 3.0 * xi_2 + 2.0 * xi_3 * xi_3 - 3.0 * xi_3 + 1.0 );
      auto scaleV1 = valueTetV1 * ( 2.0 * xi_1 * xi_1 - 1.0 * xi_1 );
      auto scaleV2 = valueTetV2 * ( 2.0 * xi_2 * xi_2 - 1.0 * xi_2 );
      auto scaleV3 = valueTetV3 * ( 2.0 * xi_3 * xi_3 - 1.0 * xi_3 );

      auto scaleE0 = valueTetE0 * ( -4.0 * xi_1 * xi_1 - 4.0 * xi_1 * xi_2 - 4.0 * xi_1 * xi_3 + 4.0 * xi_1 );
      auto scaleE1 = valueTetE1 * ( -4.0 * xi_1 * xi_2 - 4.0 * xi_2 * xi_2 - 4.0 * xi_2 * xi_3 + 4.0 * xi_2 );
      auto scaleE2 = valueTetE2 * ( 4.0 * xi_1 * xi_2 );
      auto scaleE3 = valueTetE3 * ( -4.0 * xi_1 * xi_3 - 4.0 * xi_2 * xi_3 - 4.0 * xi_3 * xi_3 + 4.0 * xi_3 );
      auto scaleE4 = valueTetE4 * ( 4.0 * xi_1 * xi_3 );
      auto scaleE5 = valueTetE5 * ( 4.0 * xi_2 * xi_3 );

      results[i] = scaleV0 + scaleV1 + scaleV2 + scaleV3 + scaleE0 + scaleE1 + scaleE2 + scaleE3 +
                   scaleE4 + scaleE5;

   }
}


void smoothSOR(
    const uint_t&                                                                                level,
    Cell&                                                                                        cell,
    const real_t&                                                                                relax,
    const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >&        vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     edgeDoFRhsId );

} // namespace macrocell
} // namespace P2
} // namespace hyteg
