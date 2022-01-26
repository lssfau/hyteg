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

#include "core/DataTypes.h"
#include "core/math/KahanSummation.h"

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/facedofspace/FaceDoFIndexing.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/types/matrix.hpp"

#ifdef DEBUG_ELEMENTWISE
#include "hyteg/format.hpp"
#endif

namespace hyteg {
namespace vertexdof {
namespace macroedge {

using indexing::Index;
using walberla::real_c;

inline indexing::Index getIndexInNeighboringMacroCell( const indexing::Index&  vertexDoFIndexInMacroEdge,
                                                       const Edge&             edge,
                                                       const uint_t&           neighborCellID,
                                                       const PrimitiveStorage& storage,
                                                       const uint_t&           level )
{
   const Cell&  neighborCell = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
   const uint_t localEdgeID  = neighborCell.getLocalEdgeID( edge.getID() );

   const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 2, 4 >(
       { neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 ),
         neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 ) } );

   const auto indexInMacroCell = indexing::basisConversion(
       vertexDoFIndexInMacroEdge, localVertexIDsAtCell, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ) );
   return indexInMacroCell;
}

inline Point3D coordinateFromIndex( const uint_t& level, const Edge& edge, const Index& index )
{
   const real_t  stepFrequency = 1.0 / levelinfo::num_microedges_per_edge( level );
   const Point3D step          = ( edge.getCoordinates()[1] - edge.getCoordinates()[0] ) * stepFrequency;
   return edge.getCoordinates()[0] + step * real_c( index.x() );
}

template < typename ValueType >
inline ValueType assembleLocal( const uint_t&                            level,
                                uint_t                                   pos,
                                const Matrix3r&                          localMatrix,
                                double*                                  src,
                                double*                                  coeff,
                                const std::array< stencilDirection, 3 >& vertices,
                                const std::array< uint_t, 3 >&           idx )
{
   ValueType meanCoeff = 1.0 / 3.0 *
                         ( coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[0] )] +
                           coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[1] )] +
                           coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[2] )] );

   ValueType tmp;
   tmp = localMatrix( idx[0], idx[0] ) * src[vertexdof::macroedge::indexFromVertex( level, pos, vertices[0] )] +
         localMatrix( idx[0], idx[1] ) * src[vertexdof::macroedge::indexFromVertex( level, pos, vertices[1] )] +
         localMatrix( idx[0], idx[2] ) * src[vertexdof::macroedge::indexFromVertex( level, pos, vertices[2] )];
   return meanCoeff * tmp;
}

template < typename ValueType >
inline void assembleLocalStencil( uint_t                                   level,
                                  uint_t                                   pos,
                                  const Matrix3r&                          localMatrix,
                                  real_t*                                  opr_data,
                                  real_t*                                  coeff,
                                  const std::array< stencilDirection, 3 >& vertices,
                                  const std::array< uint_t, 3 >&           idx )
{
   ValueType meanCoeff = 1.0 / 3.0 *
                         ( coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[0] )] +
                           coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[1] )] +
                           coeff[vertexdof::macroedge::indexFromVertex( level, pos, vertices[2] )] );

   opr_data[vertexdof::stencilIndexFromVertex( vertices[0] )] += meanCoeff * localMatrix( idx[0], idx[0] );
   opr_data[vertexdof::stencilIndexFromVertex( vertices[1] )] += meanCoeff * localMatrix( idx[0], idx[1] );
   opr_data[vertexdof::stencilIndexFromVertex( vertices[2] )] += meanCoeff * localMatrix( idx[0], idx[2] );
}

template < typename ValueType >
inline void interpolate( const uint_t&                                               level,
                         Edge&                                                       edge,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& edgeMemoryId,
                         const ValueType&                                            scalar )
{
   size_t rowsize  = levelinfo::num_microvertices_per_edge( level );
   auto   edgeData = edge.getData( edgeMemoryId )->getPointer( level );
   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      edgeData[i] = scalar;
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                                                               level,
                         Edge&                                                                                       edge,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                                 edgeMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >&                  srcIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& ) >& expr )
{
   ValueType* edgeData = edge.getData( edgeMemoryId )->getPointer( level );

   std::vector< ValueType* > srcPtr;
   for ( const auto& src : srcIds )
   {
      srcPtr.push_back( edge.getData( src )->getPointer( level ) );
   }

   std::vector< ValueType > srcVector( srcIds.size() );

   Point3D xBlend;

   for ( const auto& it : vertexdof::macroedge::Iterator( level, 1 ) )
   {
      const Point3D coordinate = coordinateFromIndex( level, edge, it );
      const uint_t  idx        = vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C );

      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
         srcVector[k] = srcPtr[k][idx];
      }
      edge.getGeometryMap()->evalF( coordinate, xBlend );
      edgeData[idx] = expr( xBlend, srcVector );
   }
}

template < typename ValueType >
inline void swap( const uint_t&                                               level,
                  Edge&                                                       edge,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstID )
{
   auto srcData = edge.getData( srcID );
   auto dstData = edge.getData( dstID );
   srcData->swap( *dstData, level );
}

template < typename ValueType >
inline void assign( const uint_t&                                                              level,
                    Edge&                                                                      edge,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      ValueType tmp =
          scalars[0] * edge.getData( srcIds[0] )
                           ->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];

      for ( size_t k = 1; k < srcIds.size(); ++k )
      {
         tmp += scalars[k] *
                edge.getData( srcIds[k] )
                    ->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      }

      edge.getData( dstId )->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] =
          tmp;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                               level,
                 const Edge&                                                 edge,
                 const ValueType&                                            scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      edge.getData( dstId )->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +=
          scalar;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                                              level,
                 Edge&                                                                      edge,
                 const std::vector< ValueType >&                                            scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      auto tmp = ValueType( 0 );

      for ( size_t k = 0; k < srcIds.size(); ++k )
      {
         tmp += scalars[k] *
                edge.getData( srcIds[k] )
                    ->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      }

      edge.getData( dstId )->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +=
          tmp;
   }
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                              level,
                             Edge&                                                                      edge,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );
   auto   dst     = edge.getData( dstId )->getPointer( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      const uint_t idx = vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C );
      ValueType    tmp = edge.getData( srcIds[0] )->getPointer( level )[idx];

      for ( size_t k = 1; k < srcIds.size(); ++k )
      {
         tmp *= edge.getData( srcIds[k] )->getPointer( level )[idx];
      }

      dst[idx] = tmp;
   }
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                               level,
                      Edge&                                                       edge,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& lhsMemoryId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsMemoryId )
{
   walberla::math::KahanAccumulator< ValueType > scalarProduct;
   size_t                                        rowsize = levelinfo::num_microvertices_per_edge( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      scalarProduct += edge.getData( lhsMemoryId )
                           ->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] *
                       edge.getData( rhsMemoryId )
                           ->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline ValueType sum( const uint_t&                                               level,
                      const Edge&                                                 edge,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dataID,
                      const bool&                                                 absolute )
{
   auto   sum     = ValueType( 0 );
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto data = edge.getData( dataID )->getPointer( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      if ( absolute )
      {
         sum += std::abs( data[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] );
      }
      else
      {
         sum += data[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      }
   }
   return sum;
}

template < typename ValueType >
inline void apply( const uint_t&                                               level,
                   Edge&                                                       edge,
                   const PrimitiveDataID< StencilMemory< ValueType >, Edge >&  operatorId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                   UpdateType                                                  update )
{
   typedef stencilDirection sD;
   size_t                   rowsize = levelinfo::num_microvertices_per_edge( level );

   auto opr_data = edge.getData( operatorId )->getPointer( level );
   auto src      = edge.getData( srcId )->getPointer( level );
   auto dst      = edge.getData( dstId )->getPointer( level );

   ValueType tmp;

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
      const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
      const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

      const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
      const auto dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
      const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

      tmp = opr_data[stencilIdxW] * src[dofIdxW] + opr_data[stencilIdxC] * src[dofIdxC] + opr_data[stencilIdxE] * src[dofIdxE];

      for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
      {
         const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
         const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
         const auto stencilWeightW          = opr_data[stencilIdxWNeighborFace];
         const auto stencilWeightE          = opr_data[stencilIdxENeighborFace];
         const auto dofIdxWNeighborFace =
             vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
         const auto dofIdxENeighborFace =
             vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
         tmp += stencilWeightW * src[dofIdxWNeighborFace] + stencilWeightE * src[dofIdxENeighborFace];
      }

      for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
      {
         const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
         const auto dofIdx =
             vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
         tmp += opr_data[stencilIdx] * src[dofIdx];
      }

      if ( update == Replace )
      {
         dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = tmp;
      }
      else if ( update == Add )
      {
         dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] += tmp;
      }
   }
}

template < typename ValueType >
inline void smooth_gs( const uint_t&                                               level,
                       Edge&                                                       edge,
                       const PrimitiveDataID< StencilMemory< ValueType >, Edge >&  operatorId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsId )
{
   typedef stencilDirection sD;
   size_t                   rowsize = levelinfo::num_microvertices_per_edge( level );

   auto opr_data = edge.getData( operatorId )->getPointer( level );
   auto rhs      = edge.getData( rhsId )->getPointer( level );
   auto dst      = edge.getData( dstId )->getPointer( level );

   const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
   const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
   const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

   const auto invCenterWeight = 1.0 / opr_data[stencilIdxC];

   ValueType tmp;

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
      const auto dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
      const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

      tmp = rhs[dofIdxC];

      tmp -= opr_data[stencilIdxW] * dst[dofIdxW] + opr_data[stencilIdxE] * dst[dofIdxE];

      for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
      {
         const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
         const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
         const auto stencilWeightW          = opr_data[stencilIdxWNeighborFace];
         const auto stencilWeightE          = opr_data[stencilIdxENeighborFace];
         const auto dofIdxWNeighborFace =
             vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
         const auto dofIdxENeighborFace =
             vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
         tmp -= stencilWeightW * dst[dofIdxWNeighborFace] + stencilWeightE * dst[dofIdxENeighborFace];
      }

      for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
      {
         const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
         const auto dofIdx =
             vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
         tmp -= opr_data[stencilIdx] * dst[dofIdx];
      }

      dst[dofIdxC] = invCenterWeight * tmp;
   }
}

template < typename ValueType >
inline void smooth_sor( const uint_t&                                               level,
                        Edge&                                                       edge,
                        const PrimitiveDataID< StencilMemory< ValueType >, Edge >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsId,
                        ValueType                                                   relax,
                        const bool&                                                 backwards = false )
{
   typedef stencilDirection sD;
   size_t                   rowsize = levelinfo::num_microvertices_per_edge( level );

   auto opr_data = edge.getData( operatorId )->getPointer( level );
   auto rhs      = edge.getData( rhsId )->getPointer( level );
   auto dst      = edge.getData( dstId )->getPointer( level );

   const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
   const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
   const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

   const auto invCenterWeight = 1.0 / opr_data[stencilIdxC];

   ValueType tmp;

   const int start = backwards ? (int) rowsize - 2 : 1;
   const int stop  = backwards ? 0 : (int) rowsize - 1;
   const int incr  = backwards ? -1 : 1;

   for ( int ii = start; ii != stop; ii += incr )
   {
      const uint_t i       = uint_c( ii );
      const auto   dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
      const auto   dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
      const auto   dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

      tmp = rhs[dofIdxC];

      tmp -= opr_data[stencilIdxW] * dst[dofIdxW] + opr_data[stencilIdxE] * dst[dofIdxE];

      for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
      {
         const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
         const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
         const auto stencilWeightW          = opr_data[stencilIdxWNeighborFace];
         const auto stencilWeightE          = opr_data[stencilIdxENeighborFace];
         const auto dofIdxWNeighborFace =
             vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
         const auto dofIdxENeighborFace =
             vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
         tmp -= stencilWeightW * dst[dofIdxWNeighborFace] + stencilWeightE * dst[dofIdxENeighborFace];
      }

      for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
      {
         const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
         const auto dofIdx =
             vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
         tmp -= opr_data[stencilIdx] * dst[dofIdx];
      }

      dst[dofIdxC] = ( 1.0 - relax ) * dst[dofIdxC] + relax * invCenterWeight * tmp;
   }
}

template < typename ValueType >
inline void smooth_jac( const uint_t&                                               level,
                        Edge&                                                       edge,
                        const PrimitiveDataID< StencilMemory< ValueType >, Edge >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& tmpId )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto opr_data = edge.getData( operatorId )->getPointer( level );
   auto dst      = edge.getData( dstId )->getPointer( level );
   auto rhs      = edge.getData( rhsId )->getPointer( level );
   auto tmp      = edge.getData( tmpId )->getPointer( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] =
          rhs[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];

      for ( auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
      {
         dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] -=
             opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
             tmp[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
      }

      for ( auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
      {
         dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] -=
             opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
             tmp[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
      }

      if ( edge.getNumNeighborFaces() == 2 )
      {
         for ( auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
         {
            dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] -=
                opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
                tmp[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
         }
      }

      dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] /=
          opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
   }
}

template < typename ValueType >
inline void enumerate( const uint_t&                                               level,
                       Edge&                                                       edge,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                       ValueType&                                                  num )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      edge.getData( dstId )->getPointer( level )[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] =
          num;
      num++;
   }
}

template < typename ValueType >
inline void integrateDG( const uint_t&                                               level,
                         Edge&                                                       edge,
                         const std::shared_ptr< PrimitiveStorage >&                  storage,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsId,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsP1Id,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId )
{
   typedef stencilDirection sD;

   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto rhs   = edge.getData( rhsId )->getPointer( level );
   auto rhsP1 = edge.getData( rhsP1Id )->getPointer( level );
   auto dst   = edge.getData( dstId )->getPointer( level );

   real_t tmp;

   Face* face = storage->getFace( edge.neighborFaces()[0] );

   real_t weightedFaceArea0 = std::pow( 4.0, -walberla::real_c( level ) ) * face->area / 3.0;
   real_t weightedFaceArea1 = real_c( 0 );

   if ( edge.getNumNeighborFaces() == 2 )
   {
      face              = storage->getFace( edge.neighborFaces()[1] );
      weightedFaceArea1 = std::pow( 4.0, -walberla::real_c( level ) ) * face->area / 3.0;
   }

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      tmp = weightedFaceArea0 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_SW )] *
            ( 0.5 * 0.5 *
                  ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                    rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W )] ) +
              0.5 * 0.5 *
                  ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                    rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_S )] ) );
      tmp += weightedFaceArea0 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_BLUE_SE )] *
             ( 0.5 * 0.5 *
                   ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                     rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_S )] ) +
               0.5 * 0.5 *
                   ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                     rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_SE )] ) );
      tmp += weightedFaceArea0 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_SE )] *
             ( 0.5 * 0.5 *
                   ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                     rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_SE )] ) +
               0.5 * 0.5 *
                   ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                     rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E )] ) );

      if ( edge.getNumNeighborFaces() == 2 )
      {
         tmp += weightedFaceArea1 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_NW )] *
                ( 0.5 * 0.5 *
                      ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                        rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W )] ) +
                  0.5 * 0.5 *
                      ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                        rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_NW )] ) );
         tmp += weightedFaceArea1 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_BLUE_NW )] *
                ( 0.5 * 0.5 *
                      ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                        rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_NW )] ) +
                  0.5 * 0.5 *
                      ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                        rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_N )] ) );
         tmp += weightedFaceArea1 * rhs[facedof::macroedge::indexFaceFromVertex( level, i, sD::CELL_GRAY_NE )] *
                ( 0.5 * 0.5 *
                      ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                        rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_N )] ) +
                  0.5 * 0.5 *
                      ( rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] +
                        rhsP1[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E )] ) );
      }

      dst[vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C )] = ValueType( tmp );
   }
}

template < typename ValueType >
inline void printFunctionMemory( const uint_t&                                               level,
                                 const Edge&                                                 edge,
                                 const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId )
{
   ValueType* edgeMemory = edge.getData( dstId )->getPointer( level );
   using namespace std;
   cout << setfill( '=' ) << setw( 100 ) << "" << endl;
   cout << edge << std::left << setprecision( 1 ) << fixed << setfill( ' ' ) << endl;
   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( uint_t i = 0; i < ( rowsize - 1 ); ++i )
      {
         cout << setw( 5 ) << edgeMemory[hyteg::vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_N )]
              << "|";
      }
      cout << endl;
   }
   for ( uint_t i = 0; i < rowsize; ++i )
   {
      cout << setw( 5 ) << edgeMemory[hyteg::vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )]
           << "|";
   }
   cout << endl << "     |";
   for ( uint_t i = 1; i < rowsize; ++i )
   {
      cout << setw( 5 ) << edgeMemory[hyteg::vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_S )]
           << "|";
   }
   cout << endl << setfill( '=' ) << setw( 100 ) << "" << endl << setfill( ' ' );
}

template < typename ValueType >
inline ValueType getMaxValue( const uint_t& level, Edge& edge, const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId )
{
   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto src      = edge.getData( srcId )->getPointer( level );
   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      localMax = std::max( localMax, src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] );
   }

   return localMax;
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Edge& edge, const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId )
{
   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto src      = edge.getData( srcId )->getPointer( level );
   auto localMax = ValueType( 0.0 );

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      localMax =
          std::max( localMax, std::abs( src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] ) );
   }

   return localMax;
}

template < typename ValueType >
inline ValueType getMinValue( const uint_t& level, Edge& edge, const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId )
{
   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto src      = edge.getData( srcId )->getPointer( level );
   auto localMin = std::numeric_limits< ValueType >::max();

   for ( size_t i = 1; i < rowsize - 1; ++i )
   {
      localMin = std::min( localMin, src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] );
   }

   return localMin;
}

inline void saveOperator( const uint_t&                                           level,
                          Edge&                                                   edge,
                          const PrimitiveStorage&                                 storage,
                          const PrimitiveDataID< StencilMemory< real_t >, Edge >& operatorId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& srcId,
                          const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& dstId,
                          const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto opr_data = edge.getData( operatorId )->getPointer( level );
   auto src      = edge.getData( srcId )->getPointer( level );
   auto dst      = edge.getData( dstId )->getPointer( level );

   for ( uint_t i = 1; i < rowsize - 1; ++i )
   {
      idx_t dstint = dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      idx_t srcint = src[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      mat->addValue(
          uint_c( dstint ), uint_c( srcint ), opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] );

      for ( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
      {
         srcint = src[vertexdof::macroedge::indexFromVertex( level, i, neighbor )];
         mat->addValue( uint_c( dstint ), uint_c( srcint ), opr_data[vertexdof::stencilIndexFromVertex( neighbor )] );
      }

      for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
      {
         srcint = src[vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, stencilDirection::VERTEX_W )];
         mat->addValue( uint_c( dstint ),
                        uint_c( srcint ),
                        opr_data[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_W, neighborFace )] );

         srcint = src[vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, stencilDirection::VERTEX_E )];
         mat->addValue( uint_c( dstint ),
                        uint_c( srcint ),
                        opr_data[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_E, neighborFace )] );
      }

      for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
      {
         const auto& neighborCell              = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
         const auto  localEdgeIDOnNeighborCell = neighborCell.getLocalEdgeID( edge.getID() );
         if ( localEdgeIDOnNeighborCell == 1 || localEdgeIDOnNeighborCell == 4 )
         {
            // Since the functions we access here carry the petsc vector indices, we cannot simply also loop over
            // ghost layer DoFs that do not exist. In the apply kernel this is okay, as we only add zeros in that case.
            // Therefore we check if there are inner vertices - this only applies for macro-edge IDs 1 and 4.
            srcint =
                src[vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCellID, edge.getNumNeighborFaces() )];
            mat->addValue(
                uint_c( dstint ),
                uint_c( srcint ),
                opr_data[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCellID, edge.getNumNeighborFaces() )] );
         }
      }
   }
}

inline void saveIdentityOperator( const uint_t&                                           level,
                                  Edge&                                                   edge,
                                  const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& dstId,
                                  const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   auto dst = edge.getData( dstId )->getPointer( level );

   for ( uint_t i = 1; i < rowsize - 1; ++i )
   {
      idx_t dstint = dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )];
      mat->addValue( uint_c( dstint ), uint_c( dstint ), 1.0 );
   }
}

template < typename ValueType >
inline void createVectorFromFunction( const uint_t&                                               level,
                                      Edge&                                                       edge,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   idx_t rowsize = (idx_t) levelinfo::num_microvertices_per_edge( level );

   auto src       = edge.getData( srcId )->getPointer( level );
   auto numerator = edge.getData( numeratorId )->getPointer( level );

   for ( uint_t i = 1; i < uint_c( rowsize - 1 ); i++ )
   {
      vec->setValue( numerator[i], src[i] );
   }
}

template < typename ValueType >
inline void createFunctionFromVector( const uint_t&                                               level,
                                      Edge&                                                       edge,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   idx_t rowsize = (idx_t) levelinfo::num_microvertices_per_edge( level );

   auto numerator = edge.getData( numeratorId )->getPointer( level );

   for ( uint_t i = 1; i < uint_c( rowsize - 1 ); i++ )
   {
      edge.getData( srcId )->getPointer( level )[i] = vec->getValue( numerator[i] );
   }
}

inline void applyDirichletBC( const uint_t&                                           level,
                              Edge&                                                   edge,
                              std::vector< idx_t >&                                   mat,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& numeratorId )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   for ( uint_t i = 1; i < rowsize - 1; i++ )
   {
      mat.push_back( edge.getData( numeratorId )->getPointer( level )[i] );
   }
}

} // namespace macroedge
} // namespace vertexdof
} // namespace hyteg
