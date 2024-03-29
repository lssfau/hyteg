/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/n1e1functionspace/N1E1MacroCell.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Edge.hpp"

namespace hyteg {
namespace n1e1 {

namespace macroedge {

using walberla::uint_t;
template < typename ValueType >
using VectorType = PointND< ValueType, 3>;

inline void add( const uint_t&                                            level,
                 Edge&                                                    edge,
                 const VectorType< real_t >&                              vector,
                 const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstId )
{
   using ValueType = real_t;

   const VectorType< ValueType > microEdgeDirection = edge.getDirection() / real_c( levelinfo::num_microedges_per_edge( level ) );

   // x ↦ ∫ₑ x·t dΓ, direction = tangent·length
   const ValueType dofScalar = vector.dot( microEdgeDirection );

   auto dstData = edge.getData( dstId )->getPointer( level );

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const uint_t idx = edgedof::macroedge::index( level, it.x() );
      dstData[idx] += dofScalar;
   }
}

inline void
    interpolate( const uint_t&                                                           level,
                 Edge&                                                                   edge,
                 Cell&                                                                   neighborCell,
                 const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                edgeMemoryId,
                 const std::vector< PrimitiveDataID< FunctionMemory< real_t >, Cell > >& srcIds,
                 const std::function< VectorType< real_t >( const Point3D&, const std::vector< VectorType< real_t > >& ) >& expr )
{
   using ValueType = real_t;

   auto                                   edgeData = edge.getData( edgeMemoryId )->getPointer( level );
   std::vector< VectorType< ValueType > > srcVector( srcIds.size() );

   const Point3D leftCoords  = edge.getCoordinates()[0];
   const Point3D rightCoords = edge.getCoordinates()[1];

   const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( level ) );
   const VectorType< ValueType > microEdgeDirection = edge.getDirection() / real_c( levelinfo::num_microedges_per_edge( level ) );

   Point3D  xBlend;
   Matrix3r DF;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 * it.x() ) * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xBlend );
      edge.getGeometryMap()->evalDF( currentCoordinates, DF );

      for ( uint_t k = 0; k < srcIds.size(); ++k )
      {
         srcVector[k] = macrocell::evaluate( level, neighborCell, currentCoordinates, srcIds[k] );
      }

      const VectorType< ValueType > vector    = DF.transpose() * expr( xBlend, srcVector );
      const ValueType               dofScalar = vector.dot( microEdgeDirection );

      edgeData[edgedof::macroedge::index( level, it.x() )] = dofScalar;
   }
}

inline void interpolate( const uint_t&                                            level,
                         Edge&                                                    edge,
                         Cell&                                                    neighborCell,
                         const PrimitiveDataID< FunctionMemory< real_t >, Edge >& edgeMemoryId,
                         const VectorType< real_t >&                              constant )
{
   using ValueType = real_t;

   const auto geometryMap = edge.getGeometryMap();
   if ( !geometryMap->isAffine() )
   {
      // If the blending map is not affine, the vector field is not constant in computational space.
      // In that case, we delegate to the non-constant interpolation routine.
      interpolate(
          level, edge, neighborCell, edgeMemoryId, {}, [&]( const Point3D&, const std::vector< VectorType< real_t > >& ) {
             return constant;
          } );
      return;
   }

   // Geometry map is affine ⇒ DF is constant
   Matrix3r DF;
   geometryMap->evalDF( {}, DF );
   const VectorType< real_t > valComp = DF.transpose() * constant;

   const VectorType< ValueType > microEdgeDirection = edge.getDirection() / real_c( levelinfo::num_microedges_per_edge( level ) );

   // x ↦ ∫ₑ x·t dΓ, direction = tangent·length
   const ValueType dofScalar = valComp.dot( microEdgeDirection );

   auto edgeData = edge.getData( edgeMemoryId )->getPointer( level );

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const uint_t idx = edgedof::macroedge::index( level, it.x() );
      edgeData[idx]    = dofScalar;
   }
}

} // namespace macroedge
} // namespace n1e1
} // namespace hyteg
