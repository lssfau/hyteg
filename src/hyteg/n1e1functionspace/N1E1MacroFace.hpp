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
#include "hyteg/primitives/Face.hpp"

namespace hyteg {
namespace n1e1 {

template < typename ValueType >
class N1E1VectorFunction;

namespace macroface {

using walberla::uint_t;
template < typename ValueType >
using VectorType = typename N1E1VectorFunction< ValueType >::VectorType;

inline Point3D microEdgeDirection( const uint_t& level, const Face& face, const edgedof::EdgeDoFOrientation& orientation )
{
   const real_t stepFrequency = 1.0 / real_c( levelinfo::num_microedges_per_edge( level ) );

   switch ( orientation )
   {
   case edgedof::EdgeDoFOrientation::X:
      return ( face.getCoordinates()[1] - face.getCoordinates()[0] ) * stepFrequency;
   case edgedof::EdgeDoFOrientation::Y:
      return ( face.getCoordinates()[2] - face.getCoordinates()[0] ) * stepFrequency;
   case edgedof::EdgeDoFOrientation::XY:
      return ( face.getCoordinates()[2] - face.getCoordinates()[1] ) * stepFrequency;
   default:
      WALBERLA_ABORT( "wrong orienation" )
   }
}

inline void add( const uint_t&                                            level,
                 Face&                                                    face,
                 const VectorType< real_t >&                              vector,
                 const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId )
{
   using ValueType = real_t;

   const Point3D horizontalMicroEdgeDirection = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::X );
   const Point3D verticalMicroEdgeDirection   = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::Y );
   const Point3D diagonalMicroEdgeDirection   = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::XY );

   // x ↦ ∫ₑ x·t dΓ, direction = tangent·length
   const ValueType dofScalarHorizontal = vector.dot( horizontalMicroEdgeDirection.vector_ );
   const ValueType dofScalarVertical   = vector.dot( verticalMicroEdgeDirection.vector_ );
   const ValueType dofScalarDiagonal   = vector.dot( diagonalMicroEdgeDirection.vector_ );

   auto dstData = face.getData( dstId )->getPointer( level );

   for ( const auto& it : edgedof::macroface::Iterator( level ) )
   {
      const uint_t idxHorizontal = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
      const uint_t idxVertical   = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
      const uint_t idxDiagonal   = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         dstData[idxHorizontal] += dofScalarHorizontal;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         dstData[idxVertical] += dofScalarVertical;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( ( it.col() + it.row() ) != idx_t( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         dstData[idxDiagonal] += dofScalarDiagonal;
      }
   }
}

inline void interpolate( const uint_t&                                            level,
                         Face&                                                    face,
                         const PrimitiveDataID< FunctionMemory< real_t >, Face >& faceMemoryId,
                         const VectorType< real_t >&                              constant )
{
   using ValueType = real_t;

   auto faceData = face.getData( faceMemoryId )->getPointer( level );

   const Point3D horizontalMicroEdgeDirection = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::X );
   const Point3D verticalMicroEdgeDirection   = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::Y );
   const Point3D diagonalMicroEdgeDirection   = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::XY );

   // x ↦ ∫ₑ x·t dΓ, direction = tangent·length
   const ValueType dofScalarHorizontal = constant.dot( horizontalMicroEdgeDirection.vector_ );
   const ValueType dofScalarVertical   = constant.dot( verticalMicroEdgeDirection.vector_ );
   const ValueType dofScalarDiagonal   = constant.dot( diagonalMicroEdgeDirection.vector_ );

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
   {
      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
         faceData[edgedof::macroface::horizontalIndex( level, it.col(), it.row() )] = dofScalarHorizontal;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         faceData[edgedof::macroface::verticalIndex( level, it.col(), it.row() )] = dofScalarVertical;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( ( it.col() + it.row() ) != idx_t( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         faceData[edgedof::macroface::diagonalIndex( level, it.col(), it.row() )] = dofScalarDiagonal;
      }
   }
}

inline void
    interpolate( const uint_t&                                                                      level,
                 Face&                                                                              face,
                 const PrimitiveDataID< FunctionMemory< real_t >, Face >&                           faceMemoryId,
                 const std::vector< std::reference_wrapper< const N1E1VectorFunction< real_t > > >& srcFunctions,
                 const std::function< VectorType< real_t >( const Point3D&, const std::vector< VectorType< real_t > >& ) >& expr )
{
   using ValueType = real_t;

   auto                                   faceData = face.getData( faceMemoryId )->getPointer( level );
   std::vector< VectorType< ValueType > > srcVector( srcFunctions.size() );

   const Point3D faceBottomLeftCoords = face.getCoordinates()[0];

   const Point3D horizontalMicroEdgeDirection = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::X );
   const Point3D verticalMicroEdgeDirection   = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::Y );
   const Point3D diagonalMicroEdgeDirection   = microEdgeDirection( level, face, edgedof::EdgeDoFOrientation::XY );

   const Point3D horizontalMicroEdgeOffset = horizontalMicroEdgeDirection * 0.5;
   const Point3D verticalMicroEdgeOffset   = verticalMicroEdgeDirection * 0.5;

   Point3D xBlend;

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
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
         face.getGeometryMap()->evalF( horizontalMicroEdgePosition, xBlend );

         for ( uint_t k = 0; k < srcFunctions.size(); ++k )
         {
            srcFunctions[k].get().evaluate( xBlend, level, srcVector[k] );
         }

         const VectorType< ValueType > vector    = expr( xBlend, srcVector );
         const ValueType               dofScalar = vector.dot( horizontalMicroEdgeDirection.vector_ );

         faceData[edgedof::macroface::horizontalIndex( level, it.col(), it.row() )] = dofScalar;
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
         face.getGeometryMap()->evalF( verticalMicroEdgePosition, xBlend );

         for ( uint_t k = 0; k < srcFunctions.size(); ++k )
         {
            srcFunctions[k].get().evaluate( xBlend, level, srcVector[k] );
         }

         const VectorType< ValueType > vector    = expr( xBlend, srcVector );
         const ValueType               dofScalar = vector.dot( verticalMicroEdgeDirection.vector_ );

         faceData[edgedof::macroface::verticalIndex( level, it.col(), it.row() )] = dofScalar;
      }

      // Do not update diagonal DoFs at diagonal border
      if ( ( it.col() + it.row() ) != idx_t( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         face.getGeometryMap()->evalF( diagonalMicroEdgePosition, xBlend );

         for ( uint_t k = 0; k < srcFunctions.size(); ++k )
         {
            srcFunctions[k].get().evaluate( xBlend, level, srcVector[k] );
         }

         const VectorType< ValueType > vector    = expr( xBlend, srcVector );
         const ValueType               dofScalar = vector.dot( diagonalMicroEdgeDirection.vector_ );

         faceData[edgedof::macroface::diagonalIndex( level, it.col(), it.row() )] = dofScalar;
      }
   }
}

} // namespace macroface
} // namespace n1e1
} // namespace hyteg
