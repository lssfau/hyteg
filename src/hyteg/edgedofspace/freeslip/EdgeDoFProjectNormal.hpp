/*
 * Copyright (c) 2020 Daniel Drzisga.
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

#include "hyteg/Levelinfo.hpp"
#include "hyteg/Macros.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/polynomial/PolynomialEvaluator.hpp"
#include "hyteg/primitives/Face.hpp"

using walberla::real_c;
using walberla::uint_t;

namespace hyteg {
namespace edgedof {

void projectionMatrix2D( const Point3D& normal, Matrix2r& projection )
{
   projection( 0, 0 ) = 1.0 - normal[0] * normal[0];
   projection( 0, 1 ) = -normal[0] * normal[1];
   projection( 1, 0 ) = projection( 0, 1 );
   projection( 1, 1 ) = 1.0 - normal[1] * normal[1];
}

void projectionMatrix3D( const Point3D& normal, Matrix3r& projection )
{
   projection( 0, 0 ) = 1.0 - normal[0] * normal[0];
   projection( 0, 1 ) = -normal[0] * normal[1];
   projection( 0, 2 ) = -normal[0] * normal[2];
   projection( 1, 0 ) = projection( 0, 1 );
   projection( 1, 1 ) = 1.0 - normal[1] * normal[1];
   projection( 1, 2 ) = -normal[1] * normal[2];
   projection( 2, 0 ) = projection( 0, 2 );
   projection( 2, 1 ) = projection( 1, 2 );
   projection( 2, 2 ) = 1.0 - normal[2] * normal[2];
}

namespace macroface {

template < typename ValueType >
inline void projectNormal3D( uint_t                                                      level,
                             const Face&                                                 face,
                             const std::shared_ptr< PrimitiveStorage >&                  storage,
                             const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdV,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdW )
{
   auto dstU = face.getData( dstIdU )->getPointer( level );
   auto dstV = face.getData( dstIdV )->getPointer( level );
   auto dstW = face.getData( dstIdW )->getPointer( level );

   Point3D  normal;
   Matrix3r projection;
   Point3D  in;
   Point3D  out;

   Point3D x;
   Point3D xBlend;

   const Point3D faceBottomLeftCoords  = face.getCoordinates()[0];
   const Point3D faceBottomRightCoords = face.getCoordinates()[1];
   const Point3D faceTopLeftCoords     = face.getCoordinates()[2];

   const Point3D horizontalMicroEdgeOffset =
       ( ( faceBottomRightCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;
   const Point3D verticalMicroEdgeOffset =
       ( ( faceTopLeftCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;

   for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
   {
      const Point3D horizontalMicroEdgePosition =
          faceBottomLeftCoords +
          ( ( real_c( it.x() ) * 2 + 1 ) * horizontalMicroEdgeOffset + ( real_c( it.y() ) * 2 ) * verticalMicroEdgeOffset );
      const Point3D verticalMicroEdgePosition =
          faceBottomLeftCoords +
          ( ( real_c( it.x() ) * 2 ) * horizontalMicroEdgeOffset + ( real_c( it.y() ) * 2 + 1 ) * verticalMicroEdgeOffset );
      const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

      // Do not update horizontal DoFs at bottom
      if ( it.y() != 0 )
      {
         face.getGeometryMap()->evalF( horizontalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );
         projectionMatrix3D( normal, projection );

         const uint_t idx = edgedof::macroface::horizontalIndex( level, it.x(), it.y() );

         in[0] = dstU[idx];
         in[1] = dstV[idx];
         in[2] = dstW[idx];

         out = projection * in;

         dstU[idx] = out[0];
         dstV[idx] = out[1];
         dstW[idx] = out[2];
      }

      // Do not update vertical DoFs at left border
      if ( it.x() != 0 )
      {
         face.getGeometryMap()->evalF( verticalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );
         projectionMatrix3D( normal, projection );

         const uint_t idx = edgedof::macroface::verticalIndex( level, it.x(), it.y() );

         in[0] = dstU[idx];
         in[1] = dstV[idx];
         in[2] = dstW[idx];

         out = projection * in;

         dstU[idx] = out[0];
         dstV[idx] = out[1];
         dstW[idx] = out[2];
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.x() + it.y() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         face.getGeometryMap()->evalF( diagonalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );
         projectionMatrix3D( normal, projection );

         const uint_t idx = edgedof::macroface::diagonalIndex( level, it.x(), it.y() );

         in[0] = dstU[idx];
         in[1] = dstV[idx];
         in[2] = dstW[idx];

         out = projection * in;

         dstU[idx] = out[0];
         dstV[idx] = out[1];
         dstW[idx] = out[2];
      }
   }
}

} // namespace macroface

namespace macroedge {

template < typename ValueType >
inline void projectNormal2D( uint_t                                                      level,
                             const Edge&                                                 edge,
                             const std::shared_ptr< PrimitiveStorage >&                  storage,
                             const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdV )
{
   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );

   const Point3D leftCoords  = edge.getCoordinates()[0];
   const Point3D rightCoords = edge.getCoordinates()[1];

   const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( level ) );

   Point3D  normal;
   Matrix2r projection;
   Point2D  in;
   Point2D  out;

   Point3D xPhy;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 ) * it.x() * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xPhy );

      normal_function( xPhy, normal );
      projectionMatrix2D( normal, projection );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );

      in[0] = dstU[idx];
      in[1] = dstV[idx];

      out = projection * in;

      dstU[idx] = out[0];
      dstV[idx] = out[1];
   }
}

template < typename ValueType >
inline void projectNormal3D( uint_t                                                      level,
                             const Edge&                                                 edge,
                             const std::shared_ptr< PrimitiveStorage >&                  storage,
                             const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdU,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdV,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdW )
{
   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );
   auto dstW = edge.getData( dstIdW )->getPointer( level );

   const Point3D leftCoords  = edge.getCoordinates()[0];
   const Point3D rightCoords = edge.getCoordinates()[1];

   const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( level ) );

   Point3D  normal;
   Matrix3r projection;
   Point3D  in;
   Point3D  out;

   Point3D xPhy;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 ) * it.x() * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xPhy );

      normal_function( xPhy, normal );
      projectionMatrix3D( normal, projection );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );

      in[0] = dstU[idx];
      in[1] = dstV[idx];
      in[2] = dstW[idx];

      out = projection * in;

      dstU[idx] = out[0];
      dstV[idx] = out[1];
      dstW[idx] = out[2];
   }
}

inline void saveProjectNormalOperator2D( uint_t                                                   level,
                                         const Edge&                                              edge,
                                         const std::shared_ptr< PrimitiveStorage >&               storage,
                                         const std::function< void( const Point3D&, Point3D& ) >& normal_function,
                                         const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&  dstIdU,
                                         const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&  dstIdV,
                                         const std::shared_ptr< SparseMatrixProxy >&              mat )
{
   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );

   const Point3D leftCoords  = edge.getCoordinates()[0];
   const Point3D rightCoords = edge.getCoordinates()[1];

   const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( level ) );

   Point3D  normal;
   Matrix2r projection;
   Point2D  in;
   Point2D  out;

   Point3D xPhy;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 ) * it.x() * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xPhy );

      normal_function( xPhy, normal );
      projectionMatrix2D( normal, projection );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );

      const auto idxU = dstU[idx];
      const auto idxV = dstV[idx];

      mat->addValue( uint_c( idxU ), uint_c( idxU ), projection(0, 0) );
      mat->addValue( uint_c( idxU ), uint_c( idxV ), projection( 0, 1 ) );
      mat->addValue( uint_c( idxV ), uint_c( idxU ), projection( 1, 0 ) );
      mat->addValue( uint_c( idxV ), uint_c( idxV ), projection( 1, 1 ) );
   }
}

} // namespace macroedge

} // namespace edgedof
} // namespace hyteg
