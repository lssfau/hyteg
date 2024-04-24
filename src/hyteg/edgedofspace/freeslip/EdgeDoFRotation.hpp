/*
 * Copyright (c) 2023 Ponsuganth Ilangovan P
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

void rotationMatrix2D( const Point3D& normal, Matrix2r& rotation )
{
   rotation( 0, 0 ) = -normal[1];
   rotation( 0, 1 ) = normal[0];
   rotation( 1, 0 ) = normal[0];
   rotation( 1, 1 ) = normal[1];
}

void rotationMatrix3D( const Point3D& n, Matrix3r& rotation, bool transpose )
{
   Point3D nCross[] = {
       n.cross( Point3D( 1.0, 0.0, 0.0 ) ), n.cross( Point3D( 0.0, 1.0, 0.0 ) ), n.cross( Point3D( 0.0, 0.0, 1.0 ) ) };

   real_t nCrossComp[] = { nCross[0].norm(), nCross[1].norm(), nCross[2].norm() };

   uint_t iMax = nCrossComp[0] > nCrossComp[1] ? 0 : 1;
   iMax        = nCrossComp[iMax] > nCrossComp[2] ? iMax : 2;

   Point3D t1 = nCross[iMax].normalized();

   Point3D t2 = n.cross( t1 );

   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "n.t1 = %f, n.t2 = %f, t1.t2 = %f", n.dot( t1 ), n.dot( t2 ), t1.dot( t2 ) ) );

   rotation( 0, 0 ) = t1[0];
   rotation( 0, 1 ) = t1[1];
   rotation( 0, 2 ) = t1[2];
   rotation( 1, 0 ) = t2[0];
   rotation( 1, 1 ) = t2[1];
   rotation( 1, 2 ) = t2[2];
   rotation( 2, 0 ) = n[0];
   rotation( 2, 1 ) = n[1];
   rotation( 2, 2 ) = n[2];

   if(transpose)
   {
      rotation.transposeInPlace();
   }
}

void addMatrix3D( const uint_t                                            idx,
                  uint_t                                                  level,
                  const std::shared_ptr< SparseMatrixProxy >&             mat,
                  const Face&                                             face,
                  const PrimitiveDataID< FunctionMemory< idx_t >, Face >& dstIdU,
                  const PrimitiveDataID< FunctionMemory< idx_t >, Face >& dstIdV,
                  const PrimitiveDataID< FunctionMemory< idx_t >, Face >& dstIdW,
                  Point3D                                                 normal,
                  bool transpose )
{
   auto dstU = face.getData( dstIdU )->getPointer( level );
   auto dstV = face.getData( dstIdV )->getPointer( level );
   auto dstW = face.getData( dstIdW )->getPointer( level );

   Matrix3r rotation;
   rotationMatrix3D( normal, rotation, transpose );

   const idx_t idxUVW[] = { dstU[idx], dstV[idx], dstW[idx] };

   for ( uint_t iMat = 0U; iMat < 3U; iMat++ )
   {
      for ( uint_t jMat = 0U; jMat < 3U; jMat++ )
      {
         mat->addValue( uint_c( idxUVW[iMat] ), uint_c( idxUVW[jMat] ), rotation( iMat, jMat ) );
      }
   }
}
namespace macroface {

template < typename ValueType >
inline void rotation3D( uint_t                                                      level,
                        const Face&                                                 face,
                        const std::shared_ptr< PrimitiveStorage >&                  storage,
                        const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdU,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdV,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstIdW, bool transpose )
{
   auto dstU = face.getData( dstIdU )->getPointer( level );
   auto dstV = face.getData( dstIdV )->getPointer( level );
   auto dstW = face.getData( dstIdW )->getPointer( level );

   Point3D  normal;
   Matrix3r rotation;
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
      const Point3D verticalMicroEdgePosition = faceBottomLeftCoords + ( ( real_c( it.x() ) * 2 ) * horizontalMicroEdgeOffset +
                                                                         ( real_c( it.y() ) * 2 + 1 ) * verticalMicroEdgeOffset );
      const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

      // Do not update horizontal DoFs at bottom
      if ( it.y() != 0 )
      {
         face.getGeometryMap()->evalF( horizontalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );
         rotationMatrix3D( normal, rotation, transpose );

         const uint_t idx = edgedof::macroface::horizontalIndex( level, it.x(), it.y() );

         in[0] = dstU[idx];
         in[1] = dstV[idx];
         in[2] = dstW[idx];

         out = rotation * in;

         dstU[idx] = out[0];
         dstV[idx] = out[1];
         dstW[idx] = out[2];
      }

      // Do not update vertical DoFs at left border
      if ( it.x() != 0 )
      {
         face.getGeometryMap()->evalF( verticalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );
         rotationMatrix3D( normal, rotation, transpose );

         const uint_t idx = edgedof::macroface::verticalIndex( level, it.x(), it.y() );

         in[0] = dstU[idx];
         in[1] = dstV[idx];
         in[2] = dstW[idx];

         out = rotation * in;

         dstU[idx] = out[0];
         dstV[idx] = out[1];
         dstW[idx] = out[2];
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.x() + it.y() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         face.getGeometryMap()->evalF( diagonalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );
         rotationMatrix3D( normal, rotation, transpose );

         const uint_t idx = edgedof::macroface::diagonalIndex( level, it.x(), it.y() );

         in[0] = dstU[idx];
         in[1] = dstV[idx];
         in[2] = dstW[idx];

         out = rotation * in;

         dstU[idx] = out[0];
         dstV[idx] = out[1];
         dstW[idx] = out[2];
      }
   }
}

inline void saveRotationOperator3D( uint_t                                                   level,
                                    const Face&                                              face,
                                    const std::shared_ptr< PrimitiveStorage >&               storage,
                                    const std::function< void( const Point3D&, Point3D& ) >& normal_function,
                                    const PrimitiveDataID< FunctionMemory< idx_t >, Face >&  dstIdU,
                                    const PrimitiveDataID< FunctionMemory< idx_t >, Face >&  dstIdV,
                                    const PrimitiveDataID< FunctionMemory< idx_t >, Face >&  dstIdW,
                                    const std::shared_ptr< SparseMatrixProxy >&              mat,
                                    bool                                                     transpose )
{
   if ( face.getNumNeighborCells() == 2 )
   {
      WALBERLA_ABORT( "Cannot project normals if not a boundary face" );
   }

   Point3D normal;
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
      const Point3D verticalMicroEdgePosition = faceBottomLeftCoords + ( ( real_c( it.x() ) * 2 ) * horizontalMicroEdgeOffset +
                                                                         ( real_c( it.y() ) * 2 + 1 ) * verticalMicroEdgeOffset );
      const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

      // Do not update horizontal DoFs at bottom
      if ( it.y() != 0 )
      {
         face.getGeometryMap()->evalF( horizontalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );

         const uint_t idx = edgedof::macroface::horizontalIndex( level, it.x(), it.y() );
         addMatrix3D( idx, level, mat, face, dstIdU, dstIdV, dstIdW, normal, transpose );
      }

      // Do not update vertical DoFs at left border
      if ( it.x() != 0 )
      {
         face.getGeometryMap()->evalF( verticalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );

         const uint_t idx = edgedof::macroface::verticalIndex( level, it.x(), it.y() );
         addMatrix3D( idx, level, mat, face, dstIdU, dstIdV, dstIdW, normal, transpose );
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.x() + it.y() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
      {
         face.getGeometryMap()->evalF( diagonalMicroEdgePosition, xBlend );
         normal_function( xBlend, normal );

         const uint_t idx = edgedof::macroface::diagonalIndex( level, it.x(), it.y() );
         addMatrix3D( idx, level, mat, face, dstIdU, dstIdV, dstIdW, normal, transpose );
      }
   }
}

} // namespace macroface

namespace macroedge {

template < typename ValueType >
inline void rotation2D( uint_t                                                      level,
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
   Matrix2r rotation;
   Point2D  in;
   Point2D  out;

   Point3D xPhy;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 ) * it.x() * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xPhy );

      normal_function( xPhy, normal );
      rotationMatrix2D( normal, rotation );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );

      in[0] = dstU[idx];
      in[1] = dstV[idx];

      out = rotation * in;

      dstU[idx] = out[0];
      dstV[idx] = out[1];
   }
}

template < typename ValueType >
inline void rotation3D( uint_t                                                      level,
                        const Edge&                                                 edge,
                        const std::shared_ptr< PrimitiveStorage >&                  storage,
                        const std::function< void( const Point3D&, Point3D& ) >&    normal_function,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdU,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdV,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstIdW, bool transpose )
{
   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );
   auto dstW = edge.getData( dstIdW )->getPointer( level );

   const Point3D leftCoords  = edge.getCoordinates()[0];
   const Point3D rightCoords = edge.getCoordinates()[1];

   const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( level ) );

   Point3D  normal;
   Matrix3r rotation;
   Point3D  in;
   Point3D  out;

   Point3D xPhy;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 ) * it.x() * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xPhy );

      normal_function( xPhy, normal );
      rotationMatrix3D( normal, rotation, transpose );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );

      in[0] = dstU[idx];
      in[1] = dstV[idx];
      in[2] = dstW[idx];

      out = rotation * in;

      dstU[idx] = out[0];
      dstV[idx] = out[1];
      dstW[idx] = out[2];
   }
}

inline void saveRotationOperator2D( uint_t                                                   level,
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
   Matrix2r rotation;
   Point2D  in;
   Point2D  out;

   Point3D xPhy;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 ) * it.x() * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xPhy );

      normal_function( xPhy, normal );
      rotationMatrix2D( normal, rotation );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );

      const auto idxU = dstU[idx];
      const auto idxV = dstV[idx];

      mat->addValue( uint_c( idxU ), uint_c( idxU ), rotation( 0, 0 ) );
      mat->addValue( uint_c( idxU ), uint_c( idxV ), rotation( 0, 1 ) );
      mat->addValue( uint_c( idxV ), uint_c( idxU ), rotation( 1, 0 ) );
      mat->addValue( uint_c( idxV ), uint_c( idxV ), rotation( 1, 1 ) );
   }
}

inline void saveRotationOperator3D( uint_t                                                   level,
                                    const Edge&                                              edge,
                                    const std::shared_ptr< PrimitiveStorage >&               storage,
                                    const std::function< void( const Point3D&, Point3D& ) >& normal_function,
                                    const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&  dstIdU,
                                    const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&  dstIdV,
                                    const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&  dstIdW,
                                    const std::shared_ptr< SparseMatrixProxy >&              mat,
                                    bool                                                     transpose )
{
   auto dstU = edge.getData( dstIdU )->getPointer( level );
   auto dstV = edge.getData( dstIdV )->getPointer( level );
   auto dstW = edge.getData( dstIdW )->getPointer( level );

   const Point3D leftCoords  = edge.getCoordinates()[0];
   const Point3D rightCoords = edge.getCoordinates()[1];

   const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( level ) );

   Point3D  normal;
   Matrix3r rotation;

   Point3D xPhy;

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const Point3D currentCoordinates = leftCoords + microEdgeOffset + real_c( 2 ) * it.x() * microEdgeOffset;
      edge.getGeometryMap()->evalF( currentCoordinates, xPhy );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );

      normal_function( xPhy, normal );
      rotationMatrix3D( normal, rotation, transpose );

      const idx_t idxUVW[] = { dstU[idx], dstV[idx], dstW[idx] };

      for ( uint_t iMat = 0U; iMat < 3U; iMat++ )
      {
         for ( uint_t jMat = 0U; jMat < 3U; jMat++ )
         {
            mat->addValue( uint_c( idxUVW[iMat] ), uint_c( idxUVW[jMat] ), rotation( iMat, jMat ) );
         }
      }
   }
}

} // namespace macroedge

} // namespace edgedof
} // namespace hyteg
