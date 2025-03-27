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

#include "hyteg/gridtransferoperators/N1E1toN1E1Prolongation.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
namespace n1e1 {

void N1E1toN1E1Prolongation::prolongateAdditively( const N1E1VectorFunction< real_t >& function,
                                                   const uint_t&                       sourceLevel,
                                                   const DoFType&                      flag,
                                                   const UpdateType&                   updateType ) const
{
   WALBERLA_ASSERT( function.getStorage()->hasGlobalCells(), "N1E1 prolongation only implemented in 3D." )

   const uint_t destinationLevel = sourceLevel + 1;

   function.communicate< Edge, Face >( sourceLevel );
   function.communicate< Face, Cell >( sourceLevel );

   for ( const auto& it : function.getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         auto data = edge.getData( function.getDoFs()->getEdgeDataID() );

         if ( updateType == Replace )
         {
            data->setToZero( destinationLevel );
         }
         prolongateMacroEdge( data->getPointer( sourceLevel ), data->getPointer( destinationLevel ), sourceLevel );
      }
   }

   for ( const auto& it : function.getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         auto data = face.getData( function.getDoFs()->getFaceDataID() );

         if ( updateType == Replace )
         {
            data->setToZero( destinationLevel );
         }
         prolongateMacroFace( data->getPointer( sourceLevel ), data->getPointer( destinationLevel ), sourceLevel );
      }
   }

   for ( const auto& it : function.getStorage()->getCells() )
   {
      const Cell& cell = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         auto data = cell.getData( function.getDoFs()->getCellDataID() );

         if ( updateType == Replace )
         {
            data->setToZero( destinationLevel );
         }
         prolongateMacroCell( data->getPointer( sourceLevel ), data->getPointer( destinationLevel ), sourceLevel );
      }
   }
}

void N1E1toN1E1Prolongation::prolongateMacroEdge( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   uint_t n = levelinfo::num_microedges_per_edge( sourceLevel );

   for ( idx_t i = 0; i < idx_t( n ); ++i )
   {
      dst[edgedof::macroedge::index( sourceLevel + 1, 2 * i )] += real_c( 0.5 ) * src[edgedof::macroedge::index( sourceLevel, i )];
      dst[edgedof::macroedge::index( sourceLevel + 1, 2 * i + 1 )] += real_c( 0.5 ) * src[edgedof::macroedge::index( sourceLevel, i )];
   }
}

void N1E1toN1E1Prolongation::prolongateMacroFace( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   using sD = stencilDirection;

   for ( auto it : edgedof::macroface::Iterator( sourceLevel ) )
   {
      const idx_t row = it.y();
      const idx_t col = it.x();

      // horizontal edges
      real_t src_val = src[edgedof::macroface::horizontalIndex( sourceLevel, col, row )];

      dst[edgedof::macroface::horizontalIndex( sourceLevel + 1, 2 * col, 2 * row )] += real_c( 0.5 ) * src_val;
      dst[edgedof::macroface::horizontalIndex( sourceLevel + 1, 2 * col + 1, 2 * row )] += real_c( 0.5 ) * src_val;

      dst[edgedof::macroface::horizontalIndex( sourceLevel + 1, 2 * col, 2 * row + 1 )] += real_c( 0.25 ) * src_val;
      dst[edgedof::macroface::indexFromHorizontalEdge( sourceLevel + 1, 2 * col + 1, 2 * row, sD::EDGE_VE_NW )] += real_c( 0.25 ) * src_val;
      dst[edgedof::macroface::indexFromHorizontalEdge( sourceLevel + 1, 2 * col, 2 * row, sD::EDGE_DI_N )] += -real_c( 0.25 ) * src_val;

      if ( row > 0 )
      {
         dst[edgedof::macroface::horizontalIndex( sourceLevel + 1, 2 * col + 1, 2 * row - 1 )] += real_c( 0.25 ) * src_val;
         dst[edgedof::macroface::indexFromHorizontalEdge( sourceLevel + 1, 2 * col, 2 * row, sD::EDGE_VE_SE )] += real_c( 0.25 ) * src_val;
         dst[edgedof::macroface::indexFromHorizontalEdge( sourceLevel + 1, 2 * col + 1, 2 * row, sD::EDGE_DI_S )] +=
             -real_c( 0.25 ) * src_val;
      }

      // vertical edges
      src_val = src[edgedof::macroface::verticalIndex( sourceLevel, col, row )];

      dst[edgedof::macroface::verticalIndex( sourceLevel + 1, 2 * col, 2 * row )] += real_c( 0.5 ) * src_val;
      dst[edgedof::macroface::verticalIndex( sourceLevel + 1, 2 * col, 2 * row + 1 )] += real_c( 0.5 ) * src_val;

      dst[edgedof::macroface::verticalIndex( sourceLevel + 1, 2 * col + 1, 2 * row )] += real_c( 0.25 ) * src_val;
      dst[edgedof::macroface::indexFromVerticalEdge( sourceLevel + 1, 2 * col, 2 * row + 1, sD::EDGE_HO_SE )] += real_c( 0.25 ) * src_val;
      dst[edgedof::macroface::indexFromVerticalEdge( sourceLevel + 1, 2 * col, 2 * row, sD::EDGE_DI_E )] += real_c( 0.25 ) * src_val;

      if ( col > 0 )
      {
         dst[edgedof::macroface::verticalIndex( sourceLevel + 1, 2 * col - 1, 2 * row + 1 )] += real_c( 0.25 ) * src_val;
         dst[edgedof::macroface::indexFromVerticalEdge( sourceLevel + 1, 2 * col, 2 * row, sD::EDGE_HO_NW )] += real_c( 0.25 ) * src_val;
         dst[edgedof::macroface::indexFromVerticalEdge( sourceLevel + 1, 2 * col, 2 * row + 1, sD::EDGE_DI_W )] += real_c( 0.25 ) * src_val;
      }

      // diagonal edges
      src_val = src[edgedof::macroface::diagonalIndex( sourceLevel, col, row )];

      dst[edgedof::macroface::diagonalIndex( sourceLevel + 1, 2 * col, 2 * row + 1 )] += real_c( 0.5 ) * src_val;
      dst[edgedof::macroface::diagonalIndex( sourceLevel + 1, 2 * col + 1, 2 * row )] += real_c( 0.5 ) * src_val;

      dst[edgedof::macroface::diagonalIndex( sourceLevel + 1, 2 * col, 2 * row )] += real_c( 0.25 ) * src_val;
      dst[edgedof::macroface::indexFromDiagonalEdge( sourceLevel + 1, 2 * col, 2 * row, sD::EDGE_HO_N )] += -real_c( 0.25 ) * src_val;
      dst[edgedof::macroface::indexFromDiagonalEdge( sourceLevel + 1, 2 * col, 2 * row, sD::EDGE_VE_E )] += real_c( 0.25 ) * src_val;

      if ( row + col < idx_t( levelinfo::num_microedges_per_edge( sourceLevel ) - 1 ) )
      {
         dst[edgedof::macroface::diagonalIndex( sourceLevel + 1, 2 * col + 1, 2 * row + 1 )] += real_c( 0.25 ) * src_val;
         dst[edgedof::macroface::indexFromDiagonalEdge( sourceLevel + 1, 2 * col + 1, 2 * row, sD::EDGE_HO_N )] +=
             -real_c( 0.25 ) * src_val;
         dst[edgedof::macroface::indexFromDiagonalEdge( sourceLevel + 1, 2 * col, 2 * row + 1, sD::EDGE_VE_E )] += real_c( 0.25 ) * src_val;
      }
   }
}

void N1E1toN1E1Prolongation::prolongateMacroCell( const real_t* src, real_t* dst, const uint_t& srcLvl ) const
{
   using edgedof::macrocell::xIndex;
   using edgedof::macrocell::xyIndex;
   using edgedof::macrocell::xyzIndex;
   using edgedof::macrocell::xzIndex;
   using edgedof::macrocell::yIndex;
   using edgedof::macrocell::yzIndex;
   using edgedof::macrocell::zIndex;

   const uint_t dstLvl = srcLvl + 1;

   const auto isValid = []( const uint_t level, const idx_t x, const idx_t y, const idx_t z ) {
      return x >= 0 && y >= 0 && z >= 0 && ( x + y + z < idx_t( levelinfo::num_microedges_per_edge( level ) ) );
   };
   const auto isValidXYZ = []( const uint_t level, const idx_t x, const idx_t y, const idx_t z ) {
      return x >= 0 && y >= 0 && z >= 0 && ( x + y + z < idx_t( levelinfo::num_microedges_per_edge( level ) ) - 1 );
   };

   for ( auto it : edgedof::macrocell::Iterator( srcLvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      // clang-format off
      if ( isValid   (dstLvl, 2*x -1, 2*y +0, 2*z +1) ) {
         dst[  xIndex(dstLvl, 2*x -1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → (1, 4, 5), (2, 6) → (2, 7, 6)
         dst[ xyIndex(dstLvl, 2*x -1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → (1, 3, 5)
         dst[ xzIndex(dstLvl, 2*x -1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → (1, 4, 5), (2, 6) → (2, 7, 6)
         dst[  yIndex(dstLvl, 2*x -1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → green up (1, 3, 4, 5)
         dst[  zIndex(dstLvl, 2*x -1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → (1, 4, 5), (2, 6) → (2, 7, 6)
      }
      if ( isValid   (dstLvl, 2*x -1, 2*y +1, 2*z +0) ) {
         dst[  xIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (1, 2) → (1, 3, 2), (5, 6) → (5, 7, 6)
         dst[ xyIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (1, 2) → (1, 3, 2), (5, 6) → (5, 7, 6)
         dst[ xzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.0000 ) * src[  yIndex(srcLvl, x, y, z)]; // (1, 2) → blue up (1, 3, 2, 5)
         dst[  yIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (1, 2) → (1, 3, 2), (5, 6) → (5, 7, 6)
         dst[ xzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → blue up (1, 3, 2, 5)
         dst[ yzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (3, 2, 5)
         dst[ xzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → blue up (1, 3, 2, 5)
         dst[  zIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → (1, 3, 5)
      }
      if ( isValid   (dstLvl, 2*x -1, 2*y +1, 2*z +1) ) {
         dst[  xIndex(dstLvl, 2*x -1, 2*y +1, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (3, 2, 5)
         dst[ xyIndex(dstLvl, 2*x -1, 2*y +1, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (2, 5, 7)
         dst[ xzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (2, 5, 7)
         dst[  yIndex(dstLvl, 2*x -1, 2*y +1, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → green down (3, 2, 5, 7)
         dst[ yzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (2, 5, 7)
      }
      if ( isValid   (dstLvl, 2*x +0, 2*y -1, 2*z +1) ) {
         dst[ xyIndex(dstLvl, 2*x +0, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (3, 7) → (3, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (3, 7) → blue down (3, 4, 5, 7)
         dst[  yIndex(dstLvl, 2*x +0, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (2, 6) → (2, 5, 6), (3, 7) → (3, 4, 7)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (2, 6) → (2, 5, 6), (3, 7) → (3, 4, 7)
         dst[  zIndex(dstLvl, 2*x +0, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (2, 6) → (2, 5, 6), (3, 7) → (3, 4, 7)
      }
      if ( isValid   (dstLvl, 2*x +0, 2*y +0, 2*z +0) ) {
         dst[  xIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.5000 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1), (3, 2), (4, 5), (7, 6)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1) → (0, 1, 3), (4, 5) → (4, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1) → (0, 1, 4), (3, 2) → (3, 2, 7)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (0, 1, 3), (5, 7) → (4, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (0, 1, 4), (2, 7) → (3, 2, 7)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3) → (0, 1, 3), (4, 7) → (4, 5, 7)
         dst[  yIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.5000 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3), (1, 2), (4, 7), (5, 6)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] += -real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3) → (0, 3, 4), (1, 2) → (1, 2, 5)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (1, 2, 5), (3, 4) → (0, 3, 4)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4) → (0, 1, 4), (3, 7) → (3, 2, 7)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4) → (0, 3, 4), (1, 5) → (1, 2, 5)
         dst[  zIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.5000 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4), (1, 5), (2, 6), (3, 7)
      }
      if ( isValid   (dstLvl, 2*x +0, 2*y +0, 2*z +1) ) {
         dst[  xIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1) → (0, 1, 4), (3, 2) → (3, 2, 7)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 4)
         dst[  xIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] += -real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (0, 1, 4), (2, 7) → (3, 2, 7)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (1, 3, 4)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.5000 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4), (2, 7)
         dst[  yIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3) → (0, 3, 4), (1, 2) → (1, 2, 5)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → (1, 3, 4)
         dst[  yIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (1, 2, 5), (3, 4) → (0, 3, 4)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.5000 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5), (3, 4)
         dst[  xIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4) → (0, 1, 4), (3, 7) → (3, 2, 7)
         dst[  yIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4) → (0, 3, 4), (1, 5) → (1, 2, 5)
         dst[  zIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.5000 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4), (1, 5), (2, 6), (3, 7)
      }
      if ( isValid   (dstLvl, 2*x +0, 2*y +1, 2*z -1) ) {
         dst[  xIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → (3, 4, 5)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → blue down (3, 4, 5, 7)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → (3, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → blue down (3, 4, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] +=  real_c( 0.0000 ) * src[  yIndex(srcLvl, x, y, z)]; // (4, 7) → blue down (3, 4, 5, 7)
         dst[  yIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (4, 7) → (3, 4, 7), (5, 6) → (2, 5, 6)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] += -real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (4, 7) → (3, 4, 7), (5, 6) → (2, 5, 6)
         dst[  zIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (4, 7) → (3, 4, 7), (5, 6) → (2, 5, 6)
      }
      if ( isValid   (dstLvl, 2*x +0, 2*y +1, 2*z +0) ) {
         dst[  xIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1) → (0, 1, 3), (4, 5) → (4, 5, 7)
         dst[  xIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (0, 1, 3), (5, 7) → (4, 5, 7)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.5000 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3), (5, 7)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 4)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (1, 3, 4)
         dst[  xIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3) → (0, 1, 3), (4, 7) → (4, 5, 7)
         dst[  yIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.5000 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3), (1, 2), (4, 7), (5, 6)
         dst[  zIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3) → (0, 3, 4), (1, 2) → (1, 2, 5)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → (1, 3, 4)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.5000 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5), (3, 4)
         dst[  zIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (1, 2, 5), (3, 4) → (0, 3, 4)
         dst[  zIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4) → (0, 3, 4), (1, 5) → (1, 2, 5)
      }
      if ( isValid   (dstLvl, 2*x +0, 2*y +1, 2*z +1) ) {
         dst[  xIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → (3, 4, 5)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → blue down (3, 4, 5, 7)
         dst[  yIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (2, 5, 6), (3, 4) → (3, 4, 7)
         dst[ yzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (2, 5, 6), (3, 4) → (3, 4, 7)
         dst[  zIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (2, 5, 6), (3, 4) → (3, 4, 7)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y -1, 2*z +0) ) {
         dst[  xIndex(dstLvl, 2*x +1, 2*y -1, 2*z +0)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → (1, 3, 2), (7, 6) → (5, 7, 6)
         dst[ xyIndex(dstLvl, 2*x +1, 2*y -1, 2*z +0)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → (1, 3, 2), (7, 6) → (5, 7, 6)
         dst[ xzIndex(dstLvl, 2*x +1, 2*y -1, 2*z +0)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → blue up (1, 3, 2, 5)
         dst[  yIndex(dstLvl, 2*x +1, 2*y -1, 2*z +0)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → (1, 3, 2), (7, 6) → (5, 7, 6)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y -1, 2*z +0)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → (3, 2, 5)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y -1, 2*z +1) ) {
         dst[  xIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → (3, 2, 5)
         dst[  yIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → green down (3, 2, 5, 7)
         dst[ xyIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (2, 7) → (2, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (2, 7) → (2, 5, 7)
         dst[  yIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.0000 ) * src[ xzIndex(srcLvl, x, y, z)]; // (2, 7) → green down (3, 2, 5, 7)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (2, 7) → (2, 5, 7)
         dst[  yIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (3, 7) → green down (3, 2, 5, 7)
         dst[  zIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (3, 7) → (3, 5, 7)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +0, 2*z -1) ) {
         dst[  xIndex(dstLvl, 2*x +1, 2*y +0, 2*z -1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → (1, 4, 5), (7, 6) → (2, 7, 6)
         dst[ xzIndex(dstLvl, 2*x +1, 2*y +0, 2*z -1)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → (1, 4, 5), (7, 6) → (2, 7, 6)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z -1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → green up (1, 3, 4, 5)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +0, 2*z -1)] += -real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → (3, 4, 5)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +0, 2*z -1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → (1, 4, 5), (7, 6) → (2, 7, 6)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +0, 2*z +0) ) {
         dst[  xIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.5000 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1), (3, 2), (4, 5), (7, 6)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1) → (0, 1, 3), (4, 5) → (4, 5, 7)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1) → (0, 1, 4), (3, 2) → (3, 2, 7)
         dst[ xyIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.5000 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3), (5, 7)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (0, 1, 3), (5, 7) → (4, 5, 7)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 4)
         dst[ xzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.5000 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4), (2, 7)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (1, 3, 4)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (0, 1, 4), (2, 7) → (3, 2, 7)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3) → (0, 1, 3), (4, 7) → (4, 5, 7)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → (1, 3, 4)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4) → (0, 1, 4), (3, 7) → (3, 2, 7)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +0, 2*z +1) ) {
         dst[ xyIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 5)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → green up (1, 3, 4, 5)
         dst[  xIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] += -real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (1, 4, 5), (2, 7) → (2, 7, 6)
         dst[ xzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (1, 4, 5), (2, 7) → (2, 7, 6)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.0000 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → green up (1, 3, 4, 5)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → (1, 4, 5), (2, 7) → (2, 7, 6)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] += -real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → green up (1, 3, 4, 5)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → (3, 4, 5)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +1, 2*z -1) ) {
         dst[ xyIndex(dstLvl, 2*x +1, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → (2, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +1, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → (2, 5, 7)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → green down (3, 2, 5, 7)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +1, 2*z -1)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → (2, 5, 7)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +1, 2*z -1)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → (3, 5, 7)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +1, 2*z +0) ) {
         dst[  xIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 2), (5, 7) → (5, 7, 6)
         dst[ xyIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 2), (5, 7) → (5, 7, 6)
         dst[ xzIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → blue up (1, 3, 2, 5)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 2), (5, 7) → (5, 7, 6)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x -1, 2*y -1, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x -1, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (2, 6) → white down (2, 5, 7, 6)
      }
      if ( isValidXYZ(dstLvl, 2*x -1, 2*y +0, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x -1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (1, 5) → (1, 3, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x -1, 2*y +1, 2*z -1) ) {
         dst[xyzIndex(dstLvl, 2*x -1, 2*y +1, 2*z -1)] +=  real_c( 0.0000 ) * src[  yIndex(srcLvl, x, y, z)]; // (5, 6) → white down (2, 5, 7, 6)
      }
      if ( isValidXYZ(dstLvl, 2*x -1, 2*y +1, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → (3, 2, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x -1, 2*y +1, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x -1, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (2, 5) → white down (2, 5, 7, 6)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y -1, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y -1, 2*z +1)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (3, 7) → (3, 5, 7)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y +0, 2*z -1) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z -1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (4, 5) → (3, 4, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y +0, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (0, 1) → white up (0, 1, 3, 4)
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → white up (0, 1, 3, 4)
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.0000 ) * src[ xzIndex(srcLvl, x, y, z)]; // (1, 4) → white up (0, 1, 3, 4)
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.0000 ) * src[  yIndex(srcLvl, x, y, z)]; // (0, 3) → white up (0, 1, 3, 4)
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → white up (0, 1, 3, 4)
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[  zIndex(srcLvl, x, y, z)]; // (0, 4) → white up (0, 1, 3, 4)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y +0, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[ yzIndex(srcLvl, x, y, z)]; // (3, 4) → (3, 4, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y +1, 2*z -1) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +1, 2*z -1)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → (3, 5, 7)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y -1, 2*z -1) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y -1, 2*z -1)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (7, 6) → white down (2, 5, 7, 6)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y -1, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y -1, 2*z +0)] +=  real_c( 0.2500 ) * src[  xIndex(srcLvl, x, y, z)]; // (3, 2) → (3, 2, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y -1, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y -1, 2*z +1)] +=  real_c( 0.0000 ) * src[ xzIndex(srcLvl, x, y, z)]; // (2, 7) → white down (2, 5, 7, 6)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y +0, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (1, 3) → (1, 3, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y +1, 2*z -1) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y +1, 2*z -1)] += -real_c( 0.2500 ) * src[ xyIndex(srcLvl, x, y, z)]; // (5, 7) → white down (2, 5, 7, 6)
      }
      // clang-format on
   }

   for ( auto it : edgedof::macrocell::IteratorXYZ( srcLvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      // clang-format off
      if ( isValid   (dstLvl, 2*x +0, 2*y +1, 2*z +1) ) {
         dst[  xIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 4, 5)
         dst[ xyIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] += -real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 5, 7)
         dst[ xzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] +=  real_c( 0.0000 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → blue down (3, 4, 5, 7)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +0, 2*z +1) ) {
         dst[ xyIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] += -real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (1, 3, 5)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.0000 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → green up (1, 3, 4, 5)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 4, 5)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +1, 2*z +0) ) {
         dst[ xzIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.0000 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → blue up (1, 3, 2, 5)
         dst[ yzIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 2, 5)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (1, 3, 5)
      }
      if ( isValid   (dstLvl, 2*x +1, 2*y +1, 2*z +1) ) {
         dst[  xIndex(dstLvl, 2*x +1, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 2, 5)
         dst[  yIndex(dstLvl, 2*x +1, 2*y +1, 2*z +1)] +=  real_c( 0.0000 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → green down (3, 2, 5, 7)
         dst[  zIndex(dstLvl, 2*x +1, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 5, 7)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y +0, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +0, 2*z +1)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 4, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y +1, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +0)] +=  real_c( 0.5000 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +0, 2*y +1, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x +0, 2*y +1, 2*z +1)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 5, 7)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y +0, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +0)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (1, 3, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y +0, 2*z +1) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y +0, 2*z +1)] +=  real_c( 0.5000 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5)
      }
      if ( isValidXYZ(dstLvl, 2*x +1, 2*y +1, 2*z +0) ) {
         dst[xyzIndex(dstLvl, 2*x +1, 2*y +1, 2*z +0)] +=  real_c( 0.2500 ) * src[xyzIndex(srcLvl, x, y, z)]; // (3, 5) → (3, 2, 5)
      }
      // clang-format on
   }
}

} // namespace n1e1
} // namespace hyteg
