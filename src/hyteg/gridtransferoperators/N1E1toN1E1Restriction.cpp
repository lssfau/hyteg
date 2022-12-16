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

#include "hyteg/gridtransferoperators/N1E1toN1E1Restriction.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/memory/FunctionMemory.hpp"

namespace hyteg {
namespace n1e1 {

void N1E1toN1E1Restriction::restrictAdditively( const N1E1VectorFunction< real_t >& function,
                                                const uint_t&                       sourceLevel,
                                                const DoFType&                      flag,
                                                const UpdateType&                   updateType ) const
{
   WALBERLA_ASSERT( function.getStorage()->hasGlobalCells(), "N1E1 restriction only implemented in 3D." )

   const uint_t destinationLevel = sourceLevel - 1;

   for ( const auto& it : function.getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         auto data = edge.getData( function.getDoFs()->getEdgeDataID() );

         setToZeroIfNeededMacroEdge( destinationLevel, edge, function.getDoFs()->getEdgeDataID(), updateType );
         restrictMacroEdge( data->getPointer( sourceLevel ), data->getPointer( destinationLevel ), sourceLevel );
      }
   }

   for ( const auto& it : function.getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         auto data = face.getData( function.getDoFs()->getFaceDataID() );

         setToZeroIfNeededMacroFace( destinationLevel, face, function.getDoFs()->getFaceDataID(), updateType );
         restrictMacroFace( data->getPointer( sourceLevel ), data->getPointer( destinationLevel ), sourceLevel );
      }
   }

   for ( const auto& it : function.getStorage()->getCells() )
   {
      const Cell& cell = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         auto data = cell.getData( function.getDoFs()->getCellDataID() );

         setToZeroIfNeededMacroCell( destinationLevel, cell, function.getDoFs()->getCellDataID(), updateType );
         restrictMacroCell( data->getPointer( sourceLevel ), data->getPointer( destinationLevel ), sourceLevel );
      }
   }

   // Face → Edge comm must come before Cell → Face comm
   // otherwise the Cell → Edge contributions are propagated twice, via the two faces adjacent to the cell and edge
   function.communicateAdditively< Face, Edge >( destinationLevel, false );
   function.communicateAdditively< Cell, Face >( destinationLevel, false );
   function.communicateAdditively< Cell, Edge >( destinationLevel, false );
}

void N1E1toN1E1Restriction::setToZeroIfNeededMacroEdge( const uint_t&                                            level,
                                                        const Edge&                                              edge,
                                                        const PrimitiveDataID< FunctionMemory< real_t >, Edge >& edgeDataID,
                                                        const UpdateType& updateType ) const
{
   if ( updateType == Replace )
   {
      edge.getData( edgeDataID )->setToZero( level );
   }
   // updateType == Add: edges have no boundary ⇒ nothing to be done
}

void N1E1toN1E1Restriction::setToZeroIfNeededMacroFace( const uint_t&                                            level,
                                                        const Face&                                              face,
                                                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& faceDataID,
                                                        const UpdateType& updateType ) const
{
   switch ( updateType )
   {
   case Replace:
      face.getData( faceDataID )->setToZero( level );
      break;
   case Add:
      edgedof::macroface::setBoundaryToZero( level, face, faceDataID );
      break;
   }
}

void N1E1toN1E1Restriction::setToZeroIfNeededMacroCell( const uint_t&                                            level,
                                                        const Cell&                                              cell,
                                                        const PrimitiveDataID< FunctionMemory< real_t >, Cell >& cellDataID,
                                                        const UpdateType& updateType ) const
{
   switch ( updateType )
   {
   case Replace:
      cell.getData( cellDataID )->setToZero( level );
      break;
   case Add:
      edgedof::macrocell::setBoundaryToZero( level, cell, cellDataID );
      break;
   }
}

void N1E1toN1E1Restriction::restrictMacroEdge( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   const uint_t destinationLevel = sourceLevel - 1;
   const uint_t n                = levelinfo::num_microedges_per_edge( destinationLevel );

   for ( idx_t i = 0; i < idx_t( n ); ++i )
   {
      dst[edgedof::macroedge::index( destinationLevel, i )] += 0.5 * src[edgedof::macroedge::index( sourceLevel, 2 * i )];
      dst[edgedof::macroedge::index( destinationLevel, i )] += 0.5 * src[edgedof::macroedge::index( sourceLevel, 2 * i + 1 )];
   }
}

void N1E1toN1E1Restriction::restrictMacroFace( const real_t* src, real_t* dst, const uint_t& srcLvl ) const
{
   using edgedof::macroface::diagonalIndex;
   using edgedof::macroface::horizontalIndex;
   using edgedof::macroface::indexFromDiagonalEdge;
   using edgedof::macroface::indexFromHorizontalEdge;
   using edgedof::macroface::indexFromVerticalEdge;
   using edgedof::macroface::verticalIndex;
   using sD = stencilDirection;

   const uint_t dstLvl = srcLvl - 1;

   for ( auto it : edgedof::macroface::Iterator( dstLvl ) )
   {
      const idx_t row = it.row();
      const idx_t col = it.col();

      // horizontal edges
      dst[horizontalIndex( dstLvl, col, row )] += 0.25 * src[horizontalIndex( srcLvl, 2 * col, 2 * row + 1 )];
      dst[horizontalIndex( dstLvl, col, row )] +=
          0.25 * src[indexFromHorizontalEdge( srcLvl, 2 * col + 1, 2 * row, sD::EDGE_VE_NW )];
      dst[horizontalIndex( dstLvl, col, row )] += -0.25 * src[indexFromHorizontalEdge( srcLvl, 2 * col, 2 * row, sD::EDGE_DI_N )];

      if ( row > 0 )
      {
         dst[horizontalIndex( dstLvl, col, row )] += 0.5 * src[horizontalIndex( srcLvl, 2 * col, 2 * row )];
         dst[horizontalIndex( dstLvl, col, row )] += 0.5 * src[horizontalIndex( srcLvl, 2 * col + 1, 2 * row )];

         dst[horizontalIndex( dstLvl, col, row )] += 0.25 * src[horizontalIndex( srcLvl, 2 * col + 1, 2 * row - 1 )];
         dst[horizontalIndex( dstLvl, col, row )] +=
             0.25 * src[indexFromHorizontalEdge( srcLvl, 2 * col, 2 * row, sD::EDGE_VE_SE )];
         dst[horizontalIndex( dstLvl, col, row )] +=
             -0.25 * src[indexFromHorizontalEdge( srcLvl, 2 * col + 1, 2 * row, sD::EDGE_DI_S )];
      }

      // vertical edges
      dst[verticalIndex( dstLvl, col, row )] += 0.25 * src[verticalIndex( srcLvl, 2 * col + 1, 2 * row )];
      dst[verticalIndex( dstLvl, col, row )] += 0.25 * src[indexFromVerticalEdge( srcLvl, 2 * col, 2 * row + 1, sD::EDGE_HO_SE )];
      dst[verticalIndex( dstLvl, col, row )] += 0.25 * src[indexFromVerticalEdge( srcLvl, 2 * col, 2 * row, sD::EDGE_DI_E )];

      if ( col > 0 )
      {
         dst[verticalIndex( dstLvl, col, row )] += 0.5 * src[verticalIndex( srcLvl, 2 * col, 2 * row )];
         dst[verticalIndex( dstLvl, col, row )] += 0.5 * src[verticalIndex( srcLvl, 2 * col, 2 * row + 1 )];

         dst[verticalIndex( dstLvl, col, row )] += 0.25 * src[verticalIndex( srcLvl, 2 * col - 1, 2 * row + 1 )];
         dst[verticalIndex( dstLvl, col, row )] += 0.25 * src[indexFromVerticalEdge( srcLvl, 2 * col, 2 * row, sD::EDGE_HO_NW )];
         dst[verticalIndex( dstLvl, col, row )] +=
             0.25 * src[indexFromVerticalEdge( srcLvl, 2 * col, 2 * row + 1, sD::EDGE_DI_W )];
      }

      // diagonal edges
      dst[diagonalIndex( dstLvl, col, row )] += 0.25 * src[diagonalIndex( srcLvl, 2 * col, 2 * row )];
      dst[diagonalIndex( dstLvl, col, row )] += -0.25 * src[indexFromDiagonalEdge( srcLvl, 2 * col, 2 * row, sD::EDGE_HO_N )];
      dst[diagonalIndex( dstLvl, col, row )] += 0.25 * src[indexFromDiagonalEdge( srcLvl, 2 * col, 2 * row, sD::EDGE_VE_E )];

      if ( row + col < idx_t( levelinfo::num_microedges_per_edge( dstLvl ) - 1 ) )
      {
         dst[diagonalIndex( dstLvl, col, row )] += 0.5 * src[diagonalIndex( srcLvl, 2 * col, 2 * row + 1 )];
         dst[diagonalIndex( dstLvl, col, row )] += 0.5 * src[diagonalIndex( srcLvl, 2 * col + 1, 2 * row )];

         dst[diagonalIndex( dstLvl, col, row )] += 0.25 * src[diagonalIndex( srcLvl, 2 * col + 1, 2 * row + 1 )];
         dst[diagonalIndex( dstLvl, col, row )] +=
             -0.25 * src[indexFromDiagonalEdge( srcLvl, 2 * col + 1, 2 * row, sD::EDGE_HO_N )];
         dst[diagonalIndex( dstLvl, col, row )] +=
             0.25 * src[indexFromDiagonalEdge( srcLvl, 2 * col, 2 * row + 1, sD::EDGE_VE_E )];
      }
   }
}

void N1E1toN1E1Restriction::restrictMacroCell( const real_t* src, real_t* dst, const uint_t& srcLvl ) const
{
   using edgedof::EdgeDoFOrientation;
   using edgedof::macrocell::isInnerEdgeDoF;
   using edgedof::macrocell::xIndex;
   using edgedof::macrocell::xyIndex;
   using edgedof::macrocell::xyzIndex;
   using edgedof::macrocell::xzIndex;
   using edgedof::macrocell::yIndex;
   using edgedof::macrocell::yzIndex;
   using edgedof::macrocell::zIndex;

   const uint_t dstLvl = srcLvl - 1;

   const auto isValid =
       []( const uint_t level, const idx_t x, const idx_t y, const idx_t z, const EdgeDoFOrientation& orientation ) {
          const idx_t maxSum =
              idx_t( levelinfo::num_microedges_per_edge( level ) ) - ( orientation == EdgeDoFOrientation::XYZ ? 1 : 0 );
          return x >= 0 && y >= 0 && z >= 0 && ( x + y + z < maxSum ) &&
                 isInnerEdgeDoF( level, indexing::Index( x, y, z ), orientation );
       };

   for ( auto it : edgedof::macrocell::Iterator( dstLvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      // clang-format off
      // 158 stencil entries
      if ( isValid(srcLvl, 2*x -1, 2*y -1, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x -1, 2*y -1, 2*z +1)]; // (2, 6) → white down (2, 5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +0, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x -1, 2*y +0, 2*z +0)]; // (1, 5) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +0, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x -1, 2*y +0, 2*z +1)]; // (1, 5) → (1, 4, 5), (2, 6) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +0, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x -1, 2*y +0, 2*z +1)]; // (1, 5) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +0, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x -1, 2*y +0, 2*z +1)]; // (1, 5) → (1, 4, 5), (2, 6) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +0, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x -1, 2*y +0, 2*z +1)]; // (1, 5) → green up (1, 3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +0, 2*z +1, EdgeDoFOrientation::Z) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x -1, 2*y +0, 2*z +1)]; // (1, 5) → (1, 4, 5), (2, 6) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z -1, EdgeDoFOrientation::XYZ) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.0000 * src[xyzIndex(srcLvl, 2*x -1, 2*y +1, 2*z -1)]; // (5, 6) → white down (2, 5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +0, EdgeDoFOrientation::X) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (1, 2) → (1, 3, 2), (5, 6) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +0, EdgeDoFOrientation::XY) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (1, 2) → (1, 3, 2), (5, 6) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (2, 5) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +0, EdgeDoFOrientation::XZ) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.0000 * src[ xzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (1, 2) → blue up (1, 3, 2, 5)
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (2, 5) → blue up (1, 3, 2, 5)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (1, 5) → blue up (1, 3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +0, EdgeDoFOrientation::Y) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (1, 2) → (1, 3, 2), (5, 6) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +0, EdgeDoFOrientation::YZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (2, 5) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +0, EdgeDoFOrientation::Z) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x -1, 2*y +1, 2*z +0)]; // (1, 5) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[  xIndex(srcLvl, 2*x -1, 2*y +1, 2*z +1)]; // (2, 5) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[ xyIndex(srcLvl, 2*x -1, 2*y +1, 2*z +1)]; // (2, 5) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +1)]; // (2, 5) → white down (2, 5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +1)]; // (2, 5) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[  yIndex(srcLvl, 2*x -1, 2*y +1, 2*z +1)]; // (2, 5) → green down (3, 2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x -1, 2*y +1, 2*z +1, EdgeDoFOrientation::YZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x -1, 2*y +1, 2*z +1)]; // (2, 5) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y -1, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y -1, 2*z +1)]; // (3, 7) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y -1, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y -1, 2*z +1)]; // (3, 7) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y -1, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y -1, 2*z +1)]; // (3, 7) → blue down (3, 4, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y -1, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +0, 2*y -1, 2*z +1)]; // (2, 6) → (2, 5, 6), (3, 7) → (3, 4, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y -1, 2*z +1, EdgeDoFOrientation::YZ) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +0, 2*y -1, 2*z +1)]; // (2, 6) → (2, 5, 6), (3, 7) → (3, 4, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y -1, 2*z +1, EdgeDoFOrientation::Z) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +0, 2*y -1, 2*z +1)]; // (2, 6) → (2, 5, 6), (3, 7) → (3, 4, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z -1, EdgeDoFOrientation::XYZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z -1)]; // (4, 5) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +0, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.5000 * src[  xIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 1), (3, 2), (4, 5), (7, 6)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +0, EdgeDoFOrientation::XY) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 1) → (0, 1, 3), (4, 5) → (4, 5, 7)
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (1, 3) → (0, 1, 3), (5, 7) → (4, 5, 7)
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 3) → (0, 1, 3), (4, 7) → (4, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 1) → white up (0, 1, 3, 4)
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (1, 3) → white up (0, 1, 3, 4)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.0000 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (1, 4) → white up (0, 1, 3, 4)
         dst[  yIndex(dstLvl, x, y, z)] +=  0.0000 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 3) → white up (0, 1, 3, 4)
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (3, 4) → white up (0, 1, 3, 4)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 4) → white up (0, 1, 3, 4)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +0, EdgeDoFOrientation::XZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 1) → (0, 1, 4), (3, 2) → (3, 2, 7)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (1, 4) → (0, 1, 4), (2, 7) → (3, 2, 7)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 4) → (0, 1, 4), (3, 7) → (3, 2, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +0, EdgeDoFOrientation::Y) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.5000 * src[  yIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 3), (1, 2), (4, 7), (5, 6)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +0, EdgeDoFOrientation::YZ) ) {
         dst[  yIndex(dstLvl, x, y, z)] += -0.2500 * src[ yzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 3) → (0, 3, 4), (1, 2) → (1, 2, 5)
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (2, 5) → (1, 2, 5), (3, 4) → (0, 3, 4)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 4) → (0, 3, 4), (1, 5) → (1, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +0, EdgeDoFOrientation::Z) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.5000 * src[  zIndex(srcLvl, 2*x +0, 2*y +0, 2*z +0)]; // (0, 4), (1, 5), (2, 6), (3, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (0, 1) → (0, 1, 4), (3, 2) → (3, 2, 7)
         dst[ xzIndex(dstLvl, x, y, z)] += -0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (1, 4) → (0, 1, 4), (2, 7) → (3, 2, 7)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (0, 4) → (0, 1, 4), (3, 7) → (3, 2, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (1, 3) → (1, 3, 4)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (1, 4) → (1, 3, 4)
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (3, 4) → (1, 3, 4)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (3, 4) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.5000 * src[ xzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (1, 4), (2, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (0, 3) → (0, 3, 4), (1, 2) → (1, 2, 5)
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[  yIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (2, 5) → (1, 2, 5), (3, 4) → (0, 3, 4)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (0, 4) → (0, 3, 4), (1, 5) → (1, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::YZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.5000 * src[ yzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (2, 5), (3, 4)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::Z) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.5000 * src[  zIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (0, 4), (1, 5), (2, 6), (3, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z -1, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (4, 5) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z -1, EdgeDoFOrientation::XY) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (5, 7) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z -1, EdgeDoFOrientation::XYZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (5, 7) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z -1, EdgeDoFOrientation::XZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (4, 5) → blue down (3, 4, 5, 7)
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (5, 7) → blue down (3, 4, 5, 7)
         dst[  yIndex(dstLvl, x, y, z)] +=  0.0000 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (4, 7) → blue down (3, 4, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z -1, EdgeDoFOrientation::Y) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (4, 7) → (3, 4, 7), (5, 6) → (2, 5, 6)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z -1, EdgeDoFOrientation::YZ) ) {
         dst[  yIndex(dstLvl, x, y, z)] += -0.2500 * src[ yzIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (4, 7) → (3, 4, 7), (5, 6) → (2, 5, 6)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z -1, EdgeDoFOrientation::Z) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +0, 2*y +1, 2*z -1)]; // (4, 7) → (3, 4, 7), (5, 6) → (2, 5, 6)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +0, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (0, 1) → (0, 1, 3), (4, 5) → (4, 5, 7)
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (1, 3) → (0, 1, 3), (5, 7) → (4, 5, 7)
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (0, 3) → (0, 1, 3), (4, 7) → (4, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +0, EdgeDoFOrientation::XY) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.5000 * src[ xyIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (1, 3), (5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +0, EdgeDoFOrientation::XZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (1, 3) → (1, 3, 4)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (1, 4) → (1, 3, 4)
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (3, 4) → (1, 3, 4)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +0, EdgeDoFOrientation::Y) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.5000 * src[  yIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (0, 3), (1, 2), (4, 7), (5, 6)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +0, EdgeDoFOrientation::YZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.5000 * src[ yzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (2, 5), (3, 4)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +0, EdgeDoFOrientation::Z) ) {
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (0, 3) → (0, 3, 4), (1, 2) → (1, 2, 5)
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (2, 5) → (1, 2, 5), (3, 4) → (0, 3, 4)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (0, 4) → (0, 3, 4), (1, 5) → (1, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (3, 4) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (3, 4) → blue down (3, 4, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[  yIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (2, 5) → (2, 5, 6), (3, 4) → (3, 4, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::YZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (2, 5) → (2, 5, 6), (3, 4) → (3, 4, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::Z) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (2, 5) → (2, 5, 6), (3, 4) → (3, 4, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z -1, EdgeDoFOrientation::XYZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +1, 2*y -1, 2*z -1)]; // (7, 6) → white down (2, 5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +0, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +1, 2*y -1, 2*z +0)]; // (3, 2) → (1, 3, 2), (7, 6) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +0, EdgeDoFOrientation::XY) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ xyIndex(srcLvl, 2*x +1, 2*y -1, 2*z +0)]; // (3, 2) → (1, 3, 2), (7, 6) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +1, 2*y -1, 2*z +0)]; // (3, 2) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +0, EdgeDoFOrientation::XZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ xzIndex(srcLvl, 2*x +1, 2*y -1, 2*z +0)]; // (3, 2) → blue up (1, 3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +0, EdgeDoFOrientation::Y) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y -1, 2*z +0)]; // (3, 2) → (1, 3, 2), (7, 6) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +0, EdgeDoFOrientation::YZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y -1, 2*z +0)]; // (3, 2) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (3, 2) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (2, 7) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.0000 * src[xyzIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (2, 7) → white down (2, 5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (2, 7) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (3, 2) → green down (3, 2, 5, 7)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.0000 * src[  yIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (2, 7) → green down (3, 2, 5, 7)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (3, 7) → green down (3, 2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +1, EdgeDoFOrientation::YZ) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (2, 7) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y -1, 2*z +1, EdgeDoFOrientation::Z) ) {
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y -1, 2*z +1)]; // (3, 7) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z -1, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +1, 2*y +0, 2*z -1)]; // (4, 5) → (1, 4, 5), (7, 6) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z -1, EdgeDoFOrientation::XZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ xzIndex(srcLvl, 2*x +1, 2*y +0, 2*z -1)]; // (4, 5) → (1, 4, 5), (7, 6) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z -1, EdgeDoFOrientation::Y) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z -1)]; // (4, 5) → green up (1, 3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z -1, EdgeDoFOrientation::YZ) ) {
         dst[  xIndex(dstLvl, x, y, z)] += -0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +0, 2*z -1)]; // (4, 5) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z -1, EdgeDoFOrientation::Z) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +0, 2*z -1)]; // (4, 5) → (1, 4, 5), (7, 6) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::X) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.5000 * src[  xIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (0, 1), (3, 2), (4, 5), (7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::XY) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.5000 * src[ xyIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (1, 3), (5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[xyzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (1, 3) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::XZ) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.5000 * src[ xzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (1, 4), (2, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::Y) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (0, 1) → (0, 1, 3), (4, 5) → (4, 5, 7)
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (1, 3) → (0, 1, 3), (5, 7) → (4, 5, 7)
         dst[  yIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (0, 3) → (0, 1, 3), (4, 7) → (4, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::YZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (1, 3) → (1, 3, 4)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (1, 4) → (1, 3, 4)
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (3, 4) → (1, 3, 4)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::Z) ) {
         dst[  xIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (0, 1) → (0, 1, 4), (3, 2) → (3, 2, 7)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (1, 4) → (0, 1, 4), (2, 7) → (3, 2, 7)
         dst[  zIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (0, 4) → (0, 1, 4), (3, 7) → (3, 2, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[ xzIndex(dstLvl, x, y, z)] += -0.2500 * src[  xIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (1, 4) → (1, 4, 5), (2, 7) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (1, 3) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (1, 4) → (1, 4, 5), (2, 7) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (1, 3) → green up (1, 3, 4, 5)
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.0000 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (1, 4) → green up (1, 3, 4, 5)
         dst[ yzIndex(dstLvl, x, y, z)] += -0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (3, 4) → green up (1, 3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::YZ) ) {
         dst[ yzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (3, 4) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::Z) ) {
         dst[ xzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (1, 4) → (1, 4, 5), (2, 7) → (2, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z -1, EdgeDoFOrientation::XY) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +1, 2*y +1, 2*z -1)]; // (5, 7) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z -1, EdgeDoFOrientation::XYZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[xyzIndex(srcLvl, 2*x +1, 2*y +1, 2*z -1)]; // (5, 7) → white down (2, 5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z -1, EdgeDoFOrientation::XZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +1, 2*y +1, 2*z -1)]; // (5, 7) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z -1, EdgeDoFOrientation::Y) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +1, 2*z -1)]; // (5, 7) → green down (3, 2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z -1, EdgeDoFOrientation::YZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +1, 2*z -1)]; // (5, 7) → (2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z -1, EdgeDoFOrientation::Z) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +1, 2*z -1)]; // (5, 7) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::X) ) {
         dst[ xyIndex(dstLvl, x, y, z)] += -0.2500 * src[  xIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (1, 3) → (1, 3, 2), (5, 7) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::XY) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xyIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (1, 3) → (1, 3, 2), (5, 7) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::XZ) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[ xzIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (1, 3) → blue up (1, 3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::Y) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[  yIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (1, 3) → (1, 3, 2), (5, 7) → (5, 7, 6)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::Z) ) {
         dst[ xyIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (1, 3) → (1, 3, 5)
      }
      // clang-format on
   }

   for ( auto it : edgedof::macrocell::IteratorXYZ( dstLvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      // clang-format off
      if ( isValid(srcLvl, 2*x +0, 2*y +0, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +0, 2*z +1)]; // (3, 5) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.5000 * src[xyzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +0)]; // (3, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (3, 5) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[xyzIndex(dstLvl, x, y, z)] += -0.2500 * src[ xyIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (3, 5) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (3, 5) → (3, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +0, 2*y +1, 2*z +1, EdgeDoFOrientation::XZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.0000 * src[ xzIndex(srcLvl, 2*x +0, 2*y +1, 2*z +1)]; // (3, 5) → blue down (3, 4, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +0)]; // (3, 5) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::XY) ) {
         dst[xyzIndex(dstLvl, x, y, z)] += -0.2500 * src[ xyIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (3, 5) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::XYZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.5000 * src[xyzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (3, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.0000 * src[  yIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (3, 5) → green up (1, 3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +0, 2*z +1, EdgeDoFOrientation::YZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +0, 2*z +1)]; // (3, 5) → (3, 4, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::XYZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[xyzIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (3, 5) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::XZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.0000 * src[ xzIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (3, 5) → blue up (1, 3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::YZ) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[ yzIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (3, 5) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +0, EdgeDoFOrientation::Z) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +1, 2*z +0)]; // (3, 5) → (1, 3, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +1, EdgeDoFOrientation::X) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  xIndex(srcLvl, 2*x +1, 2*y +1, 2*z +1)]; // (3, 5) → (3, 2, 5)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +1, EdgeDoFOrientation::Y) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.0000 * src[  yIndex(srcLvl, 2*x +1, 2*y +1, 2*z +1)]; // (3, 5) → green down (3, 2, 5, 7)
      }
      if ( isValid(srcLvl, 2*x +1, 2*y +1, 2*z +1, EdgeDoFOrientation::Z) ) {
         dst[xyzIndex(dstLvl, x, y, z)] +=  0.2500 * src[  zIndex(srcLvl, 2*x +1, 2*y +1, 2*z +1)]; // (3, 5) → (3, 5, 7)
      }
      // clang-format on
   }
}

} // namespace n1e1
} // namespace hyteg
