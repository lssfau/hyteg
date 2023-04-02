/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/indexing/MacroCellIndexing.hpp"

#include <cassert>

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"

namespace hyteg {
namespace indexing {

using walberla::int_c;
using walberla::uint_c;
using walberla::uint_t;

std::set< uint_t > isOnCellFace( const indexing::Index& index, const uint_t& width )
{
   std::set< uint_t > cellFaceIndices;
   const auto         dstIndex = toDistanceIndex( index, { { 0, 1, 2, 3 } }, width );
   const uint_t       maxDist  = width - 1;
   if ( dstIndex.d0() == maxDist )
      cellFaceIndices.insert( 3 ); // face with vertices 1, 2, 3
   if ( dstIndex.d1() == maxDist )
      cellFaceIndices.insert( 2 ); // face with vertices 0, 2, 3
   if ( dstIndex.d2() == maxDist )
      cellFaceIndices.insert( 1 ); // face with vertices 0, 1, 3
   if ( dstIndex.d3() == maxDist )
      cellFaceIndices.insert( 0 ); // face with vertices 0, 1, 2
   return cellFaceIndices;
}

std::set< uint_t > isOnCellEdge( const indexing::Index& index, const uint_t& width )
{
   std::set< uint_t > cellEdgeIndices;
   const auto         onFaces = isOnCellFace( index, width );
   if ( onFaces.size() <= 1 ) // index on edge <=> index on >= 2 faces
      return cellEdgeIndices;
   if ( onFaces.count( 0 ) == 1 && onFaces.count( 1 ) == 1 )
      cellEdgeIndices.insert( 0 );
   if ( onFaces.count( 0 ) == 1 && onFaces.count( 2 ) == 1 )
      cellEdgeIndices.insert( 1 );
   if ( onFaces.count( 0 ) == 1 && onFaces.count( 3 ) == 1 )
      cellEdgeIndices.insert( 2 );
   if ( onFaces.count( 1 ) == 1 && onFaces.count( 2 ) == 1 )
      cellEdgeIndices.insert( 3 );
   if ( onFaces.count( 1 ) == 1 && onFaces.count( 3 ) == 1 )
      cellEdgeIndices.insert( 4 );
   if ( onFaces.count( 2 ) == 1 && onFaces.count( 3 ) == 1 )
      cellEdgeIndices.insert( 5 );
   return cellEdgeIndices;
}

std::set< uint_t > isOnCellVertex( const indexing::Index& index, const uint_t& width )
{
   std::set< uint_t > cellVertexIndices;
   const auto         dstIndex = toDistanceIndex( index, { { 0, 1, 2, 3 } }, width );
   if ( dstIndex.d0() == 0 )
      cellVertexIndices.insert( 0 );
   if ( dstIndex.d1() == 0 )
      cellVertexIndices.insert( 1 );
   if ( dstIndex.d2() == 0 )
      cellVertexIndices.insert( 2 );
   if ( dstIndex.d3() == 0 )
      cellVertexIndices.insert( 3 );
   return cellVertexIndices;
}

CellIterator::CellIterator( const uint_t& width, const uint_t& offsetToCenter, const bool& end )
: width_( width )
, internalWidth_( width - ( offsetToCenter == 0 ? 0 : ( 2 + 2 * offsetToCenter ) ) )
, offsetToCenter_( offsetToCenter )
,
// Number of vertices in a tetrahedron with edge length n:
// T(n) = ( (n+2) * (n+1) * n ) / 6
// Number of _inner_ vertices of a tetrahedron with edge length n:
// T(n-4) = ( (n-2) * (n-3) * (n-4) ) / 6
totalNumberOfDoFs_( ( ( internalWidth_ + 2 ) * ( internalWidth_ + 1 ) * ( internalWidth_ ) ) / 6 )
, step_( 0 )
{
   WALBERLA_ASSERT_LESS_EQUAL( offsetToCenter, width, "Offset to center is beyond cell width!" );

   coordinates_.x() = idx_t( offsetToCenter_ );
   coordinates_.y() = idx_t( offsetToCenter_ );
   coordinates_.z() = idx_t( offsetToCenter_ );

   internalCoordinates_.x() = 0;
   internalCoordinates_.y() = 0;
   internalCoordinates_.z() = 0;

   if ( end )
   {
      step_ = totalNumberOfDoFs_;
   }
}

CellIterator& CellIterator::operator++() // prefix
{
   WALBERLA_ASSERT_LESS_EQUAL( step_, totalNumberOfDoFs_, "Incrementing iterator beyond end!" );

   step_++;

   const idx_t currentDep = internalCoordinates_.dep();
   const idx_t currentRow = internalCoordinates_.row();
   const idx_t currentCol = internalCoordinates_.col();

   const idx_t lengthOfCurrentRow   = idx_t( internalWidth_ ) - currentRow - currentDep;
   const idx_t heightOfCurrentSlice = idx_t( internalWidth_ ) - currentDep;

   if ( currentCol < lengthOfCurrentRow - 1 )
   {
      internalCoordinates_.col()++;
   }
   else if ( currentRow < heightOfCurrentSlice - 1 )
   {
      internalCoordinates_.row()++;
      internalCoordinates_.col() = 0;
   }
   else
   {
      internalCoordinates_.dep()++;
      internalCoordinates_.row() = 0;
      internalCoordinates_.col() = 0;
   }

   coordinates_ = internalCoordinates_ + IndexIncrement( (int) offsetToCenter_, (int) offsetToCenter_, (int) offsetToCenter_ );

   return *this;
}

CellIterator CellIterator::operator++( int ) // postfix
{
   const CellIterator tmp( *this );
   ++*this;
   return tmp;
}

CellBoundaryIterator::CellBoundaryIterator( const uint_t&                  width,
                                            const std::array< uint_t, 3 >& vertices,
                                            const uint_t&                  offsetToCenter,
                                            const bool&                    end )
: width_( width )
, vertices_( vertices )
, offsetToCenter_( offsetToCenter )
, totalNumberOfSteps_( levelinfo::num_microvertices_per_face_from_width( width - offsetToCenter ) )
, firstDirIncrement_( calculateIncrement( vertices_[0], vertices_[1] ) )
, secondDirIncrement_( calculateIncrement( vertices_[0], vertices_[2] ) )
, wrapAroundStep_( width - offsetToCenter )
, wrapArounds_( uint_c( 0 ) )
{
   WALBERLA_ASSERT_LESS_EQUAL( vertices_[0], 3 );
   WALBERLA_ASSERT_LESS_EQUAL( vertices_[1], 3 );
   WALBERLA_ASSERT_LESS_EQUAL( vertices_[2], 3 );

   WALBERLA_ASSERT_NOT_IDENTICAL( vertices_[0], vertices_[1] );
   WALBERLA_ASSERT_NOT_IDENTICAL( vertices_[1], vertices_[2] );
   WALBERLA_ASSERT_NOT_IDENTICAL( vertices_[2], vertices_[0] );

   if ( end )
   {
      step_ = totalNumberOfSteps_;
   }
   else
   {
      step_ = 0;
   }

   switch ( vertices_[0] )
   {
   case 0:
      coordinates_ = Index( 0, 0, 0 );
      break;
   case 1:
      coordinates_ = Index( idx_t( width_ - 1 ), 0, 0 );
      break;
   case 2:
      coordinates_ = Index( 0, idx_t( width_ - 1 ), 0 );
      break;
   case 3:
      coordinates_ = Index( 0, 0, idx_t( width_ - 1 ) );
      break;
   default:
      WALBERLA_ASSERT( false, "Invalid coordinates in CellBoundaryIterator!" );
      break;
   }

   switch ( tup4( vertices_[0], vertices_[1], vertices_[2] ) )
   {
   // Front face
   case tup4( 0, 1, 2 ):
   case tup4( 0, 2, 1 ):
      coordinates_ += IndexIncrement( 0, 0, (int) offsetToCenter );
      break;
   case tup4( 1, 0, 2 ):
   case tup4( 1, 2, 0 ):
      coordinates_ += IndexIncrement( -(int) offsetToCenter, 0, (int) offsetToCenter );
      break;
   case tup4( 2, 0, 1 ):
   case tup4( 2, 1, 0 ):
      coordinates_ += IndexIncrement( 0, -(int) offsetToCenter, (int) offsetToCenter );
      break;

   // Bottom face
   case tup4( 0, 1, 3 ):
   case tup4( 0, 3, 1 ):
      coordinates_ += IndexIncrement( 0, (int) offsetToCenter, 0 );
      break;
   case tup4( 1, 0, 3 ):
   case tup4( 1, 3, 0 ):
      coordinates_ += IndexIncrement( -(int) offsetToCenter, (int) offsetToCenter, 0 );
      break;
   case tup4( 3, 0, 1 ):
   case tup4( 3, 1, 0 ):
      coordinates_ += IndexIncrement( 0, (int) offsetToCenter, -(int) offsetToCenter );
      break;

   // Left face
   case tup4( 0, 2, 3 ):
   case tup4( 0, 3, 2 ):
      coordinates_ += IndexIncrement( (int) offsetToCenter, 0, 0 );
      break;
   case tup4( 2, 0, 3 ):
   case tup4( 2, 3, 0 ):
      coordinates_ += IndexIncrement( (int) offsetToCenter, -(int) offsetToCenter, 0 );
      break;
   case tup4( 3, 0, 2 ):
   case tup4( 3, 2, 0 ):
      coordinates_ += IndexIncrement( (int) offsetToCenter, 0, -(int) offsetToCenter );
      break;

   // Back/diagonal face
   case tup4( 3, 1, 2 ):
   case tup4( 3, 2, 1 ):
      coordinates_ += IndexIncrement( 0, 0, -(int) offsetToCenter );
      break;
   case tup4( 1, 3, 2 ):
   case tup4( 1, 2, 3 ):
      coordinates_ += IndexIncrement( -(int) offsetToCenter, 0, 0 );
      break;
   case tup4( 2, 3, 1 ):
   case tup4( 2, 1, 3 ):
      coordinates_ += IndexIncrement( 0, -(int) offsetToCenter, 0 );
      break;
   }

   wrapAroundCoordinates_ = coordinates_;
}

CellBoundaryIterator::CellBoundaryIterator( const uint_t& width,
                                            const uint_t& vertex0,
                                            const uint_t& vertex1,
                                            const uint_t& vertex2,
                                            const uint_t& offsetToCenter,
                                            const bool&   end )
: CellBoundaryIterator( width, { { vertex0, vertex1, vertex2 } }, offsetToCenter, end )
{}

CellBoundaryIterator& CellBoundaryIterator::operator++() // prefix
{
   WALBERLA_ASSERT_LESS_EQUAL( step_, totalNumberOfSteps_, "Incrementing iterator beyond end!" );

   step_++;

   if ( step_ == totalNumberOfSteps_ )
   {
      return *this;
   }

   if ( step_ == wrapAroundStep_ )
   {
      wrapAroundCoordinates_ += secondDirIncrement_;
      coordinates_ = wrapAroundCoordinates_;

      WALBERLA_ASSERT_GREATER_EQUAL( width_, wrapArounds_ + 1 );
      wrapAroundStep_ += ( width_ - offsetToCenter_ ) - ( wrapArounds_ + 1 );
      wrapArounds_++;
   }
   else
   {
      coordinates_ += firstDirIncrement_;
   }

   return *this;
}

CellBoundaryIterator CellBoundaryIterator::operator++( int ) // postfix
{
   const CellBoundaryIterator tmp( *this );
   ++*this;
   return tmp;
}

IndexIncrement CellBoundaryIterator::calculateIncrement( const uint_t& vertex0, const uint_t& vertex1 ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( vertex0, vertex1 );
   WALBERLA_ASSERT_LESS_EQUAL( vertex0, 3 );
   WALBERLA_ASSERT_LESS_EQUAL( vertex1, 3 );

   switch ( tup4( vertex0, vertex1 ) )
   {
   case tup4( 0, 1 ):
      return IndexIncrement( 1, 0, 0 );
   case tup4( 1, 0 ):
      return IndexIncrement( -1, 0, 0 );
   case tup4( 0, 2 ):
      return IndexIncrement( 0, 1, 0 );
   case tup4( 2, 0 ):
      return IndexIncrement( 0, -1, 0 );
   case tup4( 1, 2 ):
      return IndexIncrement( -1, 1, 0 );
   case tup4( 2, 1 ):
      return IndexIncrement( 1, -1, 0 );
   case tup4( 0, 3 ):
      return IndexIncrement( 0, 0, 1 );
   case tup4( 3, 0 ):
      return IndexIncrement( 0, 0, -1 );
   case tup4( 1, 3 ):
      return IndexIncrement( -1, 0, 1 );
   case tup4( 3, 1 ):
      return IndexIncrement( 1, 0, -1 );
   case tup4( 2, 3 ):
      return IndexIncrement( 0, -1, 1 );
   case tup4( 3, 2 ):
      return IndexIncrement( 0, 1, -1 );
   default:
      WALBERLA_ASSERT( false, "Invalid tuple in increment calculation!" );
      break;
   }
   return IndexIncrement( 0, 0, 0 );
}

} // namespace indexing
} // namespace hyteg
