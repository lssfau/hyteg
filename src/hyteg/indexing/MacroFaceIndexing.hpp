/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl.
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

#include "core/Abort.h"
#include "core/debug/Debug.h"

#include "hyteg/indexing/Common.hpp"

using walberla::uint_t;
using walberla::uint_c;

namespace hyteg {
namespace indexing {

/// Array access layouts - wrapped by general access function.
/// The general memory layout can thus be switched globally by setting
/// the called functions in the general size and access function.
namespace layout {

/// Required memory for the linear macro face layout
inline constexpr uint_t linearMacroFaceSize( const uint_t& width )
{
   return ( ( width + 1 ) * width ) / 2;
}

/// General linear memory layout indexing function for macro faces
inline constexpr uint_t linearMacroFaceIndex( const uint_t& width, const idx_t& x, const idx_t& y )
{
   const idx_t rowOffset = y * ( idx_t( width ) + 1 ) - ( ( ( y + 1 ) * ( y ) ) / 2 );
   return uint_t( rowOffset + x );
}

} // namespace layout

inline constexpr uint_t macroFaceSize( const uint_t& width )
{
   return layout::linearMacroFaceSize( width );
}

inline constexpr uint_t macroFaceIndex( const uint_t& width, const idx_t& x, const idx_t& y )
{
   return layout::linearMacroFaceIndex( width, x, y );
}

/// Iterator over a face.
/// Iterates from bottom to top in a row wise fashion.
/// It is possible to parameterize the iterator to only iterate over a inner part of the face.
/// This is done by setting the offset parameter to the distance to the edge.
/// If set to zero, the iterator iterates over the whole face (including the border).
class FaceIterator
{
 public:
   using iterator_category = std::input_iterator_tag;
   using value_type        = Index;
   using reference         = value_type const&;
   using pointer           = value_type const*;
   using difference_type   = ptrdiff_t;

   FaceIterator( const uint_t& width, const uint_t& offsetToCenter = 0, const bool& end = false )
   : width_( width )
   , offsetToCenter_( offsetToCenter )
   , totalNumberOfDoFs_( ( ( width - 3 * offsetToCenter + 1 ) * ( width - 3 * offsetToCenter ) ) / 2 )
   , step_( 0 )
   {
      WALBERLA_ASSERT_LESS_EQUAL( offsetToCenter, width, "Offset to center is beyond face width!" );

      coordinates_.z() = 0;

      coordinates_.x() = static_cast< idx_t >( offsetToCenter );
      coordinates_.y() = static_cast< idx_t >( offsetToCenter );

      if ( end )
      {
         step_ = totalNumberOfDoFs_;
      }
   }

   FaceIterator begin() { return FaceIterator( width_, offsetToCenter_ ); }
   FaceIterator end() { return FaceIterator( width_, offsetToCenter_, true ); }

   bool operator==( const FaceIterator& other ) const { return other.step_ == step_; }

   bool operator!=( const FaceIterator& other ) const { return other.step_ != step_; }

   reference operator*() const { return coordinates_; };
   pointer   operator->() const { return &coordinates_; };

   FaceIterator& operator++() // prefix
   {
      WALBERLA_ASSERT_LESS_EQUAL( step_, totalNumberOfDoFs_, "Incrementing iterator beyond end!" );

      step_++;

      const uint_t currentRow = uint_c( coordinates_.y() );
      const uint_t currentCol = uint_c( coordinates_.x() );

      const uint_t lengthOfCurrentRowWithoutOffset = width_ - currentRow;

      if ( currentCol < lengthOfCurrentRowWithoutOffset - offsetToCenter_ - 1 )
      {
         coordinates_.x()++;
      }
      else
      {
         coordinates_.y()++;
         coordinates_.x() = static_cast< idx_t >( offsetToCenter_ );
      }

      return *this;
   }

   FaceIterator operator++( int ) // postfix
   {
      const FaceIterator tmp( *this );
      ++*this;
      return tmp;
   }

 private:
   uint_t width_;
   uint_t offsetToCenter_;
   uint_t totalNumberOfDoFs_;
   uint_t step_;
   Index  coordinates_;
};

enum class FaceBoundaryDirection
{
   BOTTOM_LEFT_TO_RIGHT,
   BOTTOM_RIGHT_TO_LEFT,
   LEFT_BOTTOM_TO_TOP,
   LEFT_TOP_TO_BOTTOM,
   DIAGONAL_BOTTOM_TO_TOP,
   DIAGONAL_TOP_TO_BOTTOM,
};

/// \param localEdgeId local Id of the edge on the face
/// \param orientation orientation of the edge; 1 for same as face; -1 for opposing
inline FaceBoundaryDirection getFaceBoundaryDirection( uint_t localEdgeId, int orientation )
{
   if ( localEdgeId == 0 )
   {
      if ( orientation == 1 )
      {
         return FaceBoundaryDirection::BOTTOM_LEFT_TO_RIGHT;
      }
      else if ( orientation == -1 )
      {
         return FaceBoundaryDirection::BOTTOM_RIGHT_TO_LEFT;
      }
   }
   else if ( localEdgeId == 2 )
   {
      if ( orientation == 1 )
      {
         return FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP;
      }
      else if ( orientation == -1 )
      {
         return FaceBoundaryDirection::DIAGONAL_TOP_TO_BOTTOM;
      }
   }
   else if ( localEdgeId == 1 )
   {
      if ( orientation == 1 )
      {
         return FaceBoundaryDirection::LEFT_BOTTOM_TO_TOP;
      }
      else if ( orientation == -1 )
      {
         return FaceBoundaryDirection::LEFT_TOP_TO_BOTTOM;
      }
   }
   else
   {
      WALBERLA_ABORT( "wrong EdgeId or orientation" );
   }
   return FaceBoundaryDirection::BOTTOM_LEFT_TO_RIGHT;
}

inline FaceBoundaryDirection getFaceBoundaryDirectionFromLocalVertexPair( uint_t v0, uint_t v1 )
{
   if ( v0 == 0 && v1 == 1 )
      return FaceBoundaryDirection::BOTTOM_LEFT_TO_RIGHT;
   if ( v0 == 1 && v1 == 0 )
      return FaceBoundaryDirection::BOTTOM_RIGHT_TO_LEFT;
   if ( v0 == 0 && v1 == 2 )
      return FaceBoundaryDirection::LEFT_BOTTOM_TO_TOP;
   if ( v0 == 2 && v1 == 0 )
      return FaceBoundaryDirection::LEFT_TOP_TO_BOTTOM;
   if ( v0 == 1 && v1 == 2 )
      return FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP;
   if ( v0 == 2 && v1 == 1 )
      return FaceBoundaryDirection::DIAGONAL_TOP_TO_BOTTOM;

   WALBERLA_ABORT( "Invalid vertex pair." );
}

/// Iterator over the boundaries of a face.
/// Decoupled from the indexing function, it returns the logical coordinates
/// which can then be inserted into the respective indexing function.
/// Since it implements all necessary methods, you can do:
///
///   for ( const auto & it : FaceBorderIterator( 9, FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP ) )
///   {
///     WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it.x() << ", row = " << it.y() );
///   }
///
/// \param width width of one edge of the face
/// \param direction FaceBorderDirection indicating the face border of interest
/// \param offsetToCenter if > 0, the iterator iterates parallel to the specified border, shifted to the center of the face by offsetToCenter
/// \param offsetFromVertices the iterator skips the first and last offsetFromVertices points of the border
class FaceBoundaryIterator
{
 public:
   using iterator_category = std::input_iterator_tag;
   using value_type        = Index;
   using reference         = value_type const&;
   using pointer           = value_type const*;
   using difference_type   = ptrdiff_t;

   FaceBoundaryIterator( const uint_t& width,
                         const uint_t& vertex0,
                         const uint_t& vertex1,
                         const uint_t& offsetToCenter     = 0,
                         const uint_t& offsetFromVertices = 0,
                         const bool&   end                = false )
   : FaceBoundaryIterator( width,
                           getFaceBoundaryDirectionFromLocalVertexPair( vertex0, vertex1 ),
                           offsetToCenter,
                           offsetFromVertices,
                           end )
   {}

   FaceBoundaryIterator( const uint_t&                  width,
                         const std::array< uint_t, 2 >& vertices,
                         const uint_t&                  offsetToCenter     = 0,
                         const uint_t&                  offsetFromVertices = 0,
                         const bool&                    end                = false )
   : FaceBoundaryIterator( width,
                           getFaceBoundaryDirectionFromLocalVertexPair( vertices[0], vertices[1] ),
                           offsetToCenter,
                           offsetFromVertices,
                           end )
   {}

   FaceBoundaryIterator( const uint_t&                width,
                         const FaceBoundaryDirection& direction,
                         const uint_t&                offsetToCenter     = 0,
                         const uint_t&                offsetFromVertices = 0,
                         const bool&                  end                = false )
   : width_( width )
   , direction_( direction )
   , offsetToCenter_( offsetToCenter )
   , offsetFromVertices_( offsetFromVertices )
   , step_( 0 )
   {
      WALBERLA_ASSERT_GREATER( width, 0, "Size of face must be larger than zero!" );

      coordinates_.z() = 0;

      switch ( direction )
      {
      case FaceBoundaryDirection::BOTTOM_LEFT_TO_RIGHT:
         coordinates_.x() = static_cast< idx_t >( 0 + offsetFromVertices );
         coordinates_.y() = static_cast< idx_t >( 0 + offsetToCenter );
         break;
      case FaceBoundaryDirection::BOTTOM_RIGHT_TO_LEFT:
         coordinates_.x() = static_cast< idx_t >( width - 1 - offsetToCenter - offsetFromVertices );
         coordinates_.y() = static_cast< idx_t >( 0 + offsetToCenter );
         break;
      case FaceBoundaryDirection::LEFT_BOTTOM_TO_TOP:
         coordinates_.x() = static_cast< idx_t >( 0 + offsetToCenter );
         coordinates_.y() = static_cast< idx_t >( 0 + offsetFromVertices );
         break;
      case FaceBoundaryDirection::LEFT_TOP_TO_BOTTOM:
         coordinates_.x() = static_cast< idx_t >( 0 + offsetToCenter );
         coordinates_.y() = static_cast< idx_t >( width - 1 - offsetToCenter - offsetFromVertices );
         break;
      case FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP:
         coordinates_.x() = static_cast< idx_t >( width - 1 - offsetToCenter - offsetFromVertices );
         coordinates_.y() = static_cast< idx_t >( 0 + offsetFromVertices );
         break;
      case FaceBoundaryDirection::DIAGONAL_TOP_TO_BOTTOM:
         coordinates_.x() = static_cast< idx_t >( 0 + offsetFromVertices );
         coordinates_.y() = static_cast< idx_t >( width - 1 - offsetToCenter - offsetFromVertices );
         break;
      default:
         WALBERLA_ASSERT( false, "Invalid direction in face border iterator" );
         break;
      }

      if ( end )
      {
         step_ = width - offsetToCenter - 2 * offsetFromVertices;
      }
   }

   FaceBoundaryIterator begin() { return FaceBoundaryIterator( width_, direction_, offsetToCenter_, offsetFromVertices_ ); }
   FaceBoundaryIterator end() { return FaceBoundaryIterator( width_, direction_, offsetToCenter_, offsetFromVertices_, true ); }

   bool operator==( const FaceBoundaryIterator& other ) const
   {
      WALBERLA_ASSERT_EQUAL( direction_, other.direction_, "Comparing two iterators of different type!" );
      return other.step_ == step_;
   }

   bool operator!=( const FaceBoundaryIterator& other ) const
   {
      WALBERLA_ASSERT_EQUAL( direction_, other.direction_, "Comparing two iterators of different type!" );
      return other.step_ != step_;
   }

   reference operator*() const { return coordinates_; };
   pointer   operator->() const { return &coordinates_; };

   FaceBoundaryIterator& operator++() // prefix
   {
      WALBERLA_ASSERT_LESS_EQUAL( step_, width_, "Incrementing iterator beyond end!" );

      step_++;

      switch ( direction_ )
      {
      case FaceBoundaryDirection::BOTTOM_LEFT_TO_RIGHT:
         coordinates_.x()++;
         break;
      case FaceBoundaryDirection::BOTTOM_RIGHT_TO_LEFT:
         coordinates_.x()--;
         break;
      case FaceBoundaryDirection::LEFT_BOTTOM_TO_TOP:
         coordinates_.y()++;
         break;
      case FaceBoundaryDirection::LEFT_TOP_TO_BOTTOM:
         coordinates_.y()--;
         break;
      case FaceBoundaryDirection::DIAGONAL_BOTTOM_TO_TOP:
         coordinates_.x()--;
         coordinates_.y()++;
         break;
      case FaceBoundaryDirection::DIAGONAL_TOP_TO_BOTTOM:
         coordinates_.x()++;
         coordinates_.y()--;
         break;
      default:
         WALBERLA_ASSERT( false, "Invalid direction in face border iterator" );
         break;
      }

      return *this;
   }

   FaceBoundaryIterator operator++( int ) // postfix
   {
      const FaceBoundaryIterator tmp( *this );
      ++*this;
      return tmp;
   }

 private:
   const uint_t                width_;
   const FaceBoundaryDirection direction_;
   const uint_t                offsetToCenter_;
   const uint_t                offsetFromVertices_;
   uint_t                      step_;
   Index                       coordinates_;
};

template < uint_t width >
inline Index getFaceBottomLeftCorner()
{
   Index idx;
   idx.x() = 0;
   idx.y() = 0;
   return idx;
}

template < uint_t width >
inline Index getFaceBottomRightCorner()
{
   Index idx;
   idx.x() = width - 1;
   idx.y() = 0;
   return idx;
}

template < uint_t width >
inline Index getFaceTopLeftCorner()
{
   Index idx;
   idx.x() = 0;
   idx.y() = width - 1;
   return idx;
}

} // namespace indexing
} // namespace hyteg
