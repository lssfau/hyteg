
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/indexing/Common.hpp"

namespace hhg {
namespace indexing {

using walberla::uint_t;

/// Required memory for the linear macro face layout
template< uint_t width >
inline constexpr uint_t linearMacroFaceSize()
{
  return ( ( width + 1 ) * width ) / 2;
}

/// General linear memory layout indexing function for macro faces
template< uint_t width >
inline constexpr uint_t linearMacroFaceIndex( const uint_t & col, const uint_t & row )
{
  const uint_t rowOffset = row * ( width + 1 ) - ( ( ( row + 1 ) * ( row ) ) / 2 );
  return rowOffset + col;
}

enum class FaceBorderDirection
{
  BOTTOM_LEFT_TO_RIGHT,
  BOTTOM_RIGHT_TO_LEFT,
  LEFT_BOTTOM_TO_TOP,
  LEFT_TOP_TO_BOTTOM,
  DIAGONAL_BOTTOM_TO_TOP,
  DIAGONAL_TOP_TO_BOTTOM,
};

/// Iterator over the borders of a face.
/// Decoupled from the indexing function, it returns the logical coordinates
/// which can then be inserted into the respective indexing function.
/// Since it implements all necessary methods, you can do:
///
///   for ( const auto & it : FaceBorderIterator< 9 >( FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP ) )
///   {
///     WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it[0] << ", row = " << it[1] );
///   }
///
template< uint_t width >
class FaceBorderIterator
{
public:

  using iterator_category = std::input_iterator_tag;
  using value_type        = Index;
  using reference         = value_type const&;
  using pointer           = value_type const*;
  using difference_type   = ptrdiff_t;

  FaceBorderIterator( const FaceBorderDirection & direction, const uint_t & offsetToCenter = 0, const bool & end = false ) :
    direction_( direction ), offsetToCenter_( offsetToCenter ), step_( 0 )
  {
    WALBERLA_ASSERT_LESS( offsetToCenter, width, "Offset to center is beyond face width!" );

    coordinates_.dep() = 0;

    switch( direction )
    {
    case FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT:
      coordinates_.col() = 0;
      coordinates_.row() = 0 + offsetToCenter;
      break;
    case FaceBorderDirection::BOTTOM_RIGHT_TO_LEFT:
      coordinates_.col() = width - 1;
      coordinates_.row() = 0 + offsetToCenter;
      break;
    case FaceBorderDirection::LEFT_BOTTOM_TO_TOP:
      coordinates_.col() = 0 + offsetToCenter;
      coordinates_.row() = 0;
      break;
    case FaceBorderDirection::LEFT_TOP_TO_BOTTOM:
      coordinates_.col() = 0 + offsetToCenter;
      coordinates_.row() = width - 1;
      break;
    case FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP:
      coordinates_.col() = width - 1 - offsetToCenter;
      coordinates_.row() = 0;
      break;
    case FaceBorderDirection::DIAGONAL_TOP_TO_BOTTOM:
      coordinates_.col() = 0;
      coordinates_.row() = width - 1 - offsetToCenter;
      break;
    default:
      WALBERLA_ASSERT( false, "Invalid direction in face border iterator" );
      break;
    }

    if ( end )
    {
      step_ = width - offsetToCenter;
    }
  }

  FaceBorderIterator begin() { return FaceBorderIterator( direction_, offsetToCenter_ ); }
  FaceBorderIterator end()   { return FaceBorderIterator( direction_, offsetToCenter_, true ); }

  bool operator==( const FaceBorderIterator & other ) const
  {
    WALBERLA_ASSERT_EQUAL( direction_, other.direction_, "Comparing two iterators of different type!" );
    return other.step_ == step_;
  }

  bool operator!=( const FaceBorderIterator & other ) const
  {
    WALBERLA_ASSERT_EQUAL( direction_, other.direction_, "Comparing two iterators of different type!" );
    return other.step_ != step_;
  }

  reference operator*()  const { return  coordinates_; };
  pointer   operator->() const { return &coordinates_; };

  FaceBorderIterator & operator++() // prefix
  {
    WALBERLA_ASSERT_LESS_EQUAL( step_, width, "Incrementing iterator beyond end!" );

    step_++;

    switch( direction_ )
    {
    case FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT:
      coordinates_.col()++;
      break;
    case FaceBorderDirection::BOTTOM_RIGHT_TO_LEFT:
      coordinates_.col()--;
      break;
    case FaceBorderDirection::LEFT_BOTTOM_TO_TOP:
      coordinates_.row()++;
      break;
    case FaceBorderDirection::LEFT_TOP_TO_BOTTOM:
      coordinates_.row()--;
      break;
    case FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP:
      coordinates_.col()--;
      coordinates_.row()++;
      break;
    case FaceBorderDirection::DIAGONAL_TOP_TO_BOTTOM:
      coordinates_.col()++;
      coordinates_.row()--;
      break;
    default:
      WALBERLA_ASSERT( false, "Invalid direction in face border iterator" );
      break;
    }

    return *this;
  }

  FaceBorderIterator operator++( int ) // postfix
  {
    const FaceBorderIterator tmp( *this );
    ++*this;
    return tmp;
  }

private:

  const FaceBorderDirection direction_;
  const uint_t    offsetToCenter_;
        uint_t    step_;
        Index     coordinates_;

};

}
}
