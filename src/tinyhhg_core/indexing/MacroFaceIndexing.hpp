
#pragma once

#include "core/debug/Debug.h"
#include "core/Abort.h"
#include "tinyhhg_core/indexing/Common.hpp"

namespace hhg {
namespace indexing {

using walberla::uint_t;

/// Array access layouts - wrapped by general access function.
/// The general memory layout can thus be switched globally by setting
/// the called functions in the general size and access function.
namespace layout {

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

}


template< uint_t width >
inline constexpr uint_t macroFaceSize()
{
  return layout::linearMacroFaceSize< width >();
}

template< uint_t width >
inline constexpr uint_t macroFaceIndex( const uint_t & col, const uint_t & row )
{
  return layout::linearMacroFaceIndex< width >( col, row );
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

  FaceIterator( const uint_t & width, const uint_t & offsetToCenter = 0, const bool & end = false ) :
    width_( width ), offsetToCenter_( offsetToCenter ),
    totalNumberOfDoFs_( ( ( width - 3 * offsetToCenter + 1 ) * ( width - 3 * offsetToCenter ) ) / 2 ), step_( 0 )
  {
    WALBERLA_ASSERT_LESS( offsetToCenter, width, "Offset to center is beyond face width!" );

    coordinates_.dep() = 0;

    coordinates_.col() = offsetToCenter;
    coordinates_.row() = offsetToCenter;

    if ( end )
    {
      step_ = totalNumberOfDoFs_;
    }
  }

  FaceIterator begin() { return FaceIterator( width_, offsetToCenter_ ); }
  FaceIterator end()   { return FaceIterator( width_, offsetToCenter_, true ); }

  bool operator==( const FaceIterator & other ) const
  {
    return other.step_ == step_;
  }

  bool operator!=( const FaceIterator & other ) const
  {
    return other.step_ != step_;
  }

  reference operator*()  const { return  coordinates_; };
  pointer   operator->() const { return &coordinates_; };

  FaceIterator & operator++() // prefix
  {
    WALBERLA_ASSERT_LESS_EQUAL( step_, totalNumberOfDoFs_, "Incrementing iterator beyond end!" );

    step_++;

    const uint_t currentRow = coordinates_.row();
    const uint_t currentCol = coordinates_.col();

    const uint_t lengthOfCurrentRowWithoutOffset = width_ - currentRow;

    if ( currentCol < lengthOfCurrentRowWithoutOffset - offsetToCenter_ - 1 )
    {
      coordinates_.col()++;
    }
    else
    {
      coordinates_.row()++;
      coordinates_.col() = offsetToCenter_;
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

  const uint_t    width_;
  const uint_t    offsetToCenter_;
  const uint_t    totalNumberOfDoFs_;
        uint_t    step_;
        Index     coordinates_;

};


enum class FaceBorderDirection
{
  BOTTOM_LEFT_TO_RIGHT,
  BOTTOM_RIGHT_TO_LEFT,
  LEFT_BOTTOM_TO_TOP,
  LEFT_TOP_TO_BOTTOM,
  DIAGONAL_BOTTOM_TO_TOP,
  DIAGONAL_TOP_TO_BOTTOM,
};

/// return
/// \param localEdgeId local Id of the edge on the face
/// \param orientation orientation of the edge; 1 for same as face; -1 for opposing
inline FaceBorderDirection getFaceBorderDirection(uint_t localEdgeId, int orientation){
  if(localEdgeId == 0) {
    if (orientation == 1) {
      return FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT;
    } else if(orientation == -1)
      return FaceBorderDirection::BOTTOM_RIGHT_TO_LEFT;
  } else if(localEdgeId == 1){
    if (orientation == 1) {
      return FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP;
    } else if(orientation == -1)
      return FaceBorderDirection::DIAGONAL_TOP_TO_BOTTOM;
  } else if(localEdgeId == 2){
    if (orientation == 1) {
      return FaceBorderDirection::LEFT_TOP_TO_BOTTOM;
    } else if(orientation == -1)
      return FaceBorderDirection::LEFT_BOTTOM_TO_TOP;
  } else {
    WALBERLA_ABORT("wrong EdgeId or orientation");
  }
  return FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT;
}


/// Iterator over the borders of a face.
/// Decoupled from the indexing function, it returns the logical coordinates
/// which can then be inserted into the respective indexing function.
/// Since it implements all necessary methods, you can do:
///
///   for ( const auto & it : FaceBorderIterator( 9, FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP ) )
///   {
///     WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it.col() << ", row = " << it.row() );
///   }
///
/// \param width width of one edge of the face
/// \param direction FaceBorderDirection indicating the face border of interest
/// \param offsetToCenter if > 0, the iterator iterates parallel to the specified border, shifted to the center of the face by offsetToCenter
/// \param offsetFromVertices the iterator skips the first and last offsetFromVertices points of the border
class FaceBorderIterator
{
public:

  using iterator_category = std::input_iterator_tag;
  using value_type        = Index;
  using reference         = value_type const&;
  using pointer           = value_type const*;
  using difference_type   = ptrdiff_t;

  FaceBorderIterator( const uint_t & width, const FaceBorderDirection & direction,
                      const uint_t & offsetToCenter = 0,
                      const uint_t & offsetFromVertices = 0, const bool & end = false ) :
    width_( width ), direction_( direction ), offsetToCenter_( offsetToCenter ),
    offsetFromVertices_( offsetFromVertices ), step_( 0 )
  {
    WALBERLA_ASSERT_LESS( offsetToCenter, width, "Offset to center is beyond face width!" );

    coordinates_.dep() = 0;

    switch( direction )
    {
    case FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT:
      coordinates_.col() = 0 + offsetFromVertices;
      coordinates_.row() = 0 + offsetToCenter;
      break;
    case FaceBorderDirection::BOTTOM_RIGHT_TO_LEFT:
      coordinates_.col() = width - 1 - offsetFromVertices;
      coordinates_.row() = 0 + offsetToCenter;
      break;
    case FaceBorderDirection::LEFT_BOTTOM_TO_TOP:
      coordinates_.col() = 0 + offsetToCenter;
      coordinates_.row() = 0 + offsetFromVertices;
      break;
    case FaceBorderDirection::LEFT_TOP_TO_BOTTOM:
      coordinates_.col() = 0 + offsetToCenter;
      coordinates_.row() = width - 1 - offsetFromVertices;
      break;
    case FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP:
      coordinates_.col() = width - 1 - offsetToCenter - offsetFromVertices;
      coordinates_.row() = 0 + offsetFromVertices;
      break;
    case FaceBorderDirection::DIAGONAL_TOP_TO_BOTTOM:
      coordinates_.col() = 0 + offsetFromVertices;
      coordinates_.row() = width - 1 - offsetToCenter - offsetFromVertices;
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

  FaceBorderIterator begin() { return FaceBorderIterator( width_, direction_, offsetToCenter_, offsetFromVertices_ ); }
  FaceBorderIterator end()   { return FaceBorderIterator( width_, direction_, offsetToCenter_, offsetFromVertices_, true ); }

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
    WALBERLA_ASSERT_LESS_EQUAL( step_, width_, "Incrementing iterator beyond end!" );

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

  const uint_t              width_;
  const FaceBorderDirection direction_;
  const uint_t              offsetToCenter_;
  const uint_t              offsetFromVertices_;
        uint_t              step_;
        Index               coordinates_;

};


template< uint_t width >
inline Index getFaceBottomLeftCorner() { Index idx; idx.col() = 0; idx.row() = 0; return idx; }

template< uint_t width >
inline Index getFaceBottomRightCorner() { Index idx; idx.col() = width - 1; idx.row() = 0; return idx; }

template< uint_t width >
inline Index getFaceTopLeftCorner() { Index idx; idx.col() = 0; idx.row() = width - 1; return idx; }

}
}
