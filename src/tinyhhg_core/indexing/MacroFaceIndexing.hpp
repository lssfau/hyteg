#include <iostream>

#include "core/debug/Debug.h"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {
namespace indexing {

typedef unsigned uint_t;

/// Templated constexpr evaluating number of DoFs per row for different types (not sure if levelinfo:: stuff works here :/)

template< uint_t level >
constexpr uint_t levelToWidthP1 = ( 1u << level ) + 1u;

template< uint_t level >
constexpr uint_t levelToWidthEdgeDoF = ( 1u << level );


/// General linear memory layout indexing function for macro faces
/// For 2D we never have ghost layers
/// For 3D we have at most one or two, their layout is also triangular but with one less column and row
/// \param layerIndex 1 (= 0b001) for owned DoFs, 2 (= 0b010) for first, 4 (= 0b100) for second ghost layer
template< uint_t width >
constexpr uint_t linearMacroFaceIndex( const uint_t & col, const uint_t & row, const uint_t & layerIndex )
{
  const uint_t ghostLayerWidth       = width - 1;
  const uint_t numEntriesMiddleLayer = ( (           width + 1 ) * (           width ) ) / 2;
  const uint_t numEntriesGhostLayer  = ( ( ghostLayerWidth + 1 ) * ( ghostLayerWidth ) ) / 2;

  const uint_t rowWidth    = ( layerIndex & (1 << 0) ) * width                  // row width on middle layer
                           + ( layerIndex & (1 << 1) ) * ghostLayerWidth        // row width on first ghost layer
                           + ( layerIndex & (1 << 2) ) * ghostLayerWidth;       // row width on second ghost layer

  const uint_t rowOffset   = row * ( rowWidth + 1 ) - ( ( ( row + 1 ) * ( row ) ) / 2 );

  const uint_t layerOffset = ( layerIndex & (1 << 0) ) * 0                      // offset for middle layer
                           + ( layerIndex & (1 << 1) ) * numEntriesMiddleLayer  // offset for first ghost layer
                           + ( layerIndex & (1 << 2) ) * numEntriesGhostLayer;  // offset for second ghost layer

  return layerOffset + rowOffset + col;
}

enum class Direction
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
/// Returns PointND< int, 2 > objects indicating ( col, row )
/// Since it implements all necessary methods, you can do:
///
///   for ( const auto & it : P1FaceBorderIterator< 3 >( Direction::DIAGONAL_BOTTOM_TO_TOP ) )
///   {
///     WALBERLA_LOG_INFO_ON_ROOT( "FaceBorderIterator: col = " << it[0] << ", row = " << it[1] << " ( idx = " << P1FaceIndexFromVertex< 3 >( it[0], it[1], VERTEX_C ) << " ) " );
///   }
///
template< uint_t width >
class FaceBorderIterator
{
public:

  using iterator_category = std::input_iterator_tag;
  using value_type        = PointND< int, 2 >;
  using reference         = value_type const&;
  using pointer           = value_type const*;
  using difference_type   = ptrdiff_t;

  FaceBorderIterator( const Direction & direction, const uint_t & offsetToCenter = 0, const bool & end = false ) :
    direction_( direction ), offsetToCenter_( offsetToCenter ), step_( 0 )
  {
    WALBERLA_ASSERT_LESS( offsetToCenter, width, "Offset to center is beyond face width!" );

    switch( direction )
    {
    case Direction::BOTTOM_LEFT_TO_RIGHT:
      col() = 0;
      row() = 0 + offsetToCenter;
      break;
    case Direction::BOTTOM_RIGHT_TO_LEFT:
      col() = width - 1;
      row() = 0 + offsetToCenter;
      break;
    case Direction::LEFT_BOTTOM_TO_TOP:
      col() = 0 + offsetToCenter;
      row() = 0;
      break;
    case Direction::LEFT_TOP_TO_BOTTOM:
      col() = 0 + offsetToCenter;
      row() = width - 1;
      break;
    case Direction::DIAGONAL_BOTTOM_TO_TOP:
      col() = width - 1 - offsetToCenter;
      row() = 0;
      break;
    case Direction::DIAGONAL_TOP_TO_BOTTOM:
      col() = 0;
      row() = width - 1 - offsetToCenter;
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
    case Direction::BOTTOM_LEFT_TO_RIGHT:
      col()++;
      break;
    case Direction::BOTTOM_RIGHT_TO_LEFT:
      col()--;
      break;
    case Direction::LEFT_BOTTOM_TO_TOP:
      row()++;
      break;
    case Direction::LEFT_TOP_TO_BOTTOM:
      row()--;
      break;
    case Direction::DIAGONAL_BOTTOM_TO_TOP:
      col()--;
      row()++;
      break;
    case Direction::DIAGONAL_TOP_TO_BOTTOM:
      col()++;
      row()--;
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

  int & col() { return coordinates_[0]; }
  int & row() { return coordinates_[1]; }

  const Direction         direction_;
  const uint_t            offsetToCenter_;
        uint_t            step_;
        PointND< int, 2 > coordinates_;

};



/// Access functions
/// Alias for different DoF types, achieving minimum code duplication.
/// Also we could switch layouts e.g. via define!

template< uint_t level >
constexpr uint_t P1FaceMacroIndex( const uint_t & col, const uint_t & row )
{
  return linearMacroFaceIndex< levelToWidthP1< level > >( col, row, 1 );
};

template< uint_t level >
constexpr uint_t EdgeDoFFaceMacroIndex( const uint_t & col, const uint_t & row )
{
  return linearMacroFaceIndex< levelToWidthEdgeDoF< level > >( col, row, 1 );
};

// Iterators
template< uint_t level >
using P1FaceBorderIterator = FaceBorderIterator< levelToWidthP1< level > >;



/// Stencil functions - as usual (could be generated?)

enum stencilDirection
{
  VERTEX_C,
  VERTEX_E,

  EDGE_HO_C
};

template< uint_t level >
constexpr uint_t P1FaceIndexFromVertex( const uint_t & col, const uint_t & row,
                                        const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_C:
    return P1FaceMacroIndex< level >( col, row );
  case sD::VERTEX_E:
    return P1FaceMacroIndex< level >( col + 1, row );

    // ...

  default:
    return 0;
    break;
  }
}

template< uint_t level >
constexpr uint_t EdgeDoFFaceIndexFromVertex( const uint_t & col, const uint_t & row,
                                             const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::EDGE_HO_C:
    return EdgeDoFFaceMacroIndex< level >( col, row );

    // ...

  default:
    return 0;
    break;
  }
}

}
}
