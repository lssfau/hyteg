
#include "tinyhhg_core/indexing/MacroCellIndexing.hpp"

#include "tinyhhg_core/levelinfo.hpp"
#include "core/debug/Debug.h"
#include "core/DataTypes.h"

#include <cassert>

namespace hhg {
namespace indexing {

using walberla::uint_t;
using walberla::uint_c;

/// Helper function to create a 'tuple' from two integers taking less than 5 bits of space.
/// Allows for more readable switch statements in some cases.
static constexpr uint_t tup( const uint_t & a, const uint_t & b )
{
  assert( a < (1 << 8) );
  assert( b < (1 << 8) );
  return a << 4 | b;
}

CellIterator::CellIterator( const uint_t & width, const uint_t & offsetToCenter, const bool & end ) :
  width_( width ), offsetToCenter_( offsetToCenter ),
  // Number of vertices in a tetrahedron with edge length n:
  // T(n) = ( (n+2) * (n+1) * n ) / 6
  // Number of _inner_ vertices of a tetrahedron with edge length n:
  // T(n-4) = ( (n-2) * (n-3) * (n-4) ) / 6
  totalNumberOfDoFs_( ( ( width - 4 * offsetToCenter + 2 ) * ( width - 4 * offsetToCenter + 1 ) * ( width - 4 * offsetToCenter ) ) / 6 ), step_( 0 )
{
  WALBERLA_ASSERT_GREATER( width, 0, "Size of cell must be larger than zero!" );
  WALBERLA_ASSERT_LESS( offsetToCenter, width, "Offset to center is beyond cell width!" );

  coordinates_.col() = offsetToCenter;
  coordinates_.row() = offsetToCenter;
  coordinates_.dep() = offsetToCenter;

  if ( end )
  {
    step_ = totalNumberOfDoFs_;
  }
}

CellIterator & CellIterator::operator++() // prefix
{
  WALBERLA_ASSERT_LESS_EQUAL( step_, totalNumberOfDoFs_, "Incrementing iterator beyond end!" );

  step_++;

  const uint_t currentDep = coordinates_.dep();
  const uint_t currentRow = coordinates_.row();
  const uint_t currentCol = coordinates_.col();

  const uint_t lengthOfCurrentRowWithoutOffset   = width_ - currentRow - currentDep;
  const uint_t heightOfCurrentSliceWithoutOffset = width_ - currentDep;

  if ( currentCol < lengthOfCurrentRowWithoutOffset - offsetToCenter_ - 1 )
  {
    coordinates_.col()++;
  }
  else if ( currentRow < heightOfCurrentSliceWithoutOffset - offsetToCenter_ - 1 )
  {
    coordinates_.row()++;
    coordinates_.col() = offsetToCenter_;
  }
  else
  {
    coordinates_.dep()++;
    coordinates_.row() = offsetToCenter_;
    coordinates_.col() = offsetToCenter_;
  }

  return *this;
}

CellIterator CellIterator::operator++( int ) // postfix
{
  const CellIterator tmp( *this );
  ++*this;
  return tmp;
}

CellBorderIterator::CellBorderIterator( const uint_t & width, const std::array< uint_t, 3 > & vertices,
                                        const bool & end ) :
    width_( width ), vertices_( vertices ), totalNumberOfSteps_( levelinfo::num_microvertices_per_face_from_width( width ) ), 
    wrapAroundStep_( width ), wrapArounds_( uint_c( 0 ) )
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
    coordinates_ = Index( width_ - 1, 0, 0 );
    break;
  case 2:
    coordinates_ = Index( 0, width_ - 1, 0 );
    break;
  case 3:
    coordinates_ = Index( 0, 0, width_ - 1 );
    break;
  default:
    WALBERLA_ASSERT( false, "Invalid coordinates in CellBorderIterator!" );
    break;
  }

  wrapAroundCoordinates_ = coordinates_;
}


CellBorderIterator::CellBorderIterator( const uint_t & width, const uint_t & vertex0,
                                        const uint_t & vertex1, const uint_t & vertex2,
                                        const bool & end ) :
    CellBorderIterator( width, {{ vertex0, vertex1, vertex2 }}, end )
{}

CellBorderIterator & CellBorderIterator::operator++() // prefix
{
  WALBERLA_ASSERT_LESS_EQUAL( step_, totalNumberOfSteps_, "Incrementing iterator beyond end!" );

  step_++;

  if ( step_ == totalNumberOfSteps_ )
  {
    return *this;
  }

  if ( step_ == wrapAroundStep_ )
  {
    const IndexIncrement secondDirIncrement = calculateIncrement( vertices_[0], vertices_[2] );
    wrapAroundCoordinates_ += secondDirIncrement;
    coordinates_ = wrapAroundCoordinates_;

    WALBERLA_ASSERT_GREATER_EQUAL( width_, wrapArounds_ + 1 );
    wrapAroundStep_ += width_ - ( wrapArounds_ + 1 );
    wrapArounds_++;
  }
  else
  {
    const IndexIncrement firstDirIncrement = calculateIncrement( vertices_[0], vertices_[1] );
    coordinates_ += firstDirIncrement;
  }

  return *this;
}

CellBorderIterator CellBorderIterator::operator++( int ) // postfix
{
  const CellBorderIterator tmp( *this );
  ++*this;
  return tmp;
}

IndexIncrement CellBorderIterator::calculateIncrement( const uint_t & vertex0, const uint_t & vertex1 ) const
{
  WALBERLA_ASSERT_NOT_IDENTICAL( vertex0, vertex1 );
  WALBERLA_ASSERT_LESS_EQUAL( vertex0, 3 );
  WALBERLA_ASSERT_LESS_EQUAL( vertex1, 3 );

  switch ( tup( vertex0, vertex1 ) )
  {
  case tup( 0, 1 ):
      return IndexIncrement(  1,  0,  0 );
  case tup( 1, 0 ):
      return IndexIncrement( -1,  0,  0 );
  case tup( 0, 2 ):
      return IndexIncrement(  0,  1,  0 );
  case tup( 2, 0 ):
      return IndexIncrement(  0, -1,  0 );
  case tup( 1, 2 ):
      return IndexIncrement( -1,  1,  0 );
  case tup( 2, 1 ):
      return IndexIncrement(  1, -1,  0 );
  case tup( 0, 3 ):
      return IndexIncrement(  0,  0,  1 );
  case tup( 3, 0 ):
      return IndexIncrement(  0,  0, -1 );
  case tup( 1, 3 ):
      return IndexIncrement( -1,  0,  1 );
  case tup( 3, 1 ):
      return IndexIncrement(  1,  0, -1 );
  case tup( 2, 3 ):
      return IndexIncrement(  0, -1,  1 );
  case tup( 3, 2 ):
      return IndexIncrement(  0,  1, -1 );
  default:
      WALBERLA_ASSERT( false, "Invalid tuple in increment calculation!" );
      break;
  }
  return IndexIncrement( 0, 0, 0 );
}



}
}
