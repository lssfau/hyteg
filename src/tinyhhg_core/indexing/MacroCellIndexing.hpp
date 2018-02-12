
#pragma once

namespace hhg {
namespace indexing {

using walberla::uint_t;

/// Array access layouts - wrapped by general access function.
/// The general memory layout can thus be switched globally by setting
/// the called functions in the general size and access function.
namespace layout {

/// Required memory for the linear macro cell layout
template< uint_t width >
inline constexpr uint_t linearMacroCellSize()
{
  return ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
}

/// General linear memory layout indexing function for macro cells
template< uint_t width >
inline constexpr uint_t linearMacroCellIndex( const uint_t & x, const uint_t & y, const uint_t & z )
{
  const uint_t widthMinusSlice = width - z;
  const uint_t sliceOffset = linearMacroCellSize< width >() - ( ( widthMinusSlice + 2 ) * ( widthMinusSlice + 1 ) * widthMinusSlice ) / 6;
  const uint_t rowOffset   = y * ( widthMinusSlice + 1 ) - ( ( ( y + 1 ) * ( y ) ) / 2 );
  return sliceOffset + rowOffset + x;
}

} // namespace layout


template< uint_t width >
inline constexpr uint_t macroCellSize()
{
  return layout::linearMacroCellSize< width >();
}

template< uint_t width >
inline constexpr uint_t macroCellIndex( const uint_t & x, const uint_t & y, const uint_t & z )
{
  return layout::linearMacroCellIndex< width >( x, y, z );
}


/// Iterator over a cell.
/// Iterates in x-direction first, then increases in y-direction, then in z-direction.
/// To be used as follows:
///
/// \code{.cpp}
/// for ( const auto & it : CellIterator( 5, 0 ) )
/// {
///   WALBERLA_LOG_INFO_ON_ROOT( "col = " << it.col() << ", row = " << it.row() << ", depth = " << it.dep() );
/// }
/// \endcode
///
/// \param width width of an edge of the cell
/// \param offsetToCenter if > 0 the iterator iterates only over the inner points of the cell (shifted by offsetToCenter)
/// \param end if true creates the end iterator (not needed in general since the iterator itself implements the end() method)
class CellIterator
{
public:

  using iterator_category = std::input_iterator_tag;
  using value_type        = Index;
  using reference         = value_type const&;
  using pointer           = value_type const*;
  using difference_type   = ptrdiff_t;

  CellIterator( const uint_t & width, const uint_t & offsetToCenter = 0, const bool & end = false ) :
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

  CellIterator begin() { return CellIterator( width_, offsetToCenter_ ); }
  CellIterator end()   { return CellIterator( width_, offsetToCenter_, true ); }

  bool operator==( const CellIterator & other ) const
  {
    return other.step_ == step_;
  }

  bool operator!=( const CellIterator & other ) const
  {
    return other.step_ != step_;
  }

  reference operator*()  const { return  coordinates_; };
  pointer   operator->() const { return &coordinates_; };

  CellIterator & operator++() // prefix
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

  CellIterator operator++( int ) // postfix
  {
    const CellIterator tmp( *this );
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

} // namespace indexing
} // namespace hhg
