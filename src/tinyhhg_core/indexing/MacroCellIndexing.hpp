
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

}


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

}
}
