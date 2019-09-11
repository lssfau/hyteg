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

#pragma once

#include "core/DataTypes.h"

namespace hyteg {
namespace indexing {
namespace optimization {

using walberla::uint_t;

/// @name Macro face indexing optimization utility functions.
/// By recalculation of the logical indices from a rectangle to the
/// triangular layout, it is possible to loop over a rectangle to access all
/// face DoFs.
///
/// This technique might help the compiler optimizing loops.
///
/// The layout is transformed as follows, using row wise access patterns (example for level 2):
///
/// 14
/// 12 13
/// 09 10 11        <-->   09 10 11 12 13
/// 05 06 07 08            05 06 07 08 14
/// 00 01 02 03 04         00 01 02 03 04
///
///
/// Example:
/// \code
///
/// // loop over rectangular domain
/// // size is given by helper functions
///
/// const uint_t width = 9 // == level 3
/// const uint_t rectWidth = unwrapCols< width >();
/// const uint_t rectHeight = unwrapRows< width >();
///
/// for ( uint_t row = 0; row < rectHeight; row++ )
/// {
///   for ( uint_t col = 0; col < rectwidth; col++ )
///   {
///     // convert to actual logical coordinates
///     const uint_t actualCol = unwrapCol< width >( col, row );
///     const uint_t actualRow = unwrapRow< width >( col, row );
///
///     // call index function
///     uint_t idx = index( actualCol, actualRow );
///   }
/// }
///
/// \code
///@{
template< uint_t width >
inline constexpr uint_t unwrapNumCols()
{
  return width % 2 == 0 ? width + 1 : width;
}

template< uint_t width >
inline constexpr uint_t unwrapNumRows()
{
  return width % 2 == 0 ? width / 2 : ( width + 1 ) / 2;
}

template< uint_t width >
inline constexpr uint_t unwrapCol( const uint_t & col, const uint_t & row )
{
  return col >= width - row ? col - ( width - row ) : col;
}

template< uint_t width >
inline constexpr uint_t unwrapRow( const uint_t & col, const uint_t & row )
{
  if ( width % 2 == 0 )
  {
    return col >= width - row ? ( width - 1 ) - row : row;
  }
  else
  {
    return col >= width - row ? ( width - 1 ) - ( row - 1 ) : row;
  }
}
///@}

}
}
}
