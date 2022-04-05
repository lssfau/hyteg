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
#include "hyteg/indexing/Common.hpp"

#include <set>

namespace hyteg {
namespace indexing {

using walberla::uint_t;
using walberla::uint16_t;

/// Array access layouts - wrapped by general access function.
/// The general memory layout can thus be switched globally by setting
/// the called functions in the general size and access function.
namespace layout {

/// Required memory for the linear macro cell layout
inline constexpr uint_t linearMacroCellSize( const uint_t & width )
{
  return ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
}

/// General linear memory layout indexing function for macro cells
inline constexpr uint_t linearMacroCellIndex( const uint_t & width, const idx_t & x, const idx_t & y, const idx_t & z )
{
  const uint_t widthMinusSlice = width - uint_t( z );
  const uint_t sliceOffset = linearMacroCellSize( width ) - (( widthMinusSlice + 2 ) * ( widthMinusSlice + 1 ) * widthMinusSlice ) / 6;
  const uint_t rowOffset   = uint_t( y ) * ( widthMinusSlice + 1 ) - ( uint_t( ( ( y + 1 ) * ( y ) ) / 2 ) );
  return sliceOffset + rowOffset + uint_t( x );
}

} // namespace layout


inline constexpr uint_t macroCellSize( const uint_t & width )
{
  return layout::linearMacroCellSize( width );
}

inline constexpr uint_t macroCellIndex( const uint_t & width, const idx_t & x, const idx_t & y, const idx_t & z )
{
  return layout::linearMacroCellIndex( width, x, y, z );
}


/// Returns the local face indices of the cell if the index is located on a face of the cell.
/// Note that an index can be located on multiple faces (e.g. if it lies on an edge).
/// If it is not located on any face, the returned set is empty.
std::set< uint_t > isOnCellFace( const indexing::Index & index, const uint_t & width );


/// Returns the local edge indices of the cell if the index is located on an edge of the cell.
/// Note that an index can be located on multiple edges (e.g. if it lies on a vertex).
/// If it is not located on any edge, the returned set is empty.
std::set< uint_t > isOnCellEdge( const indexing::Index & index, const uint_t & width );


/// Returns a set with the vertex index of the cell if the index is located on a vertex of the cell
/// and empty set otherwise.
std::set< uint_t > isOnCellVertex( const indexing::Index & index, const uint_t & width );


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

  CellIterator( const uint_t & width, const uint_t & offsetToCenter = 0, const bool & end = false );

  CellIterator begin() { return CellIterator( width_, offsetToCenter_ ); }
  CellIterator end()   { return CellIterator( width_, offsetToCenter_, true ); }

  bool operator==( const CellIterator & other ) const { return other.step_ == step_; }
  bool operator!=( const CellIterator & other ) const { return other.step_ != step_; }

  reference operator*()  const { return  coordinates_; };
  pointer   operator->() const { return &coordinates_; };

  CellIterator & operator++(); // prefix
  CellIterator operator++( int );

private:

  uint_t    width_;
  uint_t    internalWidth_;
  uint_t    offsetToCenter_;
  uint_t    totalNumberOfDoFs_;
  uint_t    step_;
  Index     coordinates_;
  Index     internalCoordinates_;

};


/// \brief Iterator to iterate over a face-boundary of a macro-cell.
///
/// Using this iterator, each face with corner-coordinates a, b and c can be iterated over in six ways:
///
/// c
/// +
/// |\
/// | \
/// |  \
/// +---+
/// a    b
///
/// - a -> b -> c (bottom to top in row-wise fashion)
/// - b -> a -> c (bottom to tip in row-wise fashion (backwards))
/// - a -> c -> b (left to right in column-wise fashion (columns bottom to top))
/// - c -> a -> b (left to right in column-wise fashion (columns top to bottom))
/// - b -> c -> a (top-right to bottom-left (diagonal rows from bottom right to top left))
/// - c -> b -> a (top-right to bottom-left (diagonal rows from top left to bottom right))
///
/// a, b and c can be set to three distinct cell-vertices (combinations of [0, 1, 2, 3], no repetition)
///
/// The cell vertices are at the following logical indices (also refer to the documentation):
///
/// 0: (      0,       0,       0)
/// 1: (width-1,       0,       0)
/// 2: (      0, width-1,       0)
/// 3: (      0,       0, width-1)
///
class CellBoundaryIterator
{
public:

  using iterator_category = std::input_iterator_tag;
  using value_type        = Index;
  using reference         = value_type const&;
  using pointer           = value_type const*;
  using difference_type   = ptrdiff_t;

  CellBoundaryIterator( const uint_t & width, const uint_t & vertex0,
                      const uint_t & vertex1, const uint_t & vertex2,
                      const uint_t & offsetToCenter = 0, const bool & end = false );

  CellBoundaryIterator( const uint_t & width, const std::array< uint_t, 3 > & vertices,
                      const uint_t & offsetToCenter = 0, const bool & end = false );

  CellBoundaryIterator begin() { return CellBoundaryIterator( width_, vertices_, offsetToCenter_       ); }
  CellBoundaryIterator end()   { return CellBoundaryIterator( width_, vertices_, offsetToCenter_, true ); }

  bool operator==( const CellBoundaryIterator & other ) const { return other.step_ == step_; }
  bool operator!=( const CellBoundaryIterator & other ) const { return other.step_ != step_; }

  reference operator*()  const { return   coordinates_; };
  pointer   operator->() const { return & coordinates_; };

  CellBoundaryIterator & operator++(); // prefix
  CellBoundaryIterator operator++( int );

private:

  IndexIncrement calculateIncrement( const uint_t & vertex0, const uint_t & vertex1 ) const;

  const uint_t                  width_;
  const std::array< uint_t, 3 > vertices_;
  const uint_t                  offsetToCenter_;
  const uint_t                  totalNumberOfSteps_;
  const IndexIncrement          firstDirIncrement_;
  const IndexIncrement          secondDirIncrement_;
        uint_t                  step_;
        uint_t                  wrapAroundStep_;
        Index                   coordinates_;
        Index                   wrapAroundCoordinates_;
        uint_t                  wrapArounds_;
};



} // namespace indexing
} // namespace hyteg
