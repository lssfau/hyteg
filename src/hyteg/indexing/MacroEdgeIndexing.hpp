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

#include "core/debug/Debug.h"
#include "hyteg/indexing/Common.hpp"

namespace hyteg {
namespace indexing {

using walberla::uint_t;

namespace layout {

/// Required memory for the linear macro edge layout
constexpr uint_t linearMacroEdgeSize( const uint_t & width )
{
  return width;
}

/// General linear memory layout indexing function for macro edges
constexpr uint_t linearMacroEdgeIndex( const uint_t & width, const idx_t & x )
{
  WALBERLA_UNUSED( width );
  return uint_t( x );
}

}


constexpr uint_t macroEdgeSize( const uint_t & width )
{
  return layout::linearMacroEdgeSize( width );
}


constexpr uint_t macroEdgeIndex( const uint_t & width, const idx_t & x )
{
  return layout::linearMacroEdgeIndex( width, x );
}

/// Iterator over an edge.
/// Does not include ghost layers!
/// It is possible to parameterize the iterator to only iterate over a inner part of the edge.
/// This is done by setting the offset parameter to the distance to the vertices.
/// If set to zero, the iterator iterates over the whole edge (including both adjacent vertices).
class EdgeIterator
{
public:

  using iterator_category = std::input_iterator_tag;
  using value_type        = Index;
  using reference         = value_type const&;
  using pointer           = value_type const*;
  using difference_type   = ptrdiff_t;

  EdgeIterator( const uint_t & width, const uint_t & offsetToCenter = 0, const bool & backwards = false, const bool & end = false ) :
    width_( width ), offsetToCenter_( offsetToCenter ), backwards_( backwards ),
    totalNumberOfDoFs_( width - 2 * offsetToCenter ), step_( 0 )
  {
    WALBERLA_ASSERT_GREATER( width, 0, "Size of edge must be larger than zero!" );
    WALBERLA_ASSERT_LESS( offsetToCenter, width, "Offset to center is beyond edge width!" );

    coordinates_.z() = 0;
    coordinates_.y() = 0;

    if ( backwards )
    {
       coordinates_.x() = static_cast< idx_t >( width - 1 - offsetToCenter );
    }
    else
    {
       coordinates_.x() = static_cast< idx_t >( offsetToCenter );
    }


    if ( end )
    {
      step_ = totalNumberOfDoFs_;
    }
  }

  EdgeIterator begin() { return EdgeIterator( width_, offsetToCenter_, backwards_ ); }
  EdgeIterator end()   { return EdgeIterator( width_, offsetToCenter_, backwards_, true ); }

  bool operator==( const EdgeIterator & other ) const
  {
    return other.step_ == step_;
  }

  bool operator!=( const EdgeIterator & other ) const
  {
    return other.step_ != step_;
  }

  reference operator*()  const { return  coordinates_; };
  pointer   operator->() const { return &coordinates_; };

  EdgeIterator & operator++() // prefix
  {
    WALBERLA_ASSERT_LESS_EQUAL( step_, totalNumberOfDoFs_, "Incrementing iterator beyond end!" );

    step_++;
    if ( backwards_ )
    {
      coordinates_.x()--;
    }
    else
    {
      coordinates_.x()++;
    }

    return *this;
  }

  EdgeIterator operator++( int ) // postfix
  {
    const EdgeIterator tmp( *this );
    ++*this;
    return tmp;
  }


private:

  uint_t    width_;
  uint_t    offsetToCenter_;
  bool      backwards_;
  uint_t    totalNumberOfDoFs_;
  uint_t    step_;
  Index     coordinates_;

};
}
}
