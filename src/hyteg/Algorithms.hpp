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
#include "core/debug/Debug.h"

#include <algorithm>
#include <numeric>
#include <set>
#include <vector>

namespace hyteg {
namespace algorithms {

using walberla::uint_t;

/// \brief Convenience function to reduce code clutter when interested if container contains an element but there is no 'std::contains'...
template< typename Container, typename ElementType >
inline bool contains( const Container & container, const ElementType & element )
{
  return std::find( container.begin(), container.end(), element ) != container.end();
}

/// Returns an array concatenates all missing natural numbers of [0, OutputArraySize) in the input array to the input array in ascending order.
/// Example: input array: [2, 1], input size: 4 => output array [2, 1, 0, 3]
template< uint_t InputArraySize, uint_t OutputArraySize >
inline std::array< uint_t, OutputArraySize > getMissingIntegersAscending( const std::array< uint_t, InputArraySize > & inputArray )
{
  static_assert( InputArraySize <= OutputArraySize, "Input array must be a subset of output array!" );
  WALBERLA_ASSERT( *std::max_element( inputArray.begin(), inputArray.end() ) < OutputArraySize );
  WALBERLA_ASSERT_EQUAL( std::set< uint_t >( inputArray.begin(), inputArray.end() ).size(), InputArraySize );

  std::array< uint_t, OutputArraySize > outputArray;
  std::copy( inputArray.begin(), inputArray.end(), outputArray.begin() );

  uint_t pos = 0;
  for ( uint_t i = 0; i < OutputArraySize; i++ )
  {
    if ( std::find( inputArray.begin(), inputArray.end(), i ) == inputArray.end() )
    {
      outputArray[pos+InputArraySize] = i;
      pos++;
    }
  }
  return outputArray;
}

/// Creates a permutation vector that is being sorted as if the passed vector would have been sorted instead.
///
/// Thanks to https://stackoverflow.com/a/17074810.
///
template < typename T, typename Compare >
std::vector< std::size_t > sortPermutation( const std::vector< T >& vec, Compare compare )
{
   std::vector< std::size_t > p( vec.size() );
   std::iota( p.begin(), p.end(), 0 );
   std::sort( p.begin(), p.end(), [&]( std::size_t i, std::size_t j ) { return compare( vec[i], vec[j] ); } );
   return p;
}

/// Sorts the given vector according to the passed permutation.
///
/// Use this combined with sortPermutation() to sort two or more vectors the same way:
///
/// // pseudocode
/// auto permutation = sortPermutation( vec1, compare );
/// vec1Sorted = applyPermutation( vec1, permutation ):
/// vec2Sorted = applyPermutation( vec2, permutation ):
///
/// Thanks to https://stackoverflow.com/a/17074810.
///
template < typename T >
std::vector< T > applyPermutation( const std::vector< T >& vec, const std::vector< std::size_t >& p )
{
   std::vector< T > sorted_vec( vec.size() );
   std::transform( p.begin(), p.end(), sorted_vec.begin(), [&]( std::size_t i ) { return vec[i]; } );
   return sorted_vec;
}

}
}