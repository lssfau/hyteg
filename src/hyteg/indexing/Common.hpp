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

#include "hyteg/types/PointND.hpp"
#include "hyteg/types/types.hpp"

using walberla::uint_t;

namespace hyteg {
namespace indexing {

/// Helper function to create a 'tuple' from two integers taking less than 5 bits of space.
/// Allows for more readable switch statements in some cases.
inline constexpr uint_t tup4( const uint_t& a, const uint_t& b )
{
#if 0
   assert( a < ( 1 << 4 ) );
   assert( b < ( 1 << 4 ) );
#endif
   return a << 4 | b;
}

inline constexpr uint_t tup4( const uint_t& a, const uint_t& b, const uint_t& c )
{
#if 0
   assert( a < ( 1 << 4 ) );
   assert( b < ( 1 << 4 ) );
   assert( c < ( 1 << 4 ) );
#endif
   return a << 8 | b << 4 | c;
}

/// Wrapper around Point3D for convenient access to logical indices.
class Index : public PointND< idx_t, 3 >
{
   using PointND< idx_t, 3 >::PointND;

 public:
   static Index max()
   {
      return { std::numeric_limits< idx_t >::max(), std::numeric_limits< idx_t >::max(), std::numeric_limits< idx_t >::max() };
   }
};

inline bool operator<( const Index& lhs, const Index& rhs )
{
   return lhs.z() < rhs.z() || ( lhs.z() == rhs.z() && lhs.y() < rhs.y() ) ||
          ( lhs.z() == rhs.z() && lhs.y() == rhs.y() && lhs.x() < rhs.x() );
}

inline std::ostream& operator<<( std::ostream& os, const Index& index )
{
   os << "( " << index.x() << ", " << index.y() << ", " << index.z() << " )";
   return os;
}

} // namespace indexing
} // namespace hyteg

// this can be changed to Index once Index is same as Index
namespace std {
template <>
struct less< hyteg::PointND< hyteg::idx_t, 3 > >
{
   bool operator()( const hyteg::PointND< hyteg::idx_t, 3 >& lhs, const hyteg::PointND< hyteg::idx_t, 3 >& rhs ) const
   {
      return lhs.z() < rhs.z() || ( lhs.z() == rhs.z() && lhs.y() < rhs.y() ) ||
             ( lhs.z() == rhs.z() && lhs.y() == rhs.y() && lhs.x() < rhs.x() );
   }
};
} // namespace std