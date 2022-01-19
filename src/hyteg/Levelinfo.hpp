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

namespace hyteg {
namespace levelinfo {

using walberla::uint_c;
using walberla::uint_t;

constexpr inline uint_t num_microvertices_per_vertex( uint_t /*level*/ )
{
   return 1;
}

constexpr inline uint_t num_microvertices_per_edge( uint_t level )
{
   return ( ( uint_t( 1 ) << level ) + uint_t( 1 ) );
   //return (uint_t) std::pow(2, level) + 1;
}

constexpr inline uint_t num_microvertices_per_edge_from_width( uint_t width )
{
   return width;
}

constexpr inline uint_t num_microedges_per_edge( uint_t level )
{
   return num_microvertices_per_edge( level ) - 1;
}

constexpr inline uint_t num_microedges_per_edge_from_width( const uint_t& width )
{
   return num_microvertices_per_edge_from_width( width ) - 1;
}

constexpr inline uint_t num_microedges_per_face_from_width( const uint_t& width )
{
   return 3 * ( ( ( num_microedges_per_edge_from_width( width ) + 1 ) * num_microedges_per_edge_from_width( width ) ) / 2 );
}

constexpr inline uint_t num_microedges_per_face( uint_t level )
{
   return 3 * ( ( ( num_microedges_per_edge( level ) + 1 ) * num_microedges_per_edge( level ) ) / 2 );
}

constexpr inline uint_t num_microvertices_per_face_from_width( const uint_t& width )
{
   return ( ( width * ( width + 1 ) ) >> uint_t( 1 ) );
}

constexpr inline uint_t num_microvertices_per_face( uint_t level )
{
   const uint_t width = ( uint_t( 1 ) << level ) + 1;
   return num_microvertices_per_face_from_width( width );
}

constexpr inline uint_t num_microfaces_per_face( uint_t level )
{
   return ( uint_t( 1 ) << ( 2 * level ) );
}

constexpr inline uint_t num_microfaces_per_face_from_width( uint_t width )
{
   return ( width * width + width ) / 2 + ( ( width - 1 ) * ( width - 1 ) + ( width - 1 ) ) / 2;
}

constexpr inline uint_t num_microvertices_per_cell_from_width( const uint_t& width )
{
   return ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
}

constexpr inline uint_t num_microedges_per_cell_from_width( const uint_t& width )
{
   return 6 * num_microvertices_per_cell_from_width( width - 1 ) + num_microvertices_per_cell_from_width( width - 2 );
}

constexpr inline uint_t num_microvertices_per_cell( const uint_t& level )
{
   const uint_t width = ( uint_t( 1 ) << level ) + 1;
   return ( ( width + 2 ) * ( width + 1 ) * width ) / 6;
}

constexpr inline uint_t num_microedges_per_cell( const uint_t& level )
{
   const uint_t width = ( uint_t( 1 ) << level ) + 1;
   return num_microedges_per_cell_from_width( width );
}

constexpr inline uint_t num_microcells_per_cell( const uint_t& level )
{
   // num_microcells_per_cell = 8 ^ level
   return ( uint_t( 1 ) << ( 3 * level ) );
}

inline uint_t num_microcells_per_cell_from_width( const uint_t& width )
{
   WALBERLA_ASSERT_GREATER_EQUAL( width, 2 )

   if ( width == 2 )
   {
      return 1;
   }

   const uint_t whiteUp   = num_microvertices_per_cell_from_width( width - 1 );
   const uint_t whiteDown = num_microvertices_per_cell_from_width( width - 3 );
   const uint_t others    = num_microvertices_per_cell_from_width( width - 2 );
   return whiteUp + whiteDown + 4 * others;
}

} // namespace levelinfo

} // namespace hyteg
