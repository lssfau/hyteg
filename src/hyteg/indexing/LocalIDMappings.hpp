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

#include <map>
#include <set>

namespace hyteg {
namespace indexing {

using walberla::uint_t;

const std::map< uint_t, std::set< uint_t > > cellLocalEdgeIDsToCellLocalNeighborFaceIDs = {
  { 0, std::set< uint_t >( { 0, 1 } ) },
  { 1, std::set< uint_t >( { 0, 2 } ) },
  { 2, std::set< uint_t >( { 0, 3 } ) },
  { 3, std::set< uint_t >( { 1, 2 } ) },
  { 4, std::set< uint_t >( { 1, 3 } ) },
  { 5, std::set< uint_t >( { 2, 3 } ) }
};

const std::map< uint_t, std::set< uint_t > > faceLocalEdgeIDsToSpanningVertexIDs = {
    { 0, std::set< uint_t >( { 0, 1 } ) },
    { 1, std::set< uint_t >( { 0, 2 } ) },
    { 2, std::set< uint_t >( { 1, 2 } ) },
};

const std::map< uint_t, std::set< uint_t > > cellLocalFaceIDsToSpanningVertexIDs = {
  { 0, std::set< uint_t >( { 0, 1, 2 } ) },
  { 1, std::set< uint_t >( { 0, 1, 3 } ) },
  { 2, std::set< uint_t >( { 0, 2, 3 } ) },
  { 3, std::set< uint_t >( { 1, 2, 3 } ) }
};

inline uint_t getCellLocalFaceIDFromCellLocalVertexIDs( const uint_t v0, const uint_t v1, const uint_t v2 )
{
	const auto v3 = 6 - (v0 + v1 + v2);
	return 3 - v3;
}

inline uint_t getCellLocalEdgeIDFromCellLocalVertexIDs( const uint_t v0, const uint_t v1 )
{
  WALBERLA_ASSERT_UNEQUAL( v0, v1 );
  WALBERLA_ASSERT_LESS_EQUAL( v0, 3 );
  WALBERLA_ASSERT_LESS_EQUAL( v1, 3 );
  switch ( v0 )
  {
    case 0:
      switch ( v1 )
      {
        case 1: return 0;
        case 2: return 1;
        case 3: return 3;
        default: WALBERLA_ABORT("Invalid vertex IDs");
      }
    case 1:
      switch ( v1 )
      {
        case 0: return 0;
        case 2: return 2;
        case 3: return 4;
        default: WALBERLA_ABORT("Invalid vertex IDs");
      }
    case 2:
      switch ( v1 )
      {
        case 0: return 1;
        case 1: return 2;
        case 3: return 5;
        default: WALBERLA_ABORT("Invalid vertex IDs");
      }
    case 3:
      switch ( v1 )
      {
        case 0: return 3;
        case 1: return 4;
        case 2: return 5;
        default: WALBERLA_ABORT("Invalid vertex IDs");
      }
    default: WALBERLA_ABORT("Invalid vertex IDs");
  }
}

inline uint_t getCellLocalOppositeEdgeID( const uint_t & cellLocalEdgeID )
{
   WALBERLA_ASSERT_LESS_EQUAL( cellLocalEdgeID, 5 );
   return 5 - cellLocalEdgeID;
}

}
}
