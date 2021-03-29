/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

/// Convert a P1 to a P2 function
///
/// The function takes a P1Function on level (P2Level-1) and converts it to a P2Function on level P2Level by
/// assigning the DoFs of the P1Function to those of the P2Function on the coarser mesh level.
template < typename ValueType >
void P1toP2Conversion( const P1Function< ValueType >& src,
                       const P2Function< ValueType >& dst,
                       const uint_t&                  P2Level,
                       const DoFType&                 flag = All );
}
