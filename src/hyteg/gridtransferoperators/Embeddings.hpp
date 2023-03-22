/*
 * Copyright (c) 2023 Marcus Mohr.
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
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

/// Embed a P1 function into the P2 function space
///
/// The free-function takes a P1Function and converts it to a P2Function by
/// embedding it into the corresponding P2 space.
template < typename ValueType >
void embedP1IntoP2( const P1Function< ValueType >& src,
                    const P2Function< ValueType >& dst,
                    const uint_t&                  level,
                    const DoFType&                 flag = All );

/// Embed a P1 vector function into the P2 vector function space
///
/// The free-function takes a P1VectorFunction and converts it to a
/// P2VectorFunction by embedding it into the corresponding P2^dim space.
template < typename ValueType >
void embedP1IntoP2( const P1VectorFunction< ValueType >& src,
                    const P2VectorFunction< ValueType >& dst,
                    const uint_t&                        level,
                    const DoFType&                       flag = All );

} // namespace hyteg
