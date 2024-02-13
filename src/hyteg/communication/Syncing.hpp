/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

// Some forward declarations
template < typename funcType >
class P2Function;

class FEFunctionRegistry;

namespace communication {

/// Describe direction of data synchronisation between primitives
enum class syncDirection_t
{
   LOW2HIGH,     ///< synchronise halos from low-dimensional to high-dimensional primitives, i.e. vertex -> edge -> face (-> cell)
   BIDIRECTIONAL ///< completely synchronise halos between all primitives
};

using walberla::uint_t;

template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function,
                                    const uint_t&   level,
                                    syncDirection_t direction = syncDirection_t::BIDIRECTIONAL );

// template < template< class > class funcType, typename vType >
template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< vType >& function,
                                          const uint_t&                    level,
                                          syncDirection_t                  direction = syncDirection_t::BIDIRECTIONAL );

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< vType >& function,
                                          const uint_t&                    level,
                                          syncDirection_t                  direction = syncDirection_t::BIDIRECTIONAL );

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const EGFunction< vType >& function,
                                          const uint_t&              level,
                                          syncDirection_t            direction = syncDirection_t::BIDIRECTIONAL );

/// Sync all functions registered with the passed FEFunctionRegistry object
///
/// \note The function currently only syncs registered functions with the following value types:
/// - double
/// - int32_t
/// - int64_t
void syncRegisteredFunctions( const FEFunctionRegistry& feFunctionRegistry,
                              uint_t                    level,
                              syncDirection_t           direction = syncDirection_t::BIDIRECTIONAL );

} // namespace communication
} // namespace hyteg
