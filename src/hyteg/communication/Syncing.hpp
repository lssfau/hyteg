/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

// Some forward declarations
template < typename funcType >
class P2Function;

class FEFunctionRegistry;

template < typename vType >
class P1VectorFunction;

template < typename vType >
class P2VectorFunction;

template < typename vType >
class EGFunction;

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
/// - float
/// - double
/// - int32_t
/// - int64_t
///
/// It will count the number of functions it syncs and compare this to the number of functions in
/// the registry to detect, if other types are present. In case of an inconsistency it will abort.
///
/// \param feFunctionRegistry      an FEFunctionRegistry object
/// \param level                   provides the level on which halo synchronisation shall be performed
/// \param excludeDGTypeFunctions  if true, functions of DG-type will be excluded from the synching (this
///                                can e.g. be done for data-export with VTKOutput to avoid unneccessary
///                                communication)
/// \param direction               specifies direction of communication (LOW2HIGH or BIDIRECTIONAL)
void syncRegisteredFunctions( const FEFunctionRegistry& feFunctionRegistry,
                              uint_t                    level,
                              bool                      excludeDGTypeFunctions,
                              syncDirection_t           direction = syncDirection_t::BIDIRECTIONAL );

} // namespace communication
} // namespace hyteg
