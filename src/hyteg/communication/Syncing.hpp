/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/p1dgefunctionspace/P1DGEFunction.hpp"

namespace hyteg {

template < typename funcType >
class P2Function;

namespace communication {

using walberla::uint_t;

template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function, const uint_t& level );

// template < template< class > class funcType, typename vType >
template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< vType >& function, const uint_t& level );

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< vType >& function, const uint_t& level );

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P1DGEFunction< vType >& function, const uint_t& level );

template < typename ValueType >
void syncP2FunctionBetweenPrimitives( const P2Function< ValueType >& function, const uint_t& level );

} // namespace communication
} // namespace hyteg
