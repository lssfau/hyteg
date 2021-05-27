/*
 * Copyright (c) 2020 Marcus Mohr.
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

using walberla::uint_t;

#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {
namespace indexing {

template < typename srcType, typename dstType >
uint_t countLocalDoFCouplings( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   WALBERLA_UNUSED( storage );
   WALBERLA_ABORT( "countDoFCouplings() not specialised for requested pair of source/destination function" );
};

template <>
uint_t countLocalDoFCouplings< P1FunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level );

template <>
uint_t countLocalDoFCouplings< P2FunctionTag, P2FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level );

template <>
uint_t countLocalDoFCouplings< P1VectorFunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                     uint_t                                     level );

template <>
uint_t countLocalDoFCouplings< P2VectorFunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                     uint_t                                     level );
template <>
uint_t countLocalDoFCouplings< P1FunctionTag, P2VectorFunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                     uint_t                                     level );

template <>
uint_t countLocalDoFCouplings< EdgeDoFFunctionTag, EdgeDoFFunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                         uint_t                                     level );
template <>
uint_t countLocalDoFCouplings< EdgeDoFFunctionTag, VertexDoFFunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                           uint_t                                     level );

} // namespace indexing
} // namespace hyteg
