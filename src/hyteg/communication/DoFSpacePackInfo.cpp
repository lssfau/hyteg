/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include "DoFSpacePackInfo.hpp"

#include "hyteg/communication/PackInfo.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/primitivedata/PrimitiveDataID.hpp"

namespace hyteg {
namespace communication {

template < typename ValueType >
DoFSpacePackInfo< ValueType >::DoFSpacePackInfo( uint_t                                                 level,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                                                 std::weak_ptr< PrimitiveStorage >                      storage )
: level_( level )
, dataIDVertex_( dataIDVertex )
, dataIDEdge_( dataIDEdge )
, dataIDFace_( dataIDFace )
, dataIDCell_( storage.lock()->generateInvalidPrimitiveDataID< FunctionMemory< ValueType >, Cell >() )
, storage_( storage )
{}

template < typename ValueType >
DoFSpacePackInfo< ValueType >::DoFSpacePackInfo( uint_t                                                 level,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                                                 std::weak_ptr< PrimitiveStorage >                      storage )
: level_( level )
, dataIDVertex_( dataIDVertex )
, dataIDEdge_( dataIDEdge )
, dataIDFace_( dataIDFace )
, dataIDCell_( dataIDCell )
, storage_( storage )
{}

template class DoFSpacePackInfo< double >;
template class DoFSpacePackInfo< float >;
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template class DoFSpacePackInfo< walberla::float16 >;
#endif
template class DoFSpacePackInfo< int >;
template class DoFSpacePackInfo< long >;
template class DoFSpacePackInfo< long long >;
template class DoFSpacePackInfo< uint_t >;

} // namespace communication
} // namespace hyteg
