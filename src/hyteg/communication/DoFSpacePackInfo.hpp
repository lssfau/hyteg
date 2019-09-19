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

#include "hyteg/communication/PackInfo.hpp"

namespace hyteg {

class Vertex;
class Edge;
class Face;
class Cell;
class PrimitiveStorage;

template < typename ValueType >
class FunctionMemory;
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

namespace communication {

using walberla::uint_t;

template < typename ValueType >
class DoFSpacePackInfo : public communication::PackInfo
{
 public:
   DoFSpacePackInfo( uint_t                                                 level,
                     PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                     PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                     PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                     std::weak_ptr< PrimitiveStorage >                      storage );

   DoFSpacePackInfo( uint_t                                                 level,
                     PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                     PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                     PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                     PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                     std::weak_ptr< PrimitiveStorage >                      storage );

 protected:
   uint_t                                                 level_;
   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell_;
   std::weak_ptr< hyteg::PrimitiveStorage >                 storage_;
};

} // namespace communication
} // namespace hyteg
