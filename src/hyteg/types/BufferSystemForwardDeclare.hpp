/*
* Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl,
 * Marcus Mohr.
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


namespace walberla {
namespace mpi {
template< typename T >
class GenericRecvBuffer;
using RecvBuffer = GenericRecvBuffer<unsigned char>;

struct OptimalGrowth;
template< typename T
   , typename G >
class GenericSendBuffer;
using SendBuffer = GenericSendBuffer<unsigned char, OptimalGrowth>;

template< typename RecvBuffer_T, typename SendBuffer_T >
class GenericOpenMPBufferSystem;
using OpenMPBufferSystem = GenericOpenMPBufferSystem<RecvBuffer, SendBuffer>;
}
}