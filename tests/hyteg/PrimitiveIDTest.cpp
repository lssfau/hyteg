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

#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "hyteg/PrimitiveID.hpp"

namespace hyteg {

// to be removed when moving to walberla namespace
using walberla::mpi::SendBuffer;
using walberla::mpi::RecvBuffer;

static void testPrimitiveID()
{

  PrimitiveID id0;
  PrimitiveID id1 = PrimitiveID::create( 42 );
  PrimitiveID id2( id1 );

  WALBERLA_ABORT( "Extend this test for new pid implementation." );
#if 0
  WALBERLA_CHECK_EQUAL( id0.getID(), uint_t(0) );
  WALBERLA_CHECK_EQUAL( id1.getID(), uint_t(42) );
  WALBERLA_CHECK_EQUAL( id2.getID(), uint_t(42) );

  SendBuffer sendBuffer;

  PrimitiveID id3, id4, id5;

  sendBuffer << id0 << id1 << id2;

  RecvBuffer recvBuffer( sendBuffer );

  recvBuffer >> id3 >> id4 >> id5;

  WALBERLA_CHECK_EQUAL( id3.getID(), uint_t(0) );
  WALBERLA_CHECK_EQUAL( id4.getID(), uint_t(42) );
  WALBERLA_CHECK_EQUAL( id5.getID(), uint_t(42) );

  PrimitiveID id6( uint_t(257) );

  WALBERLA_CHECK_EQUAL( id3.getUsedBits(), 0 );
  WALBERLA_CHECK_EQUAL( id4.getUsedBits(), 6 );
  WALBERLA_CHECK_EQUAL( id6.getUsedBits(), 9 );

  WALBERLA_CHECK_EQUAL( id3.getUsedBytes(), 0 );
  WALBERLA_CHECK_EQUAL( id4.getUsedBytes(), 1 );
  WALBERLA_CHECK_EQUAL( id6.getUsedBytes(), 2 );

#endif
}

} // namespace hyteg


int main()
{
   walberla::debug::enterTestMode();

   hyteg::testPrimitiveID();

   return EXIT_SUCCESS;
}
