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

#include "hyteg/PrimitiveID.hpp"

#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

namespace hyteg {

// to be removed when moving to walberla namespace
using walberla::mpi::RecvBuffer;
using walberla::mpi::SendBuffer;

static void testPrimitiveID()
{
   PrimitiveID id0;
   PrimitiveID id1 = PrimitiveID::create( 42 );
   PrimitiveID id2( id1 );

   SendBuffer sendBuffer;

   PrimitiveID id3, id4, id5;

   sendBuffer << id0 << id1 << id2;

   RecvBuffer recvBuffer( sendBuffer );

   recvBuffer >> id3 >> id4 >> id5;

   WALBERLA_CHECK_EQUAL( id3, id0 );
   WALBERLA_CHECK_EQUAL( id4, id1 );
   WALBERLA_CHECK_EQUAL( id5, id2 );

   auto parent = id1;
   for ( uint_t k = 1; k <= PrimitiveID::maxRefinement(); ++k )
   {
      auto children = parent.createChildren();

      SendBuffer snd;
      for ( auto& child : children )
      {
         WALBERLA_CHECK_EQUAL( child.numAncestors(), k );
         WALBERLA_CHECK_EQUAL( child.getParent(), parent );
         snd << child;
      }

      RecvBuffer rcv( snd );
      for ( auto& child : children )
      {
         rcv >> id5;
         WALBERLA_CHECK_EQUAL( child, id5 );
      }

      parent = children[0];
   }
}

} // namespace hyteg

int main()
{
   walberla::debug::enterTestMode();

   hyteg::testPrimitiveID();

   return EXIT_SUCCESS;
}
