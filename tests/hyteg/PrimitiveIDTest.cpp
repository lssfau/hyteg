
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
  PrimitiveID id1( 42 );
  PrimitiveID id2( id1 );

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
}

} // namespace hyteg


int main()
{
   walberla::debug::enterTestMode();

   hyteg::testPrimitiveID();

   return EXIT_SUCCESS;
}
