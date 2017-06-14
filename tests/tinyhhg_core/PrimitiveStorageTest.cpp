
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

static void testPrimitiveStorage()
{



}

} // namespace hhg


int main()
{
   walberla::debug::enterTestMode();

   hhg::testPrimitiveStorage();

   return EXIT_SUCCESS;
}
