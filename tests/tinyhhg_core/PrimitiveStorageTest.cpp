
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace walberla {
namespace hhg {

static void testPrimitiveStorage()
{

  PrimitiveStorage storage;

}

} // namespace hhg
} // namespace walberla


int main()
{
   walberla::debug::enterTestMode();

   walberla::hhg::testPrimitiveStorage();

   return EXIT_SUCCESS;
}
