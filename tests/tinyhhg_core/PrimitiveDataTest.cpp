
#include "tinyhhg_core/primitives/Primitive.hpp"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"

namespace walberla {
namespace hhg {

class TestData
{
public:
  bool a = false;
  int i = 100;
  std::vector<bool> aa;
};

class OtherTestData
{
public:
  bool a = false;
  int i = 200;
  std::vector<bool> aa;
};

static void testPrimitiveData()
{

  Primitive primitive;

  PrimitiveDataID< TestData >      testDataID( 0 );
  PrimitiveDataID< OtherTestData > otherTestDataID( 42 );

  TestData testData;
  OtherTestData otherTestData;

  internal::PrimitiveData testDataWrapper( &testData );
  internal::PrimitiveData otherTestDataWrapper( &otherTestData );

  primitive.addData( testDataID, &testDataWrapper );
  primitive.addData( otherTestDataID, &otherTestDataWrapper );

  TestData *p_testData = primitive.getData( testDataID );
  OtherTestData *p_otherTestData = primitive.getData( otherTestDataID );

  WALBERLA_CHECK_EQUAL( p_testData->i, 100 );
  WALBERLA_CHECK_EQUAL( p_otherTestData->i, 200 );

}

} // namespace hhg
} // namespace walberla


int main()
{
   walberla::debug::enterTestMode();

   walberla::hhg::testPrimitiveData();

   return EXIT_SUCCESS;
}
