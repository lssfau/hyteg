
#include "tinyhhg_core/tinyhhg.hpp"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

namespace hhg {

class TestData
{
public:
  bool a = false;
  int i = 100;
  std::vector<bool> aa;
};

class VertexTestData
{
public:
  bool a = false;
  int i = 200;
  std::vector<bool> aa;
};

class TestDataHandling : public NoSerializePrimitiveDataHandling< TestData >
{
public:

  TestData * initialize( const Primitive * const block ) const
  {
    TestData * testData = new TestData();
    testData->i = 7777;
    return testData;
  }

};

class VertexTestDataHandling : public NoSerializePrimitiveDataHandling< VertexTestData >
{
public:

  VertexTestData * initialize( const Primitive * const block ) const
  {
    VertexTestData * testData = new VertexTestData();
    testData->i = 8888;
    return testData;
  }

};

static void testPrimitiveData()
{


  PrimitiveStorage storage("");

  PrimitiveID primitiveID = storage.addVertex();
  Primitive *primitive = storage.getVertex( primitiveID );


  TestDataHandling testDataHandling;
  VertexTestDataHandling vertexTestDataHandling;

  PrimitiveDataID< TestData > testDataID = storage.addPrimitiveData( testDataHandling, "primitive data" );
  TestData * testData = primitive->getData( testDataID );
  WALBERLA_CHECK_EQUAL( testData->i, 7777 );

  PrimitiveDataID< VertexTestData > vertexTestDataID = storage.addVertexData( vertexTestDataHandling, "vertex data" );
  VertexTestData * vertexTestData = primitive->getData( vertexTestDataID );
  WALBERLA_CHECK_EQUAL( vertexTestData->i, 8888 );


#if 0
  TestData *p_testData = primitive.getData( testDataID );
  OtherTestData *p_otherTestData = primitive.getData( otherTestDataID );

  WALBERLA_CHECK_EQUAL( p_testData->i, 100 );
  WALBERLA_CHECK_EQUAL( p_otherTestData->i, 200 );
#endif
}

} // namespace hhg


int main(int argc, char* argv[])
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testPrimitiveData();

   return EXIT_SUCCESS;
}
