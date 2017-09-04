
#include "tinyhhg_core/tinyhhg.hpp"
#include "core/Environment.h"
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

  ~VertexTestData() { delete[] f; }

  bool a = false;
  int i = 200;
  std::vector<bool> aa;
  float * f;
};

class EdgeTestData
{
public:
  bool a = false;
  int i = 300;
  std::vector<bool> aa;
};

class CellTestData
{
public:
  bool a = false;
  int i = 400;
  std::vector<bool> aa;
};


class TestDataHandling : public OnlyInitializeDataHandling< TestData, Primitive >
{
public:

  std::shared_ptr< TestData > initialize( const Primitive * const ) const
  {
    auto testData = std::make_shared< TestData >();
    testData->i = 5555;
    return testData;
  }

};

class VertexTestDataHandling : public OnlyInitializeDataHandling< VertexTestData, Vertex >
{
public:

  std::shared_ptr< VertexTestData > initialize( const Vertex * const ) const
  {
    auto testData = std::make_shared< VertexTestData >();
    testData->i = 6666;
    testData->f = new float[10000];
    return testData;
  }

};

class EdgeTestDataHandling : public OnlyInitializeDataHandling< EdgeTestData, Edge >
{
public:

  std::shared_ptr< EdgeTestData > initialize( const Edge * const ) const
  {
    auto testData = std::make_shared< EdgeTestData >();
    testData->i = 7777;
    return testData;
  }

};

class CellTestDataHandling : public OnlyInitializeDataHandling< CellTestData, Cell >
{
public:

  std::shared_ptr< CellTestData > initialize( const Cell * const ) const
  {
    auto testData = std::make_shared< CellTestData >();
    testData->i = 9999;
    return testData;
  }

};

static void testPrimitiveData()
{

  const std::string meshFileName = "../../data/meshes/3D/cube_24el.msh";
  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  WALBERLA_LOG_INFO( setupStorage );

  PrimitiveStorage storage( setupStorage );

  auto testDataHandling       = std::make_shared< TestDataHandling >();
  auto vertexTestDataHandling = std::make_shared< VertexTestDataHandling >();
  auto edgeTestDataHandling   = std::make_shared< EdgeTestDataHandling >();
  auto cellTestDataHandling   = std::make_shared< CellTestDataHandling >();

  // Adding data to all primitives
  PrimitiveDataID< TestData, Primitive > testDataID;
  storage.addPrimitiveData( testDataID, testDataHandling, "primitive data" );
  // Adding data only to vertices
  PrimitiveDataID< VertexTestData, Vertex > vertexTestDataID;
  storage.addVertexData( vertexTestDataID, vertexTestDataHandling, "vertex data" );
  // Adding data only to cells
  PrimitiveDataID< CellTestData, Cell > cellTestDataID;
  storage.addCellData( cellTestDataID, cellTestDataHandling, "cell data" );

  std::vector< PrimitiveID > primitiveIDs;
  storage.getPrimitiveIDs( primitiveIDs );
  for ( const auto & id : primitiveIDs )
  {
    WALBERLA_LOG_PROGRESS( "Checking content of primitive with ID: " << id.getID() );
    auto primitive = storage.getPrimitive( id );
    TestData * testData = primitive->getData( testDataID );
    WALBERLA_CHECK_EQUAL( testData->i, 5555 );
  }

  // Obtaining initialized vertex data from a vertex
  for ( const auto & it : storage.getVertices() )
  {
    WALBERLA_LOG_PROGRESS( "Checking content of vertex with ID: " << it.second->getID().getID() );
    VertexTestData * vertexTestData = it.second->getData( vertexTestDataID );
    WALBERLA_CHECK_EQUAL( vertexTestData->i, 6666 );
  }

  // Obtaining initialized cell data from a cell
  for ( const auto & it : storage.getCells() )
  {
    WALBERLA_LOG_PROGRESS( "Checking content of cell with ID: " << it.second->getID().getID() );
    CellTestData * cellTestData = it.second->getData( cellTestDataID );
    WALBERLA_CHECK_EQUAL( cellTestData->i, 9999 );
  }

  WALBERLA_UNUSED( testDataID );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testPrimitiveData();

   return EXIT_SUCCESS;
}
