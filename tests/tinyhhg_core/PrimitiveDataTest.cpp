
#include "tinyhhg_core/tinyhhg.hpp"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

namespace hhg {

std::string getExampleMeshFileContent()
{
  std::string content =
  "$MeshFormat\n"
  "2.2 0 8\n"
  "$EndMeshFormat\n"
  "$Nodes\n"
  "4\n"
  "1 0 0 0\n"
  "2 1 0 0\n"
  "3 0 1 0\n"
  "4 1 1 0\n"
  "$EndNodes\n"
  "$Elements\n"
  "9\n"
  "1 15 2 1 1 1\n"
  "2 15 2 1 2 2\n"
  "3 15 2 1 3 3\n"
  "4 1 2 1 1 1 2\n"
  "5 1 2 1 2 2 4\n"
  "6 1 2 1 2 4 3\n"
  "7 1 2 1 3 3 1\n"
  "8 2 2 0 5 2 4 1\n"
  "9 2 2 0 5 1 4 3\n"
  "$EndElements\n";

  return content;
}

void writeTestMeshFile( const std::string & meshFileName )
{
  std::string meshFileContent = getExampleMeshFileContent();
  std::ofstream file( meshFileName );
  file << meshFileContent;
  file.close();
}

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

class TestDataHandling : public OnlyInitializeDataHandling< TestData, Primitive >
{
public:

  std::shared_ptr< TestData > initialize( const Primitive * const ) const
  {
    auto testData = std::make_shared< TestData >();
    testData->i = 7777;
    return testData;
  }

};

class VertexTestDataHandling : public OnlyInitializeDataHandling< VertexTestData, Vertex >
{
public:

  std::shared_ptr< VertexTestData > initialize( const Vertex * const ) const
  {
    auto testData = std::make_shared< VertexTestData >();
    testData->i = 8888;
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
    testData->i = 9999;
    return testData;
  }

};

static void testPrimitiveData()
{
  std::string meshFileName = "./tmpMeshFile.msh";
  writeTestMeshFile( meshFileName );

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  WALBERLA_LOG_INFO( setupStorage );

  PrimitiveStorage storage( setupStorage );

  auto testDataHandling = std::make_shared< TestDataHandling >();
  auto vertexTestDataHandling = std::make_shared< VertexTestDataHandling >();
  auto edgeTestDataHandling = std::make_shared< EdgeTestDataHandling >();

  // Adding data to all primitives
  PrimitiveDataID< TestData, Primitive > testDataID;
  storage.addPrimitiveData( testDataID, testDataHandling, "primitive data" );
  // Adding data only to vertices
  PrimitiveDataID< VertexTestData, Vertex > vertexTestDataID;
  storage.addVertexData( vertexTestDataID, vertexTestDataHandling, "vertex data" );

  // Obtaining initialized vertex data from a vertex
  for ( auto it = storage.beginVertices(); it != storage.endVertices(); it++ )
  {
    WALBERLA_LOG_PROGRESS( "Checking content of vertex with ID: " << it->second->getID().getID() );
    VertexTestData * vertexTestData = it->second->getData( vertexTestDataID );
    WALBERLA_CHECK_EQUAL( vertexTestData->i, 8888 );
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
