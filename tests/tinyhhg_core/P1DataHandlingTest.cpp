
#include "tinyhhg_core/tinyhhg.hpp"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

namespace hhg {

static void testP1DataHandling()
{
  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  uint_t minLevel =  1;
  uint_t maxLevel = 10;

  std::string meshFileName = "../../data/meshes/tri_2el.msh";

  MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
  PrimitiveStorage      storage( setupStorage );

  WALBERLA_LOG_INFO( setupStorage );

  VertexP1FunctionMemoryDataHandling vertexP1FunctionMemoryDataHandling( minLevel, maxLevel );
  EdgeP1FunctionMemoryDataHandling   edgeP1FunctionMemoryDataHandling  ( minLevel, maxLevel );
  FaceP1FunctionMemoryDataHandling   faceP1FunctionMemoryDataHandling  ( minLevel, maxLevel );

  PrimitiveDataID< VertexP1FunctionMemory, Vertex > vertexP1FunctionMemoryID = storage.addVertexData( vertexP1FunctionMemoryDataHandling, "P1 vertex data" );
  PrimitiveDataID< EdgeP1FunctionMemory,   Edge   > edgeP1FunctionMemoryID   = storage.addEdgeData  ( edgeP1FunctionMemoryDataHandling,   "P1 edge data"   );
  PrimitiveDataID< FaceP1FunctionMemory,   Face >   faceP1FunctionMemoryID   = storage.addFaceData  ( faceP1FunctionMemoryDataHandling,   "P1 face data"   );

  for ( const auto & it : storage.getVertices() )
  {
    PrimitiveID vertexID = it.first;
    Vertex *    vertex   = it.second;

    VertexP1FunctionMemory * mem = vertex->getData( vertexP1FunctionMemoryID );

    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      WALBERLA_LOG_INFO( "Vertex: " << vertexID.getID() << " | data size (level " << level << ") : " << mem->getSize( level ) );
    }
  }

  for ( const auto & it : storage.getEdges() )
  {
    PrimitiveID edgeID = it.first;
    Edge *      edge   = it.second;

    EdgeP1FunctionMemory * mem = edge->getData( edgeP1FunctionMemoryID );

    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      WALBERLA_LOG_INFO( "Edge:   " << edgeID.getID() << " | data size (level " << level << ") : " << mem->getSize( level ) );
    }
  }

  for ( const auto & it : storage.getFaces() )
  {
    PrimitiveID faceID = it.first;
    Face *      face   = it.second;

    FaceP1FunctionMemory * mem = face->getData( faceP1FunctionMemoryID );

    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      WALBERLA_LOG_INFO( "Face:   " << faceID.getID() << " | data size (level " << level << ") : " << mem->getSize( level ) );
    }
  }

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testP1DataHandling();

   return EXIT_SUCCESS;
}
