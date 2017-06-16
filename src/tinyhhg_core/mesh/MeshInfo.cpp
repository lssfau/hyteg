
#include "tinyhhg_core/mesh/MeshInfo.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

#include <array>

namespace hhg {

MeshInfo MeshInfo::fromGmshFile( const std::string & meshFileName )
{
  MeshInfo meshInfo;
  const std::array< DoFType, 3 > BOUNDARY_TYPE_TO_FLAG = { Inner, DirichletBoundary, NeumannBoundary };

  WALBERLA_LOG_INFO_ON_ROOT( "[Mesh] Opening mesh file: " << meshFileName );

  std::ifstream meshFile;
  meshFile.open( meshFileName.c_str() );

  WALBERLA_CHECK( !!meshFile, "[Mesh] Error opening file: " << meshFileName );

  std::string token;
  meshFile >> token; // $MeshFormat

  WALBERLA_CHECK_EQUAL( token, "$MeshFormat", "[Mesh] Missing: $MeshFormat" );

  meshFile >> token; // version

  WALBERLA_CHECK_EQUAL( token, "2.2", "[Mesh] Meshfile version should be 2.2");

  meshFile >> token; // 0
  meshFile >> token; // 8
  meshFile >> token; // $EndMeshFormat
  meshFile >> token; // $Nodes

  WALBERLA_CHECK_EQUAL( token, "$Nodes", "[Mesh] Missing: $Nodes");

  uint_t numVertices;
  meshFile >> numVertices;

  for ( uint_t vertex = 0; vertex < numVertices; ++vertex )
  {
    uint_t id;
    real_t x[3];
    meshFile >> id;
    meshFile >> x[0];
    meshFile >> x[1];
    meshFile >> x[2];

    meshInfo.vertices_[id] = Point3D( x );
  }

  WALBERLA_ASSERT_EQUAL( meshInfo.vertices_.size(), numVertices );

  meshFile >> token; // $EndNodes
  meshFile >> token; // $Elements

  WALBERLA_CHECK_EQUAL( token, "$Elements", "[Mesh] Missing: $Elements" )

  uint_t numPrimitives;
  meshFile >> numPrimitives;



  for ( uint_t primitive = 0; primitive < numPrimitives; ++primitive )
  {
    uint_t ig;
    meshFile >> ig; // ignore

    uint_t primitiveType;
    meshFile >> primitiveType;

    if (primitiveType == 1) // edge
    {
      uint_t type;
      uint_t v0, v1;
      meshFile >> ig; // ignore
      meshFile >> type;
      meshFile >> ig; // ignore
      meshFile >> v0;
      meshFile >> v1;

      WALBERLA_ASSERT_LESS( type, BOUNDARY_TYPE_TO_FLAG.size() );
      DoFType dofType = BOUNDARY_TYPE_TO_FLAG[type];

      meshInfo.addEdge( v0, v1, dofType );
    }

    else if (primitiveType == 2) // triangle
    {
      uint_t v0, v1, v2;
      meshFile >> ig; // ignore
      meshFile >> ig; // ignore
      meshFile >> ig; // ignore
      meshFile >> v0;
      meshFile >> v1;
      meshFile >> v2;

      DoFType dofType = BOUNDARY_TYPE_TO_FLAG[0];

      meshInfo.addEdge( v0, v1, dofType );
      meshInfo.addEdge( v1, v2, dofType );
      meshInfo.addEdge( v2, v0, dofType );

      meshInfo.faces_.insert( std::make_tuple( v0, v1, v2 ) );
    }

    else if (primitiveType == 15) // vertex
    {
      // do nothing
      meshFile >> ig;
      meshFile >> ig;
      meshFile >> ig;
      meshFile >> ig;
    }

    else
    {
      WALBERLA_ABORT( "[Mesh] Unknown element type: " << primitiveType );
    }
  }

  meshFile.close();

  return meshInfo;
}


void MeshInfo::addEdge( uint_t v0, uint_t v1, DoFType dofType )
{
  WALBERLA_CHECK_UNEQUAL( v0, v1, "[Mesh] File contains edge with zero size." );

  std::pair< uint_t, uint_t > edge;
  if ( v0 < v1 )
  {
    edge = std::make_pair( v0, v1 );
  }
  else if ( v1 < v0 )
  {
  edge = std::make_pair( v1, v0 );
  }

  edges_.insert( std::make_pair( edge, dofType ) );
}

} // namespace hhg
