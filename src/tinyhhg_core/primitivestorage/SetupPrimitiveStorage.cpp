
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"

namespace hhg {

MeshInfo::MeshInfo( const std::string & meshFileName )
{
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

    vertices_[id] = Point3D( x );
  }

  WALBERLA_ASSERT_EQUAL( vertices_.size(), numVertices );

  meshFile >> token; // $EndNodes
  meshFile >> token; // $Elements

  WALBERLA_CHECK_EQUAL( token, "$Elements", "[Mesh] Missing: $Elements" )

  uint_t numElements;
  meshFile >> numElements;

#if 0
  Mesh::EdgeMap map;
  std::vector<std::array<size_t, 3> > face_edge_indices;

  for (size_t i=0; i < numElements; ++i)
  {
    size_t ig;
    meshFile >> ig; // ignore

    size_t elemType;
    meshFile >> elemType;

    if (elemType == 1) // edge
    {
      size_t type;
      size_t v0, v1;
      meshFile >> ig; // ignore
      meshFile >> type;
      meshFile >> ig; // ignore
      meshFile >> v0;
      --v0;
      meshFile >> v1;
      --v1;

      addEdge(v0, v1, type, map);
    }
    else if (elemType == 2) // triangle
    {
      size_t v0, v1, v2;
      meshFile >> ig; // ignore
      meshFile >> ig; // ignore
      meshFile >> ig; // ignore
      meshFile >> v0;
      --v0;
      meshFile >> v1;
      --v1;
      meshFile >> v2;
      --v2;

      size_t e0 = addEdge(v0, v1, 0, map);
      size_t e1 = addEdge(v1, v2, 0, map);
      size_t e2 = addEdge(v2, v0, 0, map);

      face_edge_indices.push_back({{e0, e1, e2}});
    }
    else if (elemType == 15) // vertex
    {
      // do nothing
      meshFile >> ig;
      meshFile >> ig;
      meshFile >> ig;
      meshFile >> ig;
    }
    else
    {
      fmt::print("[Mesh] Unknown element type: {}\n", elemType);
      std::exit(-1);
    }
  }
#endif
  meshFile.close();
}

}
