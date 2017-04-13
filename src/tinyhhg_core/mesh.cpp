#include "mesh.hpp"

#include <fstream>
#include <fmt/format.h>
#include "types/pointnd.hpp"
#include "types/flags.hpp"
#include "comm.hpp"

namespace hhg
{

Boundary BoundaryTypeToFlag[] = { Inner, DirichletBoundary, NeumannBoundary };

Mesh::Mesh(const std::string& filename)
{
  if (Comm::get().rk == 0)
  {
    fmt::printf("[Mesh] Opening mesh file: %s\n", filename);
  }

  std::ifstream meshFile;
  meshFile.open(filename.c_str());

  if(!meshFile)
  {
    fmt::printf("[Mesh] Error opening file: %s\n", filename);
    std::exit(-1);
  }

  std::string token;
  meshFile >> token; // $MeshFormat

  if (token != "$MeshFormat")
  {
    fmt::printf("[Mesh] Missing: $MeshFormat\n");
    std::exit(-1);
  }

  meshFile >> token; // version

  if (token != "2.2")
  {
    fmt::printf("[Mesh] Meshfile version should be 2.2\n");
    std::exit(-1);
  }

  meshFile >> token; // 0
  meshFile >> token; // 8
  meshFile >> token; // $EndMeshFormat
  meshFile >> token; // $Nodes

  if (token != "$Nodes")
  {
    fmt::printf("[Mesh] Missing: $Nodes\n");
    std::exit(-1);
  }

  size_t numVertices;
  size_t numElements;

  meshFile >> numVertices;
  vertices.reserve(numVertices);

  for (size_t i=0; i < numVertices; ++i)
  {
    size_t id;
    double x[3];
    meshFile >> id; // ignore
    meshFile >> x[0];
    meshFile >> x[1];
    meshFile >> x[2];

    vertices.push_back(Vertex(i, Point3D(x)));

    // fmt::print("[Mesh] Added: {}\n", vertices[i]);
  }

  meshFile >> token; // $EndNodes
  meshFile >> token; // $Elements

  if (token != "$Elements")
  {
    fmt::printf("[Mesh] Missing: $Elements\n");
    std::exit(-1);
  }

  meshFile >> numElements;

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

      face_edge_indices.push_back({e0, e1, e2});
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

  meshFile.close();

  // setting up pointers
  faces.reserve(face_edge_indices.size());
  for (auto& fei : face_edge_indices)
  {
    Edge* e[] = { &edges[fei[0]], &edges[fei[1]], &edges[fei[2]] };
    faces.push_back(Face(faces.size(), e));

    for (size_t i = 0; i < 3; ++i)
    {
      e[i]->addFace(&faces[faces.size()-1]);
    }

    for (Vertex* v : faces[faces.size()-1].vertices)
    {
      v->addFace(&faces[faces.size()-1]);
    }

    // fmt::print("[Mesh] Added: {}\n", faces[faces.size()-1]);
  }

  for (Edge& edge : edges)
  {
    edge.v0->addEdge(&edge);
    edge.v1->addEdge(&edge);
  }

  // improve rank distribution locally
  for (Face& face : faces)
  {
    for (Edge* edge : face.edges)
    {
      edge->rank = face.rank;
    }

    for (Vertex* vertex : face.vertices)
    {
      vertex->rank = face.rank;
    }
  }
}

size_t Mesh::addEdge(size_t v0, size_t v1, size_t type, Mesh::EdgeMap& map)
{
  std::pair<size_t, size_t> key = std::make_pair(v0, v1);

  if (key.first > key.second)
  {
    std::swap(key.first, key.second);
  }

  if (map.find(key) != map.end())
  {
    return map[key];
  }

  size_t edgeId = this->edges.size();
  map[key] = edgeId;

  this->edges.push_back(Edge(edgeId, BoundaryTypeToFlag[type], &this->vertices[v0], &this->vertices[v1]));

  // TODO: this is a workaround
  if (this->edges[edgeId].v0->type != DirichletBoundary)
  {
    this->edges[edgeId].v0->type = BoundaryTypeToFlag[type];
  }

  if (this->edges[edgeId].v1->type != DirichletBoundary)
  {
    this->edges[edgeId].v1->type = BoundaryTypeToFlag[type];
  }

  // fmt::print("[Mesh] Added: {}\n", this->edges[edgeId]);

  return edgeId;
}

}