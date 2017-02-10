#ifndef MESH_HPP
#define MESH_HPP

#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

#include "primitives/vertex.hpp"
#include "primitives/edge.hpp"
#include "primitives/face.hpp"

namespace hhg
{

class Mesh
{
public:

  struct pairhash {
  public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const
    {
      return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
  };

  typedef std::unordered_map<std::pair<size_t, size_t>, size_t, pairhash> EdgeMap;

  Mesh(const std::string& filename);

  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  std::vector<Face> faces;

private:
  size_t addEdge(size_t v0, size_t v1, size_t type, EdgeMap& map);

};

}

#endif /* MESH_HPP */