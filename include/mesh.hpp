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

/// \brief  Macro mesh implementation with \a GMSH file reader
/// \author Daniel Drzisga (drzisga@ma.tum.de)
/// \date   March, 2017
///
/// The Mesh class is responsible for reading a \a GMSH mesh file and creating all
/// the appropriate primitive objects. These primitives are stored in separate std::vector
/// objects for each primitive type. Since pointers into these vectors are later used, modifications
/// of these vectors is prohibited.
class Mesh
{
public:

  /// Mesh constructor that reads the given \a GMSH file and mesh builds the mesh topology using primitives.
  /// \param filename filename of \a GMSH file.
  Mesh(const std::string& filename);

  /// Vector containing all macro vertices
  std::vector<Vertex> vertices;

  /// Vector containing all macro edges
  std::vector<Edge> edges;

  /// Vector containing all macro faces
  std::vector<Face> faces;

private:
  struct pairhash {
  public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const
    {
      return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
  };

  typedef std::unordered_map<std::pair<size_t, size_t>, size_t, pairhash> EdgeMap;

  size_t addEdge(size_t v0, size_t v1, size_t type, EdgeMap& map);

};

}

#endif /* MESH_HPP */