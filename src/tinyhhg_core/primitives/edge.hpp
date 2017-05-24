#ifndef EDGE_HPP
#define EDGE_HPP

#include "tinyhhg_core/types/pointnd.hpp"

#include <vector>

#include <tinyhhg_core/types/pointnd.hpp>
#include <tinyhhg_core/types/flags.hpp>

#include <core/DataTypes.h>

namespace hhg
{

class Vertex;
class Face;

class EdgeMemory;

class Edge
{
public:

  Edge(size_t id, DoFType type, Vertex* v0, Vertex* v1);
  void addFace(Face* face);

  size_t vertex_index(const Vertex& vertex) const;
  size_t face_index(const Face& face) const;

  Vertex* get_opposite_vertex(const Vertex& vertex) const;

  size_t id;
  walberla::uint_t rank;
  DoFType type;
  Vertex* v0;
  Vertex* v1;

  Point3D direction;
  real_t length;
  Point3D tangent;
  Point3D normal_2d;

  std::vector<Face*> faces;

  std::vector<EdgeMemory*> memory;

  friend std::ostream &operator<<(std::ostream &os, const Edge &edge);
};



class EdgeMemory
{
public:

  const MemoryType type;

  virtual void free() = 0;

protected:
  EdgeMemory(MemoryType t) : type(t) { ; }
};
}

#endif /* EDGE_HPP */
