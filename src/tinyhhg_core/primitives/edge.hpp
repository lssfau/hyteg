#ifndef EDGE_HPP
#define EDGE_HPP

#include "tinyhhg_core/types/pointnd.hpp"

#include <vector>
#include <fmt/ostream.h>

namespace hhg
{

class Vertex;
class Face;

class Edge
{
public:

  Edge(size_t id, size_t type, Vertex* v0, Vertex* v1);
  void addFace(Face* face);

  size_t vertex_index(const Vertex& vertex) const;
  size_t face_index(const Face& face) const;

  Vertex* get_opposite_vertex(const Vertex& vertex) const;

  size_t id;
  int rank;
  size_t type;
  Vertex* v0;
  Vertex* v1;

  Point3D direction;
  double length;
  Point3D tangent;
  Point3D normal_2d;

  std::vector<Face*> faces;
  std::vector<std::vector<double*> > data;
  std::vector<std::vector<double*> > opr_data;

  friend std::ostream &operator<<(std::ostream &os, const Edge &edge);
};

}

#endif /* EDGE_HPP */
