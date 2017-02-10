#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <vector>
#include <fmt/ostream.h>
#include <types/pointnd.hpp>

namespace hhg
{

class Edge;
class Face;

class Vertex
{
public:
  Vertex(size_t id, const Point3D& coords);
  void addEdge(Edge* edge);
  void addFace(Face* face);

  size_t edge_index(const Edge& edge) const;

  size_t id;
  int rank;
  size_t type;
  Point3D coords;

  std::vector<Edge*> edges;
  std::vector<Face*> faces;

  std::vector<std::vector<double*> > data;
  std::vector<std::vector<double*> > opr_data;

  friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex);
};

}

#endif /* VERTEX_HPP */