#ifndef FACE_HPP
#define FACE_HPP

#include "tinyhhg_core/types/pointnd.hpp"
#include <core/DataTypes.h>

#include <array>
#include <vector>

namespace hhg
{

class Edge;

class Face
{
public:
  Face(size_t id, Edge* edges[3]);

  size_t vertex_index(const Vertex& vertex) const;
  size_t edge_index(const Edge& edge) const;

  std::vector<Edge*> adjacent_edges(const Vertex& vertex) const;
  Vertex* get_vertex_opposite_to_edge(const Edge& edge) const;

  size_t id;
  walberla::uint_t rank;
  walberla::uint_t type;
  double area;
  Edge* edges[3];
  std::vector<Vertex*> vertices;
  std::array<int, 3> edge_orientation;
  std::array<Point3D, 3> coords;

  std::vector<std::vector<double*> > data;
  std::vector<std::vector<double*> > opr_data;

  friend std::ostream &operator<<(std::ostream &os, const Face &face);
};

}

#endif /* FACE_HPP */
