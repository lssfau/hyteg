#ifndef FACE_HPP
#define FACE_HPP

#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"

#include <array>
#include <vector>
#include <map>

namespace hhg
{

class Edge;
class FaceMemory;

class Face
{
public:
  Face(size_t id, Edge* edges[3]);

  size_t vertex_index(const Vertex& vertex) const;
  size_t edge_index(const Edge& edge) const;

  std::vector<Edge*> adjacent_edges(const Vertex& vertex) const;
  Vertex* get_vertex_opposite_to_edge(const Edge& edge) const;

  size_t id;
  int rank;
  size_t type;
  double area;
  Edge* edges[3];
  std::vector<Vertex*> vertices;
  std::array<int, 3> edge_orientation;
  std::array<Point3D, 3> coords;

  std::vector<FaceMemory*> memory;

  friend std::ostream &operator<<(std::ostream &os, const Face &face);
};

class FaceMemory
{
public:

	const MemoryType type;

	virtual void free() = 0;

protected:
	FaceMemory(MemoryType t) : type(t) { ; }
};


class FaceStencilMemory
	:public FaceMemory
{
public:
	FaceStencilMemory() : FaceMemory(Stencil) { ; }

	std::map<size_t,double*> data;

	virtual void free();

};


class FaceP1Memory
	:public FaceMemory
{
public:
	FaceP1Memory() : FaceMemory(P1) { ; }

	std::map<size_t,double*> data;

	virtual void free();
};


}

#endif /* FACE_HPP */
