#ifndef EDGE_HPP
#define EDGE_HPP

#include "tinyhhg_core/types/pointnd.hpp"

#include <vector>
#include <map>
#include <fmt/ostream.h>

namespace hhg
{

class Vertex;
class Face;

class EdgeMemory;

class Edge
{
public:

  Edge(size_t id, size_t type, Vertex* v0, Vertex* v1);
  void addFace(Face* face);

  size_t vertex_index(const Vertex& vertex) const;
  size_t face_index(const Face& face) const;

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

	std::vector<EdgeMemory*> memory;

  friend std::ostream &operator<<(std::ostream &os, const Edge &edge);
};



class EdgeMemory
{
public:
	enum EdgeMemoryType { Base, Stencil, P1 };

	const EdgeMemoryType type;

	virtual void free() = 0;

protected:
	EdgeMemory(EdgeMemoryType t) : type(t) { ; }
};


class EdgeStencilMemory
	:public EdgeMemory
{
public:
	EdgeStencilMemory() : EdgeMemory(Stencil) { ; }

	std::map<size_t, double*> data;

	virtual void free();

};


class EdgeP1Memory
	:public EdgeMemory
{
public:
	EdgeP1Memory() : EdgeMemory(P1) { ; }

	std::map<size_t, double*> data;

	virtual void free();
};


}

#endif /* EDGE_HPP */