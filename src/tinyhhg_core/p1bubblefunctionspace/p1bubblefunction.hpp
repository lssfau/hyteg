#pragma once

#include "tinyhhg_core/mesh.hpp"
#include "tinyhhg_core/OldFunction.hpp"

#include <core/mpi/all.h>

#include "p1bubblevertex.hpp"
#include "p1bubbleedge.hpp"
#include "p1bubbleface.hpp"
#include "P1BubblePackInfo.hpp"

namespace hhg
{
#if 0

//FIXME remove after we are in walberla namespace
using namespace walberla::mpistubs;

class P1BubbleFunction : public Function
{
public:

  P1BubbleFunction( const std::string& name,
                    const std::shared_ptr< PrimitiveStorage > & storage,
                    uint_t minLevel,
                    uint_t maxLevel );

  ~P1BubbleFunction();

  /// Interpolates a given expression to a P1Function
  void interpolate(std::function<real_t(const hhg::Point3D&)>& expr, uint_t level, DoFType flag = All);

  void assign(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All);

  void add(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All);

  real_t dot(P1BubbleFunction& rhs, size_t level, DoFType flag = All);

  void prolongate(size_t level, DoFType flag = All);

  void prolongateQuadratic(size_t level, DoFType flag = All);

  void restrict(size_t level, DoFType flag = All);


const PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> &getVertexDataID() const { return vertexDataID_; }

const PrimitiveDataID<EdgeP1BubbleFunctionMemory, Edge> &getEdgeDataID() const { return edgeDataID_; }

const PrimitiveDataID<FaceP1BubbleFunctionMemory, Face> &getFaceDataID() const { return faceDataID_; }

private:
  PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> vertexDataID_;
  PrimitiveDataID<EdgeP1BubbleFunctionMemory, Edge> edgeDataID_;
  PrimitiveDataID<FaceP1BubbleFunctionMemory, Face> faceDataID_;
};


class P1BubbleFunction : public OldFunction
{
public:

  P1BubbleFunction(const std::string& _name, Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : OldFunction(_name, _mesh, _minLevel, _maxLevel)
  {
    for (Vertex& v : mesh.vertices)
    {
      if (v.rank == rank)
      {
        memory_id = v.memory.size();
        break;
      }
    }

    if (memory_id == -1)
    {
      for (Edge& e : mesh.edges)
      {
        if (e.rank == rank)
        {
          memory_id = e.memory.size();
          break;
        }
      }
    }

    if (memory_id == -1)
    {
      for (Face& f : mesh.faces)
      {
        if (f.rank == rank)
        {
          memory_id = f.memory.size();
          break;
        }
      }
    }

    if (memory_id == -1)
    {
      fmt::printf("There was an error allocating P1 memory\n");
      std::exit(-1);
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank)
      {
        P1BubbleVertex::allocate(vertex, memory_id, minLevel, maxLevel);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank)
      {
        P1BubbleEdge::allocate(edge, memory_id, minLevel, maxLevel);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank)
      {
        P1BubbleFace::allocate(face, memory_id, minLevel, maxLevel);
      }
    }
  }

  ~P1BubbleFunction()
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank)
      {
        P1BubbleVertex::free(vertex, memory_id);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank)
      {
        P1BubbleEdge::free(edge, memory_id);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank)
      {
        P1BubbleFace::free(face, memory_id);
      }
    }
  }

  void interpolate(std::function<real_t(const hhg::Point3D&)>& expr, size_t level, DoFType flag = All)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1BubbleVertex::interpolate(vertex, memory_id, expr, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1BubbleEdge::interpolate(level, edge, memory_id, expr);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1BubbleFace::interpolate(level, face, memory_id, expr);
      }
    }
  }

  void assign(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All)
  {
    std::vector<size_t> src_ids(functions.size());
    for (size_t i = 0; i < functions.size(); ++i)
    {
      src_ids[i] = functions[i]->memory_id;
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank)
      {
        P1BubbleVertex::assign(vertex, scalars, src_ids, memory_id, level, flag);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1BubbleEdge::assign(edge, scalars, src_ids, memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1BubbleFace::assign(level, face, scalars, src_ids, memory_id);
      }
    }
  }

  void add(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All)
  {
    std::vector<size_t> src_ids(functions.size());
    for (size_t i = 0; i < functions.size(); ++i)
    {
      src_ids[i] = functions[i]->memory_id;
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank)
      {
        P1BubbleVertex::add(vertex, scalars, src_ids, memory_id, level, flag);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1BubbleEdge::add(edge, scalars, src_ids, memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1BubbleFace::add(level, face, scalars, src_ids, memory_id);
      }
    }
  }

  real_t dot(P1BubbleFunction& rhs, size_t level, DoFType flag = All)
  {
    real_t sp_l = 0.0;

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank)
      {
        sp_l += P1BubbleVertex::dot(vertex, memory_id, rhs.memory_id, level, flag);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        sp_l += P1BubbleEdge::dot(edge, memory_id, rhs.memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        sp_l += P1BubbleFace::dot(level, face, memory_id, rhs.memory_id);
      }
    }



#ifdef WALBERLA_BUILD_WITH_MPI
    walberla::mpi::allReduceInplace(sp_l,walberla::mpi::SUM,walberla::MPIManager::instance().get()->comm());
#endif
    return sp_l;
  }

  void enumerate(size_t level, size_t& num)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      P1BubbleVertex::enumerate(level, vertex, memory_id, num);
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::enumerate(level, edge, memory_id, num);
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::enumerate(level, face, memory_id, num);
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_halos(edge, memory_id, level);
    }

    for (Vertex& vertex : mesh.vertices)
    {
      P1BubbleVertex::pull_halos(vertex, memory_id, level);
    }
  }

  void enumerate_p1(size_t level, size_t& num)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      P1BubbleVertex::enumerate_p1(level, vertex, memory_id, num);
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      // enumerate_p1 for edges is the same as enumerate
      P1BubbleEdge::enumerate(level, edge, memory_id, num);
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::enumerate_p1(level, face, memory_id, num);
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_halos(edge, memory_id, level);
    }

    for (Vertex& vertex : mesh.vertices)
    {
      P1BubbleVertex::pull_halos(vertex, memory_id, level);
    }
  }

//  void prolongate(size_t level, DoFType flag = All)
//  {
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (vertex.rank == rank && testFlag(vertex.type, flag))
//      {
//        P1BubbleVertex::prolongate(vertex, memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      P1BubbleEdge::pull_vertices(edge, memory_id, level+1);
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      if (edge.rank == rank && testFlag(edge.type, flag))
//      {
//        P1BubbleEdge::prolongate(edge, memory_id, level);
//      }
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      P1BubbleFace::pull_edges(face, memory_id, level+1);
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      if (face.rank == rank && testFlag(face.type, flag))
//      {
//        P1BubbleFace::prolongate(level, face, memory_id);
//      }
//    }
//  }
//
//  void restrict(size_t level, DoFType flag = All)
//  {
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (testFlag(vertex.type, flag))
//      {
//        P1BubbleVertex::pull_halos(vertex, memory_id, level);
//      }
//    }
//
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (vertex.rank == rank && testFlag(vertex.type, flag))
//      {
//        P1BubbleVertex::restrict(vertex, memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      P1BubbleEdge::pull_vertices(edge, memory_id, level-1);
//      if (testFlag(edge.type, flag))
//      {
//        P1BubbleEdge::pull_halos(edge, memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      if (edge.rank == rank && testFlag(edge.type, flag))
//      {
//        P1BubbleEdge::restrict(edge, memory_id, level);
//      }
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      P1BubbleFace::pull_edges(face, memory_id, level-1);
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      if (face.rank == rank && testFlag(face.type, flag))
//      {
//        P1BubbleFace::restrict(level, face, memory_id);
//      }
//    }
//  }
};
#endif

}
