#ifndef P1FUNCTION_HPP
#define P1FUNCTION_HPP

#include "tinyhhg_core/mesh.hpp"
#include "tinyhhg_core/function.hpp"

#include <core/mpi/all.h>

#include "p1vertex.hpp"
#include "p1edge.hpp"
#include "p1face.hpp"

namespace hhg
{

class P1Function : public Function
{
public:

  P1Function(const std::string& _name, Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : Function(_name, _mesh, _minLevel, _maxLevel)
  {
    for (Vertex& v : mesh.vertices)
    {
      if (v.rank == rank)
      {
        memory_id = v.data.size();
        break;
      }
    }

    if (memory_id == -1)
    {
      for (Edge& e : mesh.edges)
      {
        if (e.rank == rank)
        {
          memory_id = e.data.size();
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
          memory_id = f.data.size();
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
        P1Vertex::allocate(vertex, memory_id, minLevel, maxLevel);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank)
      {
        P1Edge::allocate(edge, memory_id, minLevel, maxLevel);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank)
      {
        P1Face::allocate(face, memory_id, minLevel, maxLevel);
      }
    }
  }

  ~P1Function()
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank)
      {
        P1Vertex::free(vertex, memory_id, minLevel, maxLevel);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank)
      {
        P1Edge::free(edge, memory_id, minLevel, maxLevel);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank)
      {
        P1Face::free(face, memory_id, minLevel, maxLevel);
      }
    }
  }

  void interpolate(std::function<double(const hhg::Point3D&)>& expr, size_t level, size_t flag = All)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::interpolate(vertex, memory_id, expr, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::interpolate(edge, memory_id, expr, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::interpolate(face, memory_id, expr, level);
      }
    }
  }

  void assign(const std::vector<double> scalars, const std::vector<P1Function*> functions, size_t level, size_t flag = All)
  {
    std::vector<size_t> src_ids(functions.size());
    for (size_t i = 0; i < functions.size(); ++i)
    {
      src_ids[i] = functions[i]->memory_id;
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::assign(vertex, scalars, src_ids, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::assign(edge, scalars, src_ids, memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::assign(face, scalars, src_ids, memory_id, level);
      }
    }
  }

  void add(const std::vector<double> scalars, const std::vector<P1Function*> functions, size_t level, size_t flag = All)
  {
    std::vector<size_t> src_ids(functions.size());
    for (size_t i = 0; i < functions.size(); ++i)
    {
      src_ids[i] = functions[i]->memory_id;
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::add(vertex, scalars, src_ids, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::add(edge, scalars, src_ids, memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::add(face, scalars, src_ids, memory_id, level);
      }
    }
  }

  double dot(Function& rhs, size_t level, size_t flag = All)
  {
    double sp_l = 0.0;

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        sp_l += P1Vertex::dot(vertex, memory_id, rhs.memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        sp_l += P1Edge::dot(edge, memory_id, rhs.memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        sp_l += P1Face::dot(face, memory_id, rhs.memory_id, level);
      }
    }

    double sp_g = 0.0;
    MPI_Allreduce(&sp_l, &sp_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sp_g;
  }

  void apply(Operator& opr, Function& dst, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos(vertex, memory_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::apply(vertex, opr.id, memory_id, dst.memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, dst.memory_id, level);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos(edge, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::apply(edge, opr.id, memory_id, dst.memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, dst.memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::apply(face, opr.id, memory_id, dst.memory_id, level);
      }
    }
  }

  void smooth_gs(Operator& opr, Function& rhs, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos(vertex, memory_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::smooth_gs(vertex, opr.id, memory_id, rhs.memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, memory_id, level);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos(edge, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::smooth_gs(edge, opr.id, memory_id, rhs.memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::smooth_gs(face, opr.id, memory_id, rhs.memory_id, level);
      }
    }
  }

  void prolongate(size_t level, size_t flag = All)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::prolongate(vertex, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, memory_id, level+1);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::prolongate(edge, memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, memory_id, level+1);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::prolongate(face, memory_id, level);
      }
    }
  }

  void restrict(size_t level, size_t flag = All)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos(vertex, memory_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::restrict(vertex, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, memory_id, level-1);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos(edge, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::restrict(edge, memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, memory_id, level-1);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::restrict(face, memory_id, level);
      }
    }
  }

  void printmatrix(Operator& opr, size_t level, size_t flag = All)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      P1Vertex::pull_halos(vertex, memory_id, level);
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::printmatrix(vertex, opr.id, memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_halos(edge, memory_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::printmatrix(edge, opr.id, memory_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::printmatrix(face, opr.id, memory_id, level);
      }
    }
  }
};

}

#endif
