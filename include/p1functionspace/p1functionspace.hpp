#ifndef P1FUNCTIONSPACE_HPP
#define P1FUNCTIONSPACE_HPP

#include "../mesh.hpp"
#include "../comm.hpp"

#include "p1vertex.hpp"
#include "p1edge.hpp"
#include "p1face.hpp"

namespace hhg
{

class P1FunctionSpace
{
public:

  P1FunctionSpace(Mesh& _mesh)
    : mesh(_mesh), rank(Comm::get().rk)
  {
  }

  size_t allocate(size_t minLevel, size_t maxLevel)
  {
    // size_t memory_id = mesh.vertices[0].data.size();
    size_t memory_id = -1;

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

    return memory_id;
  }

  void free(size_t memory_id, size_t minLevel, size_t maxLevel)
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

  void interpolate(size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level, size_t flag)
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

  template<size_t N>
  void assign(std::array<double, N>& scalars, std::array<size_t, N>& src_ids, size_t dst_id, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::assign<N>(vertex, scalars, src_ids, dst_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, dst_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::assign<N>(edge, scalars, src_ids, dst_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, dst_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::assign<N>(face, scalars, src_ids, dst_id, level);
      } 
    }
  }

  template<size_t N>
  void add(std::array<double, N>& scalars, std::array<size_t, N>& src_ids, size_t dst_id, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::add<N>(vertex, scalars, src_ids, dst_id, level);
      }
    }
    
    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, dst_id, level);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::add<N>(edge, scalars, src_ids, dst_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, dst_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::add<N>(face, scalars, src_ids, dst_id, level);
      } 
    }
  }

  double dot(size_t lhs_id, size_t rhs_id, size_t level, size_t flag)
  {
    double sp_l = 0.0;

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        sp_l += P1Vertex::dot(vertex, lhs_id, rhs_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        sp_l += P1Edge::dot(edge, lhs_id, rhs_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        sp_l += P1Face::dot(face, lhs_id, rhs_id, level);
      } 
    }

    double sp_g = 0.0;
    MPI_Allreduce(&sp_l, &sp_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sp_g;
  }

  void apply(size_t opr_id, size_t src_id, size_t dst_id, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos(vertex, src_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::apply(vertex, opr_id, src_id, dst_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, dst_id, level);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos(edge, src_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::apply(edge, opr_id, src_id, dst_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, dst_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::apply(face, opr_id, src_id, dst_id, level);
      } 
    }
  }

  void smooth_gs(size_t opr_id, size_t f_id, size_t rhs_id, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos(vertex, f_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::smooth_gs(vertex, opr_id, f_id, rhs_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, f_id, level);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos(edge, f_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::smooth_gs(edge, opr_id, f_id, rhs_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, f_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::smooth_gs(face, opr_id, f_id, rhs_id, level);
      } 
    }
  }

  void prolongate(size_t f_id, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::prolongate(vertex, f_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, f_id, level+1);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::prolongate(edge, f_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, f_id, level+1);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::prolongate(face, f_id, level);
      } 
    }
  }

  void restrict(size_t f_id, size_t level, size_t flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos(vertex, f_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::restrict(vertex, f_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices(edge, f_id, level-1);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos(edge, f_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::restrict(edge, f_id, level);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges(face, f_id, level-1);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::restrict(face, f_id, level);
      } 
    }
  }

  Mesh& mesh;

private:
  int rank;
};

}

#endif