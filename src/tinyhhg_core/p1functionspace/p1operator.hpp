#ifndef P1OPERATOR_HPP
#define P1OPERATOR_HPP

#include <fmt/format.h>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/operator.hpp"

#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"

namespace hhg
{

enum ElementType
{
  UPWARD,
  DOWNWARD
};

void compute_micro_coords(const Face& face, size_t level, walberla::real_t coords[6], ElementType element_type)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D d0 = face.edge_orientation[0] * face.edges[0]->direction / (rowsize-1);
  Point3D d2 = -face.edge_orientation[2] * face.edges[2]->direction / (rowsize-1);

  walberla::real_t orientation = 1.0;

  if (element_type == DOWNWARD) {
    orientation = -1.0;
  }

  coords[0] = 0.0;
  coords[1] = 0.0;
  coords[2] = orientation * d0[0];
  coords[3] = orientation * d0[1];
  coords[4] = orientation * d2[0];
  coords[5] = orientation * d2[1];
}

template<class UFCOperator>
void compute_local_stiffness(const Face& face, size_t level, walberla::real_t local_stiffness[3][3], ElementType element_type)
{
  walberla::real_t A[9];
  walberla::real_t coords[6];
  compute_micro_coords(face, level, coords, element_type);
  UFCOperator gen;
  gen.tabulate_tensor(A, NULL, coords, 0);

  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      local_stiffness[i][j] = A[3 * j + i];
    }
  }
}

template<class UFCOperator>
class P1Operator : public Operator
{
public:
  P1Operator(Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : Operator(_mesh, _minLevel, _maxLevel)
  {
    memory_id = mesh.faces[0].opr_data.size();
//    fmt::printf("Creating Laplace operator with id %d\n", id);

    for (size_t level = minLevel; level <= maxLevel; ++level)
    {

      for (Face& face : mesh.faces)
      {
        walberla::real_t local_stiffness_up[3][3];
        walberla::real_t local_stiffness_down[3][3];
        compute_local_stiffness<UFCOperator>(face, level, local_stiffness_up, UPWARD);
        compute_local_stiffness<UFCOperator>(face, level, local_stiffness_down, DOWNWARD);

        walberla::real_t* face_stencil = new walberla::real_t[7]();

        face_stencil[0] = local_stiffness_down[0][2] + local_stiffness_up[2][0];
        face_stencil[1] = local_stiffness_down[1][2] + local_stiffness_up[2][1];
        face_stencil[2] = local_stiffness_down[0][1] + local_stiffness_up[1][0];

        face_stencil[4] = local_stiffness_down[1][0] + local_stiffness_up[0][1];
        face_stencil[5] = local_stiffness_down[2][1] + local_stiffness_up[1][2];
        face_stencil[6] = local_stiffness_down[2][0] + local_stiffness_up[0][2];

        face_stencil[3] = local_stiffness_up[0][0] + local_stiffness_up[1][1] + local_stiffness_up[2][2]
                            + local_stiffness_down[0][0] + local_stiffness_down[1][1] + local_stiffness_down[2][2];

        if (level == minLevel)
        {
          face.opr_data.push_back(std::vector<walberla::real_t*>());
        }
        face.opr_data[memory_id].push_back(face_stencil);

//        fmt::printf("&face = %p\n", (void*) &fs.mesh.faces[0]);
//        fmt::print("face_stencil = {}\n", PointND<walberla::real_t, 7>(face_stencil));
      }

      for (Edge& edge : mesh.edges)
      {
        if (level == minLevel)
        {
          edge.opr_data.push_back(std::vector<walberla::real_t*>());
        }

        walberla::real_t* edge_stencil = new walberla::real_t[7]();
        edge.opr_data[memory_id].push_back(edge_stencil);

        walberla::real_t local_stiffness_up[3][3];
        walberla::real_t local_stiffness_down[3][3];
        // first face
        Face* face = edge.faces[0];
        compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_up, UPWARD);
        compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_down, DOWNWARD);

        size_t start_id = face->vertex_index(*edge.v0);
        size_t end_id = face->vertex_index(*edge.v1);
        size_t opposite_id = face->vertex_index(*face->get_vertex_opposite_to_edge(edge));

        edge_stencil[0] = local_stiffness_up[end_id][opposite_id] + local_stiffness_down[opposite_id][end_id];
        edge_stencil[1] = local_stiffness_up[start_id][opposite_id] + local_stiffness_down[opposite_id][start_id];

        edge_stencil[2] = local_stiffness_up[end_id][start_id];
        edge_stencil[4] = local_stiffness_up[start_id][end_id];

        edge_stencil[3] = local_stiffness_up[start_id][start_id] + local_stiffness_up[end_id][end_id] + local_stiffness_down[opposite_id][opposite_id];

        if (edge.faces.size() == 2)
        {
          // second face
          Face* face = edge.faces[1];
          compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_up, UPWARD);
          compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_down, DOWNWARD);

          size_t start_id = face->vertex_index(*edge.v0);
          size_t end_id = face->vertex_index(*edge.v1);
          size_t opposite_id = face->vertex_index(*face->get_vertex_opposite_to_edge(edge));

          edge_stencil[5] = local_stiffness_up[end_id][opposite_id] + local_stiffness_down[opposite_id][end_id];
          edge_stencil[6] = local_stiffness_up[start_id][opposite_id] + local_stiffness_down[opposite_id][start_id];

          edge_stencil[2] += local_stiffness_up[end_id][start_id];
          edge_stencil[4] += local_stiffness_up[start_id][end_id];

          edge_stencil[3] += local_stiffness_up[start_id][start_id] + local_stiffness_up[end_id][end_id] + local_stiffness_down[opposite_id][opposite_id];
        }
      }

      for (Vertex& vertex : mesh.vertices)
      {
        // allocate new level-vector if first level
        if (level == minLevel)
        {
          vertex.opr_data.push_back(std::vector<walberla::real_t*>());
        }

        walberla::real_t* vertex_stencil = new walberla::real_t[1 + vertex.edges.size()]();
        vertex.opr_data[memory_id].push_back(vertex_stencil);

        // iterate over adjacent faces
        for (Face* face : vertex.faces)
        {
          walberla::real_t local_stiffness[3][3];
          compute_local_stiffness<UFCOperator>(*face, level, local_stiffness, UPWARD);

          size_t v_i = face->vertex_index(vertex);

          std::vector<Edge*> adj_edges = face->adjacent_edges(vertex);

          // iterate over adjacent edges
          for (Edge* edge : adj_edges)
          {
            size_t edge_idx = vertex.edge_index(*edge) + 1;
            Vertex* vertex_j = edge->get_opposite_vertex(vertex);

            size_t v_j = face->vertex_index(*vertex_j);

            vertex_stencil[edge_idx] += local_stiffness[v_i][v_j];
          }

          vertex_stencil[0] += local_stiffness[v_i][v_i];
        }
      }
    }

  }

  ~P1Operator()
  {
    for (Vertex& v : mesh.vertices)
    {
      delete[] v.opr_data[memory_id][0];
    }

    for (Edge& e : mesh.edges)
    {
      delete[] e.opr_data[memory_id][0];
    }

    for (Face& f : mesh.faces)
    {
      delete[] f.opr_data[memory_id][0];
    }
  }

  template<size_t Level>
  void apply(const P1Function& src, P1Function& dst, DoFType flag, UpdateType updateType = Replace)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos<Level>(vertex, src.memory_id);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::apply<Level>(vertex, this->memory_id, src.memory_id, dst.memory_id, updateType);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices<Level>(edge, dst.memory_id);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos<Level>(edge, src.memory_id);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::apply<Level>(edge, this->memory_id, src.memory_id, dst.memory_id, updateType);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges<Level>(face, dst.memory_id);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::apply<Level>(face, this->memory_id, src.memory_id, dst.memory_id, updateType);
      }
    }
  }

  template<size_t Level>
  void smooth_gs(P1Function& dst, const P1Function& rhs, DoFType flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1Vertex::pull_halos<Level>(vertex, dst.memory_id);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::smooth_gs<Level>(vertex, this->memory_id, dst.memory_id, rhs.memory_id);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_vertices<Level>(edge, dst.memory_id);
      if (testFlag(edge.type, flag))
      {
        P1Edge::pull_halos<Level>(edge, dst.memory_id);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::smooth_gs<Level>(edge, this->memory_id, dst.memory_id, rhs.memory_id);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1Face::pull_edges<Level>(face, dst.memory_id);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::smooth_gs<Level>(face, this->memory_id, dst.memory_id, rhs.memory_id);
      }
    }
  }

  template<size_t Level>
  void printmatrix(const P1Function& src, DoFType flag = All)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      P1Vertex::pull_halos<Level>(vertex, src.memory_id);
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1Vertex::printmatrix<Level>(vertex, this->memory_id, src.memory_id);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1Edge::pull_halos<Level>(edge, src.memory_id);
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1Edge::printmatrix<Level>(edge, this->memory_id, src.memory_id);
      }
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1Face::printmatrix<Level>(face, this->memory_id, src.memory_id);
      }
    }
  }

};

typedef P1Operator<p1_diffusion_cell_integral_0_otherwise> P1LaplaceOperator;

typedef P1Operator<p1_div_cell_integral_0_otherwise> P1DivxOperator;
typedef P1Operator<p1_div_cell_integral_1_otherwise> P1DivyOperator;

typedef P1Operator<p1_divt_cell_integral_0_otherwise> P1DivTxOperator;
typedef P1Operator<p1_divt_cell_integral_1_otherwise> P1DivTyOperator;

typedef P1Operator<p1_mass_cell_integral_0_otherwise> P1MassOperator;

typedef P1Operator<p1_pspg_cell_integral_0_otherwise> P1PSPGOperator;

}

#endif /* P1OPERATOR_HPP */
