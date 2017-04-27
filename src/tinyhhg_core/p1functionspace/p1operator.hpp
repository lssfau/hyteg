#ifndef P1OPERATOR_HPP
#define P1OPERATOR_HPP

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/operator.hpp"

#include <fmt/format.h>

#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"

namespace hhg
{

void compute_micro_coords(const std::array<Point3D, 3>& macro_coords, size_t level, double coords[6])
{
  auto& x = macro_coords;
  double fac = std::pow(2.0, -(double) level);

  coords[0] = fac*(x[0][0]-x[0][0]);
  coords[1] = fac*(x[0][1]-x[0][1]);
  coords[2] = fac*(x[1][0]-x[0][0]);
  coords[3] = fac*(x[1][1]-x[0][1]);
  coords[4] = fac*(x[2][0]-x[0][0]);
  coords[5] = fac*(x[2][1]-x[0][1]);
}

template<class UFCOperator>
void compute_local_stiffness(const Face& face, size_t level, double local_stiffness[3][3])
{
  double A[9];
  double coords[6];
  compute_micro_coords(face.coords, level, coords);
  UFCOperator gen;
  gen.tabulate_tensor(A, NULL, coords, 0);

  // double local_stiffness[3][3] = { { t[0], t[1], t[2] }, {t[3], t[4], t[5]}, {t[6], t[7], t[8]} };
  // double local_stiffness[3][3] = { { t[0], t[3], t[6] }, {t[1], t[4], t[7]}, {t[2], t[5], t[8]} };

  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      local_stiffness[i][j] = A[3 * i + j];
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
    id = mesh.faces[0].opr_data.size();
    fmt::printf("Creating Laplace operator with id %d\n", id);

    for (size_t level = minLevel; level <= maxLevel; ++level)
    {

      for (Face& face : mesh.faces)
      {
        double local_stiffness[3][3];
        compute_local_stiffness<UFCOperator>(face, level, local_stiffness);

        double* face_stencil = new double[7]();

        face_stencil[0] = local_stiffness[0][2] + local_stiffness[2][0];
        face_stencil[1] = local_stiffness[1][2] + local_stiffness[2][1];
        face_stencil[2] = local_stiffness[1][0] + local_stiffness[0][1];

        face_stencil[4] = local_stiffness[0][1] + local_stiffness[1][0];
        face_stencil[5] = local_stiffness[2][1] + local_stiffness[1][2];
        face_stencil[6] = local_stiffness[2][0] + local_stiffness[0][2];

        face_stencil[3] = 2.0 * (local_stiffness[0][0] + local_stiffness[1][1] + local_stiffness[2][2]);

        if (level == minLevel)
        {
          face.opr_data.push_back(std::vector<double*>());
        }
        face.opr_data[id].push_back(face_stencil);

        // fmt::printf("&face = %p\n", (void*) &fs.mesh.faces[0]);
        // fmt::print("face_stencil = {}\n", PointND<double, 7>(face_stencil));

        for (size_t i = 0; i < 3; ++i)
        {
          Edge& edge = *face.edges[i];
          size_t face_idx = edge.face_index(face);

          double* edge_stencil;

          if (level == minLevel && edge.opr_data.size() == id)
          {
            edge.opr_data.push_back(std::vector<double*>());
          }

          if (edge.opr_data[id].size() == level-minLevel)
          {
            edge_stencil = new double[7]();
            edge.opr_data[id].push_back(edge_stencil);
          }
          else
          {
            edge_stencil = edge.opr_data[id][level-minLevel];
          }

          std::pair<size_t, size_t> base;

          if (i == 0)
          {
            base = std::make_pair(0, 1);
          }
          else if(i == 1)
          {
            base = std::make_pair(1, 2);
          }
          else
          {
            base = std::make_pair(0, 2);
          }

          if (face_idx == 0)
          {
            if (face.edge_orientation[i] == -1)
            {
              edge_stencil[0] = 2 * local_stiffness[(i+0) % 3][(i+2) % 3];
              edge_stencil[1] = 2 * local_stiffness[(i+1) % 3][(i+2) % 3];
            }
            else
            {
              edge_stencil[1] = 2 * local_stiffness[(i+0) % 3][(i+2) % 3];
              edge_stencil[0] = 2 * local_stiffness[(i+1) % 3][(i+2) % 3];
            }

            edge_stencil[3] -= edge_stencil[0] + edge_stencil[1];
          }
          else
          {
            if (face.edge_orientation[i] == 1)
            {
              edge_stencil[5] = 2 * local_stiffness[(i+1) % 3][(i+2) % 3];
              edge_stencil[6] = 2 * local_stiffness[(i+0) % 3][(i+2) % 3];
            }
            else
            {
              edge_stencil[6] = 2 * local_stiffness[(i+1) % 3][(i+2) % 3];
              edge_stencil[5] = 2 * local_stiffness[(i+0) % 3][(i+2) % 3];
            }

            edge_stencil[3] -= edge_stencil[5] + edge_stencil[6];
          }

          double tmp = local_stiffness[base.first][base.second];

          edge_stencil[2] += tmp;
          edge_stencil[4] += tmp;

          edge_stencil[3] -= 2 * tmp;

          // fmt::print("edge_stencil = {}\n", PointND<double, 7>(edge_stencil));
        }

        // create vertex stencil
        // for (size_t v = 0; v < 3; ++v)
        // {
        //   Vertex& vertex = *face.vertices[v];

        //   double* vertex_stencil = NULL;

        //   if (level == minLevel && vertex.opr_data.size() == id)
        //   {
        //     vertex.opr_data.push_back(std::vector<double*>());
        //   }

        //   if (vertex.opr_data[id].size() == level-minLevel)
        //   {
        //     vertex_stencil = new double[1 + vertex.edges.size()];
        //     vertex.opr_data[id].push_back(vertex_stencil);
        //   }
        //   else
        //   {
        //     vertex_stencil = vertex.opr_data[id][level-minLevel];
        //   }

        //   std::vector<Edge*> local_edges;

        //   for (size_t i = 0; i < 3; ++i)
        //   {
        //     if (&vertex == face.edges[i]->v0 || &vertex == face.edges[i]->v1)
        //     {
        //       local_edges.push_back(face.edges[i]);
        //     }
        //   }

        //   size_t e1i = vertex.edge_index(*local_edges[0]);
        //   size_t e2i = vertex.edge_index(*local_edges[1]);

        //   std::vector<size_t> idx;

        //   for (size_t e = 0; e < 2; ++e)
        //   {
        //     if (&vertex != local_edges[e]->v0)
        //     {
        //       idx.push_back(face.vertex_index(*local_edges[e]->v0));
        //     }
        //     else
        //     {
        //       idx.push_back(face.vertex_index(*local_edges[e]->v1));
        //     }
        //   }

        //   double val1 = local_stiffness[v][idx[0]];
        //   double val2 = local_stiffness[v][idx[1]];

        //   vertex_stencil[e1i] += val1;
        //   vertex_stencil[e2i] += val2;

        //   vertex_stencil[0] -= val1 + val2;
        // }
      }

      for (Vertex& vertex : mesh.vertices)
      {

        // allocate new level-vector if first level
        if (level == minLevel)
        {
          vertex.opr_data.push_back(std::vector<double*>());
        }

        double* vertex_stencil = new double[1 + vertex.edges.size()]();
        vertex.opr_data[id].push_back(vertex_stencil);

        // iterate over adjacent faces
        for (Face* face : vertex.faces)
        {
          double local_stiffness[3][3];
          compute_local_stiffness<UFCOperator>(*face, level, local_stiffness);

          size_t v_i = face->vertex_index(vertex);

          std::vector<Edge*> adj_edges = face->adjacent_edges(vertex);

          // iterate over adjacent edges
          for (Edge* edge : adj_edges)
          {
            size_t edge_idx = vertex.edge_index(*edge);
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
      delete[] v.opr_data[id][0];
    }

    for (Edge& e : mesh.edges)
    {
      delete[] e.opr_data[id][0];
    }

    for (Face& f : mesh.faces)
    {
      delete[] f.opr_data[id][0];
    }
  }

};

typedef P1Operator<p1_diffusion_cell_integral_0_otherwise> P1LaplaceOperator;

typedef P1Operator<p1_div_cell_integral_0_otherwise> P1DivxOperator;
typedef P1Operator<p1_div_cell_integral_1_otherwise> P1DivyOperator;

typedef P1Operator<p1_divt_cell_integral_0_otherwise> P1DivTxOperator;
typedef P1Operator<p1_divt_cell_integral_1_otherwise> P1DivTyOperator;

typedef P1Operator<p1_mass_cell_integral_0_otherwise> P1MassOperator;

}

#endif /* P1OPERATOR_HPP */
