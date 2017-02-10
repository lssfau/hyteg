#ifndef P1OPERATOR_HPP
#define P1OPERATOR_HPP

#include <array>
#include <types/pointnd.hpp>
#include <p1functionspace/p1functionspace.hpp>

#include <fmt/format.h>

namespace hhg
{

inline std::array<Point2D, 3> grad_shape(std::array<Point3D, 3>& coords)
{
  double x_1 = coords[0][0];
  double x_2 = coords[1][0];
  double x_3 = coords[2][0];

  double y_1 = coords[0][1];
  double y_2 = coords[1][1];
  double y_3 = coords[2][1];

  Point2D a({(y_2 - y_3)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), -(x_2 - x_3)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2)});
  Point2D b({-(y_1 - y_3)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), (x_1 - x_3)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2)});
  Point2D c({(y_1 - y_2)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), -(x_1 - x_2)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2)});

  return { a, b, c };
}

class P1LaplaceOperator
{
public:
  P1LaplaceOperator(P1FunctionSpace& _fs, size_t _minLevel, size_t _maxLevel)
    : fs(_fs), minLevel(_minLevel), maxLevel(_maxLevel)
  {
    id = fs.mesh.faces[0].opr_data.size();
    fmt::printf("Creating Laplace operator with id %d\n", id);

    for (Face& face : fs.mesh.faces)
    {
      std::array<Point2D, 3> gshape = grad_shape(face.coords);

      double local_stiffness[3][3];

      for (size_t i = 0; i < 3; ++i)
      {
        for (size_t j = 0; j < 3; ++j)
        {
          local_stiffness[i][j] = face.area * gshape[i].dot(gshape[j]);
        }
      }

      double* face_stencil = new double[7];

      face_stencil[0] = 2 * local_stiffness[0][2];
      face_stencil[1] = 2 * local_stiffness[1][2];
      face_stencil[2] = 2 * local_stiffness[1][0];

      face_stencil[4] = 2 * local_stiffness[0][1];
      face_stencil[5] = 2 * local_stiffness[2][1];
      face_stencil[6] = 2 * local_stiffness[2][0];

      face_stencil[3] = 0.0;
      double sum = 0.0;
      for (size_t i = 0; i < 7; ++i)
      {
        sum += face_stencil[i];
      }
      face_stencil[3] = -sum;

      face.opr_data.push_back(std::vector<double*>());
      face.opr_data[id].push_back(face_stencil);

      // fmt::printf("&face = %p\n", (void*) &fs.mesh.faces[0]);
      // fmt::print("face_stencil = {}\n", PointND<double, 7>(face_stencil));

      for (size_t i = 0; i < 3; ++i)
      {
        Edge& edge = *face.edges[i];
        size_t face_idx = edge.face_index(face);

        double* edge_stencil;

        if (edge.opr_data.size() == id)
        {
          edge_stencil = new double[7];
          edge.opr_data.push_back(std::vector<double*>());
          edge.opr_data[id].push_back(edge_stencil);
        }
        else
        {
          edge_stencil = edge.opr_data[id][0];
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
      for (size_t v = 0; v < 3; ++v)
      {
        Vertex& vertex = *face.vertices[v];

        double* vertex_stencil = NULL;

        if (vertex.opr_data.size() == id)
        {
          vertex_stencil = new double[1 + vertex.edges.size()];
          vertex.opr_data.push_back(std::vector<double*>());
          vertex.opr_data[id].push_back(vertex_stencil);
        }
        else
        {
          vertex_stencil = vertex.opr_data[id][0];
        }

        std::vector<Edge*> local_edges;

        for (size_t i = 0; i < 3; ++i)
        {
          if (&vertex == face.edges[i]->v0 || &vertex == face.edges[i]->v1)
          {
            local_edges.push_back(face.edges[i]);
          }
        }

        size_t e1i = vertex.edge_index(*local_edges[0]);
        size_t e2i = vertex.edge_index(*local_edges[1]);

        std::vector<size_t> idx;

        for (size_t e = 0; e < 2; ++e)
        {
          if (&vertex != local_edges[e]->v0)
          {
            idx.push_back(face.vertex_index(*local_edges[e]->v0));
          }
          else
          {
            idx.push_back(face.vertex_index(*local_edges[e]->v1));
          }
        }

        double val1 = local_stiffness[v][idx[0]];
        double val2 = local_stiffness[v][idx[1]];

        vertex_stencil[e1i] += val1;
        vertex_stencil[e2i] += val2;

        vertex_stencil[0] -= val1 + val2;
      }
    }

  }

  ~P1LaplaceOperator()
  {
    for (Vertex& v : fs.mesh.vertices)
    {
      delete[] v.opr_data[id][0];
    }

    for (Edge& e : fs.mesh.edges)
    {
      delete[] e.opr_data[id][0];
    }

    for (Face& f : fs.mesh.faces)
    {
      delete[] f.opr_data[id][0];
    }
  }

  size_t id;
  P1FunctionSpace& fs;
  size_t minLevel;
  size_t maxLevel;

};

}

#endif /* P1OPERATOR_HPP */