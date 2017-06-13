#pragma once

#include <fmt/format.h>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/operator.hpp"

#include "tinyhhg_core/p1bubblefunctionspace/generated/p1bubble_diffusion.h"

#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblefaceindex.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubbleedge.hpp"

namespace hhg
{

namespace P1Bubble {

enum ElementType {
  GRAY,
  BLUE
};

void compute_micro_coords(const Face &face, size_t level, real_t coords[6], ElementType element_type) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D d0 = face.edge_orientation[0] * face.edges[0]->direction / (rowsize - 1);
  Point3D d2 = -face.edge_orientation[2] * face.edges[2]->direction / (rowsize - 1);

  real_t orientation = 1.0;

  if (element_type == BLUE) {
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
void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[4][4], ElementType element_type) {
  real_t A[16];
  real_t coords[6];
  compute_micro_coords(face, level, coords, element_type);
  UFCOperator gen;
  gen.tabulate_tensor(A, NULL, coords, 0);

  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      local_stiffness[i][j] = A[4 * j + i];
    }
  }
}
}

template<class UFCOperator>
class P1BubbleOperator : public Operator
{
public:
  P1BubbleOperator(Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : Operator(_mesh, _minLevel, _maxLevel)
  {
    for (Vertex& v : mesh.vertices)
    {
      if (v.rank == rank)
      {
        memory_id = v.memory.size();
        break;
      }
    }

    //if (id == -1)
    //{
      for (Edge& e : mesh.edges)
      {
        if (e.rank == rank)
        {
          if (memory_id == -1)
          {
            memory_id = e.memory.size();
            break;
          }
          else if (memory_id != e.memory.size())
            WALBERLA_LOGLEVEL_WARNING("ID of Vertex and Edge are not the same");

        }
      }
    //}

    //if (id == -1)
    //{
      for (Face& f : mesh.faces)
      {
        if (f.rank == rank)
        {

          if (memory_id == -1)
          {
            memory_id = f.memory.size();
            break;
          }
          else if (memory_id != f.memory.size())
            WALBERLA_LOGLEVEL_WARNING("ID of Vertex and Face are not the same");
        }
      }
    //}

    if (memory_id == -1)
    {
      WALBERLA_ABORT("Could not determine memory id of P1 operator");
    }

    WALBERLA_LOG_DEVEL("Created P1Bubble Operator with ID " + std::to_string(memory_id))


    for (size_t level = minLevel; level <= maxLevel; ++level)
    {

      for (Face& face : mesh.faces)
      {
        if (face.rank != rank)
        {
          continue;
        }

        if (level == minLevel)
        {
          face.memory.push_back(new FaceP1BubbleStencilMemory());
        }

        auto& face_stencil_stack = P1Bubble::getFaceStencilMemory(face, memory_id)->addlevel(level);

        auto& face_vertex_stencil = face_stencil_stack[0];
        auto& face_gray_stencil = face_stencil_stack[1];
        auto& face_blue_stencil = face_stencil_stack[2];

        real_t local_stiffness_up[4][4];
        real_t local_stiffness_down[4][4];
        P1Bubble::compute_local_stiffness<UFCOperator>(face, level, local_stiffness_up, P1Bubble::GRAY);
        P1Bubble::compute_local_stiffness<UFCOperator>(face, level, local_stiffness_down, P1Bubble::BLUE);

        face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_S] = local_stiffness_down[0][2] + local_stiffness_up[2][0];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_SE] = local_stiffness_down[1][2] + local_stiffness_up[2][1];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_W] = local_stiffness_down[0][1] + local_stiffness_up[1][0];

        face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_E] = local_stiffness_down[1][0] + local_stiffness_up[0][1];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_NW] = local_stiffness_down[2][1] + local_stiffness_up[1][2];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_N] = local_stiffness_down[2][0] + local_stiffness_up[0][2];

        face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_C] = local_stiffness_up[0][0] + local_stiffness_up[1][1] + local_stiffness_up[2][2]
                            + local_stiffness_down[0][0] + local_stiffness_down[1][1] + local_stiffness_down[2][2];

        face_vertex_stencil[P1BubbleFace::CoordsVertex::CELL_GRAY_SE] = local_stiffness_up[2][3];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::CELL_GRAY_NW] = local_stiffness_up[1][3];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::CELL_GRAY_NE] = local_stiffness_up[0][3];

        face_vertex_stencil[P1BubbleFace::CoordsVertex::CELL_BLUE_SW] = local_stiffness_down[0][3];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::CELL_BLUE_SE] = local_stiffness_down[1][3];
        face_vertex_stencil[P1BubbleFace::CoordsVertex::CELL_BLUE_NW] = local_stiffness_down[2][3];

        face_gray_stencil[P1BubbleFace::CoordsCellGray::VERTEX_SW] = local_stiffness_up[3][0];
        face_gray_stencil[P1BubbleFace::CoordsCellGray::VERTEX_SE] = local_stiffness_up[3][1];
        face_gray_stencil[P1BubbleFace::CoordsCellGray::VERTEX_NW] = local_stiffness_up[3][2];
        face_gray_stencil[P1BubbleFace::CoordsCellGray::CELL_GRAY_C] = local_stiffness_up[3][3];

        face_blue_stencil[P1BubbleFace::CoordsCellBlue::VERTEX_SE] = local_stiffness_down[3][2];
        face_blue_stencil[P1BubbleFace::CoordsCellBlue::VERTEX_NW] = local_stiffness_down[3][1];
        face_blue_stencil[P1BubbleFace::CoordsCellBlue::VERTEX_NE] = local_stiffness_down[3][0];
        face_blue_stencil[P1BubbleFace::CoordsCellBlue::CELL_BLUE_C] = local_stiffness_down[3][3];
      }

      for (Edge& edge : mesh.edges)
      {
        if (edge.rank != rank)
        {
          continue;
        }

        if (level == minLevel)
        {
          edge.memory.push_back(new EdgeP1BubbleStencilMemory());
        }
        //WALBERLA_LOG_DEVEL("Edge.memory.size() = " + std::to_string(edge.memory.size()));

        auto& edge_stencil = P1Bubble::getEdgeStencilMemory(edge, memory_id)->addlevel(level);

        real_t local_stiffness_gray[4][4];
        real_t local_stiffness_blue[4][4];
        // first face
        Face* face = edge.faces[0];
        P1Bubble::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_gray, P1Bubble::GRAY);
        P1Bubble::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_blue, P1Bubble::BLUE);

        size_t start_id = face->vertex_index(*edge.v0);
        size_t end_id = face->vertex_index(*edge.v1);
        size_t opposite_id = face->vertex_index(*face->get_vertex_opposite_to_edge(edge));

        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_S] = local_stiffness_gray[end_id][opposite_id] + local_stiffness_blue[opposite_id][end_id];
        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_SE] = local_stiffness_gray[start_id][opposite_id] + local_stiffness_blue[opposite_id][start_id];

        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_W] = local_stiffness_gray[end_id][start_id];
        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_E] = local_stiffness_gray[start_id][end_id];

        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_C] = local_stiffness_gray[start_id][start_id] + local_stiffness_gray[end_id][end_id] + local_stiffness_blue[opposite_id][opposite_id];

        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SW] = local_stiffness_gray[end_id][3];
        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::CELL_BLUE_SE] = local_stiffness_blue[opposite_id][3];
        edge_stencil[P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SE] = local_stiffness_gray[start_id][3];

        if (edge.faces.size() == 2)
        {
          // second face
          Face* face = edge.faces[1];
          P1Bubble::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_gray, P1Bubble::GRAY);
          P1Bubble::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_blue, P1Bubble::BLUE);

          size_t start_id = face->vertex_index(*edge.v0);
          size_t end_id = face->vertex_index(*edge.v1);
          size_t opposite_id = face->vertex_index(*face->get_vertex_opposite_to_edge(edge));

          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_NW] = local_stiffness_gray[end_id][opposite_id] + local_stiffness_blue[opposite_id][end_id];
          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_N] = local_stiffness_gray[start_id][opposite_id] + local_stiffness_blue[opposite_id][start_id];

          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_W] += local_stiffness_gray[end_id][start_id];
          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_E] += local_stiffness_gray[start_id][end_id];

          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_C] += local_stiffness_gray[start_id][start_id] + local_stiffness_gray[end_id][end_id] + local_stiffness_blue[opposite_id][opposite_id];

          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NW] = local_stiffness_gray[end_id][3];
          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::CELL_BLUE_NW] = local_stiffness_blue[opposite_id][3];
          edge_stencil[P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NE] = local_stiffness_gray[start_id][3];
        }
      }

      for (Vertex& vertex : mesh.vertices)
      {
        if (vertex.rank != rank)
        {
          continue;
        }

        // allocate new level-vector if first level
        if (level == minLevel)
        {
          vertex.memory.push_back(new VertexP1BubbleStencilMemory());
        }

        auto& vertex_stencil = P1Bubble::getVertexStencilMemory(vertex, memory_id)->addlevel(level, vertex.edges.size());

        // iterate over adjacent faces
        for (Face* face : vertex.faces)
        {
          real_t local_stiffness[4][4];
          P1Bubble::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness, P1Bubble::GRAY);

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

  ~P1BubbleOperator()
  {
    for (Vertex& v : mesh.vertices)
    {
      if (v.rank == rank) {
        delete v.memory[memory_id];
      }
    }

    for (Edge& e : mesh.edges)
    {
      if (e.rank == rank) {
        delete e.memory[memory_id];
      }
    }

    for (Face& f : mesh.faces)
    {
      if (f.rank == rank) {
        delete f.memory[memory_id];
      }
    }
  }

  void apply(const P1BubbleFunction& src, P1BubbleFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      if (testFlag(vertex.type, flag))
      {
        P1BubbleVertex::pull_halos(vertex, src.memory_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1BubbleVertex::apply(vertex, this->memory_id, src.memory_id, dst.memory_id, level, updateType);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleEdge::pull_vertices(edge, dst.memory_id, level);
      if (testFlag(edge.type, flag))
      {
        P1BubbleEdge::pull_halos(edge, src.memory_id, level);
      }
    }

    for (Edge& edge : mesh.edges)
    {
      if (edge.rank == rank && testFlag(edge.type, flag))
      {
        P1BubbleEdge::apply(edge, this->memory_id, src.memory_id, dst.memory_id, level, updateType);
      }
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleFace::pull_edges(face, dst.memory_id, level);
    }

    for (Face& face : mesh.faces)
    {
      if (face.rank == rank && testFlag(face.type, flag))
      {
        P1BubbleFace::apply(level, face, this->memory_id, src.memory_id, dst.memory_id, updateType);
      }
    }
  }

//  void smooth_gs(P1BubbleFunction& dst, const P1BubbleFunction& rhs, size_t level, DoFType flag)
//  {
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (testFlag(vertex.type, flag))
//      {
//        P1BubbleVertex::pull_halos(vertex, dst.memory_id, level);
//      }
//    }
//
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (vertex.rank == rank && testFlag(vertex.type, flag))
//      {
//        P1BubbleVertex::smooth_gs(vertex, this->memory_id, dst.memory_id, rhs.memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      P1BubbleEdge::pull_vertices(edge, dst.memory_id, level);
//      if (testFlag(edge.type, flag))
//      {
//        P1BubbleEdge::pull_halos(edge, dst.memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      if (edge.rank == rank && testFlag(edge.type, flag))
//      {
//        P1BubbleEdge::smooth_gs(edge, this->memory_id, dst.memory_id, rhs.memory_id, level);
//      }
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      P1BubbleFace::pull_edges(face, dst.memory_id, level);
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      if (face.rank == rank && testFlag(face.type, flag))
//      {
//        P1BubbleFace::smooth_gs(level, face, this->memory_id, dst.memory_id, rhs.memory_id);
//      }
//    }
//  }
//
//  void printmatrix(const P1BubbleFunction& src, size_t level, DoFType flag = All)
//  {
//    for (Vertex& vertex : mesh.vertices)
//    {
//      P1BubbleVertex::pull_halos(vertex, src.memory_id, level);
//    }
//
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (vertex.rank == rank && testFlag(vertex.type, flag))
//      {
//        P1BubbleVertex::printmatrix(vertex, this->memory_id, src.memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      P1BubbleEdge::pull_halos(edge, src.memory_id, level);
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      if (edge.rank == rank && testFlag(edge.type, flag))
//      {
//        P1BubbleEdge::printmatrix(edge, this->memory_id, src.memory_id, level);
//      }
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      if (face.rank == rank && testFlag(face.type, flag))
//      {
//        P1BubbleFace::printmatrix(face, this->memory_id, src.memory_id, level);
//      }
//    }
//  }

};

typedef P1BubbleOperator<p1bubble_diffusion_cell_integral_0_otherwise> P1BubbleLaplaceOperator;
//
//typedef P1Operator<p1_div_cell_integral_0_otherwise> P1DivxOperator;
//typedef P1Operator<p1_div_cell_integral_1_otherwise> P1DivyOperator;
//
//typedef P1Operator<p1_divt_cell_integral_0_otherwise> P1DivTxOperator;
//typedef P1Operator<p1_divt_cell_integral_1_otherwise> P1DivTyOperator;
//
//typedef P1Operator<p1_mass_cell_integral_0_otherwise> P1MassOperator;
//
//typedef P1Operator<p1_pspg_cell_integral_0_otherwise> P1PSPGOperator;

}