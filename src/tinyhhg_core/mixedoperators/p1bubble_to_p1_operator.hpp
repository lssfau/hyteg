#pragma once

#include <fmt/format.h>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/operator.hpp"

#include "tinyhhg_core/p1bubblefunctionspace/generated/p1bubble_div.h"

namespace hhg
{

namespace P1BubbleToP1Vertex
{
inline void apply(Vertex& vertex, size_t opr_id, size_t src_id, size_t dst_id, size_t level, UpdateType update)
{
  auto& stencil_stack = P1Bubble::getVertexStencilMemory(vertex, opr_id)->data[level];
  auto& src = P1Bubble::getVertexFunctionMemory(vertex, src_id)->data[level];
  auto& dst = P1Bubble::getVertexFunctionMemory(vertex, dst_id)->data[level];

  // apply first stencil to vertex dof
  auto& opr_data = stencil_stack[0];

  if (update == Replace) {
    dst[0] = opr_data[0] * src[0];
  }
  else if (update == Add) {
    dst[0] += opr_data[0] * src[0];
  }

  for (size_t i = 0; i < vertex.edges.size() + vertex.faces.size(); ++i)
  {
    dst[0] += opr_data[i+1] * src[i+1];
  }
}

inline void saveOperator(size_t level, Vertex& vertex, std::ostream& out, size_t opr_id, size_t src_id, size_t dst_id, DoFType flag)
{
  auto& stencil_stack = P1Bubble::getVertexStencilMemory(vertex, opr_id)->data[level];
  auto& src = P1Bubble::getVertexFunctionMemory(vertex, src_id)->data[level];
  auto& dst = P1Bubble::getVertexFunctionMemory(vertex, dst_id)->data[level];

  // apply first stencil to vertex dof
  auto& opr_data = stencil_stack[0];

  out << fmt::format("{}\t{}\t{}\n", dst[0], src[0], opr_data[0]);

  for (size_t i = 0; i < vertex.edges.size() + vertex.faces.size(); ++i)
  {
    out << fmt::format("{}\t{}\t{}\n", dst[0], src[i + 1], opr_data[i + 1]);
  }
}
}

namespace P1BubbleToP1Edge
{
template<size_t Level>
inline void apply_tmpl(Edge& edge, size_t opr_id, size_t src_id, size_t dst_id, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto& edge_vertex_stencil = P1Bubble::getEdgeStencilMemory(edge, opr_id)->data[Level];
  auto& src = P1Bubble::getEdgeFunctionMemory(edge, src_id)->data[Level];
  auto& dst = P1Bubble::getEdgeFunctionMemory(edge, dst_id)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    tmp = edge_vertex_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_C] * src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)];

    for (auto neighbor : P1BubbleEdge::EdgeCoordsVertex::neighbors_edge)
    {
      tmp += edge_vertex_stencil[neighbor] * src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)];
    }

    for (auto neighbor : P1BubbleEdge::EdgeCoordsVertex::neighbors_south)
    {
      tmp += edge_vertex_stencil[neighbor] * src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)];
    }

    if (edge.faces.size() == 2)
    {
      for (auto neighbor : P1BubbleEdge::EdgeCoordsVertex::neighbors_north)
      {
        tmp += edge_vertex_stencil[neighbor] * src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)];
      }
    }

    if (update == Replace) {
      dst[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)] = tmp;
    } else if (update == Add) {
      dst[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)] += tmp;
    }
  }
}

SPECIALIZE(void, apply_tmpl, apply)

template<size_t Level>
inline void saveOperator_tmpl(Edge& edge, std::ostream& out, size_t opr_id, size_t src_id, size_t dst_id, DoFType flag)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto& edge_vertex_stencil = P1Bubble::getEdgeStencilMemory(edge, opr_id)->data[Level];
  auto& src = P1Bubble::getEdgeFunctionMemory(edge, src_id)->data[Level];
  auto& dst = P1Bubble::getEdgeFunctionMemory(edge, dst_id)->data[Level];

  for (size_t i = 1; i < rowsize-1; ++i)
  {
    out << fmt::format("{}\t{}\t{}\n", dst[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)], src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)], edge_vertex_stencil[P1BubbleEdge::EdgeCoordsVertex::VERTEX_C]);

    for (auto neighbor : P1BubbleEdge::EdgeCoordsVertex::neighbors_edge)
    {
      out << fmt::format("{}\t{}\t{}\n", dst[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)], src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)], edge_vertex_stencil[neighbor]);
    }

    for (auto neighbor : P1BubbleEdge::EdgeCoordsVertex::neighbors_south)
    {
      out << fmt::format("{}\t{}\t{}\n", dst[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)], src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)], edge_vertex_stencil[neighbor]);
    }

    if (edge.faces.size() == 2)
    {
      for (auto neighbor : P1BubbleEdge::EdgeCoordsVertex::neighbors_north)
      {
        out << fmt::format("{}\t{}\t{}\n", dst[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)], src[P1BubbleEdge::EdgeCoordsVertex::index<Level>(i, neighbor)], edge_vertex_stencil[neighbor]);
      }
    }
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)
}

namespace P1BubbleToP1Face
{
template<size_t Level>
inline void apply_tmpl(Face& face, size_t opr_id, size_t src_id, size_t dst_id, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = P1Bubble::getFaceStencilMemory(face, opr_id)->data[Level];

  auto& face_vertex_stencil = opr_data[0];

  auto& src = P1Bubble::getFaceFunctionMemory(face, src_id)->data[Level];
  auto& dst = P1Bubble::getFaceFunctionMemory(face, dst_id)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      tmp = face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_C] * src[P1BubbleFace::CoordsVertex::index<Level>(i, j, P1BubbleFace::CoordsVertex::VERTEX_C)];

      for (auto neighbor : P1BubbleFace::CoordsVertex::neighbors)
      {
        tmp += face_vertex_stencil[neighbor] * src[P1BubbleFace::CoordsVertex::index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[P1BubbleFace::CoordsVertex::index<Level>(i, j, P1BubbleFace::CoordsVertex::VERTEX_C)] = tmp;
      } else if (update == Add) {
        dst[P1BubbleFace::CoordsVertex::index<Level>(i, j, P1BubbleFace::CoordsVertex::VERTEX_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)

template<size_t Level>
inline void saveOperator_tmpl(Face& face, std::ostream& out, size_t opr_id, size_t src_id, size_t dst_id, DoFType flag)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = P1Bubble::getFaceStencilMemory(face, opr_id)->data[Level];

  auto& face_vertex_stencil = opr_data[0];

  auto& src = P1Bubble::getFaceFunctionMemory(face, src_id)->data[Level];
  auto& dst = P1Bubble::getFaceFunctionMemory(face, dst_id)->data[Level];

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      out << fmt::format("{}\t{}\t{}\n", dst[P1BubbleFace::CoordsVertex::index<Level>(i, j, P1BubbleFace::CoordsVertex::VERTEX_C)], src[P1BubbleFace::CoordsVertex::index<Level>(i, j, P1BubbleFace::CoordsVertex::VERTEX_C)], face_vertex_stencil[P1BubbleFace::CoordsVertex::VERTEX_C]);

      for (auto neighbor : P1BubbleFace::CoordsVertex::neighbors)
      {
        out << fmt::format("{}\t{}\t{}\n", dst[P1BubbleFace::CoordsVertex::index<Level>(i, j, P1BubbleFace::CoordsVertex::VERTEX_C)], src[P1BubbleFace::CoordsVertex::index<Level>(i, j, neighbor)], face_vertex_stencil[neighbor]);
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)
}

template<class UFCOperator>
class P1BubbleToP1Operator : public Operator
{
public:
  P1BubbleToP1Operator(Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
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

        auto& vertex_stencil_stack = P1Bubble::getVertexStencilMemory(vertex, memory_id)->addlevel(level, vertex.edges.size(), vertex.faces.size());

        // build vertex stencil

        size_t f = 1;
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

            vertex_stencil_stack[0][edge_idx] += local_stiffness[v_i][v_j];

            vertex_stencil_stack[f][v_j+1] = local_stiffness[3][v_j];
          }

          vertex_stencil_stack[f][v_i+1] = local_stiffness[3][v_i];
          vertex_stencil_stack[f][0] = local_stiffness[3][3];

          size_t f_i = vertex.face_index(*face) + 1 + vertex.edges.size();
          vertex_stencil_stack[0][f_i] = local_stiffness[v_i][3];

          vertex_stencil_stack[0][0] += local_stiffness[v_i][v_i];
          ++f;
        }
      }
    }

  }

  ~P1BubbleToP1Operator()
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
//      if (testFlag(vertex.type, flag))
      {
        P1BubbleVertex::pull_halos(vertex, src.memory_id, level);
      }
    }

    for (Vertex& vertex : mesh.vertices)
    {
      if (vertex.rank == rank && testFlag(vertex.type, flag))
      {
        P1BubbleToP1Vertex::apply(vertex, this->memory_id, src.memory_id, dst.memory_id, level, updateType);
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
        P1BubbleToP1Edge::apply(level, edge, this->memory_id, src.memory_id, dst.memory_id, updateType);
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
        P1BubbleToP1Face::apply(level, face, this->memory_id, src.memory_id, dst.memory_id, updateType);
      }
    }
  }

  void save(const P1BubbleFunction& src, const P1BubbleFunction& dst, std::ostream& out, size_t level, DoFType flag)
  {
    for (Vertex& vertex : mesh.vertices)
    {
      P1BubbleToP1Vertex::saveOperator(level, vertex, out, this->memory_id, src.memory_id, dst.memory_id, flag);
    }

    for (Edge& edge : mesh.edges)
    {
      P1BubbleToP1Edge::saveOperator(level, edge, out, this->memory_id, src.memory_id, dst.memory_id, flag);
    }

    for (Face& face : mesh.faces)
    {
      P1BubbleToP1Face::saveOperator(level, face, out, this->memory_id, src.memory_id, dst.memory_id, flag);
    }
  }
};

typedef P1BubbleToP1Operator<p1bubble_div_cell_integral_0_otherwise> P1BubbleToP1DivxOperator;
typedef P1BubbleToP1Operator<p1bubble_div_cell_integral_1_otherwise> P1BubbleToP1DivyOperator;

}
