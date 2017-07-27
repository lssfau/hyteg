#ifndef P1OPERATOR_HPP
#define P1OPERATOR_HPP

#include <fmt/format.h>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/operator.hpp"

#include "P1DataHandling.hpp"

#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"

#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg
{

namespace P1Space {
enum ElementType {
  UPWARD,
  DOWNWARD
};

void compute_micro_coords(const Face &face, size_t level, real_t coords[6], ElementType element_type) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D d0 = (face.coords[1] - face.coords[0]) / walberla::real_c((rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0]) / walberla::real_c((rowsize - 1));

  real_t orientation = 1.0;

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
void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[3][3], ElementType element_type) {
  real_t A[9];
  real_t coords[6];
  compute_micro_coords(face, level, coords, element_type);
  UFCOperator gen;
  gen.tabulate_tensor(A, NULL, coords, 0);

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      local_stiffness[i][j] = A[3 * j + i];
    }
  }
}
}

template<class UFCOperator>
class P1Operator : public Operator
{
public:
  P1Operator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  {
    FaceP1StencilMemoryDataHandling faceP1StencilMemoryDataHandling(minLevel_, maxLevel_);
    EdgeP1StencilMemoryDataHandling edgeP1StencilMemoryDataHandling(minLevel_, maxLevel_);
    VertexP1StencilMemoryDataHandling vertexP1StencilMemoryDataHandling(minLevel_, maxLevel_);
    faceStencilID_ = storage->addFaceData(faceP1StencilMemoryDataHandling, "P1OperatorFaceStencil");
    edgeStencilID_ = storage->addEdgeData(edgeP1StencilMemoryDataHandling, "P1OperatorEdgeStencil");
    vertexStencilID_ = storage->addVertexData(vertexP1StencilMemoryDataHandling, "P1OperatorVertexStencil");

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        auto& face_stencil = face.getData(faceStencilID_)->data[level];

        real_t local_stiffness_up[3][3];
        real_t local_stiffness_down[3][3];
        P1Space::compute_local_stiffness<UFCOperator>(face, level, local_stiffness_up, P1Space::UPWARD);
        P1Space::compute_local_stiffness<UFCOperator>(face, level, local_stiffness_down, P1Space::DOWNWARD);

        face_stencil[0] = local_stiffness_down[0][2] + local_stiffness_up[2][0];
        face_stencil[1] = local_stiffness_down[1][2] + local_stiffness_up[2][1];
        face_stencil[2] = local_stiffness_down[0][1] + local_stiffness_up[1][0];

        face_stencil[4] = local_stiffness_down[1][0] + local_stiffness_up[0][1];
        face_stencil[5] = local_stiffness_down[2][1] + local_stiffness_up[1][2];
        face_stencil[6] = local_stiffness_down[2][0] + local_stiffness_up[0][2];

        face_stencil[3] = local_stiffness_up[0][0] + local_stiffness_up[1][1] + local_stiffness_up[2][2]
            + local_stiffness_down[0][0] + local_stiffness_down[1][1] + local_stiffness_down[2][2];
      }

      for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        auto& edge_stencil = edge.getData(edgeStencilID_)->data[level];

        real_t local_stiffness_up[3][3];
        real_t local_stiffness_down[3][3];
        // first face
        Face* face = storage_->getFace(edge.neighborFaces()[0]);
        P1Space::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_up, P1Space::UPWARD);
        P1Space::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_down, P1Space::DOWNWARD);

        size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
        size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
        size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

        edge_stencil[0] = local_stiffness_up[end_id][opposite_id] + local_stiffness_down[opposite_id][end_id];
        edge_stencil[1] = local_stiffness_up[start_id][opposite_id] + local_stiffness_down[opposite_id][start_id];

        edge_stencil[2] = local_stiffness_up[end_id][start_id];
        edge_stencil[4] = local_stiffness_up[start_id][end_id];

        edge_stencil[3] = local_stiffness_up[start_id][start_id] + local_stiffness_up[end_id][end_id] + local_stiffness_down[opposite_id][opposite_id];

        if (edge.getNumNeighborFaces() == 2)
        {
          // second face
          Face* face = storage_->getFace(edge.neighborFaces()[1]);
          P1Space::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_up, P1Space::UPWARD);
          P1Space::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness_down, P1Space::DOWNWARD);

          size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
          size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
          size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

          edge_stencil[5] = local_stiffness_up[end_id][opposite_id] + local_stiffness_down[opposite_id][end_id];
          edge_stencil[6] = local_stiffness_up[start_id][opposite_id] + local_stiffness_down[opposite_id][start_id];

          edge_stencil[2] += local_stiffness_up[end_id][start_id];
          edge_stencil[4] += local_stiffness_up[start_id][end_id];

          edge_stencil[3] += local_stiffness_up[start_id][start_id] + local_stiffness_up[end_id][end_id] + local_stiffness_down[opposite_id][opposite_id];
        }
      }

      for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        auto& vertex_stencil = vertex.getData(vertexStencilID_)->data[level];

        // iterate over adjacent faces
        for (auto& faceId : vertex.neighborFaces())
        {
          Face* face = storage_->getFace(faceId);

          real_t local_stiffness[3][3];
          P1Space::compute_local_stiffness<UFCOperator>(*face, level, local_stiffness, P1Space::UPWARD);

          uint_t v_i = face->vertex_index(vertex.getID());

          std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

          // iterate over adjacent edges
          for (auto& edgeId : adj_edges)
          {
            uint_t edge_idx = vertex.edge_index(edgeId) + 1;
            Edge* edge = storage_->getEdge(edgeId);
            PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

            uint_t v_j = face->vertex_index(vertex_j);

            vertex_stencil[edge_idx] += local_stiffness[v_i][v_j];
          }

          vertex_stencil[0] += local_stiffness[v_i][v_i];
        }
      }

    }

  }

  ~P1Operator()
  {
  }

  void apply(const P1Function& src, P1Function& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (testFlag(vertex.type, flag))
//      {
//        P1Vertex::pull_halos(vertex, src.memory_id, level);
//      }
//    }
//
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (vertex.rank == rank && testFlag(vertex.type, flag))
//      {
//        P1Vertex::apply(vertex, this->memory_id, src.memory_id, dst.memory_id, level, updateType);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      P1Edge::pull_vertices(edge, dst.memory_id, level);
//      if (testFlag(edge.type, flag))
//      {
//        P1Edge::pull_halos(edge, src.memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      if (edge.rank == rank && testFlag(edge.type, flag))
//      {
//        P1Edge::apply(edge, this->memory_id, src.memory_id, dst.memory_id, level, updateType);
//      }
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      P1Face::pull_edges(face, dst.memory_id, level);
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      if (face.rank == rank && testFlag(face.type, flag))
//      {
//        P1Face::apply(level, face, this->memory_id, src.memory_id, dst.memory_id, updateType);
//      }
//    }
  }

  void smooth_gs(P1Function& dst, const P1Function& rhs, size_t level, DoFType flag)
  {
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (testFlag(vertex.type, flag))
//      {
//        P1Vertex::pull_halos(vertex, dst.memory_id, level);
//      }
//    }
//
//    for (Vertex& vertex : mesh.vertices)
//    {
//      if (vertex.rank == rank && testFlag(vertex.type, flag))
//      {
//        P1Vertex::smooth_gs(vertex, this->memory_id, dst.memory_id, rhs.memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      P1Edge::pull_vertices(edge, dst.memory_id, level);
//      if (testFlag(edge.type, flag))
//      {
//        P1Edge::pull_halos(edge, dst.memory_id, level);
//      }
//    }
//
//    for (Edge& edge : mesh.edges)
//    {
//      if (edge.rank == rank && testFlag(edge.type, flag))
//      {
//        P1Edge::smooth_gs(edge, this->memory_id, dst.memory_id, rhs.memory_id, level);
//      }
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      P1Face::pull_edges(face, dst.memory_id, level);
//    }
//
//    for (Face& face : mesh.faces)
//    {
//      if (face.rank == rank && testFlag(face.type, flag))
//      {
//        P1Face::smooth_gs(level, face, this->memory_id, dst.memory_id, rhs.memory_id);
//      }
//    }
  }

 private:
  PrimitiveDataID<VertexP1StencilMemory, Vertex> vertexStencilID_;
  PrimitiveDataID<EdgeP1StencilMemory, Edge> edgeStencilID_;
  PrimitiveDataID<FaceP1StencilMemory, Face> faceStencilID_;

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
