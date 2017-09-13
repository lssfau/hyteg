#pragma once

#include <fmt/format.h>
#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "BubbleToP1Memory.hpp"
#include "BubbleToP1DataHandling.hpp"

#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics.hpp"

#include "generated/bubble_to_p1_divt.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "BubbleToP1Vertex.hpp"
#include "BubbleToP1Edge.hpp"
#include "BubbleToP1Face.hpp"

namespace hhg
{

template<class UFCOperator>
class BubbleToP1Operator : public Operator< BubbleFunction, P1Function< real_t > >
{
 public:
  BubbleToP1Operator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel)
  {
    auto faceBubbleToP1StencilMemoryDataHandling = std::make_shared< FaceBubbleToP1StencilMemoryDataHandling >(minLevel_, maxLevel_);
    auto edgeBubbleToP1StencilMemoryDataHandling = std::make_shared< EdgeBubbleToP1StencilMemoryDataHandling >(minLevel_, maxLevel_);
    auto vertexeBubbleToP1StencilMemoryDataHandling = std::make_shared< VertexBubbleToP1StencilMemoryDataHandling >(minLevel_, maxLevel_);

    storage->addFaceData(faceStencilID_, faceBubbleToP1StencilMemoryDataHandling, "BubbleToP1OperatorFaceStencil");
    storage->addEdgeData(edgeStencilID_, edgeBubbleToP1StencilMemoryDataHandling, "BubbleToP1OperatorEdgeStencil");
    storage->addVertexData(vertexStencilID_, vertexeBubbleToP1StencilMemoryDataHandling, "BubbleToP1OperatorVertexStencil");

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      // assemble face stencil
      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        auto& face_stencil = face.getData(faceStencilID_)->data[level];

        real_t local_stiffness_gray[3][1];
        real_t local_stiffness_blue[3][1];
        compute_local_stiffness(face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(face, level, local_stiffness_blue, fenics::BLUE);

        face_stencil[BubbleFace::CoordsVertex::CELL_GRAY_SE] = local_stiffness_gray[2][0];
        face_stencil[BubbleFace::CoordsVertex::CELL_GRAY_NW] = local_stiffness_gray[1][0];
        face_stencil[BubbleFace::CoordsVertex::CELL_GRAY_NE] = local_stiffness_gray[0][0];

        face_stencil[BubbleFace::CoordsVertex::CELL_BLUE_SW] = local_stiffness_blue[0][0];
        face_stencil[BubbleFace::CoordsVertex::CELL_BLUE_SE] = local_stiffness_blue[1][0];
        face_stencil[BubbleFace::CoordsVertex::CELL_BLUE_NW] = local_stiffness_blue[2][0];
      }

      // assemble edge stencil
      for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        auto& edge_stencil = edge.getData(edgeStencilID_)->data[level];

        real_t local_stiffness_gray[3][1];
        real_t local_stiffness_blue[3][1];
        // first face
        Face* face = storage_->getFace(edge.neighborFaces()[0]);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);

        size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
        size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
        size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

        edge_stencil[BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SW] = local_stiffness_gray[end_id][0];
        edge_stencil[BubbleEdge::EdgeCoordsVertex::CELL_BLUE_SE] = local_stiffness_blue[opposite_id][0];
        edge_stencil[BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SE] = local_stiffness_gray[start_id][0];

        if (edge.getNumNeighborFaces() == 2)
        {
          // second face
          Face* face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);

          size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
          size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
          size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

          edge_stencil[BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NW] = local_stiffness_gray[end_id][0];
          edge_stencil[BubbleEdge::EdgeCoordsVertex::CELL_BLUE_NW] = local_stiffness_blue[opposite_id][0];
          edge_stencil[BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NE] = local_stiffness_gray[start_id][0];
        }
      }

      // assemble vertex stencil
      for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        auto& vertex_stencil = vertex.getData(vertexStencilID_)->data[level];

        // iterate over adjacent faces
        for (auto& faceId : vertex.neighborFaces())
        {
          Face* face = storage_->getFace(faceId);

          real_t local_stiffness[3][1];
          compute_local_stiffness(*face, level, local_stiffness, fenics::GRAY);

          uint_t v_i = face->vertex_index(vertex.getID());

          size_t f_i = vertex.face_index(face->getID());
          vertex_stencil[f_i] = local_stiffness[v_i][0];
        }
      }
    }
  }

  ~BubbleToP1Operator()
  {
  }

 private:

  void apply_impl(BubbleFunction& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    // Since the Bubble dofs are in the interior, we have to pull them through the edges first
    src.getCommunicator(level)->startCommunication<Face, Edge>();
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    src.getCommunicator(level)->startCommunication<Edge, Vertex>();
    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        BubbleToP1Vertex::apply(vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), level, updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        BubbleToP1Edge::apply(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        BubbleToP1Face::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

#ifdef HHG_BUILD_WITH_PETSC
  void createMatrix_impl(BubbleFunction& src, P1Function<real_t>& dst, Mat &mat, size_t level, DoFType flag)
  {
    // Since the Bubble dofs are in the interior, we have to pull them through the edges first //TODO: Implement!
    /*src.getCommunicator(level)->startCommunication<Face, Edge>();
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    src.getCommunicator(level)->startCommunication<Edge, Vertex>();
    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        BubbleToP1Vertex::saveOperator(vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), out, level);
      }
    }

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        BubbleToP1Edge::saveOperator(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), out);
      }
    }

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        BubbleToP1Face::saveOperator(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), out);
      }
    }
     */
  }
#endif

  PrimitiveDataID<VertexBubbleToP1StencilMemory, Vertex> vertexStencilID_;
  PrimitiveDataID<EdgeBubbleToP1StencilMemory, Edge> edgeStencilID_;
  PrimitiveDataID<FaceBubbleToP1StencilMemory, Face> faceStencilID_;

  void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[3][1], fenics::ElementType element_type) {
    real_t A[3];
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(A, NULL, coords, 0);

    for (size_t i = 0; i < 3; ++i) {
      size_t j = 0;
      local_stiffness[i][j] = A[3 * j + i];
    }
  }
};

typedef BubbleToP1Operator<bubble_to_p1_divt_cell_integral_0_otherwise> BubbleToP1DivTxOperator;
typedef BubbleToP1Operator<bubble_to_p1_divt_cell_integral_1_otherwise> BubbleToP1DivTyOperator;

}
