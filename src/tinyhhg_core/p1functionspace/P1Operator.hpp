
#pragma once

#include <fmt/format.h>
#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "P1DataHandling.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics.hpp"

#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "tinyhhg_core/p1functionspace/P1Memory.hpp"

#include "P1Vertex.hpp"
#include "P1Edge.hpp"
#include "P1Face.hpp"

namespace hhg
{

template<class UFCOperator,  bool Diagonal = false>
class P1Operator : public Operator< P1Function, P1Function >
{
public:
  P1Operator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  {
    auto faceP1StencilMemoryDataHandling = std::make_shared< FaceP1StencilMemoryDataHandling >(minLevel_, maxLevel_);
    auto edgeP1StencilMemoryDataHandling = std::make_shared< EdgeP1StencilMemoryDataHandling >(minLevel_, maxLevel_);
    auto vertexP1StencilMemoryDataHandling = std::make_shared< VertexP1StencilMemoryDataHandling >(minLevel_, maxLevel_);
    storage->addFaceData(faceStencilID_, faceP1StencilMemoryDataHandling, "P1OperatorFaceStencil");
    storage->addEdgeData(edgeStencilID_, edgeP1StencilMemoryDataHandling, "P1OperatorEdgeStencil");
    storage->addVertexData(vertexStencilID_, vertexP1StencilMemoryDataHandling, "P1OperatorVertexStencil");

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        auto& face_stencil = face.getData(faceStencilID_)->data[level];

        Matrix3r local_stiffness_up;
        Matrix3r local_stiffness_down;
        compute_local_stiffness(face, level, local_stiffness_up, fenics::GRAY);
        compute_local_stiffness(face, level, local_stiffness_down, fenics::BLUE);

        face_stencil[0] = local_stiffness_down(0,2) + local_stiffness_up(2,0);
        face_stencil[1] = local_stiffness_down(1,2) + local_stiffness_up(2,1);
        face_stencil[2] = local_stiffness_down(0,1) + local_stiffness_up(1,0);

        face_stencil[4] = local_stiffness_down(1,0) + local_stiffness_up(0,1);
        face_stencil[5] = local_stiffness_down(2,1) + local_stiffness_up(1,2);
        face_stencil[6] = local_stiffness_down(2,0) + local_stiffness_up(0,2);

        face_stencil[3] = local_stiffness_up(0,0) + local_stiffness_up(1,1) + local_stiffness_up(2,2)
            + local_stiffness_down(0,0) + local_stiffness_down(1,1) + local_stiffness_down(2,2);

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("FACE.id = {}:face_stencil = {}", face.getID().getID(), PointND<real_t, 7>(&face_stencil[0])));
      }

      for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        auto& edge_stencil = edge.getData(edgeStencilID_)->data[level];

        Matrix3r local_stiffness_up;
        Matrix3r local_stiffness_down;
        // first face
        Face* face = storage_->getFace(edge.neighborFaces()[0]);
        compute_local_stiffness(*face, level, local_stiffness_up, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_down, fenics::BLUE);

        size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
        size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
        size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

        edge_stencil[0] = local_stiffness_up(end_id, opposite_id) + local_stiffness_down(opposite_id, end_id);
        edge_stencil[1] = local_stiffness_up(start_id, opposite_id) + local_stiffness_down(opposite_id, start_id);

        edge_stencil[2] = local_stiffness_up(end_id, start_id);
        edge_stencil[4] = local_stiffness_up(start_id, end_id);

        edge_stencil[3] = local_stiffness_up(start_id, start_id) + local_stiffness_up(end_id, end_id) + local_stiffness_down(opposite_id, opposite_id);

        if (edge.getNumNeighborFaces() == 2)
        {
          // second face
          Face* face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, local_stiffness_up, fenics::GRAY);
          compute_local_stiffness(*face, level, local_stiffness_down, fenics::BLUE);

          size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
          size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
          size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

          edge_stencil[5] = local_stiffness_up(end_id, opposite_id) + local_stiffness_down(opposite_id, end_id);
          edge_stencil[6] = local_stiffness_up(start_id, opposite_id) + local_stiffness_down(opposite_id, start_id);

          edge_stencil[2] += local_stiffness_up(end_id, start_id);
          edge_stencil[4] += local_stiffness_up(start_id, end_id);

          edge_stencil[3] += local_stiffness_up(start_id, start_id) + local_stiffness_up(end_id, end_id) + local_stiffness_down(opposite_id, opposite_id);
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("EDGE.id = {}:edge_stencil = {}", edge.getID().getID(), PointND<real_t, 7>(&edge_stencil[0])));
      }

      for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        auto& vertex_stencil = vertex.getData(vertexStencilID_)->data[level];

        // iterate over adjacent faces
        for (auto& faceId : vertex.neighborFaces())
        {
          Face* face = storage_->getFace(faceId);

          Matrix3r local_stiffness;
          compute_local_stiffness(*face, level, local_stiffness, fenics::GRAY);

          uint_t v_i = face->vertex_index(vertex.getID());

          std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

          // iterate over adjacent edges
          for (auto& edgeId : adj_edges)
          {
            uint_t edge_idx = vertex.edge_index(edgeId) + 1;
            Edge* edge = storage_->getEdge(edgeId);
            PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

            uint_t v_j = face->vertex_index(vertex_j);

            vertex_stencil[edge_idx] += local_stiffness(v_i, v_j);
          }

          vertex_stencil[0] += local_stiffness(v_i, v_i);

//          WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("VERTEX.id = {}:vertex_stencil = {}", vertex.getID().getID(), PointND<real_t, 3>(&vertex_stencil[0])));
        }
      }

    }

  }

  ~P1Operator()
  {
  }

private:

  void apply_impl(P1Function& src, P1Function& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    // start pulling vertex halos
    src.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    src.getCommunicator(level)->startCommunication<Face, Edge>();

    // end pulling vertex halos
    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::apply(vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), level, updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::apply(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_gs_impl(P1Function& dst, P1Function& rhs, size_t level, DoFType flag)
  {
    // start pulling vertex halos
    dst.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    dst.getCommunicator(level)->startCommunication<Face, Edge>();

    // end pulling vertex halos
    dst.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::smooth_gs(vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), level);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    dst.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::smooth_gs(level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::smooth_gs(level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_jac_impl(P1Function& dst, P1Function& rhs, P1Function& tmp, size_t level, DoFType flag)
  {
    // start pulling vertex halos
    tmp.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    tmp.getCommunicator(level)->startCommunication<Face, Edge>();

    // end pulling vertex halos
    tmp.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::smooth_jac(vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), tmp.getVertexDataID(), level);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    tmp.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::smooth_jac(level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::smooth_jac(level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

#ifdef HHG_BUILD_WITH_PETSC
  void createMatrix_impl(P1Function& src, P1Function& dst, Mat& mat, size_t level, DoFType flag)
  {
    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        P1Vertex::saveOperator(vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), mat, level);
      }
    }

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::saveOperator(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), mat);
      }
    }

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::saveOperator(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), mat);
      }
    }
  }
#endif

  PrimitiveDataID<VertexP1StencilMemory, Vertex> vertexStencilID_;
  PrimitiveDataID<EdgeP1StencilMemory, Edge> edgeStencilID_;
  PrimitiveDataID<FaceP1StencilMemory, Face> faceStencilID_;

  void compute_local_stiffness(const Face &face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);

    if (Diagonal) {
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          if (i != j) {
            local_stiffness(i,j) = real_t(0);
          }
        }
      }
    }
  }

};

typedef P1Operator<p1_diffusion_cell_integral_0_otherwise> P1LaplaceOperator;
typedef P1Operator<p1_diffusion_cell_integral_0_otherwise, true> P1DiagonalLaplaceOperator;

typedef P1Operator<p1_div_cell_integral_0_otherwise> P1DivxOperator;
typedef P1Operator<p1_div_cell_integral_1_otherwise> P1DivyOperator;

typedef P1Operator<p1_divt_cell_integral_0_otherwise> P1DivTxOperator;
typedef P1Operator<p1_divt_cell_integral_1_otherwise> P1DivTyOperator;

typedef P1Operator<p1_mass_cell_integral_0_otherwise> P1MassOperator;

typedef P1Operator<p1_pspg_cell_integral_0_otherwise> P1PSPGOperator;

}


