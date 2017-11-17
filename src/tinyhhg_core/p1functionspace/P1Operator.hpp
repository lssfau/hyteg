
#pragma once

#include <fmt/format.h>
#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "P1DataHandling.hpp"
#include "P1Elements.hpp"

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

using walberla::real_t;

template<class UFCOperator,  bool Diagonal = false, bool Lumped = false, bool InvertDiagonal = false>
class P1Operator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  P1Operator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  {
    using namespace P1Elements;

    auto faceP1StencilMemoryDataHandling = std::make_shared< FaceP1StencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);
    auto edgeP1StencilMemoryDataHandling = std::make_shared< EdgeP1StencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);
    auto vertexP1StencilMemoryDataHandling = std::make_shared< VertexP1StencilMemoryDataHandling< real_t > >(minLevel_, maxLevel_);

    storage->addFaceData(faceStencilID_, faceP1StencilMemoryDataHandling, "P1OperatorFaceStencil");
    storage->addEdgeData(edgeStencilID_, edgeP1StencilMemoryDataHandling, "P1OperatorEdgeStencil");
    storage->addVertexData(vertexStencilID_, vertexP1StencilMemoryDataHandling, "P1OperatorVertexStencil");

    Matrix3r local_stiffness_gray;
    Matrix3r local_stiffness_blue;

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        auto face_stencil = face.getData(faceStencilID_)->getPointer( level );
        compute_local_stiffness(face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(face, level, local_stiffness_blue, fenics::BLUE);

        for (uint_t i = 0; i < FaceVertexDoF::P1GrayElements.size(); ++i) {
          assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[i], FaceVertexDoF::P1GrayDoFMaps[i], local_stiffness_gray, face_stencil);
        }

        for (uint_t i = 0; i < FaceVertexDoF::P1BlueElements.size(); ++i) {
          assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[i], FaceVertexDoF::P1BlueDoFMaps[i], local_stiffness_blue, face_stencil);
        }

        if (Lumped) {
          face_stencil[3] += face_stencil[0] + face_stencil[1] + face_stencil[2] + face_stencil[4] + face_stencil[5]
                             + face_stencil[6];

          face_stencil[0] = 0;
          face_stencil[1] = 0;
          face_stencil[2] = 0;
          face_stencil[4] = 0;
          face_stencil[5] = 0;
          face_stencil[6] = 0;
        }

        if (InvertDiagonal) {
          face_stencil[3] = 1.0 / face_stencil[3];
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("FACE.id = {}:face_stencil = {}", face.getID().getID(), PointND<real_t, 7>(&face_stencil[0])));
      }

      for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        auto edge_stencil = edge.getData(edgeStencilID_)->getPointer( level );

        // first face
        Face* face = storage_->getFace(edge.neighborFaces()[0]);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);

        size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
        size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
        size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

        assembleP1LocalStencil({{3, 2, 0}}, {{end_id, start_id, opposite_id}}, local_stiffness_gray, edge_stencil);
        assembleP1LocalStencil({{3, 0, 1}}, {{opposite_id, end_id, start_id}}, local_stiffness_blue, edge_stencil);
        assembleP1LocalStencil({{3, 1, 4}}, {{start_id, opposite_id, end_id}}, local_stiffness_gray, edge_stencil);

        if (edge.getNumNeighborFaces() == 2)
        {
          // second face
          Face* face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);

          size_t start_id = face->vertex_index(edge.neighborVertices()[0]);
          size_t end_id = face->vertex_index(edge.neighborVertices()[1]);
          size_t opposite_id = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

          assembleP1LocalStencil({{3, 4, 6}}, {{start_id, end_id, opposite_id}}, local_stiffness_gray, edge_stencil);
          assembleP1LocalStencil({{3, 6, 5}}, {{opposite_id, start_id, end_id}}, local_stiffness_blue, edge_stencil);
          assembleP1LocalStencil({{3, 5, 2}}, {{end_id, opposite_id, start_id}}, local_stiffness_gray, edge_stencil);
        }

        if (Lumped) {
          edge_stencil[3] += edge_stencil[0] + edge_stencil[1] + edge_stencil[2] + edge_stencil[4];

          edge_stencil[0] = 0;
          edge_stencil[1] = 0;
          edge_stencil[2] = 0;
          edge_stencil[4] = 0;

          if (edge.getNumNeighborFaces() == 2)
          {
            edge_stencil[3] += edge_stencil[5] + edge_stencil[6];
            edge_stencil[5] = 0;
            edge_stencil[6] = 0;
          }
        }

        if (InvertDiagonal) {
          edge_stencil[3] = 1.0 / edge_stencil[3];
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("EDGE.id = {}:edge_stencil = {}", edge.getID().getID(), PointND<real_t, 7>(&edge_stencil[0])));
      }

      for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        auto vertex_stencil = vertex.getData(vertexStencilID_)->getPointer( level );

        // iterate over adjacent faces
        for (auto& faceId : vertex.neighborFaces())
        {
          Face* face = storage_->getFace(faceId);

          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);

          uint_t v_i = face->vertex_index(vertex.getID());

          std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

          std::array<uint_t, 3> stencilMap;
          stencilMap[0] = 0;

          std::array<uint_t, 3> dofMap;
          dofMap[0] = v_i;

          // iterate over adjacent edges
          for (uint_t i = 0; i < adj_edges.size(); ++i)
          {
            uint_t edge_idx = vertex.edge_index(adj_edges[i]) + 1;
            Edge* edge = storage_->getEdge(adj_edges[i]);
            PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

            stencilMap[i+1] = edge_idx;
            dofMap[i+1] = face->vertex_index(vertex_j);

          }

          assembleP1LocalStencil(stencilMap, dofMap, local_stiffness_gray, vertex_stencil);

//          WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("VERTEX.id = {}:vertex_stencil = {}", vertex.getID().getID(), PointND<real_t, 3>(&vertex_stencil[0])));
        }

        if (Lumped) {
          for (uint_t i = 1; i < vertex.getData(vertexStencilID_)->getSize(level); ++i) {
            vertex_stencil[0] += vertex_stencil[i];
            vertex_stencil[i] = 0;
          }
        }

        if (InvertDiagonal) {
          vertex_stencil[0] = 1.0 / vertex_stencil[0];
        }

      }

    }

  }

  ~P1Operator()
  {
  }

  void scale(real_t scalar) {
    for (uint_t level = minLevel_; level <= maxLevel_; ++level) {
      for (auto &it : storage_->getFaces()) {
        Face &face = *it.second;
        auto face_stencil = face.getData(faceStencilID_)->getPointer( level );
        for (uint_t i = 0; i < 7; ++i) {
          face_stencil[i] *= scalar;
        }
      }

      for (auto& it : storage_->getEdges()) {
        Edge &edge = *it.second;
        auto edge_stencil = edge.getData(edgeStencilID_)->getPointer( level );
        for (uint_t i = 0; i < 5; ++i) {
          edge_stencil[i] *= scalar;
        }

        if (edge.getNumNeighborFaces() == 2) {
          edge_stencil[5] *= scalar;
          edge_stencil[6] *= scalar;
        }
      }

      for (auto& it : storage_->getVertices()) {
        Vertex &vertex = *it.second;
        auto vertex_stencil = vertex.getData(vertexStencilID_)->getPointer( level );
        for (uint_t i = 0; i < vertex.getData(vertexStencilID_)->getSize(level); ++i) {
          vertex_stencil[i] *= scalar;
        }
      }
    }
  }

  const PrimitiveDataID<VertexP1StencilMemory< real_t >, Vertex> &getVertexStencilID() const { return vertexStencilID_; }

  const PrimitiveDataID<EdgeP1StencilMemory< real_t >, Edge> &getEdgeStencilID() const { return edgeStencilID_; }

  const PrimitiveDataID<FaceP1StencilMemory< real_t >, Face> &getFaceStencilID() const { return faceStencilID_; }

private:

  void apply_impl(P1Function< real_t > & src, P1Function< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace)
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
        P1Vertex::apply< real_t >(vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), level, updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::apply< real_t >(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::apply< real_t >(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_gs_impl(P1Function< real_t > & dst, P1Function< real_t > & rhs, size_t level, DoFType flag)
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
        P1Edge::smooth_gs< real_t >(level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::smooth_gs< real_t >(level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_sor_impl(P1Function< real_t > & dst, P1Function< real_t > & rhs, real_t relax, size_t level, DoFType flag)
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
        P1Vertex::smooth_sor(vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), level, relax);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    dst.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::smooth_sor< real_t >(level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), relax);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::smooth_sor< real_t >(level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), relax);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_jac_impl(P1Function< real_t > & dst, P1Function< real_t > & rhs, P1Function< real_t > & tmp, size_t level, DoFType flag)
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
        P1Edge::smooth_jac< real_t >(level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1Face::smooth_jac< real_t >(level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  PrimitiveDataID<VertexP1StencilMemory< real_t >, Vertex> vertexStencilID_;
  PrimitiveDataID<EdgeP1StencilMemory< real_t >, Edge> edgeStencilID_;
  PrimitiveDataID<FaceP1StencilMemory< real_t >, Face> faceStencilID_;

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
typedef P1Operator<p1_mass_cell_integral_0_otherwise, false, true, true> P1LumpedInvMassOperator;

typedef P1Operator<p1_pspg_cell_integral_0_otherwise> P1PSPGOperator;

}


