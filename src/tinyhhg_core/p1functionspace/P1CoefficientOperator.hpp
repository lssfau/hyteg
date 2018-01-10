
#pragma once

#include <fmt/format.h>
#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
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
#include "tinyhhg_core/p1functionspace/generated/p1_stokes_epsilon.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "tinyhhg_core/p1functionspace/P1Memory.hpp"

#include <tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp>
#include "P1Edge.hpp"
#include "P1Face.hpp"

namespace hhg
{

template<class UFCOperator,  bool Diagonal = false>
class P1CoefficientOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  P1CoefficientOperator(const std::shared_ptr< PrimitiveStorage > & storage, const std::shared_ptr<P1Function< real_t >>& coefficient, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel), coefficientP1_(coefficient)
  {
    init(minLevel, maxLevel);
  }

  P1CoefficientOperator(const std::shared_ptr< PrimitiveStorage > & storage, const std::shared_ptr<DGFunction< real_t >>& coefficient, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel), coefficientDG_(coefficient)
  {
    init(minLevel, maxLevel);
  }

  ~P1CoefficientOperator()
  {
  }

private:

  void init(uint_t minLevel, uint_t maxLevel)
  {
    auto faceP1LocalMatrixMemoryDataHandling = std::make_shared< FaceP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
    auto edgeP1LocalMatrixMemoryDataHandling = std::make_shared< EdgeP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
    auto vertexP1LocalMatrixMemoryDataHandling = std::make_shared< VertexP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
    storage_->addFaceData(faceLocalMatrixID_, faceP1LocalMatrixMemoryDataHandling, "P1OperatorFaceLocalMatrix");
    storage_->addEdgeData(edgeLocalMatrixID_, edgeP1LocalMatrixMemoryDataHandling, "P1OperatorEdgeLocalMatrix");
    storage_->addVertexData(vertexLocalMatrixID_, vertexP1LocalMatrixMemoryDataHandling, "P1OperatorVertexLocalMatrix");

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        auto faceLocalMatrices = face.getData(faceLocalMatrixID_);

        compute_local_stiffness(face, level, faceLocalMatrices->getGrayMatrix(level), fenics::GRAY);
        compute_local_stiffness(face, level, faceLocalMatrices->getBlueMatrix(level), fenics::BLUE);
      }

      for (auto& it : storage_->getEdges()) {
        Edge& edge = *it.second;

        auto edgeLocalMatrices = edge.getData(edgeLocalMatrixID_);

        // first face
        Face* face = storage_->getFace(edge.neighborFaces()[0]);
        compute_local_stiffness(*face, level, edgeLocalMatrices->getGrayMatrix(level, 0), fenics::GRAY);
        compute_local_stiffness(*face, level, edgeLocalMatrices->getBlueMatrix(level, 0), fenics::BLUE);



        if (edge.getNumNeighborFaces() == 2)
        {
          // second face
          Face* face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, edgeLocalMatrices->getGrayMatrix(level, 1), fenics::GRAY);
          compute_local_stiffness(*face, level, edgeLocalMatrices->getBlueMatrix(level, 1), fenics::BLUE);
        }
      }

      for (auto& it : storage_->getVertices()) {
        Vertex& vertex = *it.second;

        auto vertexLocalMatrices = vertex.getData(vertexLocalMatrixID_);

        // iterate over adjacent faces
        uint_t neighborId = 0;
        for (auto& faceId : vertex.neighborFaces())
        {
          Face* face = storage_->getFace(faceId);

          compute_local_stiffness(*face, level, vertexLocalMatrices->getGrayMatrix(level, neighborId), fenics::GRAY);
          ++neighborId;
        }
      }

    }
  }

  void apply_impl(P1Function< real_t >& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    // start pulling vertex halos
    src.getCommunicator(level)->startCommunication<Edge, Vertex>();
    if (coefficientP1_ != nullptr) {
      coefficientP1_->getCommunicator(level)->startCommunication<Edge, Vertex>();
    } else {
      coefficientDG_->getCommunicator(level)->startCommunication<Edge, Vertex>();
    }

    // start pulling edge halos
    src.getCommunicator(level)->startCommunication<Face, Edge>();
    if (coefficientP1_ != nullptr) {
      coefficientP1_->getCommunicator(level)->startCommunication<Face, Edge>();
    } else {
      coefficientDG_->getCommunicator(level)->startCommunication<Face, Edge>();
    }

    // end pulling vertex halos
    if (coefficientP1_ != nullptr) {
      coefficientP1_->getCommunicator(level)->endCommunication<Edge, Vertex>();
    } else {
      coefficientDG_->getCommunicator(level)->endCommunication<Edge, Vertex>();
    }
    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        vertexdof::macrovertex::applyCoefficient< real_t >(vertex, storage_, vertexLocalMatrixID_, src.getVertexDataID(), dst.getVertexDataID(), coefficientP1_->getVertexDataID(), level, updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    if (coefficientP1_ != nullptr) {
      coefficientP1_->getCommunicator(level)->endCommunication<Face, Edge>();
    } else {
      coefficientDG_->getCommunicator(level)->endCommunication<Face, Edge>();
    }
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        vertexdof::macroedge::applyCoefficient< real_t >(level, edge, storage_, edgeLocalMatrixID_, src.getEdgeDataID(), dst.getEdgeDataID(), coefficientP1_->getEdgeDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        if (coefficientP1_ != nullptr) {
          vertexdof::macroface::applyCoefficient< real_t >(level, face, faceLocalMatrixID_, src.getFaceDataID(), dst.getFaceDataID(), coefficientP1_->getFaceDataID(), updateType);
        } else {
          vertexdof::macroface::applyCoefficientDG< real_t >(level, face, faceLocalMatrixID_, src.getFaceDataID(), dst.getFaceDataID(), coefficientDG_->getFaceDataID(), updateType);
        }
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_gs_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag)
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
        WALBERLA_ABORT("To be implemented")
//        P1Vertex::smooth_gs(vertex, vertexLocalMatrixID_, dst.getVertexDataID(), rhs.getVertexDataID(), level);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    dst.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Edge::smooth_gs(level, edge, edgeLocalMatrixID_, dst.getEdgeDataID(), rhs.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Face::smooth_gs(level, face, faceLocalMatrixID_, dst.getFaceDataID(), rhs.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_jac_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, P1Function< real_t >& tmp, size_t level, DoFType flag)
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
        WALBERLA_ABORT("To be implemented")
//        P1Vertex::smooth_jac(vertex, vertexLocalMatrixID_, dst.getVertexDataID(), rhs.getVertexDataID(), tmp.getVertexDataID(), level);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    tmp.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Edge::smooth_jac(level, edge, edgeLocalMatrixID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Face::smooth_jac(level, face, faceLocalMatrixID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

#ifdef HHG_BUILD_WITH_PETSC
  void createMatrix_impl(P1Function< real_t >& src, P1Function< real_t >& dst, Mat& mat, size_t level, DoFType flag)
  {
    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Vertex::saveOperator(vertex, vertexLocalMatrixID_, src.getVertexDataID(), dst.getVertexDataID(), mat, level);
      }
    }

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Edge::saveOperator(level, edge, edgeLocalMatrixID_, src.getEdgeDataID(), dst.getEdgeDataID(), mat);
      }
    }

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Face::saveOperator(level, face, faceLocalMatrixID_, src.getFaceDataID(), dst.getFaceDataID(), mat);
      }
    }
  }
#endif

  PrimitiveDataID<VertexP1LocalMatrixMemory, Vertex> vertexLocalMatrixID_;
  PrimitiveDataID<EdgeP1LocalMatrixMemory, Edge> edgeLocalMatrixID_;
  PrimitiveDataID<FaceP1LocalMatrixMemory, Face> faceLocalMatrixID_;

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

private:
  std::shared_ptr<P1Function< real_t >> coefficientP1_;
  std::shared_ptr<DGFunction< real_t >> coefficientDG_;
};

typedef P1CoefficientOperator<p1_diffusion_cell_integral_0_otherwise> P1VariableCoefficientLaplaceOperator;

typedef P1CoefficientOperator<p1_stokes_epsilon_cell_integral_0_otherwise> P1CoefficientEpsilonOperator_uu;
typedef P1CoefficientOperator<p1_stokes_epsilon_cell_integral_1_otherwise> P1CoefficientEpsilonOperator_uv;
typedef P1CoefficientOperator<p1_stokes_epsilon_cell_integral_2_otherwise> P1CoefficientEpsilonOperator_vu;
typedef P1CoefficientOperator<p1_stokes_epsilon_cell_integral_3_otherwise> P1CoefficientEpsilonOperator_vv;

}


