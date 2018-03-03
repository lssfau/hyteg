
#pragma once

#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "P1DataHandling.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics/fenics.hpp"

#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"
#include "tinyhhg_core/p1functionspace/generated/p1_stokes_epsilon.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div_K_grad.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"

#include <tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFBlending.hpp>

namespace hhg
{

class P1BlendingOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  P1BlendingOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                        size_t minLevel,
                        size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  {
  }

  ~P1BlendingOperator()
  {
  }

  void addOperatorAndCoefficient(const fenics::TabulateTensor& fenicsTabulateTensor) {
    operators_.push_back(fenicsTabulateTensor);
  }

  void init()
  {
    WALBERLA_ASSERT(operators_.size() > 0, "At least one operator needs to be added before calling init()");

    vertexLocalMatrixIDs_.resize(operators_.size());
    edgeLocalMatrixIDs_.resize(operators_.size());
    faceLocalMatrixIDs_.resize(operators_.size());

    for(uint_t oprIdx = 0; oprIdx < operators_.size(); ++oprIdx) {
      auto faceP1LocalMatrixMemoryDataHandling = std::make_shared<FaceP1LocalMatrixMemoryDataHandling>(minLevel_,
                                                                                                       maxLevel_);
      auto edgeP1LocalMatrixMemoryDataHandling = std::make_shared<EdgeP1LocalMatrixMemoryDataHandling>(minLevel_,
                                                                                                       maxLevel_);
      auto vertexP1LocalMatrixMemoryDataHandling = std::make_shared<VertexP1LocalMatrixMemoryDataHandling>(minLevel_,
                                                                                                           maxLevel_);

      auto faceLocalMatrixID_ = PrimitiveID();

      storage_->addFaceData(faceLocalMatrixIDs_[oprIdx], faceP1LocalMatrixMemoryDataHandling, "P1OperatorFaceLocalMatrix");
      storage_->addEdgeData(edgeLocalMatrixIDs_[oprIdx], edgeP1LocalMatrixMemoryDataHandling, "P1OperatorEdgeLocalMatrix");
      storage_->addVertexData(vertexLocalMatrixIDs_[oprIdx], vertexP1LocalMatrixMemoryDataHandling,
                              "P1OperatorVertexLocalMatrix");

      for (uint_t level = minLevel_; level <= maxLevel_; ++level) {

        for (auto &it : storage_->getFaces()) {
          Face &face = *it.second;

          auto faceLocalMatrices = face.getData(faceLocalMatrixIDs_[oprIdx]);

          compute_local_stiffness(operators_[oprIdx], face, level, faceLocalMatrices->getGrayMatrix(level), fenics::GRAY);
          compute_local_stiffness(operators_[oprIdx], face, level, faceLocalMatrices->getBlueMatrix(level), fenics::BLUE);
        }

        for (auto &it : storage_->getEdges()) {
          Edge &edge = *it.second;

          auto edgeLocalMatrices = edge.getData(edgeLocalMatrixIDs_[oprIdx]);

          // first face
          Face *face = storage_->getFace(edge.neighborFaces()[0]);
          compute_local_stiffness(operators_[oprIdx], *face, level, edgeLocalMatrices->getGrayMatrix(level, 0), fenics::GRAY);
          compute_local_stiffness(operators_[oprIdx], *face, level, edgeLocalMatrices->getBlueMatrix(level, 0), fenics::BLUE);


          if (edge.getNumNeighborFaces() == 2) {
            // second face
            Face *face = storage_->getFace(edge.neighborFaces()[1]);
            compute_local_stiffness(operators_[oprIdx], *face, level, edgeLocalMatrices->getGrayMatrix(level, 1), fenics::GRAY);
            compute_local_stiffness(operators_[oprIdx], *face, level, edgeLocalMatrices->getBlueMatrix(level, 1), fenics::BLUE);
          }
        }

        for (auto &it : storage_->getVertices()) {
          Vertex &vertex = *it.second;

          auto vertexLocalMatrices = vertex.getData(vertexLocalMatrixIDs_[oprIdx]);

          // iterate over adjacent faces
          uint_t neighborId = 0;
          for (auto &faceId : vertex.neighborFaces()) {
            Face *face = storage_->getFace(faceId);

            compute_local_stiffness(operators_[oprIdx], *face, level, vertexLocalMatrices->getGrayMatrix(level, neighborId), fenics::GRAY);
            ++neighborId;
          }
        }

      }
    }
  }

  void applyPartial(P1Function< real_t >& src, P1Function< real_t >& dst, uint_t memoryLevel, uint_t coarseLevel, DoFType flag, UpdateType updateType = Replace)
  {
    // start pulling vertex halos
    src.getCommunicator(memoryLevel)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    src.getCommunicator(memoryLevel)->startCommunication<Face, Edge>();

    src.getCommunicator(memoryLevel)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        vertexdof::blending::macrovertex::applyBlending< real_t >(memoryLevel, vertex, storage_, vertexLocalMatrixIDs_, src.getVertexDataID(), dst.getVertexDataID(), updateType);
      }
    }

    dst.getCommunicator(memoryLevel)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    src.getCommunicator(memoryLevel)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        vertexdof::blending::macroedge::applyBlending< real_t >(memoryLevel, edge, storage_, edgeLocalMatrixIDs_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
    }

    dst.getCommunicator(memoryLevel)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(memoryLevel)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        vertexdof::blending::macroface::applyPartialBlending< real_t >(memoryLevel, coarseLevel, face, faceLocalMatrixIDs_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }

    dst.getCommunicator(memoryLevel)->endCommunication<Edge, Face>();
  }

private:
  void apply_impl(P1Function< real_t >& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    // start pulling vertex halos
    src.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    src.getCommunicator(level)->startCommunication<Face, Edge>();

    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        vertexdof::blending::macrovertex::applyBlending< real_t >(level, vertex, storage_, vertexLocalMatrixIDs_, src.getVertexDataID(), dst.getVertexDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        vertexdof::blending::macroedge::applyBlending< real_t >(level, edge, storage_, edgeLocalMatrixIDs_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        vertexdof::blending::macroface::applyBlending< real_t >(level, face, faceLocalMatrixIDs_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
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
        vertexdof::blending::macrovertex::smooth_gs_blending(level, vertex, storage_, vertexLocalMatrixIDs_, dst.getVertexDataID(), rhs.getVertexDataID());
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    dst.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        vertexdof::blending::macroedge::smoothGSBlending<real_t>(level, edge, storage_, edgeLocalMatrixIDs_, dst.getEdgeDataID(), rhs.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        vertexdof::blending::macroface::smoothGSBlending<real_t>(level, face, faceLocalMatrixIDs_, dst.getFaceDataID(), rhs.getFaceDataID());
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

  std::vector<PrimitiveDataID<VertexP1LocalMatrixMemory, Vertex>> vertexLocalMatrixIDs_;
  std::vector<PrimitiveDataID<EdgeP1LocalMatrixMemory, Edge>> edgeLocalMatrixIDs_;
  std::vector<PrimitiveDataID<FaceP1LocalMatrixMemory, Face>> faceLocalMatrixIDs_;

  void compute_local_stiffness(const fenics::TabulateTensor& opr, const Face &face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    opr(local_stiffness.data(), NULL, coords, 0);
  }

private:
  std::vector<fenics::TabulateTensor> operators_;
};

template<class FenicsOperator>
class P1ScalarBlendingOperator : public P1BlendingOperator {
public:
   P1ScalarBlendingOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                           size_t minLevel,
                           size_t maxLevel)
      : P1BlendingOperator(storage, minLevel, maxLevel) {
    fenicsOperator = std::make_shared<FenicsOperator>();

    this->addOperatorAndCoefficient(std::bind(&FenicsOperator::tabulate_tensor, fenicsOperator.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));

    this->init();
  }

 private:
  std::shared_ptr<FenicsOperator> fenicsOperator;
};

typedef P1ScalarBlendingOperator<p1_mass_cell_integral_0_otherwise> P1BlendingMassOperator;

template<class FenicsOperator1, class FenicsOperator2, class FenicsOperator3>
class P1TensorBlendingOperator : public P1BlendingOperator {
public:
  P1TensorBlendingOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                              size_t minLevel,
                              size_t maxLevel)
      : P1BlendingOperator(storage, minLevel, maxLevel) {
    fenicsOperator1 = std::make_shared<FenicsOperator1>();
    fenicsOperator2 = std::make_shared<FenicsOperator2>();
    fenicsOperator3 = std::make_shared<FenicsOperator3>();

    this->addOperatorAndCoefficient(std::bind(&FenicsOperator1::tabulate_tensor, fenicsOperator1.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));

    this->addOperatorAndCoefficient(std::bind(&FenicsOperator2::tabulate_tensor, fenicsOperator2.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));

    this->addOperatorAndCoefficient(std::bind(&FenicsOperator3::tabulate_tensor, fenicsOperator3.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));

    this->init();
  }

private:
  std::shared_ptr<FenicsOperator1> fenicsOperator1;
  std::shared_ptr<FenicsOperator2> fenicsOperator2;
  std::shared_ptr<FenicsOperator3> fenicsOperator3;
};

typedef P1TensorBlendingOperator<p1_div_k_grad_cell_integral_0_otherwise,
                                    p1_div_k_grad_cell_integral_1_otherwise,
                                    p1_div_k_grad_cell_integral_2_otherwise> P1BlendingLaplaceOperator;

}


