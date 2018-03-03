
#pragma once

#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "P1DataHandling.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "VertexDoFBlending.hpp"
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

#include "VertexDoFMacroVertex.hpp"
#include "VertexDoFMacroEdge.hpp"
#include "VertexDoFMacroFace.hpp"

#include "tinyhhg_core/polynomial/LSQPInterpolator.hpp"

namespace hhg
{

class P1PolynomialBlendingOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  typedef LSQPInterpolator<MonomialBasis2D> Interpolator;

  P1PolynomialBlendingOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                       uint_t minLevel,
                       uint_t maxLevel,
                       uint_t interpolationLevel)
      : Operator(storage, minLevel, maxLevel),
        interpolationLevel_(interpolationLevel)
  {
  }

  ~P1PolynomialBlendingOperator()
  {
  }

  void addOperator(const fenics::TabulateTensor& fenicsTabulateTensor) {
    operators_.push_back(fenicsTabulateTensor);
  }

  void init() {
    WALBERLA_ASSERT(operators_.size() > 0, "At least one operator needs to be added before calling init()");
    initLocalStiffnessMatrices();
  }

private:

  void initLocalStiffnessMatrices()
  {
    auto faceP1PolynomialMemoryDataHandling = std::make_shared< FaceP1PolynomialMemoryDataHandling >(polyDegree_);
    storage_->addFaceData(facePolynomialID_, faceP1PolynomialMemoryDataHandling, "P1OperatorFacePolynomial");

    vertexLocalMatrixIDs_.resize(operators_.size());
    edgeLocalMatrixIDs_.resize(operators_.size());
    faceLocalMatrixIDs_.resize(operators_.size());

    for(uint_t oprIdx = 0; oprIdx < operators_.size(); ++oprIdx) {

      auto faceP1LocalMatrixMemoryDataHandling = std::make_shared< FaceP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
      auto edgeP1LocalMatrixMemoryDataHandling = std::make_shared< EdgeP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
      auto vertexP1LocalMatrixMemoryDataHandling = std::make_shared< VertexP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);

      storage_->addFaceData(faceLocalMatrixIDs_[oprIdx], faceP1LocalMatrixMemoryDataHandling, "P1OperatorFaceLocalMatrix");
      storage_->addEdgeData(edgeLocalMatrixIDs_[oprIdx], edgeP1LocalMatrixMemoryDataHandling, "P1OperatorEdgeLocalMatrix");
      storage_->addVertexData(vertexLocalMatrixIDs_[oprIdx], vertexP1LocalMatrixMemoryDataHandling, "P1OperatorVertexLocalMatrix");

      // Initialize other local stiffness matrices on lower dimensional primitives as usual
      for (uint_t level = minLevel_; level <= maxLevel_; ++level) {

        for (auto& it : storage_->getFaces()) {
          Face& face = *it.second;

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

public:
  void interpolateStencils(uint_t fineLevel, uint_t polyDegree) {
    typedef stencilDirection SD;
    using namespace P1Elements;

    std::vector<real_t> faceStencil(7);

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      auto facePolynomials = face.getData(facePolynomialID_);
      facePolynomials->addDegree(polyDegree);

      uint_t rowsize = levelinfo::num_microvertices_per_edge(interpolationLevel_);
      uint_t rowsizeFine = levelinfo::num_microvertices_per_edge(fineLevel);
      uint_t inner_rowsize = rowsize;
      real_t coeffWeight;

      Interpolator horiInterpolator(polyDegree, interpolationLevel_);
      Interpolator vertInterpolator(polyDegree, interpolationLevel_);
      Interpolator diagInterpolator(polyDegree, interpolationLevel_);

      Point3D x, x0;
      x0 = face.coords[0];

      Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
      Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

      // fine directions
      Point3D d0f = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsizeFine - 1));
      Point3D d2f = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsizeFine - 1));

      real_t ref_H = walberla::real_c(1)/(walberla::real_c(rowsize - 1));
      real_t ref_h = walberla::real_c(1)/(walberla::real_c(rowsizeFine - 1));
      Point2D ref_x;

      Point3D offsetSW = -1.0/3.0 * d0f - 1.0/3.0 * d2f;
      Point3D offsetS = 1.0/3.0 * d0f - 2.0/3.0 * d2f;
      Point3D offsetSE = 2.0/3.0 * d0f - 1.0/3.0 * d2f;

      Point3D offsetNE = 1.0/3.0 * d0f + 1.0/3.0 * d2f;
      Point3D offsetN = -1.0/3.0 * d0f + 2.0/3.0 * d2f;
      Point3D offsetNW = -2.0/3.0 * d0f + 1.0/3.0 * d2f;

      Matrix2r mappingTensor;
      real_t* coeffs[] = { &mappingTensor(0,0), &mappingTensor(0,1), &mappingTensor(1,1) };
      std::vector<FaceP1LocalMatrixMemory *> localMatricesVector;
      for (auto operatorId : faceLocalMatrixIDs_) {
        localMatricesVector.push_back(face.getData(operatorId));
      }

      for (uint_t j = 1; j < rowsize - 2; ++j) {

        ref_x[1] = j * ref_H;

        x = x0;
        x += real_c(j)*d2 + d0;

        uint_t i;
        for (i = 1; i < inner_rowsize - 2; ++i) {

          ref_x[0] = i * ref_H;

          vertexdof::blending::macroface::assembleTensorStencil(fineLevel, face, faceStencil, x, localMatricesVector, mappingTensor, coeffs, offsetS, offsetSE, offsetSW, offsetNW, offsetN, offsetNE);

          horiInterpolator.addInterpolationPoint(ref_x + Point2D{{ -0.5 * ref_h, 0.0 }}, faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_W)]);

          vertInterpolator.addInterpolationPoint(ref_x + Point2D{{ 0.0, -0.5 * ref_h }}, faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_S)]);

          diagInterpolator.addInterpolationPoint(ref_x + Point2D{{ 0.5 * ref_h, -0.5 * ref_h }}, faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_SE)]);

          if (i == 1) {
            diagInterpolator.addInterpolationPoint(ref_x + Point2D{{ -0.5 * ref_h, 0.5 * ref_h }}, faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_NW)]);
          }

          if (i == inner_rowsize - 2 - 1) {
            horiInterpolator.addInterpolationPoint(ref_x + Point2D{{ 0.5 * ref_h, 0.0 }}, faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_E)]);
          }

          x += d0;
        }

        vertInterpolator.addInterpolationPoint(ref_x + Point2D{{ 0.0, 0.5 * ref_h }}, faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_N)]);

        --inner_rowsize;
      }

      horiInterpolator.interpolate(facePolynomials->getHoriPolynomial(polyDegree));
      vertInterpolator.interpolate(facePolynomials->getVertPolynomial(polyDegree));
      diagInterpolator.interpolate(facePolynomials->getDiagPolynomial(polyDegree));
    }
  }

  void useDegree(uint_t degree) {
    polyDegree_ = degree;
  }

  real_t lInfinityError(uint_t degreeA, uint_t degreeB, uint_t level) {

    real_t error = std::numeric_limits<real_t>::min();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;
      auto polynomials = face.getData(facePolynomialID_);

      auto horiPolyA = polynomials->getHoriPolynomial(degreeA);
      auto vertPolyA = polynomials->getVertPolynomial(degreeA);
      auto diagPolyA = polynomials->getDiagPolynomial(degreeA);

      auto horiPolyB = polynomials->getHoriPolynomial(degreeB);
      auto vertPolyB = polynomials->getVertPolynomial(degreeB);
      auto diagPolyB = polynomials->getDiagPolynomial(degreeB);

//      error = std::max(error, horiPolyA.lInfinityError(horiPolyB));
//      error = std::max(error, vertPolyA.lInfinityError(vertPolyB));
//      error = std::max(error, diagPolyA.lInfinityError(diagPolyB));

      uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
      uint_t inner_rowsize = rowsize;
      real_t ref_H = walberla::real_c(1)/(walberla::real_c(rowsize - 1));

      Point2D ref_x;

      for (uint_t j = 1; j < rowsize - 2; ++j) {
        ref_x[1] = j * ref_H;
        for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
          ref_x[0] = i * ref_H;

          error = std::max(error, std::abs(horiPolyA.eval(ref_x) - horiPolyB.eval(ref_x)));
          error = std::max(error, std::abs(vertPolyA.eval(ref_x) - vertPolyB.eval(ref_x)));
          error = std::max(error, std::abs(diagPolyA.eval(ref_x) - diagPolyB.eval(ref_x)));

        }

        --inner_rowsize;
      }
    }

    walberla::mpi::allReduceInplace(error, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm());

    return error;
  }

private:

  void apply_impl(P1Function< real_t >& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    std::vector<PrimitiveDataID<FunctionMemory< real_t >, Vertex>> vertexCoeffIds;
    std::vector<PrimitiveDataID<FunctionMemory< real_t >, Edge>> edgeCoeffIds;

    for (auto coefficient : coefficients_) {
      vertexCoeffIds.push_back(coefficient->getVertexDataID());
      edgeCoeffIds.push_back(coefficient->getEdgeDataID());
    }

    // start pulling vertex halos
    src.getCommunicator(level)->startCommunication<Edge, Vertex>();
    for (auto coefficient : coefficients_) {
      coefficient->getCommunicator(level)->startCommunication<Edge, Vertex>();
    }

    // start pulling edge halos
    src.getCommunicator(level)->startCommunication<Face, Edge>();
    for (auto coefficient : coefficients_) {
      coefficient->getCommunicator(level)->startCommunication<Face, Edge>();
    }

    // end pulling vertex halos
    for (auto coefficient : coefficients_) {
      coefficient->getCommunicator(level)->endCommunication<Edge, Vertex>();
    }
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
    for (auto coefficient : coefficients_) {
      coefficient->getCommunicator(level)->endCommunication<Face, Edge>();
    }
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
        vertexdof::macroface::applyPolynomial< real_t >(level, polyDegree_, face, facePolynomialID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_gs_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag)
  {
    std::vector<PrimitiveDataID<FunctionMemory< real_t >, Vertex>> vertexCoeffIds;
    std::vector<PrimitiveDataID<FunctionMemory< real_t >, Edge>> edgeCoeffIds;

    for (auto coefficient : coefficients_) {
      vertexCoeffIds.push_back(coefficient->getVertexDataID());
      edgeCoeffIds.push_back(coefficient->getEdgeDataID());
    }

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
        vertexdof::macroface::smooth_gs_polynomial< real_t >(level, polyDegree_, face, facePolynomialID_, dst.getFaceDataID(), rhs.getFaceDataID());
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
  PrimitiveDataID<FaceP1PolynomialMemory, Face> facePolynomialID_;

  void compute_local_stiffness(const fenics::TabulateTensor& opr, const Face &face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    opr(local_stiffness.data(), NULL, coords, 0);
  }

private:
  uint_t polyDegree_;
  uint_t interpolationLevel_;
  std::vector<std::shared_ptr<P1Function< real_t >>> coefficients_;
  std::vector<std::function<real_t(const hhg::Point3D&)>> analyticCoefficients_;
  std::vector<fenics::TabulateTensor> operators_;
};

template<class FenicsOperator1, class FenicsOperator2, class FenicsOperator3>
class P1PolynomialTensorBlendingOperator : public P1PolynomialBlendingOperator {
public:
  P1PolynomialTensorBlendingOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                                        size_t minLevel,
                                        size_t maxLevel,
                                        uint_t interpolationLevel)
      : P1PolynomialBlendingOperator(storage, minLevel, maxLevel, interpolationLevel)
  {
    fenicsOperator1 = std::make_shared<FenicsOperator1>();
    fenicsOperator2 = std::make_shared<FenicsOperator2>();
    fenicsOperator3 = std::make_shared<FenicsOperator3>();
    this->addOperator(std::bind(&FenicsOperator1::tabulate_tensor, fenicsOperator1.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
    this->addOperator(std::bind(&FenicsOperator2::tabulate_tensor, fenicsOperator2.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
    this->addOperator(std::bind(&FenicsOperator3::tabulate_tensor, fenicsOperator3.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));

    this->init();
  }
private:
  std::shared_ptr<FenicsOperator1> fenicsOperator1;
  std::shared_ptr<FenicsOperator2> fenicsOperator2;
  std::shared_ptr<FenicsOperator3> fenicsOperator3;
};

using P1PolynomialBlendingLaplaceOperator = P1PolynomialTensorBlendingOperator<p1_div_k_grad_cell_integral_0_otherwise,
    p1_div_k_grad_cell_integral_1_otherwise,
    p1_div_k_grad_cell_integral_2_otherwise>;

}


