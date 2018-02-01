
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

#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"

#include "VertexDoFMacroVertex.hpp"
#include "VertexDoFMacroEdge.hpp"
#include "VertexDoFMacroFace.hpp"

#include "tinyhhg_core/polynomial/LSQPInterpolator.hpp"

namespace hhg
{

class P1PolynomialOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  typedef LSQPInterpolator<MonomialBasis2D> Interpolator;

  P1PolynomialOperator(const std::vector<fenics::TabulateTensor>& operators,
                       const std::shared_ptr< PrimitiveStorage > & storage,
                       const std::vector<std::shared_ptr<P1Function< real_t >>>& coefficients,
                       const std::vector<std::function<real_t(const hhg::Point3D&)>>& analyticCoefficients,
                       uint_t minLevel,
                       uint_t maxLevel,
                       uint_t polyDegree,
                       uint_t interpolationLevel)
      : Operator(storage, minLevel, maxLevel),
        polyDegree_(polyDegree),
        interpolationLevel_(interpolationLevel),
        coefficients_(coefficients),
        analyticCoefficients_(analyticCoefficients)
  {
    WALBERLA_ASSERT_EQUAL(operators.size(), coefficients.size());
    WALBERLA_ASSERT_EQUAL(operators.size(), analyticCoefficients.size());
    initLocalStiffnessMatrices(operators);
    interpolateStencils(operators);
  }

  ~P1PolynomialOperator()
  {
  }

private:

  void initLocalStiffnessMatrices(const std::vector<fenics::TabulateTensor>& operators)
  {
    auto faceP1PolynomialMemoryDataHandling = std::make_shared< FaceP1PolynomialMemoryDataHandling >(polyDegree_);
    storage_->addFaceData(facePolynomialID_, faceP1PolynomialMemoryDataHandling, "P1OperatorFacePolynomial");

    vertexLocalMatrixIDs_.resize(operators.size());
    edgeLocalMatrixIDs_.resize(operators.size());
    faceLocalMatrixIDs_.resize(operators.size());

    for(uint_t oprIdx = 0; oprIdx < operators.size(); ++oprIdx) {

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

          compute_local_stiffness(operators[oprIdx], face, level, faceLocalMatrices->getGrayMatrix(interpolationLevel_), fenics::GRAY);
          compute_local_stiffness(operators[oprIdx], face, level, faceLocalMatrices->getBlueMatrix(interpolationLevel_), fenics::BLUE);
        }

        for (auto &it : storage_->getEdges()) {
          Edge &edge = *it.second;

          auto edgeLocalMatrices = edge.getData(edgeLocalMatrixIDs_[oprIdx]);

          // first face
          Face *face = storage_->getFace(edge.neighborFaces()[0]);
          compute_local_stiffness(operators[oprIdx], *face, level, edgeLocalMatrices->getGrayMatrix(level, 0), fenics::GRAY);
          compute_local_stiffness(operators[oprIdx], *face, level, edgeLocalMatrices->getBlueMatrix(level, 0), fenics::BLUE);


          if (edge.getNumNeighborFaces() == 2) {
            // second face
            Face *face = storage_->getFace(edge.neighborFaces()[1]);
            compute_local_stiffness(operators[oprIdx], *face, level, edgeLocalMatrices->getGrayMatrix(level, 1), fenics::GRAY);
            compute_local_stiffness(operators[oprIdx], *face, level, edgeLocalMatrices->getBlueMatrix(level, 1), fenics::BLUE);
          }
        }

        for (auto &it : storage_->getVertices()) {
          Vertex &vertex = *it.second;

          auto vertexLocalMatrices = vertex.getData(vertexLocalMatrixIDs_[oprIdx]);

          // iterate over adjacent faces
          uint_t neighborId = 0;
          for (auto &faceId : vertex.neighborFaces()) {
            Face *face = storage_->getFace(faceId);

            compute_local_stiffness(operators[oprIdx], *face, level, vertexLocalMatrices->getGrayMatrix(level, neighborId), fenics::GRAY);
            ++neighborId;
          }
        }
      }
    }
  }

  void interpolateStencils(const std::vector<fenics::TabulateTensor>& operators) {
    typedef stencilDirection SD;
    using namespace P1Elements;

    std::vector<real_t> faceStencil(7);

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      auto facePolynomials = face.getData(facePolynomialID_);

      uint_t rowsize = levelinfo::num_microvertices_per_edge(interpolationLevel_);
      uint_t rowsizeFine = levelinfo::num_microvertices_per_edge(maxLevel_);
      uint_t inner_rowsize = rowsize;
      real_t coeffWeight;

      Interpolator horiInterpolator(polyDegree_, interpolationLevel_);
      Interpolator vertInterpolator(polyDegree_, interpolationLevel_);
      Interpolator diagInterpolator(polyDegree_, interpolationLevel_);

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

      for (uint_t j = 1; j < rowsize - 2; ++j) {

        ref_x[1] = j * ref_H;

        x = x0;
        x += real_c(j)*d2 + d0;

        uint_t i;
        for (i = 1; i < inner_rowsize - 2; ++i) {

          ref_x[0] = i * ref_H;

          std::fill(faceStencil.begin(), faceStencil.end(), walberla::real_c(0.0));

          for(uint_t oprIdx = 0; oprIdx < operators.size(); ++oprIdx) {

            auto faceLocalMatrices = face.getData(faceLocalMatrixIDs_[oprIdx]);

            // elementS
            coeffWeight = 1.0 / 3.0 * (analyticCoefficients_[oprIdx](x) + analyticCoefficients_[oprIdx](x - d2f) +
                                       analyticCoefficients_[oprIdx](x + d0f - d2f));
            assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[0], FaceVertexDoF::P1GrayDoFMaps[0],
                                   faceLocalMatrices->getGrayMatrix(interpolationLevel_), faceStencil, coeffWeight);

            // elementNE
            coeffWeight = 1.0 / 3.0 * (analyticCoefficients_[oprIdx](x) + analyticCoefficients_[oprIdx](x + d0f) +
                                       analyticCoefficients_[oprIdx](x + d2f));
            assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[1], FaceVertexDoF::P1GrayDoFMaps[1],
                                   faceLocalMatrices->getGrayMatrix(interpolationLevel_), faceStencil, coeffWeight);

            // elementNW
            coeffWeight = 1.0 / 3.0 * (analyticCoefficients_[oprIdx](x) + analyticCoefficients_[oprIdx](x - d0f + d2f) +
                                       analyticCoefficients_[oprIdx](x - d0f));
            assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[2], FaceVertexDoF::P1GrayDoFMaps[2],
                                   faceLocalMatrices->getGrayMatrix(interpolationLevel_), faceStencil, coeffWeight);

            // elementSW
            coeffWeight = 1.0 / 3.0 * (analyticCoefficients_[oprIdx](x) + analyticCoefficients_[oprIdx](x - d0f) +
                                       analyticCoefficients_[oprIdx](x - d2f));
            assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[0], FaceVertexDoF::P1BlueDoFMaps[0],
                                   faceLocalMatrices->getBlueMatrix(interpolationLevel_), faceStencil, coeffWeight);

            // elementSE
            coeffWeight = 1.0 / 3.0 * (analyticCoefficients_[oprIdx](x) + analyticCoefficients_[oprIdx](x + d0f - d2f) +
                                       analyticCoefficients_[oprIdx](x + d0f));
            assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[1], FaceVertexDoF::P1BlueDoFMaps[1],
                                   faceLocalMatrices->getBlueMatrix(interpolationLevel_), faceStencil, coeffWeight);

            // elementN
            coeffWeight = 1.0 / 3.0 * (analyticCoefficients_[oprIdx](x) + analyticCoefficients_[oprIdx](x + d2f) +
                                       analyticCoefficients_[oprIdx](x - d0f + d2f));
            assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[2], FaceVertexDoF::P1BlueDoFMaps[2],
                                   faceLocalMatrices->getBlueMatrix(interpolationLevel_), faceStencil, coeffWeight);

          }

//          if (i == 1 && j == 1) {
//            WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("nodal = {}", PointND<real_t, 7>(&faceStencil[0])));
//          }

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

      horiInterpolator.interpolate(facePolynomials->getHoriPolynomial(polyDegree_));
      vertInterpolator.interpolate(facePolynomials->getVertPolynomial(polyDegree_));
      diagInterpolator.interpolate(facePolynomials->getDiagPolynomial(polyDegree_));

//      std::fill(faceStencil.begin(), faceStencil.end(), walberla::real_c(0.0));
//      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_W)] = facePolynomials->getHoriPolynomial().eval(Point2D({ 1 * ref_H - 0.5 * ref_h, 1*ref_H  }));
//      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_E)] = facePolynomials->getHoriPolynomial().eval(Point2D({ 1 * ref_H + 0.5 * ref_h, 1*ref_H  }));
//      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_S)] = facePolynomials->getVertPolynomial().eval(Point2D({ 1 * ref_H, 1*ref_H - 0.5 * ref_h  }));
//      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_N)] = facePolynomials->getVertPolynomial().eval(Point2D({ 1 * ref_H, 1*ref_H + 0.5 * ref_h  }));
//      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_SE)] = facePolynomials->getDiagPolynomial().eval(Point2D({ 1 * ref_H + 0.5 * ref_h, 1*ref_H - 0.5 * ref_h  }));
//      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_NW)] = facePolynomials->getDiagPolynomial().eval(Point2D({ 1 * ref_H - 0.5 * ref_h, 1*ref_H + 0.5 * ref_h  }));
//      faceStencil[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] = - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]
//                                                                                - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]
//                                                                                - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]
//                                                                                - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]
//                                                                                - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]
//                                                                                - faceStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])];
//      WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("poly  = {}", PointND<real_t, 7>(&faceStencil[0])));

//      WALBERLA_LOG_DEVEL("polynomials[0] = " << facePolynomials->getHoriPolynomial());
//      WALBERLA_LOG_DEVEL("polynomials[1] = " << facePolynomials->getVertPolynomial());
//      WALBERLA_LOG_DEVEL("polynomials[2] = " << facePolynomials->getDiagPolynomial());
    }
  }

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
        vertexdof::macrovertex::applyCoefficient< real_t >(vertex, storage_, vertexLocalMatrixIDs_, src.getVertexDataID(), dst.getVertexDataID(), vertexCoeffIds, level, updateType);
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
        vertexdof::macroedge::applyCoefficient< real_t >(level, edge, storage_, edgeLocalMatrixIDs_, src.getEdgeDataID(), dst.getEdgeDataID(), edgeCoeffIds, updateType);
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
        vertexdof::macrovertex::smooth_gs_coefficient(vertex, storage_, vertexLocalMatrixIDs_, dst.getVertexDataID(), rhs.getVertexDataID(), vertexCoeffIds, level);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    dst.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        vertexdof::macroedge::smooth_gs_coefficient<real_t>(level, edge, storage_, edgeLocalMatrixIDs_, dst.getEdgeDataID(), rhs.getEdgeDataID(), edgeCoeffIds);
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
  const std::vector<std::function<real_t(const hhg::Point3D&)>>& analyticCoefficients_;
};

template<class FenicsOperator>
class P1PolynomialScalarCoefficientOperator : public P1PolynomialOperator {
public:
  P1PolynomialScalarCoefficientOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                                        const std::shared_ptr<P1Function< real_t >>& coefficient,
                                        const std::function<real_t(const hhg::Point3D&)>& analyticCoefficient,
                                        size_t minLevel,
                                        size_t maxLevel,
                                        uint_t polyDegree,
                                        uint_t interpolationLevel)
      : P1PolynomialOperator(fillOperators(),
                             storage,
                             std::vector<std::shared_ptr<P1Function< real_t >>>{coefficient},
                             std::vector<std::function<real_t(const hhg::Point3D&)>>{analyticCoefficient},
                             minLevel, maxLevel, polyDegree, interpolationLevel)
  { }

  std::vector<fenics::TabulateTensor> fillOperators() {
    fenicsOperator = std::make_shared<FenicsOperator>();
    std::vector<fenics::TabulateTensor> operators;
    operators.push_back(std::bind(&FenicsOperator::tabulate_tensor, fenicsOperator.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
    return operators;
  }
private:
  std::shared_ptr<FenicsOperator> fenicsOperator;
};

template<class FenicsOperator1, class FenicsOperator2, class FenicsOperator3>
class P1PolynomialTensorCoefficientOperator : public P1PolynomialOperator {
public:
  P1PolynomialTensorCoefficientOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                                        const std::vector<std::shared_ptr<P1Function< real_t >>>& coefficients,
                                        const std::vector<std::function<real_t(const hhg::Point3D&)>>& analyticCoefficients,
                                        size_t minLevel,
                                        size_t maxLevel,
                                        uint_t polyDegree,
                                        uint_t interpolationLevel)
      : P1PolynomialOperator(fillOperators(),
                              storage, coefficients, analyticCoefficients, minLevel, maxLevel, polyDegree, interpolationLevel)
  { }

  std::vector<fenics::TabulateTensor> fillOperators() {
    fenicsOperator1 = std::make_shared<FenicsOperator1>();
    fenicsOperator2 = std::make_shared<FenicsOperator2>();
    fenicsOperator3 = std::make_shared<FenicsOperator3>();
    std::vector<fenics::TabulateTensor> operators;
    operators.push_back(std::bind(&FenicsOperator1::tabulate_tensor, fenicsOperator1.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
    operators.push_back(std::bind(&FenicsOperator2::tabulate_tensor, fenicsOperator2.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
    operators.push_back(std::bind(&FenicsOperator3::tabulate_tensor, fenicsOperator3.get(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4));
    return operators;
  }
private:
  std::shared_ptr<FenicsOperator1> fenicsOperator1;
  std::shared_ptr<FenicsOperator2> fenicsOperator2;
  std::shared_ptr<FenicsOperator3> fenicsOperator3;
};

using P1PolynomialScalarCoefficientLaplaceOperator = P1PolynomialScalarCoefficientOperator<p1_diffusion_cell_integral_0_otherwise>;

}


