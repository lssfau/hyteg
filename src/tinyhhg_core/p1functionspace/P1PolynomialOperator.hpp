
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

#include "P1Vertex.hpp"
#include "P1Edge.hpp"
#include "P1Face.hpp"

#include "tinyhhg_core/polynomial/LSQPInterpolator.hpp"

namespace hhg
{

template<class UFCOperator, uint_t PolyDegree, uint_t InterpolationLevel>
class P1PolynomialOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  typedef LSQPInterpolator<PolyDegree, InterpolationLevel, MonomialBasis2D> Interpolator;

  P1PolynomialOperator(const std::shared_ptr< PrimitiveStorage > & storage, const std::shared_ptr<P1Function< real_t >>& coefficient, const std::function<real_t(const hhg::Point3D&)>& analyticCoefficient, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel), coefficientP1_(coefficient), analyticCoefficient_(analyticCoefficient)
  {
    initLocalStiffnessMatrices();
    interpolateStencils();
  }

  ~P1PolynomialOperator()
  {
  }

private:

  void initLocalStiffnessMatrices()
  {
    auto faceP1PolynomialMemoryDataHandling = std::make_shared< FaceP1PolynomialMemoryDataHandling<PolyDegree> >();
    auto faceP1LocalMatrixMemoryDataHandling = std::make_shared< FaceP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
    auto edgeP1LocalMatrixMemoryDataHandling = std::make_shared< EdgeP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
    auto vertexP1LocalMatrixMemoryDataHandling = std::make_shared< VertexP1LocalMatrixMemoryDataHandling >(minLevel_, maxLevel_);
    storage_->addFaceData(facePolynomialID_, faceP1PolynomialMemoryDataHandling, "P1OperatorFacePolynomial");
    storage_->addFaceData(faceLocalMatrixID_, faceP1LocalMatrixMemoryDataHandling, "P1OperatorFaceLocalMatrix");
    storage_->addEdgeData(edgeLocalMatrixID_, edgeP1LocalMatrixMemoryDataHandling, "P1OperatorEdgeLocalMatrix");
    storage_->addVertexData(vertexLocalMatrixID_, vertexP1LocalMatrixMemoryDataHandling, "P1OperatorVertexLocalMatrix");

    // Initialize face stiffness matrices only on InterpolationLevel
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      auto faceLocalMatrices = face.getData(faceLocalMatrixID_);

      compute_local_stiffness(face, InterpolationLevel, faceLocalMatrices->getGrayMatrix(InterpolationLevel), fenics::GRAY);
      compute_local_stiffness(face, InterpolationLevel, faceLocalMatrices->getBlueMatrix(InterpolationLevel), fenics::BLUE);
    }

    // Initialize other local stiffness matrices on lower dimensional primitives as usual
    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {
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

  void interpolateStencils() {
    typedef stencilDirection SD;
    using namespace P1Elements;

    std::vector<real_t> horiValues(Interpolator::NumVertices);
    std::vector<real_t> vertValues(Interpolator::NumVertices);
    std::vector<real_t> diagValues(Interpolator::NumVertices);

    std::vector<real_t> faceStencil(7);

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      auto facePolynomials = face.getData(facePolynomialID_);

      // compute stiffness matrices on maxLevel
      Matrix3r local_stiffness_gray;
      Matrix3r local_stiffness_blue;

      compute_local_stiffness(face, maxLevel_, local_stiffness_gray, fenics::GRAY);
      compute_local_stiffness(face, maxLevel_, local_stiffness_blue, fenics::BLUE);

      uint_t rowsize = levelinfo::num_microvertices_per_edge(InterpolationLevel);
      uint_t rowsizeFine = levelinfo::num_microvertices_per_edge(maxLevel_);
      uint_t inner_rowsize = rowsize;
      real_t coeffWeight;

      Interpolator horiInterpolator;
      Interpolator vertInterpolator;
      Interpolator diagInterpolator;

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

          // elementS
          coeffWeight = analyticCoefficient_(x + 1./3. * d0f - 2./3. * d2f );
          assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[0], FaceVertexDoF::P1GrayDoFMaps[0], local_stiffness_gray, faceStencil, coeffWeight);

          // elementNE
          coeffWeight = analyticCoefficient_(x + 1./3. * d0f + 1./3. * d2f);
          assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[1], FaceVertexDoF::P1GrayDoFMaps[1], local_stiffness_gray, faceStencil, coeffWeight);

          // elementNW
          coeffWeight = analyticCoefficient_(x - 2./3. * d0f + 1./3. * d2f);
          assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[2], FaceVertexDoF::P1GrayDoFMaps[2], local_stiffness_gray, faceStencil, coeffWeight);

          // elementSW
          coeffWeight = analyticCoefficient_(x - 1./3. * d0f - 1./3. * d2f);
          assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[0], FaceVertexDoF::P1BlueDoFMaps[0], local_stiffness_blue, faceStencil, coeffWeight);

          // elementSE
          coeffWeight = analyticCoefficient_(x + 2./3. * d0f - 1./3. * d2f);
          assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[1], FaceVertexDoF::P1BlueDoFMaps[1], local_stiffness_blue, faceStencil, coeffWeight);

          // elementN
          coeffWeight = analyticCoefficient_(x - 1./3. * d0f + 2./3. * d2f);
          assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[2], FaceVertexDoF::P1BlueDoFMaps[2], local_stiffness_blue, faceStencil, coeffWeight);

//          WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("FACE.id = {}:face_stencil = {}", face.getID().getID(), PointND<real_t, 7>(&faceStencil[0])));

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


      horiInterpolator.interpolate(facePolynomials->getHoriPolynomial());
      vertInterpolator.interpolate(facePolynomials->getVertPolynomial());
      diagInterpolator.interpolate(facePolynomials->getDiagPolynomial());


//      WALBERLA_LOG_DEVEL("polynomials[0] = " << facePolynomials->getHoriPolynomial());
//      WALBERLA_LOG_DEVEL("polynomials[1] = " << facePolynomials->getVertPolynomial());
//      WALBERLA_LOG_DEVEL("polynomials[2] = " << facePolynomials->getDiagPolynomial());
    }
  }

  void apply_impl(P1Function< real_t >& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    // start pulling vertex halos
    src.getCommunicator(level)->startCommunication<Edge, Vertex>();
    coefficientP1_->getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    src.getCommunicator(level)->startCommunication<Face, Edge>();
    coefficientP1_->getCommunicator(level)->startCommunication<Face, Edge>();

    // end pulling vertex halos
    coefficientP1_->getCommunicator(level)->endCommunication<Edge, Vertex>();
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
    coefficientP1_->getCommunicator(level)->endCommunication<Face, Edge>();
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
        vertexdof::macroface::applyPolynomial< real_t, PolyDegree >(level, face, facePolynomialID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
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
        vertexdof::macrovertex::smooth_gs_coefficient(vertex, storage_, vertexLocalMatrixID_, dst.getVertexDataID(), rhs.getVertexDataID(), coefficientP1_->getVertexDataID(), level);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    dst.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        vertexdof::macroedge::smooth_gs_coefficient<real_t>(level, edge, storage_, edgeLocalMatrixID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), coefficientP1_->getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        vertexdof::macroface::smooth_gs_polynomial< real_t, PolyDegree >(level, face, facePolynomialID_, dst.getFaceDataID(), rhs.getFaceDataID());
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
  PrimitiveDataID<FaceP1PolynomialMemory<PolyDegree>, Face> facePolynomialID_;

  void compute_local_stiffness(const Face &face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
  }

private:
  std::shared_ptr<P1Function< real_t >> coefficientP1_;
  const std::function<real_t(const hhg::Point3D&)>& analyticCoefficient_;
};

template<uint_t PolyDegree, uint_t InterpolationLevel>
using P1PolynomialLaplaceOperator = P1PolynomialOperator<p1_diffusion_cell_integral_0_otherwise, PolyDegree, InterpolationLevel>;

}


