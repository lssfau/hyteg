
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

#include "tinyhhg_core/polynomial/lsqinterpolation.hpp"

namespace hhg
{

template<class UFCOperator, uint_t MaxPolyDegree, uint_t InterpolationLevel>
class P1PolynomialOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  typedef LSQInterpolator<MaxPolyDegree, InterpolationLevel, HorizontalEdgeBasis> HorizontalEdgeInterpolator;
  typedef LSQInterpolator<MaxPolyDegree, InterpolationLevel, VerticalEdgeBasis> VerticalEdgeInterpolator;

  P1PolynomialOperator(const std::shared_ptr< PrimitiveStorage > & storage, const std::shared_ptr<P1Function< real_t >>& coefficient, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel), coefficientP1_(coefficient)
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
    auto faceP1PolynomialMemoryDataHandling = std::make_shared< FaceP1PolynomialMemoryDataHandling<MaxPolyDegree, InterpolationLevel> >();
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

    std::array<SD,3> triangleBlueSW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_S  };
    std::array<SD,3> triangleGrayS  = { SD::VERTEX_C, SD::VERTEX_S,  SD::VERTEX_SE };
    std::array<SD,3> triangleBlueSE = { SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E  };
    std::array<SD,3> triangleGrayNW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_NW };
    std::array<SD,3> triangleBlueN  = { SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N  };
    std::array<SD,3> triangleGrayNE = { SD::VERTEX_C, SD::VERTEX_N,  SD::VERTEX_E  };

    std::vector<real_t> horiValues(HorizontalEdgeInterpolator::NumVertices);
    std::vector<real_t> vertValues(VerticalEdgeInterpolator::NumVertices);
//    std::vector<real_t> diagValues(Interpolator::NumVertices);

    std::vector<real_t> faceStencil(7);

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      auto facePolynomials = face.getData(facePolynomialID_);
      auto faceLocalMatrices = face.getData(faceLocalMatrixID_);
      auto coeff = face.getData(coefficientP1_->getFaceDataID())->getPointer(InterpolationLevel);

      uint_t rowsize = levelinfo::num_microvertices_per_edge(InterpolationLevel);
      uint_t inner_rowsize = rowsize;
      uint_t horiOffset = 0;
      uint_t vertOffset = 0;
      real_t coeffWeight;

      for (uint_t j = 1; j < rowsize - 2; ++j) {

        uint_t i;
        for (i = 1; i < inner_rowsize - 2; ++i) {

          std::fill(faceStencil.begin(), faceStencil.end(), walberla::real_c(0.0));

          for (uint_t k = 0; k < FaceVertexDoF::P1GrayElements.size(); ++k) {
            coeffWeight = 1.0/3.0 * (coeff[vertexdof::macroface::indexFromVertex<InterpolationLevel>(i, j, FaceVertexDoF::P1GrayElements[k][0])]
                            + coeff[vertexdof::macroface::indexFromVertex<InterpolationLevel>(i, j, FaceVertexDoF::P1GrayElements[k][1])]
                            + coeff[vertexdof::macroface::indexFromVertex<InterpolationLevel>(i, j, FaceVertexDoF::P1GrayElements[k][2])]);
            assembleP1LocalStencil(FaceVertexDoF::P1GrayStencilMaps[k], FaceVertexDoF::P1GrayDoFMaps[k], faceLocalMatrices->getGrayMatrix(InterpolationLevel), faceStencil, coeffWeight);
          }

          for (uint_t k = 0; k < FaceVertexDoF::P1BlueElements.size(); ++k) {
            coeffWeight = 1.0/3.0 * (coeff[vertexdof::macroface::indexFromVertex<InterpolationLevel>(i, j, FaceVertexDoF::P1BlueElements[k][0])]
                                     + coeff[vertexdof::macroface::indexFromVertex<InterpolationLevel>(i, j, FaceVertexDoF::P1BlueElements[k][1])]
                                     + coeff[vertexdof::macroface::indexFromVertex<InterpolationLevel>(i, j, FaceVertexDoF::P1BlueElements[k][2])]);
            assembleP1LocalStencil(FaceVertexDoF::P1BlueStencilMaps[k], FaceVertexDoF::P1BlueDoFMaps[k], faceLocalMatrices->getBlueMatrix(InterpolationLevel), faceStencil, coeffWeight);
          }

//          WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("FACE.id = {}:face_stencil = {}", face.getID().getID(), PointND<real_t, 7>(&faceStencil[0])));

          horiValues[horiOffset] = faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_W)];
          ++horiOffset;

          vertValues[vertOffset] = faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_S)];
          ++vertOffset;

          if (i == inner_rowsize - 2 - 1) {
            horiValues[horiOffset] = faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_E)];
            ++horiOffset;
          }
        }

        vertValues[vertOffset] = faceStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_N)];
        ++vertOffset;

        --inner_rowsize;
      }

      HorizontalEdgeInterpolator horiInterpolator;
      horiInterpolator.interpolate(horiValues, facePolynomials->getHoriPolynomial());

      VerticalEdgeInterpolator vertInterpolator;
      vertInterpolator.interpolate(vertValues, facePolynomials->getVertPolynomial());
//      interpolator.interpolate(diagValues, facePolynomials->getDiagPolynomial());

      WALBERLA_LOG_DEVEL("polynomials[0] = " << facePolynomials->getHoriPolynomial());
      WALBERLA_LOG_DEVEL("polynomials[1] = " << facePolynomials->getVertPolynomial());
//      WALBERLA_LOG_DEVEL("polynomials[2] = " << facePolynomials->getDiagPolynomial());
      WALBERLA_LOG_INFO("Warning: Diagonal polynomials not yet implemented!")
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
        vertexdof::macroface::applyPolynomial< real_t, MaxPolyDegree, InterpolationLevel >(level, face, facePolynomialID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
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
  PrimitiveDataID<FaceP1PolynomialMemory<MaxPolyDegree, InterpolationLevel>, Face> facePolynomialID_;

  void compute_local_stiffness(const Face &face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
  }

private:
  std::shared_ptr<P1Function< real_t >> coefficientP1_;
};

template<uint_t MaxPolyDegree, uint_t InterpolationLevel>
using P1PolynomialLaplaceOperator = P1PolynomialOperator<p1_diffusion_cell_integral_0_otherwise, MaxPolyDegree, InterpolationLevel>;

}


