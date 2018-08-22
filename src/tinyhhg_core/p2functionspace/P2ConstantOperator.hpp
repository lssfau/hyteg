#pragma once

#include "P2Function.hpp"
#include "P2Elements.hpp"
#include "P2Smooth.hpp"
#include "P2MacroFace.hpp"
#include "P2MacroEdge.hpp"

#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "generated/p2_diffusion.h"
#include "generated/p2_mass.h"
#include "generated/p2_divt.h"
#include "generated/p2_div.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

namespace hhg {

using walberla::real_t;

template<class UFCOperator>
class P2ConstantOperator : public Operator<P2Function < real_t>, P2Function<real_t> > {
public:

  P2ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel), vertexToVertex(storage, minLevel, maxLevel),
        edgeToVertex(storage, minLevel, maxLevel),
        vertexToEdge(storage, minLevel, maxLevel),
        edgeToEdge(storage, minLevel, maxLevel)
  {
    using namespace P2Elements;

    // Initialize memory for local 6x6 matrices
    Matrix6r local_stiffness_gray;
    Matrix6r local_stiffness_blue;

    // Assemble stencils on all levels
    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      // Assemble face stencils
      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        // Compute both local stiffness matrices
        compute_local_stiffness(face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(face, level, local_stiffness_blue, fenics::BLUE);

//        WALBERLA_LOG_DEVEL_ON_ROOT("local_stiffness_gray =\n" << local_stiffness_gray);
//        WALBERLA_LOG_DEVEL_ON_ROOT("local_stiffness_blue =\n" << local_stiffness_blue);

        // Assemble vertexToVertex stencil
        real_t* vStencil = storage_->getFace(face.getID())->getData(vertexToVertex.getFaceStencilID())->getPointer(level);
        P2Face::VertexToVertex::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Face = {}", PointND<real_t, 7>(&vStencil[0])));

        // Assemble edgeToVertex stencil
        vStencil = storage_->getFace(face.getID())->getData(edgeToVertex.getFaceStencilID())->getPointer(level);
        P2Face::EdgeToVertex::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Face = {}", PointND<real_t, 12>(&vStencil[0])));

        // Assemble vertexToEdge stencil
        vStencil = storage_->getFace(face.getID())->getData(vertexToEdge.getFaceStencilID())->getPointer(level);
        P2Face::VertexToEdge::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Face = {}", PointND<real_t, 12>(&vStencil[0])));

        // Assemble edgeToEdge stencil
        vStencil = storage_->getFace(face.getID())->getData(edgeToEdge.getFaceStencilID())->getPointer(level);
        P2Face::EdgeToEdge::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToEdge/Face = {}", PointND<real_t, 15>(&vStencil[0])));
      }

      // Assemble edge stencils
      for (auto& it : storage_->getEdges()) {
        Edge &edge = *it.second;

        // Assemble vertexToVertex stencil
        Face *face = storage_->getFace(edge.neighborFaces()[0]);
        real_t* vStencil = storage_->getEdge(edge.getID())->getData(vertexToVertex.getEdgeStencilID())->getPointer(level);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::VertexToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

        if (edge.getNumNeighborFaces() == 2) {
          face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
          P2Edge::VertexToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Edge = {}", PointND<real_t, 7>(&vStencil[0])));

        // Assemble edgeToVertex
        face = storage_->getFace(edge.neighborFaces()[0]);
        vStencil = storage_->getEdge(edge.getID())->getData(edgeToVertex.getEdgeStencilID())->getPointer(level);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::EdgeToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

        if (edge.getNumNeighborFaces() == 2) {
          face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
          P2Edge::EdgeToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Edge = {}", PointND<real_t, 7>(&vStencil[0])));

        // Assemble vertexToEdge stencil
        face = storage_->getFace(edge.neighborFaces()[0]);
        vStencil = storage_->getEdge(edge.getID())->getData(vertexToEdge.getEdgeStencilID())->getPointer(level);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::VertexToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

        if (edge.getNumNeighborFaces() == 2) {
          face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
          P2Edge::VertexToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Edge = {}", PointND<real_t, 4>(&vStencil[0])));

        // Assemble edgeToEdge stencil
        face = storage_->getFace(edge.neighborFaces()[0]);
        vStencil = storage_->getEdge(edge.getID())->getData(edgeToEdge.getEdgeStencilID())->getPointer(level);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::EdgeToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

        if (edge.getNumNeighborFaces() == 2) {
          face = storage_->getFace(edge.neighborFaces()[1]);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
          P2Edge::EdgeToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToEdge/Edge = {}", PointND<real_t, 5>(&vStencil[0])));
      }

      for (auto& it : storage_->getVertices()) {
        Vertex &vertex = *it.second;

        // Assemble VertexToVertex
        real_t* vStencil = storage_->getVertex(vertex.getID())->getData(vertexToVertex.getVertexStencilID())->getPointer(level);
        for (auto& faceId : vertex.neighborFaces())
        {
          Face* face = storage_->getFace(faceId);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          P2Vertex::VertexToVertex::assembleStencil(vertex, *face, local_stiffness_gray, vStencil, storage_);
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Vertex = {}", PointND<real_t, 5>(&vStencil[0])));

        // Assemble EdgeToVertex
        vStencil = storage_->getVertex(vertex.getID())->getData(edgeToVertex.getVertexStencilID())->getPointer(level);
        for (auto& faceId : vertex.neighborFaces())
        {
          Face* face = storage_->getFace(faceId);
          compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
          P2Vertex::EdgeToVertex::assembleStencil(vertex, *face, local_stiffness_gray, vStencil, storage_);
        }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Vertex = {}", PointND<real_t, 5>(&vStencil[0])));
      }

    }

  }

  P1ConstantOperator<fenics::NoAssemble>& getVertexToVertexOpr() {
    return vertexToVertex;
  }

  GenericEdgeDoFToVertexDoFOperator& getEdgeToVertexOpr() {
    return edgeToVertex;
  }

  GenericVertexDoFToEdgeDoFOperator& getVertexToEdgeOpr() {
    return vertexToEdge;
  }

  EdgeDoFOperator& getEdgeToEdgeOpr() {
    return edgeToEdge;
  }

private:

  void apply_impl(P2Function< real_t > & src, P2Function< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    vertexToVertex.apply(*src.getVertexDoFFunction(), *dst.getVertexDoFFunction(), level, flag, updateType);
    edgeToVertex.apply(*src.getEdgeDoFFunction(), *dst.getVertexDoFFunction(), level, flag, Add);

    edgeToEdge.apply(*src.getEdgeDoFFunction(), *dst.getEdgeDoFFunction(), level, flag, updateType);
    vertexToEdge.apply(*src.getVertexDoFFunction(), *dst.getEdgeDoFFunction(), level, flag, Add);
  }

  void smooth_gs_impl(P2Function< real_t > & dst, P2Function< real_t > & rhs, size_t level, DoFType flag) override
  {
    dst.getVertexDoFFunction()->communicate<Face, Edge>( level );
    dst.getVertexDoFFunction()->communicate<Edge, Vertex>( level );
    dst.getEdgeDoFFunction()->communicate<Face, Edge>( level );
    dst.getEdgeDoFFunction()->communicate<Edge, Vertex>( level );

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if (testFlag(vertexBC, flag))
      {
        P2::vertex::smoothGSvertexDoF(level, vertex,
                                  vertexToVertex.getVertexStencilID(), dst.getVertexDoFFunction()->getVertexDataID(),
                                  edgeToVertex.getVertexStencilID(), dst.getEdgeDoFFunction()->getVertexDataID(),
                                  rhs.getVertexDoFFunction()->getVertexDataID());
      }
    }

    dst.getVertexDoFFunction()->communicate<Vertex, Edge>( level );
    dst.getEdgeDoFFunction()->communicate<Vertex, Edge>( level );

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if (testFlag(edgeBC, flag))
      {
         P2::macroedge::smoothGaussSeidel(level,
                                          edge,
                                          vertexToVertex.getEdgeStencilID(),
                                          edgeToVertex.getEdgeStencilID(),
                                          dst.getVertexDoFFunction()->getEdgeDataID(),
                                          vertexToEdge.getEdgeStencilID(),
                                          edgeToEdge.getEdgeStencilID(),
                                          dst.getEdgeDoFFunction()->getEdgeDataID(),
                                          rhs.getVertexDoFFunction()->getEdgeDataID(),
                                          rhs.getEdgeDoFFunction()->getEdgeDataID());
      }
    }

    dst.getVertexDoFFunction()->communicate<Edge, Face>( level );
    dst.getEdgeDoFFunction()->communicate<Edge, Face>( level );

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if (testFlag(faceBC, flag))
      {
         P2::macroface::smoothGaussSeidel(level,
                                          face,
                                          vertexToVertex.getFaceStencilID(),
                                          edgeToVertex.getFaceStencilID(),
                                          dst.getVertexDoFFunction()->getFaceDataID(),
                                          vertexToEdge.getFaceStencilID(),
                                          edgeToEdge.getFaceStencilID(),
                                          dst.getEdgeDoFFunction()->getFaceDataID(),
                                          rhs.getVertexDoFFunction()->getFaceDataID(),
                                          rhs.getEdgeDoFFunction()->getFaceDataID());
      }
    }
  }

  void smooth_jac_impl(P2Function< real_t > & dst, P2Function< real_t > & rhs, P2Function< real_t > & src, size_t level, DoFType flag) override {
    ///TODO: remove unneccessary communication here
    src.getVertexDoFFunction()->communicate<Face, Edge>( level );
    src.getVertexDoFFunction()->communicate<Edge, Vertex>( level );
    src.getVertexDoFFunction()->communicate<Vertex, Edge>( level );
    src.getVertexDoFFunction()->communicate<Edge, Face>( level );
    src.getEdgeDoFFunction()->communicate<Face, Edge>( level );
    src.getEdgeDoFFunction()->communicate<Edge, Vertex>( level );
    src.getEdgeDoFFunction()->communicate<Vertex, Edge>( level );
    src.getEdgeDoFFunction()->communicate<Edge, Face>( level );

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if (testFlag(faceBC, flag))
      {
        P2::macroface::smoothJacobiVertexDoF(level,
                                             face,
                                             vertexToVertex.getFaceStencilID(),
                                             src.getVertexDoFFunction()->getFaceDataID(),
                                             dst.getVertexDoFFunction()->getFaceDataID(),
                                             edgeToVertex.getFaceStencilID(),
                                             src.getEdgeDoFFunction()->getFaceDataID(),
                                             rhs.getVertexDoFFunction()->getFaceDataID());
      }
    }
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if (testFlag(faceBC, flag))
      {
        P2::macroface::smoothJacobiEdgeDoF(level,
                                           face,
                                           vertexToEdge.getFaceStencilID(),
                                           src.getVertexDoFFunction()->getFaceDataID(),
                                           edgeToEdge.getFaceStencilID(),
                                           src.getEdgeDoFFunction()->getFaceDataID(),
                                           dst.getEdgeDoFFunction()->getFaceDataID(),
                                           rhs.getEdgeDoFFunction()->getFaceDataID());
      }
    }

  }

  P1ConstantOperator<fenics::NoAssemble> vertexToVertex;
  GenericEdgeDoFToVertexDoFOperator edgeToVertex;
  GenericVertexDoFToEdgeDoFOperator vertexToEdge;
  EdgeDoFOperator edgeToEdge;

  void compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
  }

};

typedef P2ConstantOperator<p2_diffusion_cell_integral_0_otherwise> P2ConstantLaplaceOperator;
typedef P2ConstantOperator<p2_mass_cell_integral_0_otherwise> P2ConstantMassOperator;

typedef P2ConstantOperator<p2_divt_cell_integral_0_otherwise> P2ConstantDivTxOperator;
typedef P2ConstantOperator<p2_divt_cell_integral_1_otherwise> P2ConstantDivTyOperator;
typedef P2ConstantOperator<p2_div_cell_integral_0_otherwise> P2ConstantDivxOperator;
typedef P2ConstantOperator<p2_div_cell_integral_1_otherwise> P2ConstantDivyOperator;

}
