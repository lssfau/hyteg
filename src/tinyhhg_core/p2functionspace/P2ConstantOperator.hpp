#pragma once

#include "P2Function.hpp"
#include "P2Elements.hpp"

#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOperator.hpp"

#include "generated/p2_diffusion.h"

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

  const P1Operator<NoAssemble>& getVertexToVertexOpr() const {
    return vertexToVertex;
  }

  const EdgeDoFToVertexDoFOperator& getEdgeToVertexOpr() const {
    return edgeToVertex;
  }

  const VertexDoFToEdgeDoFOperator& getVertexToEdgeOpr() const {
    return vertexToEdge;
  }

  const EdgeDoFOperator& getEdgeToEdgeOpr() const {
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

  P1Operator<NoAssemble> vertexToVertex;
  EdgeDoFToVertexDoFOperator edgeToVertex;
  VertexDoFToEdgeDoFOperator vertexToEdge;
  EdgeDoFOperator edgeToEdge;

  void compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
  }

};

typedef P2ConstantOperator<p2_diffusion_cell_integral_0_otherwise> P2ConstantLaplaceOperator;

}
