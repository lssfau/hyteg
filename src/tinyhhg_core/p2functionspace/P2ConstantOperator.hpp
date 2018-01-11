#pragma once

#include "P2Function.hpp"
#include "P2Elements.hpp"

#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"

#include "generated/p2_diffusion.h"

namespace hhg {

using walberla::real_t;

template<class UFCOperator>
class P2ConstantOperator : public Operator<P2Function < real_t>, P2Function<real_t> > {
public:

  P2ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel), vertexToVertex(storage, minLevel, maxLevel),
        edgeToVertex(storage, minLevel, maxLevel)
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
        real_t* vStencil = vertexToVertex.getFaceStencil(face.getID(), level);

        for (uint_t i = 0; i < P2Face::P2GrayElements.size(); ++i) {
          P2Face::VertexToVertex::assembleStencil(P2Face::VertexToVertex::P2GrayStencilMaps[i], P2Face::VertexToVertex::P2GrayDoFMaps[i], local_stiffness_gray, vStencil);
        }

        for (uint_t i = 0; i < P2Face::P2BlueElements.size(); ++i) {
          P2Face::VertexToVertex::assembleStencil(P2Face::VertexToVertex::P2BlueStencilMaps[i], P2Face::VertexToVertex::P2BlueDoFMaps[i], local_stiffness_blue, vStencil);
        }

        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Face = {}", PointND<real_t, 7>(&vStencil[0])));

        // Assemble edgeToVertex stencil
        vStencil = edgeToVertex.getFaceStencil(face.getID(), level);

        for (uint_t i = 0; i < P2Face::P2GrayElements.size(); ++i) {
          P2Face::EdgeToVertex::assembleStencil(P2Face::EdgeToVertex::P2GrayStencilMaps[i], P2Face::EdgeToVertex::P2GrayDoFMaps[i], local_stiffness_gray, vStencil);
        }

        for (uint_t i = 0; i < P2Face::P2BlueElements.size(); ++i) {
          P2Face::EdgeToVertex::assembleStencil(P2Face::EdgeToVertex::P2BlueStencilMaps[i], P2Face::EdgeToVertex::P2BlueDoFMaps[i], local_stiffness_blue, vStencil);
        }

        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Face = {}", PointND<real_t, 12>(&vStencil[0])));
      }

    }

  }

private:

  void apply_impl(P2Function< real_t > & src, P2Function< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    WALBERLA_ABORT("Not implemented!");
  }

  P1Operator<NoAssemble> vertexToVertex;
  EdgeDoFToVertexDoFOperator edgeToVertex;
//  VertexDoFToEdgeDoFOperator vertexToEdge;
//  EdgeDoFOperator edgeToEdge;

  void compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
  }

};

typedef P2ConstantOperator<p2_diffusion_cell_integral_0_otherwise> P2ConstantLaplaceOperator;

}