
#pragma once

#include <fmt/format.h>
#include "tinyhhg_core/Operator.hpp"
#include <tinyhhg_core/p1functionspace/P1Operator.hpp>
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/matrix.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics.hpp"

#include "tinyhhg_core/p2functionspace/generated/p2_diffusion.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

//#include "tinyhhg_core/p1functionspace/P2Memory.hpp"

//#include "P2Vertex.hpp"
//#include "P2Edge.hpp"
//#include "P2Face.hpp"

namespace hhg
{

template<class UFCOperator,  bool Diagonal = false>
class P2Operator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
public:
  P2Operator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  ,vToEOperator(storage, minLevel, maxLevel)
  {

    for (uint_t level = minLevel_; level <= maxLevel_; ++level) {
      for (auto &it : storage_->getFaces()) {
        Face &face = *it.second;

        Matrix6r local_stiffness_gray;
        Matrix6r local_stiffness_blue;

        WALBERLA_ABORT("implement me");
      }
    }


  }

  ~P2Operator()
  {
  }

  //P1Operator vToVOperator;
  VertexDoFToEdgeDoFOperator vToEOperator;
//  EdgeDoFToVertexDoFOperator eToVOperator;
//  EdgeDoFOperator eToEOperator;

private:

  void apply_impl(P2Function< real_t > & src, P2Function< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    WALBERLA_ABORT("not implemented");
  }

  void smooth_gs_impl(P2Function< real_t > & dst, P2Function< real_t > & rhs, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

  void smooth_sor_impl(P2Function< real_t > & dst, P2Function< real_t > & rhs, real_t relax, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

  void smooth_jac_impl(P2Function< real_t > & dst, P2Function< real_t > & rhs, P2Function< real_t > & tmp, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

  void compute_local_stiffness(const Face &face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);

    if (Diagonal) {
      for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
          if (i != j) {
            local_stiffness(i,j) = real_t(0);
          }
        }
      }
    }
  }

};

typedef P2Operator<p2_diffusion_cell_integral_0_otherwise> P2LaplaceOperator;

}


