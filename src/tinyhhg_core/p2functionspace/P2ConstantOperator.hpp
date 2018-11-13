#pragma once

#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"

class p2_diffusion_cell_integral_0_otherwise;
class p2_tet_diffusion_cell_integral_0_otherwise;
class p2_mass_cell_integral_0_otherwise;
class p2_divt_cell_integral_0_otherwise;
class p2_divt_cell_integral_1_otherwise;
class p2_div_cell_integral_0_otherwise;
class p2_div_cell_integral_1_otherwise;

namespace hhg {

using walberla::real_t;

template< class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P2ConstantOperator : public Operator<P2Function < real_t>, P2Function<real_t> >
{
public:

  P2ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel);

  P1ConstantOperator< fenics::NoAssemble, UFCOperator3D >& getVertexToVertexOpr() {
    return vertexToVertex;
  }

  EdgeDoFToVertexDoFOperator< fenics::NoAssemble, UFCOperator3D > & getEdgeToVertexOpr() {
    return edgeToVertex;
  }

  VertexDoFToEdgeDoFOperator< fenics::NoAssemble, UFCOperator3D > & getVertexToEdgeOpr() {
    return vertexToEdge;
  }

  EdgeDoFOperator& getEdgeToEdgeOpr() {
    return edgeToEdge;
  }

private:

  void assembleStencils();

  void assembleStencils3D();

  void apply_impl(P2Function< real_t > & src, P2Function< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace) override;

  void smooth_gs_impl(P2Function< real_t > & dst, P2Function< real_t > & rhs, size_t level, DoFType flag) override;

  void smooth_jac_impl(P2Function< real_t > & dst, P2Function< real_t > & rhs, P2Function< real_t > & src, size_t level, DoFType flag) override;

  P1ConstantOperator< fenics::NoAssemble, UFCOperator3D > vertexToVertex;
  EdgeDoFToVertexDoFOperator< fenics::NoAssemble, UFCOperator3D > edgeToVertex;
  VertexDoFToEdgeDoFOperator< fenics::NoAssemble, UFCOperator3D > vertexToEdge;
  EdgeDoFOperator edgeToEdge;

  void compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type);

};

typedef P2ConstantOperator<p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise> P2ConstantLaplaceOperator;
typedef P2ConstantOperator<p2_mass_cell_integral_0_otherwise> P2ConstantMassOperator;

typedef P2ConstantOperator<p2_divt_cell_integral_0_otherwise> P2ConstantDivTxOperator;
typedef P2ConstantOperator<p2_divt_cell_integral_1_otherwise> P2ConstantDivTyOperator;
typedef P2ConstantOperator<p2_div_cell_integral_0_otherwise> P2ConstantDivxOperator;
typedef P2ConstantOperator<p2_div_cell_integral_1_otherwise> P2ConstantDivyOperator;

}
