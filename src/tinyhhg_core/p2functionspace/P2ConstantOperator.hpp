#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFOperator.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"

class p2_diffusion_cell_integral_0_otherwise;
class p2_tet_diffusion_cell_integral_0_otherwise;
class p2_mass_cell_integral_0_otherwise;
class p2_pspg_cell_integral_0_otherwise;
class p2_divt_cell_integral_0_otherwise;
class p2_divt_cell_integral_1_otherwise;
class p2_tet_divt_tet_cell_integral_0_otherwise;
class p2_tet_divt_tet_cell_integral_1_otherwise;
class p2_tet_divt_tet_cell_integral_2_otherwise;
class p2_div_cell_integral_0_otherwise;
class p2_div_cell_integral_1_otherwise;
class p2_tet_div_tet_cell_integral_0_otherwise;
class p2_tet_div_tet_cell_integral_1_otherwise;
class p2_tet_div_tet_cell_integral_2_otherwise;

namespace hhg {

using walberla::real_t;

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P2ConstantOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   const P1ConstantOperator< fenics::NoAssemble, UFCOperator3D >& getVertexToVertexOpr() const { return vertexToVertex; }

   const EdgeDoFToVertexDoFOperator< fenics::NoAssemble, UFCOperator3D >& getEdgeToVertexOpr() const { return edgeToVertex; }

   const VertexDoFToEdgeDoFOperator< fenics::NoAssemble, UFCOperator3D >& getVertexToEdgeOpr() const { return vertexToEdge; }

   const EdgeDoFOperator& getEdgeToEdgeOpr() const { return edgeToEdge; }

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const;

   void smooth_gs( const P2Function< real_t >& dst, const P2Function< real_t >& rhs, size_t level, DoFType flag ) const;

   void smooth_gs_backwards( const P2Function< real_t >& dst, const P2Function< real_t >& rhs, size_t level, DoFType flag ) const
   {
      smooth_sor_backwards( dst, rhs, 1.0, level, flag );
   }

   void smooth_sor( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const real_t &              relax,
                    size_t                      level,
                    DoFType                     flag,
                    const bool &                backwards = false ) const;

   void smooth_sor_backwards( const P2Function< real_t >& dst,
                     const P2Function< real_t >& rhs,
                     const real_t &              relax,
                     size_t                      level,
                     DoFType                     flag ) const
   {
     smooth_sor( dst, rhs, relax, level, flag, true );
   }

   void smooth_jac( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const P2Function< real_t >& src,
                    size_t                      level,
                    DoFType                     flag ) const;

 private:
   void assembleStencils();

   void assembleStencils3D();

   P1ConstantOperator< fenics::NoAssemble, UFCOperator3D >         vertexToVertex;
   EdgeDoFToVertexDoFOperator< fenics::NoAssemble, UFCOperator3D > edgeToVertex;
   VertexDoFToEdgeDoFOperator< fenics::NoAssemble, UFCOperator3D > vertexToEdge;
   EdgeDoFOperator                                                 edgeToEdge;

   void compute_local_stiffness( const Face& face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type );
};

typedef P2ConstantOperator< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise >
                                                                                                       P2ConstantLaplaceOperator;
typedef P2ConstantOperator< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > P2ConstantMassOperator;

typedef P2ConstantOperator< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > P2ConstantPSPGOperator;

typedef P2ConstantOperator< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > P2ConstantDivTxOperator;
typedef P2ConstantOperator< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > P2ConstantDivTyOperator;
typedef P2ConstantOperator< fenics::NoAssemble,                p2_tet_divt_tet_cell_integral_2_otherwise > P2ConstantDivTzOperator;
typedef P2ConstantOperator< p2_div_cell_integral_0_otherwise,  p2_tet_div_tet_cell_integral_0_otherwise >  P2ConstantDivxOperator;
typedef P2ConstantOperator< p2_div_cell_integral_1_otherwise,  p2_tet_div_tet_cell_integral_1_otherwise >  P2ConstantDivyOperator;
typedef P2ConstantOperator< fenics::NoAssemble,                p2_tet_div_tet_cell_integral_2_otherwise >  P2ConstantDivzOperator;

} // namespace hhg
