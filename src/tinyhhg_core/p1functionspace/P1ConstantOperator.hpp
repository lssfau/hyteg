#pragma once

#include <array>

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/HHGDefinitions.hpp"
#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "P1DataHandling.hpp"

#include "P1Function.hpp"
#include "generatedKernels/generatedKernels.hpp"

class p1_diffusion_cell_integral_0_otherwise;
class p1_tet_diffusion_cell_integral_0_otherwise;
class p1_stokes_epsilon_cell_integral_0_otherwise;
class p1_stokes_epsilon_cell_integral_1_otherwise;
class p1_stokes_epsilon_cell_integral_2_otherwise;
class p1_stokes_epsilon_cell_integral_3_otherwise;
class p1_div_cell_integral_0_otherwise;
class p1_tet_div_tet_cell_integral_0_otherwise;
class p1_div_cell_integral_1_otherwise;
class p1_tet_div_tet_cell_integral_1_otherwise;
class p1_tet_div_tet_cell_integral_2_otherwise;
class p1_divt_cell_integral_0_otherwise;
class p1_tet_divt_tet_cell_integral_0_otherwise;
class p1_divt_cell_integral_1_otherwise;
class p1_tet_divt_tet_cell_integral_1_otherwise;
class p1_tet_divt_tet_cell_integral_2_otherwise;
class p1_mass_cell_integral_0_otherwise;
class p1_tet_mass_cell_integral_0_otherwise;
class p1_pspg_cell_integral_0_otherwise;
class p1_tet_pspg_tet_cell_integral_0_otherwise;

namespace hhg {

using walberla::real_t;

template < class UFCOperator2D,
           class UFCOperator3D = fenics::UndefinedAssembly,
           bool Diagonal       = false,
           bool Lumped         = false,
           bool InvertDiagonal = false >
class P1ConstantOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   ~P1ConstantOperator() override {}

   void scale( real_t scalar );

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >& getVertexStencilID() const { return vertexStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Edge >& getEdgeStencilID() const { return edgeStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Face >& getFaceStencilID() const { return faceStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Cell >& getCellStencilID() const { return cellStencilID_; }

 private:
   void assembleStencils();

   void assembleStencils3D();

 private:
   void apply_impl( P1Function< real_t >& src,
                    P1Function< real_t >& dst,
                    size_t                level,
                    DoFType               flag,
                    UpdateType            updateType = Replace ) override;

   void smooth_gs_impl( P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag ) override;

   void smooth_sor_impl( P1Function< real_t >& dst, P1Function< real_t >& rhs, real_t relax, size_t level, DoFType flag ) override;

   void smooth_jac_impl( P1Function< real_t >& dst,
                         P1Function< real_t >& rhs,
                         P1Function< real_t >& tmp,
                         size_t                level,
                         DoFType               flag ) override;

   PrimitiveDataID< StencilMemory< real_t >, Vertex > vertexStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >   edgeStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >   faceStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Cell >   cellStencilID_;

   void compute_local_stiffness( const Face& face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type );
};

typedef P1ConstantOperator< fenics::NoAssemble, fenics::NoAssemble > P1ZeroOperator;

typedef P1ConstantOperator< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >
                                                                                                      P1ConstantLaplaceOperator;
typedef P1ConstantOperator< p1_diffusion_cell_integral_0_otherwise, fenics::UndefinedAssembly, true > P1DiagonalLaplaceOperator;

typedef P1ConstantOperator< p1_stokes_epsilon_cell_integral_0_otherwise > P1ConstantEpsilonOperator_11;
typedef P1ConstantOperator< p1_stokes_epsilon_cell_integral_1_otherwise > P1ConstantEpsilonOperator_12;
typedef P1ConstantOperator< p1_stokes_epsilon_cell_integral_2_otherwise > P1ConstantEpsilonOperator_21;
typedef P1ConstantOperator< p1_stokes_epsilon_cell_integral_3_otherwise > P1ConstantEpsilonOperator_22;

typedef P1ConstantOperator< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > P1DivxOperator;
typedef P1ConstantOperator< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > P1DivyOperator;
typedef P1ConstantOperator< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise >               P1DivzOperator;

typedef P1ConstantOperator< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > P1DivTxOperator;
typedef P1ConstantOperator< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > P1DivTyOperator;
typedef P1ConstantOperator< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise >                P1DivTzOperator;

typedef P1ConstantOperator< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > P1MassOperator;
typedef P1ConstantOperator< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise, false, true, true >
    P1LumpedInvMassOperator;

typedef P1ConstantOperator< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > P1PSPGOperator;

} // namespace hhg
