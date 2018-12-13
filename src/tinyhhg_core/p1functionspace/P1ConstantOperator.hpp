#pragma once

#include <array>

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/StencilMemory.hpp"


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

class p2_to_p1_tet_div_tet_cell_integral_0_otherwise;
class p2_to_p1_tet_div_tet_cell_integral_1_otherwise;
class p2_to_p1_tet_div_tet_cell_integral_2_otherwise;

class p1_to_p2_tet_divt_tet_cell_integral_0_otherwise;
class p1_to_p2_tet_divt_tet_cell_integral_1_otherwise;
class p1_to_p2_tet_divt_tet_cell_integral_2_otherwise;


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

   ~P1ConstantOperator() override = default;

   void scale( real_t scalar );

   void apply(const P1Function< real_t >& src,
              const P1Function< real_t >& dst,
              size_t                level,
              DoFType               flag,
              UpdateType            updateType = Replace ) const;

   void smooth_gs( const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag ) const;

   void smooth_sor( const P1Function< real_t >& dst,
                         const P1Function< real_t >& rhs,
                         real_t                      relax,
                         size_t                      level,
                         DoFType                     flag ) const;

   void smooth_jac( const P1Function< real_t >& dst,
                         const P1Function< real_t >& rhs,
                         const P1Function< real_t >& tmp,
                         size_t                      level,
                         DoFType                     flag ) const;

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >& getVertexStencilID() const { return vertexStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Edge >& getEdgeStencilID() const { return edgeStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Face >& getFaceStencilID() const { return faceStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Cell >& getCellStencilID() const { return cellStencilID_; }

 private:
   void assembleStencils();

   void assembleStencils3D();

 private:





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

typedef P1ConstantOperator< fenics::NoAssemble,        p2_to_p1_tet_div_tet_cell_integral_0_otherwise > P2ToP1DivxVertexToVertexConstantOperator;
typedef P1ConstantOperator< fenics::NoAssemble,        p2_to_p1_tet_div_tet_cell_integral_1_otherwise > P2ToP1DivyVertexToVertexConstantOperator;
typedef P1ConstantOperator< fenics::UndefinedAssembly, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > P2ToP1DivzVertexToVertexConstantOperator;

typedef P1ConstantOperator< fenics::NoAssemble,        p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > P1ToP1DivTxVertexToVertexConstantOperator;
typedef P1ConstantOperator< fenics::NoAssemble,        p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > P1ToP1DivTyVertexToVertexConstantOperator;
typedef P1ConstantOperator< fenics::UndefinedAssembly, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > P1ToP1DivTzVertexToVertexConstantOperator;

} // namespace hhg
