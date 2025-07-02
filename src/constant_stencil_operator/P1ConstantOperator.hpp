/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/p1functionspace/P1Operator.hpp"

namespace hyteg {

class P1LinearCombinationForm;

using walberla::real_t;

template < class P1Form, bool Diagonal = false, bool Lumped = false, bool InvertDiagonal = false, typename ValueType = real_t >
class P1ConstantOperator : public P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >
{
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::P1Operator;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::storage_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::form_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::minLevel_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::maxLevel_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::vertexStencilID_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::edgeStencilID_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::faceStencilID_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::edgeStencil3DID_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::faceStencil3DID_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::cellStencilID_;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_edge_init;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_face_init;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_cell_init;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_edge;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_edge3D;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_face;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_face3D;
   using P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::assemble_variableStencil_cell;

 public:
   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form );

   void scale( ValueType scalar );

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< idx_t >&                  src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override;

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() override
   {
      // This callback to the parent class is need for SFINAE.hpp to detect the method. PLEASE DO NOT REMOVE!
      P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::computeInverseDiagonalOperatorValues();
   }

   std::shared_ptr< P1Function< ValueType > > getInverseDiagonalValues() const override
   {
      // This callback to the parent class is need for SFINAE.hpp to detect the method. PLEASE DO NOT REMOVE!
      return P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal, ValueType >::getInverseDiagonalValues();
   }

 protected:
   /// stencil assembly: stencils are pre-assembled -> nothing to do here! ///////////

   inline void assemble_stencil_edge_init( Edge&, const uint_t ) const override {}

   inline void assemble_stencil_edge( ValueType*, const uint_t ) const override {}

   inline void assemble_stencil_face_init( Face&, const uint_t ) const override {}

   inline void assemble_stencil_face( ValueType*, const uint_t, const uint_t ) const override {}

   inline void assemble_stencil_face3D( vertexdof::macroface::StencilMap_T&, const uint_t, const uint_t ) const override {}

   inline void assemble_stencil_cell_init( Cell&, const uint_t ) const override {}

   inline void
       assemble_stencil_cell( vertexdof::macrocell::StencilMap_T&, const uint_t, const uint_t, const uint_t ) const override
   {}

   /////////////////////////////////////////////////////////////////////////////////////////////////

   void apply_face3D_generated( Face&                                                       face,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                                const uint_t&                                               level,
                                UpdateType                                                  update ) const override;

   void apply_face_generated( Face&                                                       face,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                              const uint_t&                                               level,
                              UpdateType                                                  update ) const override;

   void apply_cell_generated( Cell&                                                       cell,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& srcId,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dstId,
                              const uint_t&                                               level,
                              UpdateType                                                  update ) const override;

   void smooth_sor_face3D_generated( Face&                                                       face,
                                     const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                                     const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                                     const uint_t&                                               level,
                                     ValueType                                                   relax,
                                     const bool& backwards = false ) const override;

   void smooth_sor_face_generated( Face&                                                       face,
                                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                                   const uint_t&                                               level,
                                   ValueType                                                   relax,
                                   const bool&                                                 backwards = false ) const override;

   void smooth_sor_cell_generated( Cell&                                                       cell,
                                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& dstId,
                                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& rhsId,
                                   const uint_t&                                               level,
                                   ValueType                                                   relax,
                                   const bool&                                                 backwards = false ) const override;

   inline bool backwards_sor_available() const override { return true; }
   inline bool variableStencil() const override { return false; }

   // assemble stencils for macro-edges, -faces and -cells
   void assembleStencils();
};

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, fenics::NoAssemble > > P1ZeroOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >
    P1ConstantLaplaceOperator;
// typedef P1ConstantOperator< forms::p1_diffusion_affine_q2 > P1ConstantLaplaceOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::UndefinedAssembly >, true >
    P1DiagonalLaplaceOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_0_otherwise > > P1ConstantEpsilonOperator_11;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_1_otherwise > > P1ConstantEpsilonOperator_12;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_2_otherwise > > P1ConstantEpsilonOperator_21;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_3_otherwise > > P1ConstantEpsilonOperator_22;

typedef P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >
    P1DivxOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >
                                                                                                           P1DivyOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > > P1DivzOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >
    P1DivTxOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >
                                                                                                            P1DivTyOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > > P1DivTzOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >
    P1ConstantMassOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                            false,
                            true,
                            false >
    P1LumpedMassOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                            false,
                            true,
                            true >
    P1LumpedInvMassOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >
    P1PSPGOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise >,
                            true,
                            false,
                            true >
    P1PSPGInvDiagOperator;

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >
    P2ToP1DivxVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >
    P2ToP1DivyVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >
    P2ToP1DivzVertexToVertexConstantOperator;

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >
    P1ToP1DivTxVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >
    P1ToP1DivTyVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >
    P1ToP1DivTzVertexToVertexConstantOperator;

typedef P1ConstantOperator< P1LinearCombinationForm > P1ConstantLinearCombinationOperator;

} // namespace hyteg
