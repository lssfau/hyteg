/*
 * Copyright (c) 2021 Benjamin Mann
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

#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_affine_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_divt_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/p1functionspace/P1Operator.hpp"

namespace hyteg {

template < class P1Form >
class P1VariableOperator : public P1Operator< P1Form >
{
   using P1Operator< P1Form >::P1Operator;
   using P1Operator< P1Form >::storage_;
   using P1Operator< P1Form >::assemble_variableStencil_edge_init;
   using P1Operator< P1Form >::assemble_variableStencil_face_init;
   using P1Operator< P1Form >::assemble_variableStencil_cell_init;
   using P1Operator< P1Form >::assemble_variableStencil_edge;
   using P1Operator< P1Form >::assemble_variableStencil_edge3D;
   using P1Operator< P1Form >::assemble_variableStencil_face;
   using P1Operator< P1Form >::assemble_variableStencil_face3D;
   using P1Operator< P1Form >::assemble_variableStencil_cell;

 public:
   P1VariableOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : P1VariableOperator( storage, minLevel, maxLevel, P1Form() )
   {}

   P1VariableOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form )
   : P1Operator< P1Form >( storage, minLevel, maxLevel, form )
   {}

 protected:
   /// stencil assembly ///////////

   /* Initialize assembly of variable edge stencil.
      Will be called before iterating over edge whenever the stencil is applied.
   */
   inline void assemble_stencil_edge_init( Edge& edge, const uint_t level ) const
   {
      assemble_variableStencil_edge_init( edge, level );
   }

   /* Assembly of edge stencil.
      Will be called before stencil is applied to a particuar edge-DoF.
   */
   inline void assemble_stencil_edge( real_t* edge_stencil, const uint_t i ) const
   {
      assemble_variableStencil_edge( edge_stencil, i );
   }

   /* Initialize assembly of face stencil.
      Will be called before iterating over face whenever the stencil is applied.
   */
   inline void assemble_stencil_face_init( Face& face, const uint_t level ) const
   {
      assemble_variableStencil_face_init( face, level );
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 2d domain.
   */
   inline void assemble_stencil_face( real_t* face_stencil, const uint_t i, const uint_t j ) const
   {
      assemble_variableStencil_face( face_stencil, i, j );
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 3D domain.
   */
   inline void assemble_stencil_face3D( vertexdof::macroface::StencilMap_T& face_stencil, const uint_t i, const uint_t j ) const
   {
      assemble_variableStencil_face3D( face_stencil, i, j );
   }

   /* Initialize assembly of cell stencil.
      Will be called before iterating over cell whenever the stencil is applied.
   */
   inline void assemble_stencil_cell_init( Cell& cell, const uint_t level ) const
   {
      assemble_variableStencil_cell_init( cell, level );
   }

   /* Assembly of cell stencil.
      Will be called before stencil is applied to a particuar cell-DoF.
   */
   inline void assemble_stencil_cell( vertexdof::macrocell::StencilMap_T& cell_stencil,
                                      const uint_t                        i,
                                      const uint_t                        j,
                                      const uint_t                        k ) const
   {
      assemble_variableStencil_cell( cell_stencil, i, j, k );
   }

   inline bool backwards_sor_available() const { return false; }
   inline bool variableStencil() const { return true; }
};

typedef P1VariableOperator< forms::p1_diffusion_blending_q1 > P1BlendingLaplaceOperator;
typedef P1VariableOperator< forms::p1_mass_blending_q4 >      P1BlendingMassOperator;

typedef P1VariableOperator< forms::p1_epsiloncc_0_0_blending_q2 > P1BlendingEpsilonOperator_11;
typedef P1VariableOperator< forms::p1_epsiloncc_0_1_blending_q2 > P1BlendingEpsilonOperator_12;
typedef P1VariableOperator< forms::p1_epsiloncc_1_0_blending_q2 > P1BlendingEpsilonOperator_21;
typedef P1VariableOperator< forms::p1_epsiloncc_1_1_blending_q2 > P1BlendingEpsilonOperator_22;

typedef P1VariableOperator< forms::p1_divt_0_blending_q1 > P1BlendingDivTOperator_1;
typedef P1VariableOperator< forms::p1_divt_1_blending_q1 > P1BlendingDivTOperator_2;

typedef P1VariableOperator< forms::p1_div_0_blending_q1 > P1BlendingDivOperator_1;
typedef P1VariableOperator< forms::p1_div_1_blending_q1 > P1BlendingDivOperator_2;

// typedef P1VariableOperator</*todo*/ > P1BlendingPSPGOperator; // todo

typedef P1VariableOperator< forms::p1_div_k_grad_blending_q3 > P1BlendingDivkGradOperator;
typedef P1VariableOperator< forms::p1_div_k_grad_affine_q3 >   P1AffineDivkGradOperator;

typedef P1VariableOperator< P1LinearCombinationForm > P1VariableLinearCombinationOperator;

} // namespace hyteg
