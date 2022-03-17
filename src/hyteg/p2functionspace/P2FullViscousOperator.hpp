/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include "hyteg/forms/form_fenics_generated/p2_stokes_full.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_stokes_full_tet.h"
#include "hyteg/forms/form_hyteg_generated/p2/p2_full_stokes_all_forms.hpp"
#include "hyteg/operators/VectorToVectorOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_t;

class P2ConstantFullViscousOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >
{
 public:
   P2ConstantFullViscousOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   {
      // clang-format off
      typedef P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise > > eps_0_0;
      typedef P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise > > eps_0_1;
      typedef P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise > > eps_0_2;

      typedef P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise > > eps_1_0;
      typedef P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise > > eps_1_1;
      typedef P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise > > eps_1_2;

      typedef P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise > > eps_2_0;
      typedef P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise > > eps_2_1;
      typedef P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise > > eps_2_2;
      // clang-format on

      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = std::make_shared< eps_0_0 >( storage, minLevel, maxLevel );
         this->subOper_[0][1] = std::make_shared< eps_0_1 >( storage, minLevel, maxLevel );
         this->subOper_[0][2] = std::make_shared< eps_0_2 >( storage, minLevel, maxLevel );

         this->subOper_[1][0] = std::make_shared< eps_1_0 >( storage, minLevel, maxLevel );
         this->subOper_[1][1] = std::make_shared< eps_1_1 >( storage, minLevel, maxLevel );
         this->subOper_[1][2] = std::make_shared< eps_1_2 >( storage, minLevel, maxLevel );

         this->subOper_[2][0] = std::make_shared< eps_2_0 >( storage, minLevel, maxLevel );
         this->subOper_[2][1] = std::make_shared< eps_2_1 >( storage, minLevel, maxLevel );
         this->subOper_[2][2] = std::make_shared< eps_2_2 >( storage, minLevel, maxLevel );
      }
      else
      {
         this->subOper_[0][0] = std::make_shared< eps_0_0 >( storage, minLevel, maxLevel );
         this->subOper_[0][1] = std::make_shared< eps_0_1 >( storage, minLevel, maxLevel );

         this->subOper_[1][0] = std::make_shared< eps_1_0 >( storage, minLevel, maxLevel );
         this->subOper_[1][1] = std::make_shared< eps_1_1 >( storage, minLevel, maxLevel );
      }
   }
};

class P2ElementwiseBlendingFullViscousOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
                                                 public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   P2ElementwiseBlendingFullViscousOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                             size_t                                     minLevel,
                                             size_t                                     maxLevel,
                                             std::function< real_t( const Point3D& ) >  viscosity )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   {
      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_0_0_blending_q3 > visc_0_0;
      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_0_1_blending_q3 > visc_0_1;
      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_0_2_blending_q3 > visc_0_2;

      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_1_0_blending_q3 > visc_1_0;
      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_1_1_blending_q3 > visc_1_1;
      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_1_2_blending_q3 > visc_1_2;

      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_2_0_blending_q3 > visc_2_0;
      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_2_1_blending_q3 > visc_2_1;
      typedef P2ElementwiseOperator< forms::p2_full_stokesvar_2_2_blending_q3 > visc_2_2;

      if ( this->dim_ == 3 )
      {
         auto form_0_0 = forms::p2_full_stokesvar_0_0_blending_q3( viscosity, viscosity );
         auto form_0_1 = forms::p2_full_stokesvar_0_1_blending_q3( viscosity, viscosity );
         auto form_0_2 = forms::p2_full_stokesvar_0_2_blending_q3( viscosity );

         auto form_1_0 = forms::p2_full_stokesvar_1_0_blending_q3( viscosity, viscosity );
         auto form_1_1 = forms::p2_full_stokesvar_1_1_blending_q3( viscosity, viscosity );
         auto form_1_2 = forms::p2_full_stokesvar_1_2_blending_q3( viscosity );

         auto form_2_0 = forms::p2_full_stokesvar_2_0_blending_q3( viscosity );
         auto form_2_1 = forms::p2_full_stokesvar_2_1_blending_q3( viscosity );
         auto form_2_2 = forms::p2_full_stokesvar_2_2_blending_q3( viscosity );

         this->subOper_[0][0] = std::make_shared< visc_0_0 >( storage, minLevel, maxLevel, form_0_0 );
         this->subOper_[0][1] = std::make_shared< visc_0_1 >( storage, minLevel, maxLevel, form_0_1 );
         this->subOper_[0][2] = std::make_shared< visc_0_2 >( storage, minLevel, maxLevel, form_0_2 );

         this->subOper_[1][0] = std::make_shared< visc_1_0 >( storage, minLevel, maxLevel, form_1_0 );
         this->subOper_[1][1] = std::make_shared< visc_1_1 >( storage, minLevel, maxLevel, form_1_1 );
         this->subOper_[1][2] = std::make_shared< visc_1_2 >( storage, minLevel, maxLevel, form_1_2 );

         this->subOper_[2][0] = std::make_shared< visc_2_0 >( storage, minLevel, maxLevel, form_2_0 );
         this->subOper_[2][1] = std::make_shared< visc_2_1 >( storage, minLevel, maxLevel, form_2_1 );
         this->subOper_[2][2] = std::make_shared< visc_2_2 >( storage, minLevel, maxLevel, form_2_2 );
      }
      else
      {
         auto form_0_0 = forms::p2_full_stokesvar_0_0_blending_q3( viscosity, viscosity );
         auto form_0_1 = forms::p2_full_stokesvar_0_1_blending_q3( viscosity, viscosity );

         auto form_1_0 = forms::p2_full_stokesvar_1_0_blending_q3( viscosity, viscosity );
         auto form_1_1 = forms::p2_full_stokesvar_1_1_blending_q3( viscosity, viscosity );

         this->subOper_[0][0] = std::make_shared< visc_0_0 >( storage, minLevel, maxLevel, form_0_0 );
         this->subOper_[0][1] = std::make_shared< visc_0_1 >( storage, minLevel, maxLevel, form_0_1 );

         this->subOper_[1][0] = std::make_shared< visc_1_0 >( storage, minLevel, maxLevel, form_1_0 );
         this->subOper_[1][1] = std::make_shared< visc_1_1 >( storage, minLevel, maxLevel, form_1_1 );
      }
   }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() override final
   {
      this->VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::computeInverseDiagonalOperatorValues();
   }
};

} // namespace hyteg
