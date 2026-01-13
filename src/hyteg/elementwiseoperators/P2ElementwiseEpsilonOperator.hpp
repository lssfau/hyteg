/*
 * Copyright (c) 2021-2025 Marcus Mohr.
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

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilon_all_forms.hpp"

#include "mixed_operator/VectorToVectorOperator.hpp"

namespace hyteg {

using walberla::real_t;

class P2ElementwiseAffineEpsilonOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
                                           public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   P2ElementwiseAffineEpsilonOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                       size_t                                     minLevel,
                                       size_t                                     maxLevel,
                                       std::function< real_t( const Point3D& ) >  viscosity )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   {
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_0_0_affine_q2 > eps_0_0;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_0_1_affine_q2 > eps_0_1;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_0_2_affine_q2 > eps_0_2;

      typedef P2ElementwiseOperator< forms::p2_epsilonvar_1_0_affine_q2 > eps_1_0;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_1_1_affine_q2 > eps_1_1;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_1_2_affine_q2 > eps_1_2;

      typedef P2ElementwiseOperator< forms::p2_epsilonvar_2_0_affine_q2 > eps_2_0;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_2_1_affine_q2 > eps_2_1;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_2_2_affine_q2 > eps_2_2;

      if ( this->dim_ == 3 )
      {
         auto form_0_0 = forms::p2_epsilonvar_0_0_affine_q2( viscosity, viscosity );
         auto form_0_1 = forms::p2_epsilonvar_0_1_affine_q2( viscosity, viscosity );
         auto form_0_2 = forms::p2_epsilonvar_0_2_affine_q2( viscosity );

         auto form_1_0 = forms::p2_epsilonvar_1_0_affine_q2( viscosity, viscosity );
         auto form_1_1 = forms::p2_epsilonvar_1_1_affine_q2( viscosity, viscosity );
         auto form_1_2 = forms::p2_epsilonvar_1_2_affine_q2( viscosity );

         auto form_2_0 = forms::p2_epsilonvar_2_0_affine_q2( viscosity );
         auto form_2_1 = forms::p2_epsilonvar_2_1_affine_q2( viscosity );
         auto form_2_2 = forms::p2_epsilonvar_2_2_affine_q2( viscosity );

         this->subOper_[0][0] = std::make_shared< eps_0_0 >( storage, minLevel, maxLevel, form_0_0 );
         this->subOper_[0][1] = std::make_shared< eps_0_1 >( storage, minLevel, maxLevel, form_0_1 );
         this->subOper_[0][2] = std::make_shared< eps_0_2 >( storage, minLevel, maxLevel, form_0_2 );

         this->subOper_[1][0] = std::make_shared< eps_1_0 >( storage, minLevel, maxLevel, form_1_0 );
         this->subOper_[1][1] = std::make_shared< eps_1_1 >( storage, minLevel, maxLevel, form_1_1 );
         this->subOper_[1][2] = std::make_shared< eps_1_2 >( storage, minLevel, maxLevel, form_1_2 );

         this->subOper_[2][0] = std::make_shared< eps_2_0 >( storage, minLevel, maxLevel, form_2_0 );
         this->subOper_[2][1] = std::make_shared< eps_2_1 >( storage, minLevel, maxLevel, form_2_1 );
         this->subOper_[2][2] = std::make_shared< eps_2_2 >( storage, minLevel, maxLevel, form_2_2 );
      }
      else
      {
         auto form_0_0 = forms::p2_epsilonvar_0_0_affine_q2( viscosity, viscosity );
         auto form_0_1 = forms::p2_epsilonvar_0_1_affine_q2( viscosity, viscosity );

         auto form_1_0 = forms::p2_epsilonvar_1_0_affine_q2( viscosity, viscosity );
         auto form_1_1 = forms::p2_epsilonvar_1_1_affine_q2( viscosity, viscosity );

         this->subOper_[0][0] = std::make_shared< eps_0_0 >( storage, minLevel, maxLevel, form_0_0 );
         this->subOper_[0][1] = std::make_shared< eps_0_1 >( storage, minLevel, maxLevel, form_0_1 );

         this->subOper_[1][0] = std::make_shared< eps_1_0 >( storage, minLevel, maxLevel, form_1_0 );
         this->subOper_[1][1] = std::make_shared< eps_1_1 >( storage, minLevel, maxLevel, form_1_1 );
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

/// Uses forms::p2_epsilonvar_*_*_blending_q2, if that is insufficient upgrade to q3
class P2ElementwiseBlendingEpsilonOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
                                             public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   P2ElementwiseBlendingEpsilonOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                         size_t                                     minLevel,
                                         size_t                                     maxLevel,
                                         std::function< real_t( const Point3D& ) >  viscosity )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   {
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_0_0_blending_q2 > eps_0_0;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_0_1_blending_q2 > eps_0_1;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_0_2_blending_q2 > eps_0_2;

      typedef P2ElementwiseOperator< forms::p2_epsilonvar_1_0_blending_q2 > eps_1_0;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_1_1_blending_q2 > eps_1_1;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_1_2_blending_q2 > eps_1_2;

      typedef P2ElementwiseOperator< forms::p2_epsilonvar_2_0_blending_q2 > eps_2_0;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_2_1_blending_q2 > eps_2_1;
      typedef P2ElementwiseOperator< forms::p2_epsilonvar_2_2_blending_q2 > eps_2_2;

      if ( this->dim_ == 3 )
      {
         auto form_0_0 = forms::p2_epsilonvar_0_0_blending_q2( viscosity, viscosity );
         auto form_0_1 = forms::p2_epsilonvar_0_1_blending_q2( viscosity, viscosity );
         auto form_0_2 = forms::p2_epsilonvar_0_2_blending_q2( viscosity );

         auto form_1_0 = forms::p2_epsilonvar_1_0_blending_q2( viscosity, viscosity );
         auto form_1_1 = forms::p2_epsilonvar_1_1_blending_q2( viscosity, viscosity );
         auto form_1_2 = forms::p2_epsilonvar_1_2_blending_q2( viscosity );

         auto form_2_0 = forms::p2_epsilonvar_2_0_blending_q2( viscosity );
         auto form_2_1 = forms::p2_epsilonvar_2_1_blending_q2( viscosity );
         auto form_2_2 = forms::p2_epsilonvar_2_2_blending_q2( viscosity );

         this->subOper_[0][0] = std::make_shared< eps_0_0 >( storage, minLevel, maxLevel, form_0_0 );
         this->subOper_[0][1] = std::make_shared< eps_0_1 >( storage, minLevel, maxLevel, form_0_1 );
         this->subOper_[0][2] = std::make_shared< eps_0_2 >( storage, minLevel, maxLevel, form_0_2 );

         this->subOper_[1][0] = std::make_shared< eps_1_0 >( storage, minLevel, maxLevel, form_1_0 );
         this->subOper_[1][1] = std::make_shared< eps_1_1 >( storage, minLevel, maxLevel, form_1_1 );
         this->subOper_[1][2] = std::make_shared< eps_1_2 >( storage, minLevel, maxLevel, form_1_2 );

         this->subOper_[2][0] = std::make_shared< eps_2_0 >( storage, minLevel, maxLevel, form_2_0 );
         this->subOper_[2][1] = std::make_shared< eps_2_1 >( storage, minLevel, maxLevel, form_2_1 );
         this->subOper_[2][2] = std::make_shared< eps_2_2 >( storage, minLevel, maxLevel, form_2_2 );
      }
      else
      {
         auto form_0_0 = forms::p2_epsilonvar_0_0_blending_q2( viscosity, viscosity );
         auto form_0_1 = forms::p2_epsilonvar_0_1_blending_q2( viscosity, viscosity );

         auto form_1_0 = forms::p2_epsilonvar_1_0_blending_q2( viscosity, viscosity );
         auto form_1_1 = forms::p2_epsilonvar_1_1_blending_q2( viscosity, viscosity );

         this->subOper_[0][0] = std::make_shared< eps_0_0 >( storage, minLevel, maxLevel, form_0_0 );
         this->subOper_[0][1] = std::make_shared< eps_0_1 >( storage, minLevel, maxLevel, form_0_1 );

         this->subOper_[1][0] = std::make_shared< eps_1_0 >( storage, minLevel, maxLevel, form_1_0 );
         this->subOper_[1][1] = std::make_shared< eps_1_1 >( storage, minLevel, maxLevel, form_1_1 );
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
