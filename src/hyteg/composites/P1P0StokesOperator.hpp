/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/composites/P1P0StokesFunction.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/dgfunctionspace/DGStokesP1P0PressureStabForm_Example.hpp"
#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_divdiv_affine_q2.hpp"
#include "hyteg/mixedoperators/P0ScalarToP1VectorOperator.hpp"
#include "hyteg/mixedoperators/P0ToP1Operator.hpp"
#include "hyteg/mixedoperators/P1ToP0Operator.hpp"
#include "hyteg/mixedoperators/P1VectorToP0ScalarOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hyteg {

/// \brief Implements the discrete operator for the form
///
///     mu ( div(u), div(v) )
///
/// with a constant mu.
///
class P1DivDivOperator : public VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction >
{
 public:
   P1DivDivOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, real_t mu )
   : VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction >( storage, minLevel, maxLevel )
   , mu_( mu )
   {
      divdiv_0_0 = std::make_shared< forms::p1_divdiv_0_0_affine_q2 >();
      divdiv_0_1 = std::make_shared< forms::p1_divdiv_0_1_affine_q2 >();
      divdiv_0_2 = std::make_shared< forms::p1_divdiv_0_2_affine_q2 >();

      divdiv_1_0 = std::make_shared< forms::p1_divdiv_1_0_affine_q2 >();
      divdiv_1_1 = std::make_shared< forms::p1_divdiv_1_1_affine_q2 >();
      divdiv_1_2 = std::make_shared< forms::p1_divdiv_1_2_affine_q2 >();

      divdiv_2_0 = std::make_shared< forms::p1_divdiv_2_0_affine_q2 >();
      divdiv_2_1 = std::make_shared< forms::p1_divdiv_2_1_affine_q2 >();
      divdiv_2_2 = std::make_shared< forms::p1_divdiv_2_2_affine_q2 >();

      P1LinearCombinationForm mudivdiv_0_0( { mu }, { divdiv_0_0.get() } );
      P1LinearCombinationForm mudivdiv_0_1( { mu }, { divdiv_0_1.get() } );
      P1LinearCombinationForm mudivdiv_0_2( { mu }, { divdiv_0_2.get() } );

      P1LinearCombinationForm mudivdiv_1_0( { mu }, { divdiv_1_0.get() } );
      P1LinearCombinationForm mudivdiv_1_1( { mu }, { divdiv_1_1.get() } );
      P1LinearCombinationForm mudivdiv_1_2( { mu }, { divdiv_1_2.get() } );

      P1LinearCombinationForm mudivdiv_2_0( { mu }, { divdiv_2_0.get() } );
      P1LinearCombinationForm mudivdiv_2_1( { mu }, { divdiv_2_1.get() } );
      P1LinearCombinationForm mudivdiv_2_2( { mu }, { divdiv_2_2.get() } );

      auto A_0_0 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_0_0 );
      auto A_0_1 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_0_1 );
      auto A_0_2 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_0_2 );

      auto A_1_0 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_1_0 );
      auto A_1_1 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_1_1 );
      auto A_1_2 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_1_2 );

      auto A_2_0 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_2_0 );
      auto A_2_1 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_2_1 );
      auto A_2_2 =
          std::make_shared< P1ElementwiseOperator< P1LinearCombinationForm > >( storage, minLevel, maxLevel, mudivdiv_2_2 );

      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = A_0_0;
         this->subOper_[0][1] = A_0_1;
         this->subOper_[0][2] = A_0_2;

         this->subOper_[1][0] = A_1_0;
         this->subOper_[1][1] = A_1_1;
         this->subOper_[1][2] = A_1_2;

         this->subOper_[2][0] = A_2_0;
         this->subOper_[2][1] = A_2_1;
         this->subOper_[2][2] = A_2_2;
      }
      else
      {
         this->subOper_[0][0] = A_0_0;
         this->subOper_[0][1] = A_0_1;

         this->subOper_[1][0] = A_1_0;
         this->subOper_[1][1] = A_1_1;
      }
   }

 private:
   std::shared_ptr< forms::p1_divdiv_0_0_affine_q2 > divdiv_0_0;
   std::shared_ptr< forms::p1_divdiv_0_1_affine_q2 > divdiv_0_1;
   std::shared_ptr< forms::p1_divdiv_0_2_affine_q2 > divdiv_0_2;

   std::shared_ptr< forms::p1_divdiv_1_0_affine_q2 > divdiv_1_0;
   std::shared_ptr< forms::p1_divdiv_1_1_affine_q2 > divdiv_1_1;
   std::shared_ptr< forms::p1_divdiv_1_2_affine_q2 > divdiv_1_2;

   std::shared_ptr< forms::p1_divdiv_2_0_affine_q2 > divdiv_2_0;
   std::shared_ptr< forms::p1_divdiv_2_1_affine_q2 > divdiv_2_1;
   std::shared_ptr< forms::p1_divdiv_2_2_affine_q2 > divdiv_2_2;

   real_t mu_;
};

class P1P0StokesOperator : public Operator< P1P0StokesFunction< real_t >, P1P0StokesFunction< real_t > >
{
 public:
   typedef P1ConstantVectorLaplaceOperator VelocityBlockOperator_T;
   typedef P1ConstantLaplaceOperator       VelocityOperator_T;
   typedef VelocityBlockOperator_T         EnergyNormOperator_T;
   // typedef P1P0StokesBlockPreconditioner BlockPreconditioner_T;

   P1P0StokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, real_t mu )
   : Operator( storage, minLevel, maxLevel )
   , Lapl( storage, minLevel, maxLevel )
   , mudivdiv( storage, minLevel, maxLevel, mu )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , stab( storage, minLevel, maxLevel, std::make_shared< DGStokesP1P0PressureStabForm_Example >() )
   , energyNormOp( Lapl )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P1P0StokesFunction< real_t >& src,
               const P1P0StokesFunction< real_t >& dst,
               const uint_t                        level,
               const DoFType                       flag ) const
   {
      Lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      mudivdiv.apply( src.uvw(), dst.uvw(), level, flag, Add );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
      stab.apply( *src.p().getDGFunction(), *dst.p().getDGFunction(), level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1P0StokesFunction< idx_t >&          src,
                  const P1P0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      Lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      mudivdiv.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
      stab.toMatrix( mat, *src.p().getDGFunction(), *dst.p().getDGFunction(), level, flag );
   }

   VelocityBlockOperator_T    Lapl;
   P1DivDivOperator           mudivdiv;
   P1ToP0ConstantDivOperator  div;
   P0ToP1ConstantDivTOperator divT;
   dg::DGOperator             stab;
   EnergyNormOperator_T&      energyNormOp;

   bool hasGlobalCells_;
};

} // namespace hyteg
