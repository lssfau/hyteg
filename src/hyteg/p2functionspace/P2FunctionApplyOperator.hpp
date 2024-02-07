/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.cpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilon_all_forms.hpp"
#include "hyteg/operators/VectorToVectorOperator.hpp"

namespace hyteg {

using walberla::real_t;

class P2FunctionApplyOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >
                             //   public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   P2FunctionApplyOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t level )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, level, level )
   {
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P2VectorFunction< real_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType = Replace ) const
   {
      //dst.assign({1}, {src}, level, flag);
      dst.multElementwise( { *diagonal_, src }, level, flag );
   }

   
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                                 mat,
                  const typename P2VectorFunction<real_t>::template FunctionType< idx_t >& src,
                  const typename P2VectorFunction<real_t>::template FunctionType< idx_t >& dst,
                  size_t                                                                      level,
                  DoFType                                                                     flag ) const override
   {
      WALBERLA_UNUSED( dst );
      //std::shared_ptr< funcType > opVals = InvertDiagonal ? oper_->getInverseDiagonalValues() : oper_->getDiagonalValues();
      diagonal_->interpolate(1,level,hyteg::DirichletBoundary);
      for ( uint_t i = 0; i < dim_; i++ )
      {
         workaround::externalDiagonalAssembly< P2Function<real_t> >( mat, diagonal_->component(i), src[i], level, flag );
      }
   }

   std::shared_ptr< P2VectorFunction< real_t > > diagonal_;
};
} // namespace hyteg