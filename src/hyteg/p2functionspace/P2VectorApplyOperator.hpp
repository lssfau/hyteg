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

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilon_all_forms.hpp"
#include "hyteg/operators/VectorToVectorOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

#include "hyteg/p2functionspace/P2VectorApplyOperator.hpp"
namespace hyteg {

using walberla::real_t;


class P2VectorApplyOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >
{
 public:
   P2VectorApplyOperator(  const std::shared_ptr< PrimitiveStorage >& storage,
                           size_t                                     level )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, level, level )
   { 
   }

   void apply( const P2VectorFunction<real_t>& src,
               const P2VectorFunction<real_t>& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const
   {
      dst.assign({1}, {src}, level, flag);
      dst.multElementwise({*diagonal_},level, flag);
   }

   std::shared_ptr<P2VectorFunction<real_t>> diagonal_;
};
}