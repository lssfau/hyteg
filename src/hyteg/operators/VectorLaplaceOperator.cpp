/*
 * Copyright (c) 2017-2021 Marcus Mohr.
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
#include "hyteg/operators/VectorLaplaceOperator.hpp"

namespace hyteg {

using walberla::real_t;

template < class VecFuncType, class SubOpType >
VectorLaplaceOperator< VecFuncType, SubOpType >::VectorLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                        size_t                                     minLevel,
                                                                        size_t                                     maxLevel )
: VectorToVectorOperator< VecFuncType, VecFuncType >( storage, minLevel, maxLevel )
{
   std::shared_ptr< SubOpType > zero( nullptr );
   std::shared_ptr< SubOpType > lapl = std::make_shared< SubOpType >( storage, minLevel, maxLevel );

   if ( this->dim_ == 3 )
   {
      this->subOper_[0][0] = lapl;
      this->subOper_[0][1] = zero;
      this->subOper_[0][2] = zero;

      this->subOper_[1][0] = zero;
      this->subOper_[1][1] = lapl;
      this->subOper_[1][2] = zero;

      this->subOper_[2][0] = zero;
      this->subOper_[2][1] = zero;
      this->subOper_[2][2] = lapl;
   }
   else
   {
      this->subOper_[0][0] = lapl;
      this->subOper_[0][1] = zero;

      this->subOper_[1][0] = zero;
      this->subOper_[1][1] = lapl;
   }
};

// P1ConstantVectorLaplaceOperator
template class VectorLaplaceOperator< P1VectorFunction< real_t >, P1ConstantLaplaceOperator >;

// P1ElementwiseVectorLaplaceOperator
template class VectorLaplaceOperator< P1VectorFunction< real_t >, P1ElementwiseLaplaceOperator >;

// P1ElementwiseBlendingVectorLaplaceOperator
template class VectorLaplaceOperator< P1VectorFunction< real_t >, P1ElementwiseBlendingLaplaceOperator >;

// P2ConstantVectorLaplaceOperator
template class VectorLaplaceOperator< P2VectorFunction< real_t >, P2ConstantLaplaceOperator >;

// P2ElementwiseVectorLaplaceOperator
template class VectorLaplaceOperator< P2VectorFunction< real_t >, P2ElementwiseLaplaceOperator >;

// P2ElementwiseBlendingVectorLaplaceOperator
template class VectorLaplaceOperator< P2VectorFunction< real_t >, P2ElementwiseBlendingLaplaceOperator >;

} // namespace hyteg
