/*
 * Copyright (c) 2017-2025 Marcus Mohr.
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

#include "hyteg/ccrfunctionspace/P2PlusBubbleVectorFunction.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg_operators/operators/mass/P2PlusBubbleElementwiseMass.hpp"

#include "VectorToVectorOperator.hpp"
#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "constant_stencil_operator/P2ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template < typename ValueType, template < typename > class VecFuncKind, class SubOpType >
class VectorMassOperator : public VectorToVectorOperator< ValueType, VecFuncKind, VecFuncKind >
{
 public:
   VectorMassOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );
};

// P1 versions
typedef VectorMassOperator< real_t, P1VectorFunction, P1ConstantMassOperator >            P1ConstantVectorMassOperator;
typedef VectorMassOperator< real_t, P1VectorFunction, P1ElementwiseMassOperator >         P1ElementwiseVectorMassOperator;
typedef VectorMassOperator< real_t, P1VectorFunction, P1ElementwiseBlendingMassOperator > P1ElementwiseBlendingVectorMassOperator;

// P2 versions
typedef VectorMassOperator< real_t, P2VectorFunction, P2ConstantMassOperator >            P2ConstantVectorMassOperator;
typedef VectorMassOperator< real_t, P2VectorFunction, P2ElementwiseMassOperator >         P2ElementwiseVectorMassOperator;
typedef VectorMassOperator< real_t, P2VectorFunction, P2ElementwiseBlendingMassOperator > P2ElementwiseBlendingVectorMassOperator;

// P2PlusBubble versions
typedef VectorMassOperator< real_t, P2PlusBubbleVectorFunction, operatorgeneration::P2PlusBubbleElementwiseMass >
    P2PlusBubbleElementwiseVectorMassOperator;

} // namespace hyteg
