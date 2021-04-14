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

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/operators/VectorToVectorOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template < class VecFuncType, class SubOpType >
class VectorLaplaceOperator : public VectorToVectorOperator< VecFuncType, VecFuncType >
{
 public:
   VectorLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );
};

// P1 versions
typedef VectorLaplaceOperator< P1VectorFunction< real_t >, P1ConstantLaplaceOperator >    P1ConstantVectorLaplaceOperator;
typedef VectorLaplaceOperator< P1VectorFunction< real_t >, P1ElementwiseLaplaceOperator > P1ElementwiseVectorLaplaceOperator;
typedef VectorLaplaceOperator< P1VectorFunction< real_t >, P1ElementwiseBlendingLaplaceOperator >
    P1ElementwiseBlendingVectorLaplaceOperator;

// P2 versions
typedef VectorLaplaceOperator< P2VectorFunction< real_t >, P2ConstantLaplaceOperator >    P2ConstantVectorLaplaceOperator;
typedef VectorLaplaceOperator< P2VectorFunction< real_t >, P2ElementwiseLaplaceOperator > P2ElementwiseVectorLaplaceOperator;
typedef VectorLaplaceOperator< P2VectorFunction< real_t >, P2ElementwiseBlendingLaplaceOperator >
    P2ElementwiseBlendingVectorLaplaceOperator;

} // namespace hyteg
