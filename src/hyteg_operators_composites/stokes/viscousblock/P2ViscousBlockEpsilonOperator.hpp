/*
* Copyright (c) 2017-2024 Nils Kohl.
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

#include "core/DataTypes.h"

#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_0_0.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_0_1.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_0_2.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_1_0.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_1_1.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_1_2.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_2_0.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_2_1.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilonIcosahedralShellMap_2_2.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_0_0.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_0_1.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_0_2.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_1_0.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_1_1.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_1_2.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_2_0.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_2_1.hpp"
#include "hyteg_operators/operators/epsilon/P2ElementwiseEpsilon_2_2.hpp"

#include "mixed_operator/VectorToVectorOperator.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Implements the block-"Epsilon" operator, which is the viscous block (or A-block) of the variable-viscosity Stokes operator
/// and corresponds to the operator (strong form)
///
///     - ∇ · ( 2μ ε(u) )
///
/// and the weak form (vector valued u and v)
///
///     ∫ 2 μ ε(u) : ε(v)
///
/// where
///
///     ε(w) := (1/2) (∇w + (∇w)ᵀ)
///
/// and μ = μ(x) is a user-supplied finite element function.
///
/// All blocks are non-zero. That is, we have the structure
///
///         /                \
///         | A_11 A_12 A_13 |
///     A = | A_21 A_22 A_23 |
///         | A_31 A_32 A_33 |
///         \                /
///
class P2ViscousBlockEpsilonOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >
{
 public:
   P2ViscousBlockEpsilonOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                  uint_t                                     minLevel,
                                  uint_t                                     maxLevel,
                                  const P2Function< real_t >&                mu )
   : VectorToVectorOperator< real_t, hyteg::P2VectorFunction, hyteg::P2VectorFunction >( storage, minLevel, maxLevel )
   {
      this->setSubOperator(
          0, 0, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_0_0 >( storage, minLevel, maxLevel, mu ) );
      this->setSubOperator(
          0, 1, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_0_1 >( storage, minLevel, maxLevel, mu ) );

      this->setSubOperator(
          1, 0, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_1_0 >( storage, minLevel, maxLevel, mu ) );
      this->setSubOperator(
          1, 1, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_1_1 >( storage, minLevel, maxLevel, mu ) );

      if ( storage->hasGlobalCells() )
      {
         this->setSubOperator(
             0, 2, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_0_2 >( storage, minLevel, maxLevel, mu ) );

         this->setSubOperator(
             1, 2, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_1_2 >( storage, minLevel, maxLevel, mu ) );

         this->setSubOperator(
             2, 0, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_2_0 >( storage, minLevel, maxLevel, mu ) );
         this->setSubOperator(
             2, 1, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_2_1 >( storage, minLevel, maxLevel, mu ) );
         this->setSubOperator(
             2, 2, std::make_shared< operatorgeneration::P2ElementwiseEpsilon_2_2 >( storage, minLevel, maxLevel, mu ) );
      }
   }
};

} // namespace operatorgeneration
} // namespace hyteg