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
#include "hyteg_operators/operators/diffusion/P2ElementwiseDiffusion.hpp"
#include "hyteg_operators/operators/diffusion/P2ElementwiseDiffusionIcosahedralShellMap.hpp"

#include "mixed_operator/VectorToVectorOperator.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Implements the block-Laplace operator, which is the viscous block (or A-block) of the vanilla Stokes operator.
/// Specifically, it corresponds to the operator (strong form)
///
///     - Δu
///
/// and the weak form (vector valued u and v)
///
///     ∫ ∇u : ∇v.
///
/// Only the diagonal blocks are non-zero and each of them corresponds to the scalar Laplace operator. That is, we have the
/// structure
///
///         /                \
///         | A_11  0    0   |
///     A = |  0   A_22  0   |
///         |  0    0   A_33 |
///         \                /
///
class P2ViscousBlockLaplaceOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >
{
 public:
   P2ViscousBlockLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : VectorToVectorOperator< real_t, hyteg::P2VectorFunction, hyteg::P2VectorFunction >( storage, minLevel, maxLevel )
   {
      this->setSubOperator( 0, 0, std::make_shared< operatorgeneration::P2ElementwiseDiffusion >( storage, minLevel, maxLevel ) );
      this->setSubOperator( 1, 1, std::make_shared< operatorgeneration::P2ElementwiseDiffusion >( storage, minLevel, maxLevel ) );
      if ( storage->hasGlobalCells() )
      {
         this->setSubOperator(
             2, 2, std::make_shared< operatorgeneration::P2ElementwiseDiffusion >( storage, minLevel, maxLevel ) );
      }
   }
};

/// P2ViscousBlockLaplaceOperator with IcosahedralShellMap blending. See documentation of P2ViscousBlockLaplaceOperator.
class P2ViscousBlockLaplaceIcosahedralShellMapOperator
: public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >
{
 public:
   P2ViscousBlockLaplaceIcosahedralShellMapOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                     uint_t                                     minLevel,
                                                     uint_t                                     maxLevel )
   : VectorToVectorOperator< real_t, hyteg::P2VectorFunction, hyteg::P2VectorFunction >( storage, minLevel, maxLevel )
   {
      this->setSubOperator(
          0,
          0,
          std::make_shared< operatorgeneration::P2ElementwiseDiffusionIcosahedralShellMap >( storage, minLevel, maxLevel ) );
      this->setSubOperator(
          1,
          1,
          std::make_shared< operatorgeneration::P2ElementwiseDiffusionIcosahedralShellMap >( storage, minLevel, maxLevel ) );
      if ( storage->hasGlobalCells() )
      {
         this->setSubOperator(
             2,
             2,
             std::make_shared< operatorgeneration::P2ElementwiseDiffusionIcosahedralShellMap >( storage, minLevel, maxLevel ) );
      }
   }
};

} // namespace operatorgeneration
} // namespace hyteg