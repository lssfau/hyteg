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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg_operators_composites/stokes/divergence/P2ToP1DivergenceOperator.hpp"
#include "hyteg_operators_composites/stokes/gradient/P1ToP2GradientOperator.hpp"
#include "hyteg_operators_composites/stokes/viscousblock/P2ViscousBlockEpsilonOperator.hpp"

namespace hyteg {
namespace operatorgeneration {

/// Implements the "epsilon" Stokes operator, i.e., the discrete operator that arises from the discretization of
///
///     - ∇ · ( 2μ ε(u) ) + ∇ p = f,
///       ∇ · u                 = g.
///
/// where
///
///     ε(w) := (1/2) (∇w + (∇w)ᵀ)
///
/// and μ = μ(x) is a user-supplied finite element function.
///
/// It combines the viscous "epsilon" operator for A with the discrete divergence B and gradient B to correspond to the block
/// matrix
///
///         /        \
///     K = |  A  Bᵀ |
///         |  B  0  |
///         \        /
///
class P2P1StokesEpsilonOperator : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   P2P1StokesEpsilonOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                              uint_t                                     minLevel,
                              uint_t                                     maxLevel,
                              const P2Function< real_t >&                mu )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel, mu )
   , BT( storage, minLevel, maxLevel )
   , B( storage, minLevel, maxLevel )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               const UpdateType                        updateType = Replace ) const
   {
      A.apply( src.uvw(), dst.uvw(), level, flag, updateType );
      BT.apply( src.p(), dst.uvw(), level, flag, Add );
      B.apply( src.uvw(), dst.p(), level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      BT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      B.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   operatorgeneration::P2ViscousBlockEpsilonOperator A;
   operatorgeneration::P1ToP2GradientOperator        BT;
   operatorgeneration::P2ToP1DivergenceOperator      B;
};

} // namespace operatorgeneration
} // namespace hyteg
