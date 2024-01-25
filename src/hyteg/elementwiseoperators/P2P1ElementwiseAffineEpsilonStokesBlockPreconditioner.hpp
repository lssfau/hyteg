/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include <hyteg/forms/form_hyteg_generated/p1/p1_invk_mass_affine_q4.hpp>

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "constant_stencil_operator/P2ConstantEpsilonOperator.hpp"
#include "mixed_operator/P2P1TaylorHoodStokesBlockPreconditioner.hpp"

namespace hyteg {

class P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     minLevel,
                                                          uint_t                                     maxLevel,
                                                          std::function< real_t( const Point3D& ) >  mu )
   : Operator( storage, minLevel, maxLevel )
   , viscOp( storage, minLevel, maxLevel, mu )
   , P( storage, minLevel, maxLevel, std::make_shared< P1RowSumForm >( std::make_shared< forms::p1_invk_mass_affine_q4 >( mu, mu ) ) )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      viscOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      P.toMatrix( mat, src.p(), dst.p(), level, flag );
   }

   P2ElementwiseAffineEpsilonOperator viscOp;
   P1BlendingLumpedDiagonalOperator P;

   bool hasGlobalCells_;
};

} // namespace hyteg
