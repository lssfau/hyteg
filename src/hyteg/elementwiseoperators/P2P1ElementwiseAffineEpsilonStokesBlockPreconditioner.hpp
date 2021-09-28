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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesBlockPreconditioner.hpp"
#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_affine_qe.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

class P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     minLevel,
                                                          uint_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A_0_0( storage, minLevel, maxLevel )
   , A_0_1( storage, minLevel, maxLevel )
   , A_0_2( storage, minLevel, maxLevel )
   , A_1_0( storage, minLevel, maxLevel )
   , A_1_1( storage, minLevel, maxLevel )
   , A_1_2( storage, minLevel, maxLevel )
   , A_2_0( storage, minLevel, maxLevel )
   , A_2_1( storage, minLevel, maxLevel )
   , A_2_2( storage, minLevel, maxLevel )
   , P( storage, minLevel, maxLevel, std::make_shared< P1RowSumForm >( std::make_shared< forms::p1_mass_affine_qe >() ) )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< matIdx_t >&   src,
                  const P2P1TaylorHoodFunction< matIdx_t >&   dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A_0_0.toMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
      A_0_1.toMatrix( mat, src.uvw[1], dst.uvw[0], level, flag );
      if ( src.getStorage()->hasGlobalCells() )
      {
         A_0_2.toMatrix( mat, src.uvw[2], dst.uvw[0], level, flag );
      }

      A_1_0.toMatrix( mat, src.uvw[0], dst.uvw[1], level, flag );
      A_1_1.toMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );
      if ( src.getStorage()->hasGlobalCells() )
      {
         A_1_2.toMatrix( mat, src.uvw[2], dst.uvw[1], level, flag );
      }

      if ( src.getStorage()->hasGlobalCells() )
      {
         A_2_0.toMatrix( mat, src.uvw[0], dst.uvw[2], level, flag );
         A_2_1.toMatrix( mat, src.uvw[1], dst.uvw[2], level, flag );
         A_2_2.toMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
      }

      P.toMatrix( mat, src.p, dst.p, level, flag );
   }

   P2ElementwiseOperator< forms::p2_epsiloncc_0_0_affine_q2 > A_0_0;
   P2ElementwiseOperator< forms::p2_epsiloncc_0_1_affine_q2 > A_0_1;
   P2ElementwiseOperator< forms::p2_epsiloncc_0_2_affine_q2 > A_0_2;

   P2ElementwiseOperator< forms::p2_epsiloncc_1_0_affine_q2 > A_1_0;
   P2ElementwiseOperator< forms::p2_epsiloncc_1_1_affine_q2 > A_1_1;
   P2ElementwiseOperator< forms::p2_epsiloncc_1_2_affine_q2 > A_1_2;

   P2ElementwiseOperator< forms::p2_epsiloncc_2_0_affine_q2 > A_2_0;
   P2ElementwiseOperator< forms::p2_epsiloncc_2_1_affine_q2 > A_2_1;
   P2ElementwiseOperator< forms::p2_epsiloncc_2_2_affine_q2 > A_2_2;

   P1BlendingLumpedDiagonalOperator P;

   bool hasGlobalCells_;
};

} // namespace hyteg
