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

#include "hyteg/composites/CCRStokesFunction.hpp"
#include "hyteg/mixedoperators/DG1ToP2PlusBubbleOperator.hpp"
#include "hyteg/mixedoperators/P2PlusBubbleToDG1Operator.hpp"

#include "mixed_operator/ScalarToVectorOperator.hpp"
#include "mixed_operator/VectorLaplaceOperator.hpp"
#include "mixed_operator/VectorToScalarOperator.hpp"

namespace hyteg {

class CCRStokesOperator : public Operator< CCRStokesFunction< real_t >, CCRStokesFunction< real_t > >
{
 public:
   typedef operatorgeneration::P2PlusBubbleElementwiseDiffusion VelocityOperator_T;

   CCRStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , Lapl( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   // , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const CCRStokesFunction< real_t >& src,
               const CCRStokesFunction< real_t >& dst,
               const size_t                       level,
               DoFType                            flag,
               UpdateType                         updateType = Replace ) const override
   {
      WALBERLA_CHECK( !hasGlobalCells_, "CCRStokesOperator does not support 3D simlations, yet!" );

      Lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   const operatorgeneration::P2PlusBubbleElementwiseDiffusion& getA() const
   {
      auto ptr = Lapl.getSubOperator( 0, 0 );
      return dynamic_cast< const operatorgeneration::P2PlusBubbleElementwiseDiffusion& >( *ptr );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const CCRStokesFunction< idx_t >&           src,
                  const CCRStokesFunction< idx_t >&           dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      Lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   P2PlusBubbleVectorLaplaceOperator Lapl;
   P2PlusBubbleToDG1DivOperator      div;
   DG1ToP2PlusBubbleDivTOperator     divT;

   bool hasGlobalCells_;
};

} // namespace hyteg
