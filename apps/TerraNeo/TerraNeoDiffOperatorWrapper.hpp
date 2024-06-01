/*
 * Copyright (c) 2024 Eugenio D'Ascoli
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

#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"
#include "hyteg_operators_composites/viscousblock/P2ViscousBlockLaplaceOperator.hpp"

#include "mixed_operator/ScalarToVectorOperator.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/helpers/TerraNeoParameters.hpp"

namespace terraneo {

using diffusionOperator = operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap;
using P2MassOperator    = hyteg::P2ElementwiseBlendingMassOperator;

class P2DiffusionOperatorWrapper : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2DiffusionOperatorWrapper( const std::shared_ptr< PrimitiveStorage >&     storage,
                               uint_t                                         minLevel,
                               uint_t                                         maxLevel,
                               real_t                                         timestepDt,
                               const std::shared_ptr< P2Function< real_t > >& diffusionFE,
                               const std::shared_ptr< P2Function< real_t > >& energyRHS )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , timestepDt_( timestepDt )
   , diffusionFE_( diffusionFE )
   , energyRHS_( energyRHS )
   {
      tmp1_ = std::make_shared< P2Function< real_t > >( "tmp1", storage_, minLevel_, maxLevel_ );
   }

   void initDifOperator()
   {
      massOperator_      = std::make_shared< hyteg::P2ElementwiseBlendingMassOperator >( storage_, minLevel_, maxLevel_ );
      diffusionOperator_ = std::make_shared< operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap >(
          storage_, minLevel_, maxLevel_, *diffusionFE_ );
   }
   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override
   {
      massOperator_->apply( src, dst, level, flag, updateType );
      diffusionOperator_->apply( src, *tmp1_, level, flag );
      dst.assign( { 1.0, timestepDt_ }, { dst, *tmp1_ }, level, flag );
   }

   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;
   real_t                              timestepDt_;
   // P2 Functions
   std::shared_ptr< P2Function< real_t > > diffusionFE_;
   std::shared_ptr< P2Function< real_t > > energyRHS_;
   std::shared_ptr< P2Function< real_t > > tmp1_;
   // Operators

   std::shared_ptr< P2MassOperator >    massOperator_;
   std::shared_ptr< diffusionOperator > diffusionOperator_;
};
} // namespace terraneo