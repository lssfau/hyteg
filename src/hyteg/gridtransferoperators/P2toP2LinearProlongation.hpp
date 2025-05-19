/*
 * Copyright (c) 2024 Ponsuganth Ilangovan
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

#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP2Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

// Highly inefficient?
class P2toP2LinearProlongation : public ProlongationOperator< P2Function< real_t > >
{
 public:
   P2toP2LinearProlongation( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   {
      temp_ = std::make_shared< P1Function< real_t > >( "temp__P2toP2LinearProlongation", storage, minLevel, maxLevel + 1 );
   }

   P2toP2LinearProlongation( const std::shared_ptr< P1Function< real_t > >& temp )
   : temp_( temp )
   {}

   inline void prolongate( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      P2toP1Conversion( function, *temp_, sourceLevel + 1 );
      p1Prolongation.prolongate(*temp_, sourceLevel + 1, flag);
      P1toP2Conversion( *temp_, function, sourceLevel + 1, flag );
   }

 private:
   std::shared_ptr< P1Function< real_t > > temp_;

   P1toP1LinearProlongation< real_t > p1Prolongation;
};

} // namespace hyteg