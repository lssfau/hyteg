/*
 * Copyright (c) 2024-2025 Ponsuganth Ilango, Andreas Burkhart.
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

#include "hyteg/gridtransferoperators/P1toP1InjectionRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP2Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

class P2toP2InjectionRestriction : public RestrictionOperator< P2Function< real_t > >
{
 public:
   void restrict( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      std::shared_ptr< P1Function< real_t > > tmp_ = getTemporaryFunction< P1Function< real_t > >(
          function.getStorage(), function.getMinLevel(), function.getMaxLevel() + 1, true );
      P2toP1Conversion( function, *tmp_, sourceLevel + 1 );
      p1Injection.restrict( *tmp_, sourceLevel + 1, flag );
      P1toP2Conversion( *tmp_, function, sourceLevel - 1, flag );
   }

 private:
   P1toP1InjectionRestriction p1Injection;
};

class P2toP2VectorInjectionRestriction : public RestrictionOperator< P2VectorFunction< real_t > >
{
 public:
   void restrict( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         scalarInjectionOperator_.restrict( function[k], sourceLevel, flag );
      }
   }

 private:
   P2toP2InjectionRestriction scalarInjectionOperator_;
};

} // namespace hyteg