/*
 * Copyright (c) 2020 Daniel Drzisga.
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

#include "hyteg/composites/P1VectorFunction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"

namespace hyteg {

class P1VectorFunctionLinearRestriction : public RestrictionOperator< P1VectorFunction< real_t > >
{
 public:

   void restrict( const P1VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      restrictionOperator_.restrict( function.u, sourceLevel, flag );
      restrictionOperator_.restrict( function.v, sourceLevel, flag );
      restrictionOperator_.restrict( function.w, sourceLevel, flag );
   }

 private:
   P1toP1LinearRestriction restrictionOperator_;
};
} // namespace hyteg