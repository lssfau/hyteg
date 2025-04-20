/*
 * Copyright (c) 2025 Ponsuganth Ilangovan.
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

#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/types/Averaging.hpp"

namespace hyteg {

class P0toP0AveragedInjection : public RestrictionOperator< P0Function< real_t > >
{
 public:
   P0toP0AveragedInjection( hyteg::AveragingType averagingType, bool volumeWeighted )
   : averagingType_( averagingType )
   , volumeWeighted_( volumeWeighted )
   {}

   void restrict( const P0Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag = All ) const override;

   void restrictToAllLowerLevels( const P0Function< real_t >& function, const uint_t& sourceLevel );

 private:
   hyteg::AveragingType averagingType_;
   bool                 volumeWeighted_;
};

} // namespace hyteg
