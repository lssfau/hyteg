/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

class P2toP2QuadraticRestriction : public RestrictionOperator< P2Function< real_t > >
{
 public:
   inline void restrict( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      if ( function.getStorage()->hasGlobalCells() )
      {
         restrictAdditively3D( function, sourceLevel, flag );
      }
      else
      {
         restrictAdditively( function, sourceLevel, flag );
      }
   }

 private:
   void restrictWithPostCommunication( const P2Function< real_t >& function,
                                       const uint_t&               sourceLevel,
                                       const DoFType&              flag ) const;
   void restrictAdditively( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
   void restrictAdditively3D( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
};

} // namespace hyteg
