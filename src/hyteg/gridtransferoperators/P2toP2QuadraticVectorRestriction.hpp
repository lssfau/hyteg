/*
 * Copyright (c) 2023 Andreas Burkhart.
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
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"

namespace hyteg {

class P2toP2QuadraticVectorRestriction : public RestrictionOperator< P2VectorFunction< real_t > >
{
 public:
   P2toP2QuadraticVectorRestriction() {}

   void restrict( const P2VectorFunction< real_t >& function,
                  const uint_t&                           sourceLevel,
                  const DoFType&                          flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         quadraticRestrictionOperator_.restrict( function[k], sourceLevel, flag );
      }
   }

 private:
   P2toP2QuadraticRestriction quadraticRestrictionOperator_;
};

class P2toP2QuadraticVectorRestrictionWithProjection : public P2toP2QuadraticVectorRestriction
{
 public:
   P2toP2QuadraticVectorRestrictionWithProjection( std::shared_ptr< P2ProjectNormalOperator > projection )
   : P2toP2QuadraticVectorRestriction()
   , projection_( projection )
   {}

   void restrict( const P2VectorFunction< real_t >& function,
                  const uint_t&                           sourceLevel,
                  const DoFType&                          flag ) const override
   {
      P2toP2QuadraticVectorRestriction::restrict(function, sourceLevel, flag);
      projection_->project( function, sourceLevel-1, FreeslipBoundary );
   }

 private:
   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

} // namespace hyteg
