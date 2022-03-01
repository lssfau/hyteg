/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"

namespace hyteg {

class P1P1StokesToP1P1StokesRestriction : public RestrictionOperator< P1StokesFunction< real_t > >
{
 public:
   P1P1StokesToP1P1StokesRestriction()
   : projectMeanAfterRestriction_( false )
   {}
   P1P1StokesToP1P1StokesRestriction( bool projectMeanAfterRestriction )
   : projectMeanAfterRestriction_( projectMeanAfterRestriction )
   {}

   typedef P1toP1LinearRestriction VelocityRestriction_T;
   typedef P1toP1LinearRestriction PressureRestriction_T;

   void restrict( const P1StokesFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      for ( uint_t k = 0; k < function.uvw().getDimension(); k++ )
      {
         restrictionOperator_.restrict( function.uvw()[k], sourceLevel, flag );
      }
      restrictionOperator_.restrict( function.p(), sourceLevel, flag );

      if ( projectMeanAfterRestriction_ )
      {
         vertexdof::projectMean( function.p(), sourceLevel - 1 );
      }
   }

 private:
   P1toP1LinearRestriction restrictionOperator_;

   bool projectMeanAfterRestriction_;
};
} // namespace hyteg
