/*
 * Copyright (c) 2024 Hamish Brown, Fatemeh Raezei, Eugenio D'Ascoli.
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

#include <cmath>
#include <core/math/Constants.h>
#include <core/math/Random.h>
#include <vector>

#include "core/DataTypes.h"
#include "core/Environment.h"

#include "hyteg/types/PointND.hpp"
#include "hyteg/types/types.hpp"

#include "typeAliases.hpp"

namespace terraneo {
template < typename FunctionType >
class TemperaturefieldConv
{
 public:
   TemperaturefieldConv( std::shared_ptr< FunctionType >& T,
                         real_t                           Tcmb,
                         real_t                           Tsurface,
                         real_t                           TsurfaceAdb,
                         real_t                           dissipationNumber,
                         real_t                           rMax,
                         real_t                           rMin,
                         uint_t                           lMax,
                         uint_t                           lMin )
   : T_( T )
   , Tcmb_( Tcmb )
   , Tsurface_( Tsurface )
   , TsurfaceAdb_( TsurfaceAdb )
   , dissipatioNumber_( dissipationNumber )
   , rMax_( rMax )
   , rMin_( rMin )
   , lMax_( lMax )
   , lMin_( lMin )
   {}

   real_t referenceTemperatureFct( const hyteg::Point3D& x )
   {
      auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      real_t temp   = TsurfaceAdb_ * std::exp( ( dissipatioNumber_ * ( rMax_ - radius ) ) );

      real_t retVal = temp / ( Tcmb_ - Tsurface_ );

      return retVal;
   }

   void initialiseTemperatureWhiteNoise( real_t noiseFactor )
   {
      std::function< real_t( const hyteg::Point3D& ) > temperatureInit = [&]( const hyteg::Point3D& x ) {
         auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         real_t retVal;

         // Boundaries
         if ( ( radius - rMin_ ) < real_c( 1e-10 ) )
         {
            return Tcmb_ / ( Tcmb_ - Tsurface_ );
         }
         else if ( ( rMax_ - radius ) < real_c( 1e-10 ) )
         {
            return Tsurface_ / ( Tcmb_ - Tsurface_ );
         }
         else
         {
            retVal = referenceTemperatureFct( x );

            // Random generator for Temperature initialisation ( Gaussian White Noise (GWN))
            retVal += noiseFactor * retVal * walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
         }
         return retVal;
      };

      for ( uint_t l = lMin_; l <= lMax_; l++ )
      {
         T_->interpolate( temperatureInit, l, hyteg::All );
      }
   }

 private:
   real_t Tcmb_;
   real_t Tsurface_;
   real_t TsurfaceAdb_;
   real_t dissipatioNumber_;
   real_t rMax_;
   real_t rMin_;
   real_t noiseFactor_;
   uint_t lMax_;
   uint_t lMin_;

   std::shared_ptr< FunctionType >& T_;
};
} // namespace terraneo