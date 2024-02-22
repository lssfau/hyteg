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
#include <core/DataTypes.h>
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/math/Random.h>
#include <vector>

#include "hyteg/types/PointND.hpp"
#include "hyteg/types/types.hpp"

#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

using terraneo::SphericalHarmonicsTool;

namespace terraneo {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

/// Initialises the temperature field for a given FunctionType (e.g. P2Function).
/// This allows to choose between two possible temperature intialisations, a random white noise temperature intialisation (Gaussian)
/// or a spherical harmonics temperature initialisation with the defined order and degree.
///
/// Usage Example:
/// \code
///
///     std::shared_ptr< FunctionType > temperature = std::make_shared< FunctionType >( "temperature", storage, minLevel, maxLevel );
///     std::shared_ptr< terraneo::TemperaturefieldConv< FunctionType > > Temperaturefield =
///     std::make_shared< terraneo::TemperaturefieldConv< FunctionType > >( temperature, TBottom, TTop, TAdiabate,
///     DissipationNumber, rMax, rMin, maxLevel, minLevel );
///     Temperaturefield->initialiseTemperatureWhiteNoise( noiseFactor );
///
/// \endcode
///
/// \tparam temperature        FunctionType of which our unknown of interest is composed of.
/// \tparam TBottom            Temperature at the bottom of the domain.
/// \tparam TTop               Temperature at the top of the domain.
/// \tparam TAdiabate          Adiabate Temperature profile throughout the domain.
/// \tparam DissipationNumber  Physical, non-dimensional dissipation Number.
/// \tparam rMax               Maximal radius.
/// \tparam rMin               Minimal radius.
/// \tparam maxLevel           Maximal refinement level.
/// \tparam minLevel           Minimal refinement level.

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
                         uint_t                           maxLevel,
                         uint_t                           minLevel )
   : T_( T )
   , Tcmb_( Tcmb )
   , Tsurface_( Tsurface )
   , TsurfaceAdb_( TsurfaceAdb )
   , dissipatioNumber_( dissipationNumber )
   , rMax_( rMax )
   , rMin_( rMin )
   , maxLevel_( maxLevel )
   , minLevel_( minLevel )
   {}

   real_t referenceTemperatureFct( const hyteg::Point3D& x )
   {
      auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      real_t temp   = TsurfaceAdb_ * std::exp( ( dissipatioNumber_ * ( rMax_ - radius ) ) );

      real_t retVal = temp / ( Tcmb_ - Tsurface_ );

      return retVal;
   }

   /// Interpolate random white noise onto the temperature field
   ///
   /// \param noiseFactor  Amount of noise [%] deviation from the background reference temperature.
   ///

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

      for ( uint_t l = minLevel_; l <= maxLevel_; l++ )
      {
         T_->interpolate( temperatureInit, l, hyteg::All );
      }
   }

   /// Interpolate spherical harmonics function as temperature deviation
   /// A spherical harmonics function of certain degree and order can be chosen
   /// as temperature deviation from the background reference temperature profile (e.b. adiabate T profile).
   ///
   /// \param tempInit                     Used to define where the anomalies are located in the domain.
   /// \param deg                          Spherical Harmonics degree.
   /// \param ord                          Spherical Harmonics order.
   /// \param lmax                         Largest spherical harmonics degree.
   /// \param lmin                         Smallest spherical harmonics degree.
   /// \param superposition                Defines a random superposition of spherical harmonics up to a specific degree lmax.
   /// \param buoyancyFactor               Factor to scale up or down the temperature anomalies.
   /// \param initialTemperatureSteepness  Factor for filtering anomalies.
   ///

   void initialiseTemperatureSPH( const uint_t& tempInit,
                                  const uint_t  deg,
                                  const int&    ord,
                                  const uint_t& lmax,
                                  const uint_t& lmin,
                                  const bool&   superposition,
                                  const real_t& buoyancyFactor,
                                  const real_t& initialTemperatureSteepness )
   {
      std::function< real_t( const hyteg::Point3D& ) > temperatureInitSPH = [&]( const hyteg::Point3D& x ) {
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

         retVal = referenceTemperatureFct( x );

         std::shared_ptr< SphericalHarmonicsTool > sphTool = std::make_shared< SphericalHarmonicsTool >( lmax );
         real_t                                    filter;

         switch ( tempInit )
         {
         // Anomalies near the CMB
         case 0:

            filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin_ ) / ( rMax_ - rMin_ ) ) );

            break;
            // Anomalies near the surface
         case 1:

            filter = std::exp( initialTemperatureSteepness * ( ( radius - rMax_ ) / ( rMax_ - rMin_ ) ) );

            break;
         // Anomalies near the surface and CMB
         case 2:

            filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin_ ) / ( rMax_ - rMin_ ) ) ) +
                     std::exp( initialTemperatureSteepness * ( ( radius - rMax_ ) / ( rMax_ - rMin_ ) ) );

            break;

         // No filtering
         default:

            filter = real_c( 1 );

            break;
         }

         if ( superposition )
         {
            uint_t                numHarmonics = ( ( lmax + 1 ) * ( lmax + 1 ) ) - ( lmin ) * ( lmin );
            std::vector< real_t > superpositionRand{};
            superpositionRand.reserve( numHarmonics );
            walberla::math::seedRandomGenerator( 42 );
            for ( uint_t i = 0; i < numHarmonics; i++ )
            {
               superpositionRand.push_back( walberla::math::realRandom( real_c( -1 ), real_c( 1 ) ) );
            }
            uint_t count = 0;

            for ( uint_t deg = lmin; deg <= lmax; ++deg )
            {
               for ( int ord = -walberla::int_c( deg ); ord <= walberla::int_c( deg ); ++ord )
               {
                  // Normalisation of 1/sqrt(4*pi) for non-dimensional temperature range [0,1]
                  retVal +=
                      ( buoyancyFactor * superpositionRand[count] * filter *
                        sphTool->shconvert_eval( deg, ord, x[0], x[1], x[2] ) / std::sqrt( real_c( 4 ) * walberla::math::pi ) );
                  ++count;
               }
            }
         }

         else
         {
            // Single spherical harmonic for the initialisation.
            retVal += ( buoyancyFactor * filter * sphTool->shconvert_eval( deg, ord, x[0], x[1], x[2] ) /
                        std::sqrt( real_c( 4 ) * walberla::math::pi ) );
         }
         return retVal;
      };
      for ( uint_t l = minLevel_; l <= maxLevel_; l++ )
      {
         T_->interpolate( temperatureInitSPH, l, hyteg::All );
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
   uint_t maxLevel_;
   uint_t minLevel_;

   std::shared_ptr< FunctionType >& T_;
};
} // namespace terraneo