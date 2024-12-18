/*
 * Copyright (c) 2024 Hamish Brown, Fatemeh Raezei, Eugenio D'Ascoli, Nils Kohl.
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

#include "terraneo/helpers/InterpolateProfile.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

using terraneo::SphericalHarmonicsTool;

namespace terraneo {

using hyteg::Point3D;
using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

/// Simple data container holding parameters that are often required for temperature initialization on the spherical shell.
struct TemperatureInitializationParameters
{
   /// Initializing the struct.
   ///
   /// \param Tcmb               Temperature at the core-mantle boundary.
   /// \param Tsurface           Temperature at the surface.
   /// \param TsurfaceAdb        Adiabatic temperature profile throughout the domain.
   /// \param dissipationNumber  Physical, non-dimensional dissipation number.
   /// \param rMin               Minimal radius.
   /// \param rMax               Maximal radius.
   /// \param devParams          Parameters for the deviation from the reference profile.
   TemperatureInitializationParameters( real_t                                        Tcmb,
                                        real_t                                        Tsurface,
                                        real_t                                        TsurfaceAdb,
                                        real_t                                        dissipationNumber,
                                        real_t                                        rMin,
                                        real_t                                        rMax,
                                        TemperatureDeivationInitialisationParameters* devParams )
   : Tcmb_( Tcmb )
   , Tsurface_( Tsurface )
   , TsurfaceAdb_( TsurfaceAdb )
   , dissipationNumber_( dissipationNumber )
   , rMin_( rMin )
   , rMax_( rMax )
   , deviationParameters_( devParams )
   {
      WALBERLA_CHECK_LESS_EQUAL( rMin, rMax );
   }

   real_t                                        Tcmb() const { return Tcmb_; }
   real_t                                        Tsurface() const { return Tsurface_; }
   real_t                                        TsurfaceAdb() const { return TsurfaceAdb_; }
   real_t                                        dissipationNumber() const { return dissipationNumber_; }
   real_t                                        rMin() const { return rMin_; }
   real_t                                        rMax() const { return rMax_; }
   TemperatureDeivationInitialisationParameters* deviationParameters() const { return deviationParameters_; }

 private:
   real_t                                        Tcmb_;
   real_t                                        Tsurface_;
   real_t                                        TsurfaceAdb_;
   real_t                                        dissipationNumber_;
   real_t                                        rMin_;
   real_t                                        rMax_;
   TemperatureDeivationInitialisationParameters* deviationParameters_;
};

/// Constructs a std::function that initializes the exponential reference profile
///
///   ( TsurfaceAdb * exp( dissipationNumber * ( rMax - r ) ) ) / ( Tcmb - Tsurface )
///
/// with parameters given as TemperatureInitializationParameters struct.
///
/// The returned std::function can be used for interpolation of a HyTeG FE function:
///
/// \code
///
///     P2Function< real_t > temperature( "temperature", storage, minLevel, maxLevel );
///
///     TemperatureInitializationParameters tempInitParams( /* ... */ );
///     auto referenceTemp = temperatureReferenceExponential( tempInitParams );
///
///     temperature.interpolate( referenceTemp, level );
///
/// \endcode
///
/// \param tempInitParams TemperatureInitializationParameters struct initialized with suitable parameters
/// \return std::function that can be passed to HyTeG's FE function interpolate() method
///
inline std::function< real_t( const Point3D& ) >
    temperatureReferenceExponential( const TemperatureInitializationParameters& tempInitParams )
{
   return [=]( const Point3D& x ) {
      auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      real_t temp =
          tempInitParams.TsurfaceAdb() * std::exp( ( tempInitParams.dissipationNumber() * ( tempInitParams.rMax() - radius ) ) );

      real_t retVal = ( temp ) / ( tempInitParams.Tcmb() - tempInitParams.Tsurface() );

      return retVal;
   };
}

inline std::function< real_t( const Point3D& ) >
    temperatureReferenceProfile( const TemperatureInitializationParameters& tempInitParams,
                                 const std::vector< real_t >                radius,
                                 const std::vector< real_t >                temperatureProfile )
{
   return [=]( const Point3D& x ) {
      real_t temp = interpolateDataValues( x, radius, temperatureProfile, tempInitParams.rMin(), tempInitParams.rMax() );

      real_t retVal = ( temp ) / ( tempInitParams.Tcmb() - tempInitParams.Tsurface() );

      return retVal;
   };
}

/// Constructs a std::function that adds white noise to a reference profile as
///
///   refProfile + noiseFactor * rand( -1, 1 )
///
/// with parameters given as TemperatureInitializationParameters struct, a noise factor and a reference temperature.
///
/// The returned std::function can be used for interpolation of a HyTeG FE function:
///
/// \code
///
///     P2Function< real_t > temperature( "temperature", storage, minLevel, maxLevel );
///
///     TemperatureInitializationParameters tempInitParams( /* ... */ );
///     auto referenceTemp  = temperatureReferenceExponential( tempInitParams );
///     auto whiteNoiseTemp = temperatureWhiteNoise( tempInitParams, referenceTemp );
///
///     temperature.interpolate( whiteNoiseTemp, level );
///
/// \endcode
///
/// \param tempInitParams TemperatureInitializationParameters struct initialized with suitable parameters
/// \param referenceTemp  a std::function that represents a some reference temperature profile
/// \return std::function that can be passed to HyTeG's FE function interpolate() method
///
inline std::function< real_t( const Point3D& ) >
    temperatureWhiteNoise( const TemperatureInitializationParameters&       tempInitParams,
                           const std::function< real_t( const Point3D& ) >& referenceTemp )
{
   return [=]( const hyteg::Point3D& x ) {
      auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      real_t retVal;

      const auto rMin = tempInitParams.rMin();
      const auto rMax = tempInitParams.rMax();

      const auto Tcmb     = tempInitParams.Tcmb();
      const auto Tsurface = tempInitParams.Tsurface();

      const auto tempDevInitParams = tempInitParams.deviationParameters();
      const auto buoyancyFactor    = tempDevInitParams->buoyancyFactor;

      // Boundaries
      if ( ( radius - rMin ) < real_c( 1e-10 ) )
      {
         return ( tempInitParams.Tcmb() ) / ( tempInitParams.Tcmb() - tempInitParams.Tsurface() );
      }
      else if ( ( rMax - radius ) < real_c( 1e-10 ) )
      {
         return ( tempInitParams.Tsurface() ) / ( tempInitParams.Tcmb() - tempInitParams.Tsurface() );
      }
      else
      {
         retVal = referenceTemp( x );

         // Random generator for Temperature initialisation ( Gaussian White Noise (GWN))
         retVal += buoyancyFactor * retVal * walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
      }
      return retVal;
   };
}

/// Constructs a std::function that initializes a temperature profile with deviation of a single spherical harmonic.
///
/// The returned std::function can be used for interpolation of a HyTeG FE function. See documentation of
/// terraneo::temperatureWhiteNoise for a similar usage example.
///
/// \param tempInitParams               TemperatureInitializationParameters struct initialized with suitable parameters
/// \param referenceTemp                A std::function that represents a some reference temperature profile
/// \return std::function that can be passed to HyTeG's FE function interpolate() method
///
inline std::function< real_t( const Point3D& ) >
    temperatureSingleSPH( const TemperatureInitializationParameters&       tempInitParams,
                          const std::function< real_t( const Point3D& ) >& referenceTemp )
{
   return [=]( const hyteg::Point3D& x ) {
      const auto rMin     = tempInitParams.rMin();
      const auto rMax     = tempInitParams.rMax();
      const auto Tcmb     = tempInitParams.Tcmb();
      const auto Tsurface = tempInitParams.Tsurface();

      auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      const auto tempDevInitParams           = tempInitParams.deviationParameters();
      const auto initialTemperatureSteepness = tempDevInitParams->initialTemperatureSteepness;
      const auto tempInit                    = tempDevInitParams->tempInit;

      const auto  deg            = tempDevInitParams->deg;
      const auto  ord            = tempDevInitParams->ord;
      const auto  buoyancyFactor = tempDevInitParams->buoyancyFactor;
      const auto& sphTool        = tempDevInitParams->sphTool;

      real_t retVal = referenceTemp( x );

      // Boundaries
      if ( ( radius - rMin ) < real_c( 1e-10 ) )
      {
         return ( Tsurface ) / ( Tcmb - Tsurface );
      }
      else if ( ( rMax - radius ) < real_c( 1e-10 ) )
      {
         return ( Tsurface ) / ( Tcmb - Tsurface );
      }

      real_t filter;

      switch ( tempInit )
      {
      case 0: // Anomalies near the CMB
         filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin ) / ( rMax - rMin ) ) );
         break;
      case 1: // Anomalies near the surface
         filter = std::exp( initialTemperatureSteepness * ( ( radius - rMax ) / ( rMax - rMin ) ) );
         break;
      case 2: // Anomalies near the surface and CMB
         filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin ) / ( rMax - rMin ) ) ) +
                  std::exp( initialTemperatureSteepness * ( ( radius - rMax ) / ( rMax - rMin ) ) );
         break;
      default: // No filtering
         filter = real_c( 1 );
         break;
      }

      retVal += ( buoyancyFactor * filter * std::sin( walberla::math::pi * ( radius - rMin ) / ( rMax - rMin ) ) *
                  sphTool->shconvert_eval( deg, ord, x[0], x[1], x[2] ) / std::sqrt( real_c( 4 ) * walberla::math::pi ) );

      return retVal;
   };
}

/// Constructs a std::function that initializes a temperature profile with superposition of spherical harmonics.
///
/// The returned std::function can be used for interpolation of a HyTeG FE function. See documentation of
/// terraneo::temperatureWhiteNoise for a similar usage example.
///
/// \param tempInitParams               TemperatureInitializationParameters struct initialized with suitable parameters
/// \param referenceTemp                A std::function that represents a some reference temperature profile
/// \return std::function that can be passed to HyTeG's FE function interpolate() method
inline std::function< real_t( const Point3D& ) >
    temperatureRandomSuperpositioneSPH( const TemperatureInitializationParameters&       tempInitParams,
                                       const std::function< real_t( const Point3D& ) >& referenceTemp )
{
   return [=]( const hyteg::Point3D& x ) {
      const auto rMin     = tempInitParams.rMin();
      const auto rMax     = tempInitParams.rMax();
      const auto Tcmb     = tempInitParams.Tcmb();
      const auto Tsurface = tempInitParams.Tsurface();

      auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      const auto tempDevInitParams           = tempInitParams.deviationParameters();
      const auto initialTemperatureSteepness = tempDevInitParams->initialTemperatureSteepness;
      const auto tempInit                    = tempDevInitParams->tempInit;

      const auto  lmax              = tempDevInitParams->lmax;
      const auto  lmin              = tempDevInitParams->lmin;
      const auto  buoyancyFactor    = tempDevInitParams->buoyancyFactor;
      const auto  superpositionRand = tempDevInitParams->superpositionRand;
      const auto& sphTool           = tempDevInitParams->sphTool;

      real_t retVal = referenceTemp( x );

      // Boundaries
      if ( ( radius - rMin ) < real_c( 1e-10 ) )
      {
         return ( Tsurface ) / ( Tcmb - Tsurface );
      }
      else if ( ( rMax - radius ) < real_c( 1e-10 ) )
      {
         return ( Tsurface ) / ( Tcmb - Tsurface );
      }

      real_t filter;

      switch ( tempInit )
      {
      case 0: // Anomalies near the CMB
         filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin ) / ( rMax - rMin ) ) );
         break;
      case 1: // Anomalies near the surface
         filter = std::exp( initialTemperatureSteepness * ( ( radius - rMax ) / ( rMax - rMin ) ) );
         break;
      case 2: // Anomalies near the surface and CMB
         filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin ) / ( rMax - rMin ) ) ) +
                  std::exp( initialTemperatureSteepness * ( ( radius - rMax ) / ( rMax - rMin ) ) );
         break;
      default: // No filtering
         filter = real_c( 1 );
         break;
      }

      uint_t count = 0;
      for ( uint_t deg = lmin; deg <= lmax; ++deg )
      {
         for ( int ord = -walberla::int_c( deg ); ord <= walberla::int_c( deg ); ++ord )
         {
            // Normalisation of 1/sqrt(4*pi) for non-dimensional temperature range [0,1]
            retVal += ( buoyancyFactor * superpositionRand[count] * filter *
                        sphTool->shconvert_eval( deg, ord, x[0], x[1], x[2] ) / std::sqrt( real_c( 4 ) * walberla::math::pi ) );
            ++count;
         }
      }

      return retVal;
   };
}
} // namespace terraneo