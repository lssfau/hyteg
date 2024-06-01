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
   TemperatureInitializationParameters( real_t Tcmb,
                                        real_t Tsurface,
                                        real_t TsurfaceAdb,
                                        real_t dissipationNumber,
                                        real_t rMin,
                                        real_t rMax )
   : Tcmb_( Tcmb )
   , Tsurface_( Tsurface )
   , TsurfaceAdb_( TsurfaceAdb )
   , dissipationNumber_( dissipationNumber )
   , rMin_( rMin )
   , rMax_( rMax )
   {
      WALBERLA_CHECK_LESS_EQUAL( rMin, rMax );
   }

   real_t Tcmb() const { return Tcmb_; }
   real_t Tsurface() const { return Tsurface_; }
   real_t TsurfaceAdb() const { return TsurfaceAdb_; }
   real_t dissipationNumber() const { return dissipationNumber_; }
   real_t rMin() const { return rMin_; }
   real_t rMax() const { return rMax_; }

 private:
   real_t Tcmb_;
   real_t Tsurface_;
   real_t TsurfaceAdb_;
   real_t dissipationNumber_;
   real_t rMin_;
   real_t rMax_;
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

      real_t retVal = temp / ( tempInitParams.Tcmb() - tempInitParams.Tsurface() );

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
///     auto whiteNoiseTemp = temperatureWhiteNoise( tempInitParams, referenceTemp, noiseFactor );
///
///     temperature.interpolate( whiteNoiseTemp, level );
///
/// \endcode
///
/// \param tempInitParams TemperatureInitializationParameters struct initialized with suitable parameters
/// \param referenceTemp  a std::function that represents a some reference temperature profile
/// \return std::function that can be passed to HyTeG's FE function interpolate() method
///
inline std::function< real_t( const Point3D& ) > temperatureWhiteNoise( const TemperatureInitializationParameters&       tempInitParams,
                                                                 const std::function< real_t( const Point3D& ) >& referenceTemp,
                                                                 real_t                                           noiseFactor )
{
   return [=]( const hyteg::Point3D& x ) {
      auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      real_t retVal;

      const auto rMin = tempInitParams.rMin();
      const auto rMax = tempInitParams.rMax();

      const auto Tcmb     = tempInitParams.Tcmb();
      const auto Tsurface = tempInitParams.Tsurface();

      // Boundaries
      if ( ( radius - rMin ) < real_c( 1e-10 ) )
      {
         return Tcmb / ( Tcmb - Tsurface );
      }
      else if ( ( rMax - radius ) < real_c( 1e-10 ) )
      {
         return Tsurface / ( Tcmb - Tsurface );
      }
      else
      {
         retVal = referenceTemp( x );

         // Random generator for Temperature initialisation ( Gaussian White Noise (GWN))
         retVal += noiseFactor * retVal * walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
      }
      return retVal;
   };
}

/// Constructs a std::function that initializes a temperature profile based on spherical harmonics.
///
/// The returned std::function can be used for interpolation of a HyTeG FE function. See documentation of
/// terraneo::temperatureWhiteNoise for a similar usage example.
///
/// \param tempInitParams               TemperatureInitializationParameters struct initialized with suitable parameters
/// \param referenceTemp                A std::function that represents a some reference temperature profile
/// \param tempInit                     Used to define where the anomalies are located in the domain.
/// \param deg                          Spherical Harmonics degree.
/// \param ord                          Spherical Harmonics order.
/// \param lmax                         Largest spherical harmonics degree.
/// \param lmin                         Smallest spherical harmonics degree.
/// \param superposition                Defines a random superposition of spherical harmonics up to a specific degree lmax.
/// \param buoyancyFactor               Factor to scale up or down the temperature anomalies.
/// \param initialTemperatureSteepness  Factor for filtering anomalies.
/// \return std::function that can be passed to HyTeG's FE function interpolate() method
///
inline std::function< real_t( const Point3D& ) > temperatureSPH( const TemperatureInitializationParameters&       tempInitParams,
                                                          const std::function< real_t( const Point3D& ) >& referenceTemp,
                                                          const uint_t                                     tempInit,
                                                          const uint_t                                     deg,
                                                          const int                                        ord,
                                                          const uint_t                                     lmax,
                                                          const uint_t                                     lmin,
                                                          const bool                                       superposition,
                                                          const real_t                                     buoyancyFactor,
                                                          const real_t initialTemperatureSteepness )
{
   return [=]( const hyteg::Point3D& x ) {
      const auto rMin     = tempInitParams.rMin();
      const auto rMax     = tempInitParams.rMax();
      const auto Tcmb     = tempInitParams.Tcmb();
      const auto Tsurface = tempInitParams.Tsurface();

      auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      real_t retVal = referenceTemp( x );

      // Boundaries
      if ( ( radius - rMin ) < real_c( 1e-10 ) )
      {
         return Tcmb / ( Tcmb - Tsurface );
      }
      else if ( ( rMax - radius ) < real_c( 1e-10 ) )
      {
         return Tsurface / ( Tcmb - Tsurface );
      }


      std::shared_ptr< SphericalHarmonicsTool > sphTool = std::make_shared< SphericalHarmonicsTool >( lmax );
      real_t                                    filter;

      switch ( tempInit )
      {
      // Anomalies near the CMB
      case 0:

         filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin ) / ( rMax - rMin ) ) );

         break;
         // Anomalies near the surface
      case 1:

         filter = std::exp( initialTemperatureSteepness * ( ( radius - rMax ) / ( rMax - rMin ) ) );

         break;
      // Anomalies near the surface and CMB
      case 2:

         filter = std::exp( -initialTemperatureSteepness * ( ( radius - rMin ) / ( rMax - rMin ) ) ) +
                  std::exp( initialTemperatureSteepness * ( ( radius - rMax ) / ( rMax - rMin ) ) );

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
                   ( buoyancyFactor * superpositionRand[count] * filter * sphTool->shconvert_eval( deg, ord, x[0], x[1], x[2] ) /
                     std::sqrt( real_c( 4 ) * walberla::math::pi ) );
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
}

} // namespace terraneo