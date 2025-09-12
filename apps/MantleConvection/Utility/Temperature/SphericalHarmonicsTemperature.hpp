/*
 * Copyright (c) 2024-2025 Andreas Burkhart, Eugenio D'Ascoli, Hamish Brown, Marcus Mohr and Markus Wiedemann.
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

#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "TemperatureModel.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

using walberla::int_c;
using walberla::real_t;
using walberla::math::pi;

namespace MantleConvection {

class SphericalHarmonicsTemperature : public TemperatureModel< real_t >
{
 public:
   using TemperatureModel< real_t >::prefix_;

   SphericalHarmonicsTemperature( NondimensionalisationParameters&                     nondimensionalisation,
                                  walberla::Config::BlockHandle&                       parameters,
                                  const std::shared_ptr< TemperatureModel< real_t > >& ModifiedModel,
                                  const uint_fast32_t                                  seed   = 392,
                                  std::string                                          prefix = "" )
   : TemperatureModel< real_t >( prefix )
   , ModifiedModel_( ModifiedModel )
   , randomSeed_( seed )
   {
      temperatureSurface_ = nondimensionalisation.temperatureSurface_;
      temperatureCMB_     = nondimensionalisation.temperatureCMB_;
      radiusCMB_          = nondimensionalisation.radiusCMB_;
      radiusSurface_      = nondimensionalisation.radiusSurface_;

      initialTemperatureFilter_    = parameters.getParameter< uint_t >( prefix + std::string( "initialTemperatureFilter" ) );
      initialTemperatureSteepness_ = parameters.getParameter< real_t >( prefix + std::string( "initialTemperatureSteepness" ) );
      degreeMinSH_                 = parameters.getParameter< uint_t >( prefix + std::string( "degreeMinSH" ) );
      degreeMaxSH_                 = parameters.getParameter< uint_t >( prefix + std::string( "degreeMaxSH" ) );
      buoyancyFactor_              = parameters.getParameter< real_t >( prefix + std::string( "buoyancyFactor" ) );

      // determine number of spherical harmonics
      uint_t numHarmonics_ = ( ( degreeMaxSH_ + 1 ) * ( degreeMaxSH_ + 1 ) ) - ( degreeMinSH_ ) * ( degreeMinSH_ );

      superpositionRand_.reserve( numHarmonics_ );

      // random generator
      std::mt19937 RndGenerator;
      RndGenerator.seed( seed );
      std::uniform_real_distribution< real_t > RndDistribution( real_c( -1 ), real_c( 1 ) );

      // generate spherical harmonics
      for ( uint_t i = 0; i < numHarmonics_; i++ )
      {
         superpositionRand_.push_back( RndDistribution( RndGenerator ) );
      }

      sphTool_ = std::make_shared< terraneo::SphericalHarmonicsTool >( degreeMaxSH_ );
   }

   real_t evaluateNoCheck( const hyteg::Point3D& x ) override
   {
      real_t radius = x.norm();
      real_t retVal = ModifiedModel_->evaluateNoCheck( x );

      real_t filter;

      switch ( initialTemperatureFilter_ )
      {
      case 0: //keeps anomalies near the CMB

         filter = std::exp( -initialTemperatureSteepness_ * ( ( radius - radiusCMB_ ) / ( radiusSurface_ - radiusCMB_ ) ) );

         break;
      case 1: //keeps anomalies near the surface

         filter = std::exp( initialTemperatureSteepness_ * ( ( radius - radiusSurface_ ) / ( radiusSurface_ - radiusCMB_ ) ) );

         break;
      case 2: //keeps anomalies near the surface and CMB

         filter = std::exp( -initialTemperatureSteepness_ * ( ( radius - radiusCMB_ ) / ( radiusSurface_ - radiusCMB_ ) ) ) +
                  std::exp( initialTemperatureSteepness_ * ( ( radius - radiusSurface_ ) / ( radiusSurface_ - radiusCMB_ ) ) );

         break;
      default: //default, no filtering

         filter = real_c( 1.0 );

         break;
      }

      // 1/sqrt(4*pi) below changes the normalisation of the spherical harmonics, more suitable for our temp range of 0, 1
      uint_t count = 0;

      for ( uint_t deg = degreeMinSH_; deg <= degreeMaxSH_; ++deg )
      {
         for ( int ord = -int_c( deg ); ord <= int_c( deg ); ++ord )
         {
            // the "buoyancy factor" allows the user to scale up/down the temperature anomalies in the initial state from the parameter file
            // the required factor, however, will need to be tuned depending on the superposition in question
            // a more user-friendly solution is in the works
            retVal += ( buoyancyFactor_ * superpositionRand_[count] * filter *
                        sphTool_->shconvert_eval( deg, ord, x[0], x[1], x[2] ) / std::sqrt( real_c( 4 ) * pi ) );
            ++count;
         }
      }

      return std::max( std::min( retVal, temperatureCMB_ ), temperatureSurface_ );
   }

   real_t evaluate( const hyteg::Point3D& x ) override
   {
      if ( checkSurface( x ) )
      {
         return temperatureSurface_;
      }
      else if ( checkCMB( x ) )
      {
         return temperatureCMB_;
      }

      return evaluateNoCheck( x );
   }

   bool checkSurface( const hyteg::Point3D& x ) override { return ModifiedModel_->checkSurface( x ); }
   bool checkCMB( const hyteg::Point3D& x ) override { return ModifiedModel_->checkCMB( x ); }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#####################################################"                                                                  << "\n";
      os << std::string( offset, ' ') << "########## Spherical Harmonics Temperature ##########"                                                                  << "\n";
      os << std::string( offset, ' ') << "#####################################################"                                                                  << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "prefix_: "                         << prefix_                               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "temperatureSurface_: "             << temperatureSurface_                   << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "temperatureCMB_: "                 << temperatureCMB_                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "radiusSurface_: "                  << radiusSurface_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "radiusCMB_: "                      << radiusCMB_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "initialTemperatureFilter_: "       << initialTemperatureFilter_             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "initialTemperatureSteepness_: "    << initialTemperatureSteepness_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "degreeMinSH_: "                    << degreeMinSH_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "degreeMaxSH_: "                    << degreeMaxSH_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "buoyancyFactor_: "                 << buoyancyFactor_                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 31 ) << std::left << "randomSeed_: "                     << randomSeed_                           << "\n";
      // clang-format on

      return ModifiedModel_->print( os, offset + 3 );
   }

 public:
   std::shared_ptr< TemperatureModel< real_t > > ModifiedModel_;

   real_t temperatureSurface_; // nondimensional temperature surface
   real_t temperatureCMB_;     // nondimensional temperature CMB
   real_t radiusSurface_;      // nondimensional radius surface, rMax
   real_t radiusCMB_;          // nondimensional radius cmb, rMin

   uint_t initialTemperatureFilter_;    // filter mode for the spherical harmonics temperature initialisation
   real_t initialTemperatureSteepness_; // initial temperature steepness for the spherical harmonics temperature initialisation
   uint_t degreeMinSH_;                 // minimum spherical harmoics degree for the temperature initialisation
   uint_t degreeMaxSH_;                 // maximum spherical harmoics degree for the temperature initialisation
   real_t buoyancyFactor_;              // buoyancy factor for the spherical harmonics temperature initialisation

   std::vector< real_t > superpositionRand_;                     // random value vector for the spherical harmoics superposition
   std::shared_ptr< terraneo::SphericalHarmonicsTool > sphTool_; // terraneo spherical harmonics tool
   uint_fast32_t                                       randomSeed_; // random seed;
};

inline std::ostream& operator<<( std::ostream& os, const SphericalHarmonicsTemperature& sht )
{
   return sht.print( os );
}

} // namespace MantleConvection