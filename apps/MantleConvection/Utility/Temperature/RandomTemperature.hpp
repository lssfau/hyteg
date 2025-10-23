/*
 * Copyright (c) 2024-2025 Andreas Burkhart, Fatemeh Rezaei.
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

#include "core/math/Random.h"

#include "TemperatureModel.hpp"

using walberla::real_t;

namespace MantleConvection {

class RandomTemperature : public TemperatureModel< real_t >
{
 public:
   using TemperatureModel< real_t >::prefix_;

   RandomTemperature( NondimensionalisationParameters&                 nondimensionalisation,
                      walberla::Config::BlockHandle&                   parameters,
                      std::function< bool( const hyteg::Point3D& x ) > surfaceFct,
                      std::function< bool( const hyteg::Point3D& x ) > CMBFct,
                      const uint_fast32_t                              seed   = 393,
                      std::string                                      prefix = "" )
   : TemperatureModel< real_t >( prefix )
   , surfaceFct_( surfaceFct )
   , CMBFct_( CMBFct )
   , RndGenerator_( seed )
   , randomSeed_( seed )
   {
      temperatureSurface_ = nondimensionalisation.temperatureSurface_;
      temperatureCMB_     = nondimensionalisation.temperatureCMB_;

      RndDistribution_ = std::uniform_real_distribution< real_t >( temperatureSurface_, temperatureCMB_ );
   }

   real_t evaluateNoCheck( const hyteg::Point3D& x ) override { return RndDistribution_( RndGenerator_ ); }

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

   bool checkSurface( const hyteg::Point3D& x ) override { return surfaceFct_( x ); }
   bool checkCMB( const hyteg::Point3D& x ) override { return CMBFct_( x ); }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "######################################################"                                    << "\n";
      os << std::string( offset, ' ') << "################# Random Temperature #################"                                    << "\n";
      os << std::string( offset, ' ') << "######################################################"                                    << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "prefix_: "             << prefix_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "temperatureSurface_: " << temperatureSurface_  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "temperatureCMB_: "     << temperatureCMB_      << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 22 ) << std::left << "randomSeed_: "         << randomSeed_                 ;
      // clang-format on

      return os;
   }

 public:
   std::function< bool( const hyteg::Point3D& x ) > surfaceFct_;         // returns true if x is on the surface
   std::function< bool( const hyteg::Point3D& x ) > CMBFct_;             // returns true if x is on the cmb
   real_t                                           temperatureSurface_; // nondimensional temperature surface
   real_t                                           temperatureCMB_;     // nondimensional temperature CMB

   std::mt19937                             RndGenerator_;    // random generator
   std::uniform_real_distribution< real_t > RndDistribution_; // random distribution
   uint_fast32_t                            randomSeed_;      // random seed;
};

inline std::ostream& operator<<( std::ostream& os, const RandomTemperature& rt )
{
   return rt.print( os );
}

} // namespace MantleConvection