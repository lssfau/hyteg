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

using walberla::int_c;
using walberla::real_t;

namespace MantleConvection {

class RelativeRandomTemperature : public TemperatureModel< real_t >
{
 public:
   using TemperatureModel< real_t >::prefix_;

   RelativeRandomTemperature( NondimensionalisationParameters&                     nondimensionalisation,
                              walberla::Config::BlockHandle&                       parameters,
                              const std::shared_ptr< TemperatureModel< real_t > >& ModifiedModel,
                              const uint_fast32_t                                  seed   = 394,
                              std::string                                          prefix = "" )
   : TemperatureModel< real_t >( prefix )
   , ModifiedModel_( ModifiedModel )
   , RndGenerator_( seed )
   , RndDistribution_( real_c( -1 ), real_c( 1 ) )
   , randomSeed_( seed )
   {
      temperatureSurface_ = nondimensionalisation.temperatureSurface_;
      temperatureCMB_     = nondimensionalisation.temperatureCMB_;

      relativeTemperatureNoiseFactor_ =
          parameters.getParameter< real_t >( prefix + std::string( "relativeTemperatureNoiseFactor" ) );
   }

   real_t evaluateNoCheck( const hyteg::Point3D& x ) override
   {
      real_t retVal = ModifiedModel_->evaluateNoCheck( x );
      retVal += relativeTemperatureNoiseFactor_ * retVal * RndDistribution_( RndGenerator_ );

      // bound the result
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
      os << std::string( offset, ' ') << "#######################################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "############# Relative Random Temperature #############"                                                                 << "\n";
      os << std::string( offset, ' ') << "#######################################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 34 ) << std::left << "prefix_: "                          << prefix_                               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 34 ) << std::left << "temperatureSurface_: "              << temperatureSurface_                   << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 34 ) << std::left << "temperatureCMB_: "                  << temperatureCMB_                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 34 ) << std::left << "relativeTemperatureNoiseFactor_: "  << relativeTemperatureNoiseFactor_       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 34 ) << std::left << "randomSeed_: "                      << randomSeed_                           << "\n";
      // clang-format on

      return ModifiedModel_->print( os, offset + 3 );
   }

 public:
   std::shared_ptr< TemperatureModel< real_t > > ModifiedModel_;

   real_t temperatureSurface_; // nondimensional temperature surface
   real_t temperatureCMB_;     // nondimensional temperature CMB

   real_t       relativeTemperatureNoiseFactor_;              // relative noise factor for the random temperature initialisation
   std::mt19937 RndGenerator_;                                // random generator
   std::uniform_real_distribution< real_t > RndDistribution_; // random distribution
   uint_fast32_t                            randomSeed_;      // random seed;
};

inline std::ostream& operator<<( std::ostream& os, const RelativeRandomTemperature& rrt )
{
   return rrt.print( os );
}

} // namespace MantleConvection