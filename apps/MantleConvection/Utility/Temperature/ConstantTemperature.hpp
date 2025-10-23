/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "TemperatureModel.hpp"

using walberla::real_c;
using walberla::real_t;

namespace MantleConvection {

class ConstantTemperature : public TemperatureModel< real_t >
{
 public:
   using TemperatureModel< real_t >::prefix_;

   ConstantTemperature( NondimensionalisationParameters&                 nondimensionalisation,
                        walberla::Config::BlockHandle&                   parameters,
                        std::function< bool( const hyteg::Point3D& x ) > surfaceFct,
                        std::function< bool( const hyteg::Point3D& x ) > CMBFct,
                        real_t                                           scaling = real_c( 1.0 ),
                        std::string                                      prefix  = "" )
   : TemperatureModel< real_t >( prefix )
   , surfaceFct_( surfaceFct )
   , CMBFct_( CMBFct )
   , scaling_( scaling )
   {
      WALBERLA_UNUSED(parameters);
      temperatureSurface_ = nondimensionalisation.temperatureSurface_;
      temperatureCMB_     = nondimensionalisation.temperatureCMB_;
   }

   real_t evaluateNoCheck( const hyteg::Point3D& x ) override
   {
      WALBERLA_UNUSED( x );
      return temperatureSurface_ + scaling_;
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

   bool checkSurface( const hyteg::Point3D& x ) override { return surfaceFct_( x ); }
   bool checkCMB( const hyteg::Point3D& x ) override { return CMBFct_( x ); }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "######################################################"             << "\n";
      os << std::string( offset, ' ') << "################ Constant Temperature ################"             << "\n";
      os << std::string( offset, ' ') << "######################################################"             << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "prefix_: "  << prefix_  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 11 ) << std::left << "scaling_: " << scaling_        ;
      // clang-format on

      return os;
   }

 public:
   std::function< bool( const hyteg::Point3D& x ) > surfaceFct_;         // returns true if x is on the surface
   std::function< bool( const hyteg::Point3D& x ) > CMBFct_;             // returns true if x is on the cmb
   real_t                                           temperatureSurface_; // nondimensional temperature surface
   real_t                                           temperatureCMB_;     // nondimensional temperature CMB
   real_t                                           scaling_;            // constant scaling
};

inline std::ostream& operator<<( std::ostream& os, const ConstantTemperature& ct )
{
   return ct.print( os );
}

} // namespace MantleConvection