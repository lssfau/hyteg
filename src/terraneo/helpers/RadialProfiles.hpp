/*
 * Copyright (c) 2023 Hamish Brown, Nils Kohl.
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
#include <vector>

#include "core/DataTypes.h"

using namespace hyteg;

namespace terraneo {

/// Just a simple struct holding std::vectors for convenient computation of and access to radial profiles.
struct ScalarRadialProfile
{
   std::vector< real_t > shellRadii;
   std::vector< real_t > mean;
   std::vector< real_t > max;
   std::vector< real_t > min;
   std::vector< uint_t > count;

   /// Simple logging of all vectors to file.
   void logToFile( std::string fileName, std::string fieldName )
   {
      WALBERLA_CHECK_EQUAL( shellRadii.size(), mean.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), max.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), min.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), count.size() );

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream outFile( fileName );

         if ( outFile.fail() )
         {
            WALBERLA_ABORT( "Failed to open file \"" << fileName << "\" for logging radial profiles. " )
         }

         outFile << std::string( "Radius  " ) << fieldName << std::string( "_Mean " ) << fieldName << std::string( "_Max " )
                 << fieldName << std::string( "_Min \n" );
         for ( uint_t shell = 0; shell < shellRadii.size(); ++shell )
         {
            outFile << walberla::format(
                "%6.4f  %7.4f  %6.4f  %6.4f \n", shellRadii.at( shell ), mean.at( shell ), max.at( shell ), min.at( shell ) );
         }
         outFile.close();

         if ( outFile.is_open() )
         {
            WALBERLA_ABORT( "Failed to close file \"" << fileName << "\" when logging radial profiles. " )
         }
      }
   }; // namespace terraneo
};

/// Computes radial profiles from a scalar FE function.
///
/// Assumes that the underlying domain is the spherical shell.
///
/// Iterating over the coefficients of the passed FE function, this function computes and returns the min, max, and mean of the
/// coefficients in radial layers.
///
/// Involves global communication!
///
///! Note: Currently only implemented for P2 Scalar functions on spherical meshes with evenly distributed radial layers
///
/// \tparam ScalarFunctionType type of a scalar FE function (e.g. P2Function< real_t >)
/// \param u                   scalar FE function that is evaluated
/// \param rMin                radius of innermost shell
/// \param rMax                radius of outermost shell
/// \param nRad                number of radial layers
/// \param level               FE function refinement level
/// \return a filled ScalarRadialProfile struct
template < typename ScalarFunctionType >
ScalarRadialProfile computeScalarRadialProfile( const ScalarFunctionType& u, real_t rMin, real_t rMax, uint_t nRad, uint_t level )
{
   WALBERLA_CHECK_LESS_EQUAL( rMin, rMax );

   ScalarRadialProfile profile;

   const auto numberOfLayers = 2 * ( nRad - 1 ) * ( levelinfo::num_microvertices_per_edge( level ) - 1 );

   auto getRadius = [&]( const uint_t& shell ) -> real_t {
      return rMin + real_c( shell ) * ( rMax - rMin ) / real_c( numberOfLayers );
   };

   auto getShell = [&]( const real_t& radius ) -> uint_t {
      return static_cast< uint_t >( std::round( real_c( numberOfLayers ) * ( ( radius - rMin ) / ( rMax - rMin ) ) ) );
   };

   for ( uint_t shell = 0; shell < numberOfLayers + 1; ++shell )
   {
      profile.shellRadii.push_back( getRadius( shell ) );
   }

   profile.min.resize( numberOfLayers + 1 );
   profile.max.resize( numberOfLayers + 1 );
   profile.mean.resize( numberOfLayers + 1 );
   profile.count.resize( numberOfLayers + 1 );

   std::fill( profile.min.begin(), profile.min.end(), std::numeric_limits< real_t >::max() );
   std::fill( profile.max.begin(), profile.max.end(), std::numeric_limits< real_t >::lowest() );
   std::fill( profile.mean.begin(), profile.mean.end(), real_c( 0 ) );
   std::fill( profile.count.begin(), profile.count.end(), uint_c( 0 ) );

   //interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with total temperature and number of DoFs

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherValues =
       [&]( const Point3D& x, const std::vector< real_t >& values ) {
          WALBERLA_ASSERT_EQUAL( values.size(), 1, "Why pass more than one function?" );
          real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

          uint_t shell = getShell( radius );

          // Manual bounds checking.
          WALBERLA_ASSERT_LESS( shell, numberOfLayers + 1 );

          // Add value to corresponding row in profile (using shell number)
          profile.min[shell] = std::min( profile.min[shell], values[0] );
          profile.max[shell] = std::max( profile.max[shell], values[0] );
          profile.mean[shell] += values[0];
          profile.count[shell] += 1;

          return values[0];
       };

   u.interpolate( gatherValues, { u }, level, All );

   // Reduce values on each shell over all processes
   walberla::mpi::allReduceInplace( profile.min, walberla::mpi::MIN );
   walberla::mpi::allReduceInplace( profile.max, walberla::mpi::MAX );
   walberla::mpi::allReduceInplace( profile.mean, walberla::mpi::SUM );
   walberla::mpi::allReduceInplace( profile.count, walberla::mpi::SUM );

   // Now compute mean with total / counter
   for ( uint_t shell = 0; shell < numberOfLayers + 1; ++shell )
   {
      profile.mean[shell] /= real_c( profile.count[shell] );
   }

   return profile;
}

} //namespace terraneo