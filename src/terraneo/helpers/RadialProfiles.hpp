/*
 * Copyright (c) 2023-2024 Hamish Brown, Nils Kohl.
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

/// This file contains functions and classes to operate on radial layers and shells of a thick spherical shell.
///
/// Includes features to gather data for exporting or to reduce data over shells/layers.

#pragma once

#include <cmath>
#include <vector>

#include "core/DataTypes.h"
#include "core/Format.hpp"

using namespace hyteg;

namespace terraneo {

/// Computes the number of shells for a spherical shell mesh.
/// On the corresponding blended micro-mesh, every node lies on exactly one shell.
/// Note that this really corresponds to the VERTICES of the mesh.
/// Layers have volume, shells do not.
uint_t numberOfShells( uint_t nRad, uint_t level )
{
   return ( nRad - 1 ) * ( levelinfo::num_microvertices_per_edge( level ) - 1 ) + 1;
}

/// Computes the number of layers for a spherical shell mesh.
/// Layers have volume, shells do not.
uint_t numberOfLayers( uint_t nRad, uint_t level )
{
   return numberOfShells( nRad, level ) - 1;
}

/// Computes the radius of a specific shell.
real_t radiusOfShell( uint_t shell, real_t rMin, real_t rMax, uint_t nRad, uint_t level )
{
   return rMin + real_c( shell ) * ( rMax - rMin ) / real_c( numberOfLayers( nRad, level ) );
}

/// Computes the index of the nearest shell from a given radius.
uint_t nearestShellFromRadius( real_t radius, real_t rMin, real_t rMax, uint_t nRad, uint_t level )
{
   return static_cast< uint_t >(
       std::round( real_c( numberOfLayers( nRad, level ) ) * ( ( radius - rMin ) / ( rMax - rMin ) ) ) );
}

/// Simple struct to store data organized by radial shells.
///
/// Can (and should) be used to store data for multiple functions if the positions are the same. This way, the positions are only
/// stored once.
///
/// So if you have multiple P2 functions for instance, you should add their data to the same instance of this struct.
///
/// For vector functions, the component is signalled as part of the value data structure.
/// Scalar functions are stored with component set to 0.
///
/// (To be really space-efficient, you can store Px functions and components of PxVectorFunctions in the same instance.
/// This avoids storing the points twice. The components must then be signalled through the key, all component indices with be
/// zero.)
///
/// No mesh information is stored - you get shell-wise point clouds.
///
/// Note that this stores only process-local data. Communication has to be performed explicitly via respective functions.
template < typename FunctionType >
class RadialShellData
{
 public:
   void addDataFromFunction( std::string key, const FunctionType& u, real_t rMin, real_t rMax, uint_t nRad, uint_t level )
   {
      WALBERLA_CHECK_LESS_EQUAL( rMin, rMax );

      WALBERLA_CHECK( u.getStorage()->hasGlobalCells(),
                      "The radial shell data can only be gathered in 3D on the spherical shell." )

      WALBERLA_CHECK_EQUAL( values_.count( key ), 0, "There already is data stored for the passed key." )

      WALBERLA_CHECK_GREATER( nRad, 0, "No layers?" );

      auto arePointsInitialized = points_.size() > 0;

      const auto numLayers = numberOfLayers( nRad, level );
      const auto numShells = numberOfShells( nRad, level );

      uint_t numComponents = 1;

      if constexpr ( std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > ||
                     std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag > )
      {
         numComponents = 3;
      }

      if ( !arePointsInitialized )
      {
         points_.resize( numShells );
      }

      // Initialize/resize arrays
      values_[key].resize( numComponents );
      for ( uint_t component = 0; component < numComponents; component++ )
      {
         values_[key][component].resize( numShells );
      }
      // TODO: Can we easily precompute the number of DoFs per shell to resize the vectors? Might be hard in parallel.

      // Interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with total temperature and
      // number of DoFs.

      for ( uint_t component = 0; component < numComponents; component++ )
      {
         std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherValues =
             [&]( const Point3D& x, const std::vector< real_t >& values ) {
                real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

                real_t scalarValue = values[0];

                uint_t shell = nearestShellFromRadius( radius, rMin, rMax, nRad, level );

                // Manual bounds checking.
                WALBERLA_ASSERT_LESS( shell, numShells );

                // Set point
                if ( !arePointsInitialized )
                {
                   points_[shell].push_back( x );
                }

                // Add value to corresponding shell in data array.
                values_[key][component][shell].push_back( scalarValue );

                // Returning the value to ensure that the values are not altered.
                return scalarValue;
             };

         if constexpr ( std::is_same_v< typename FunctionType::Tag, P1FunctionTag > ||
                        std::is_same_v< typename FunctionType::Tag, P2FunctionTag > )
         {
            u.interpolate( gatherValues, { u }, level, All );
         }
         else if constexpr ( std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > ||
                             std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag > )
         {
            u[component].interpolate( gatherValues, { u[component] }, level, All );
         }
         else
         {
            WALBERLA_ABORT( "Radial data cannot be collected for the selected function type." );
         }

         auto arePointsInitialized = true;
      }
   }

   const std::vector< Point3D >& points( uint_t shellId ) const { return points_.at( shellId ); }

   const std::vector< real_t >& values( std::string key, uint_t component, uint_t shellId ) const
   {
      return values_.at( key ).at( component ).at( shellId );
   }

 private:
   /// Stores the point positions for each shell (local to the process).
   /// Layout: positions[shell][pointID] = pos3D
   std::vector< std::vector< Point3D > > points_;

   /// Stores scalar data for each point. Multiple data arrays can be added.
   /// Layout: values["my_function"][component][shellID][pointID] = scalar
   std::map< std::string, std::vector< std::vector< std::vector< real_t > > > > values_;
};

/// Just a simple struct holding std::vectors for convenient computation of and access to radial profiles.
struct RadialProfile
{
   std::vector< real_t > shellRadii;
   std::vector< real_t > mean;
   std::vector< real_t > max;
   std::vector< real_t > min;
   std::vector< uint_t > numDoFsPerShell;

   /// Simple logging of all vectors to file.
   void logToFile( std::string fileName, std::string fieldName )
   {
      WALBERLA_CHECK_EQUAL( shellRadii.size(), mean.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), max.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), min.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), numDoFsPerShell.size() );

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
      }
   };
};

/// Computes radial profiles from a scalar or vector-valued FE function.
///
/// Assumes that the underlying domain is the spherical shell.
///
/// Iterating over the coefficients of the passed FE function, this function computes and returns the min, max, and mean of the
/// coefficients (for scalar functions) or of the magnitude (for vector-valued functions) in radial layers.
///
/// Involves global communication!
///
///! Note: Currently only implemented for P1Functions, P2Functions, P1VectorFunctions, and P2VectorFunctions - but can easily be
/// extended.
///
/// \tparam FunctionType       type of a FE function (e.g. P2Function< real_t > or P2VectorFunction< real_t >)
/// \param u                   FE function that is evaluated
/// \param rMin                radius of innermost shell
/// \param rMax                radius of outermost shell
/// \param nRad                number of radial layers
/// \param level               FE function refinement level
/// \return a filled RadialProfile struct
template < typename FunctionType >
RadialProfile computeRadialProfile( const FunctionType& u, real_t rMin, real_t rMax, uint_t nRad, uint_t level )
{
   WALBERLA_CHECK_LESS_EQUAL( rMin, rMax );

   WALBERLA_CHECK( u.getStorage()->hasGlobalCells(), "The radial profile can only be computed in 3D on the spherical shell." )

   RadialProfile profile;

   const auto numLayers = numberOfLayers( nRad, level );
   const auto numShells = numberOfShells( nRad, level );

   profile.shellRadii.resize( numShells );
   profile.min.resize( numShells, std::numeric_limits< real_t >::max() );
   profile.max.resize( numShells, std::numeric_limits< real_t >::lowest() );
   profile.mean.resize( numShells );
   profile.numDoFsPerShell.resize( numShells );

   for ( uint_t shell = 0; shell < numShells; ++shell )
   {
      profile.shellRadii[shell] = radiusOfShell( shell, rMin, rMax, nRad, level );
   }

   // Interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with total temperature and
   // number of DoFs.

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherValues =
       [&]( const Point3D& x, const std::vector< real_t >& values ) {
          real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

          auto scalarValue = real_c( 0.0 );

          if constexpr ( std::is_same_v< typename FunctionType::Tag, P1FunctionTag > ||
                         std::is_same_v< typename FunctionType::Tag, P2FunctionTag > )
          {
             scalarValue = values[0];
          }
          else if constexpr ( std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > ||
                              std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag > )
          {
             scalarValue = std::sqrt( values[0] * values[0] + values[1] * values[1] + values[2] * values[2] );
          }
          else
          {
             WALBERLA_ABORT( "Radial profile cannot be computed for the selected function type." );
          }

          uint_t shell = nearestShellFromRadius( radius, rMin, rMax, nRad, level );

          // Manual bounds checking.
          WALBERLA_ASSERT_LESS( shell, numShells );

          // Add value to corresponding row in profile (using the shell number)
          profile.min[shell] = std::min( profile.min[shell], scalarValue );
          profile.max[shell] = std::max( profile.max[shell], scalarValue );
          profile.mean[shell] += scalarValue;
          profile.numDoFsPerShell[shell] += 1;

          // Returning the value of the first function to ensure that the values are not altered.
          // This should be called on the first component of a CSFVectorFunction until we have a better implementation that
          // does not abuse interpolate().
          return values[0];
       };

   if constexpr ( std::is_same_v< typename FunctionType::Tag, P1FunctionTag > ||
                  std::is_same_v< typename FunctionType::Tag, P2FunctionTag > )
   {
      u.interpolate( gatherValues, { u }, level, All );
   }
   else if constexpr ( std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > ||
                       std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag > )
   {
      u[0].interpolate( gatherValues, { u[0], u[1], u[2] }, level, All );
   }
   else
   {
      WALBERLA_ABORT( "Radial profile cannot be computed for the selected function type." );
   }

   // Reduce values on each shell over all processes
   walberla::mpi::allReduceInplace( profile.min, walberla::mpi::MIN );
   walberla::mpi::allReduceInplace( profile.max, walberla::mpi::MAX );
   walberla::mpi::allReduceInplace( profile.mean, walberla::mpi::SUM );
   walberla::mpi::allReduceInplace( profile.numDoFsPerShell, walberla::mpi::SUM );

   // Now compute mean with total / counter
   for ( uint_t shell = 0; shell < numShells; ++shell )
   {
      profile.mean[shell] /= real_c( profile.numDoFsPerShell[shell] );
   }

   return profile;
}

} //namespace terraneo