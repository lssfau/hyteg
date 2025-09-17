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
#include "core/debug/CheckFunctions.h"

using namespace hyteg;

namespace terraneo {

/// Computes the number of shells for a spherical shell mesh.
/// On the corresponding blended micro-mesh, every node lies on exactly one shell.
///
/// Note that for P2 discretizations there are more shells, since there is one for every edge unknown along the radial direction, too.
/// This is treated by exploiting the overlap of edge dofs with new vertex dofs after refinement.
///
/// Layers have volume, shells do not.
///
/// \param nRad                           number of radial shells of the coarse mesh
/// \param level                          refinement level
/// \param polynomialOrderOfLagrangeDiscr either 1 or 2 - defines the degree of the Lagrangian polynomials (only P1 and P2
///                                       functions supported)
/// \return number of shells
inline uint_t numberOfShells( uint_t nRad, uint_t level, uint_t polynomialOrderOfLagrangeDiscr )
{
   WALBERLA_CHECK_GREATER_EQUAL( polynomialOrderOfLagrangeDiscr, 1 );
   WALBERLA_CHECK_LESS_EQUAL( polynomialOrderOfLagrangeDiscr, 2 );
   switch ( polynomialOrderOfLagrangeDiscr )
   {
   case 1:
      return ( nRad - 1 ) * ( levelinfo::num_microvertices_per_edge( level ) - 1 ) + 1;
   case 2:
      return numberOfShells( nRad, level + 1, 1 );
   default:
      return 0;
   }
}

/// Computes the number of layers for a spherical shell mesh.
/// Layers have volume, shells do not.
inline uint_t numberOfLayers( uint_t nRad, uint_t level, uint_t polynomialOrderOfLagrangeDiscr )
{
   return numberOfShells( nRad, level, polynomialOrderOfLagrangeDiscr ) - 1;
}

/// Computes the radius of a specific shell.
inline real_t
    radiusOfShell( uint_t shell, real_t rMin, real_t rMax, uint_t nRad, uint_t level, uint_t polynomialOrderOfLagrangeDiscr )
{
   return rMin + real_c( shell ) * ( rMax - rMin ) / real_c( numberOfLayers( nRad, level, polynomialOrderOfLagrangeDiscr ) );
}

/// Computes the radius of a specific shell of the macro mesh.
inline real_t radiusOfMacroMeshShell( uint_t layer, real_t rMin, real_t rMax, uint_t nRad )
{
   return rMin + real_c( layer ) * ( rMax - rMin ) / real_c( nRad - 1 );
}

std::vector< real_t >
    computeShellRadii( const std::vector< real_t >& layers, uint_t level, uint_t polynomialOrderOfLagrangeDiscr )
{
   uint_t effectiveLevel = level + polynomialOrderOfLagrangeDiscr - 1;

   uint_t                nRad      = layers.size();
   uint_t                numShells = numberOfShells( nRad, level, polynomialOrderOfLagrangeDiscr );
   std::vector< real_t > shellRadii( numShells, 0.0 );

   for ( int iShell = 0U; iShell < nRad - 1; iShell++ )
   {
      shellRadii[( 1 << effectiveLevel ) * iShell] = layers[iShell];
      //uint_t shell                                 = ( 1 << effectiveLevel ) * iShell;
      // WALBERLA_LOG_INFO_ON_ROOT("Shell Layer " << shell << " = " << profile.shellRadii[shell]);
      for ( int jShell = 0U; jShell < ( 1 << effectiveLevel ) - 1; jShell++ )
      {
         shellRadii[( 1 << effectiveLevel ) * iShell + jShell + 1] =
             layers[iShell] + ( jShell + 1 ) * ( layers[iShell + 1] - layers[iShell] ) / ( 1 << effectiveLevel );

         //shell = ( 1 << effectiveLevel ) * iShell + jShell + 1;
         // WALBERLA_LOG_INFO_ON_ROOT("Shell Layer " << shell << " = " << profile.shellRadii[shell]);
      }
   }

   shellRadii[numShells - 1] = layers[nRad - 1];

   return shellRadii;
}

/// Computes the index of the nearest shell from a given radius.
inline uint_t nearestShellFromRadius( real_t radius,
                                      real_t rMin,
                                      real_t rMax,
                                      uint_t nRad,
                                      uint_t level,
                                      uint_t polynomialOrderOfLagrangeDiscr )
{
   return static_cast< uint_t >( std::round( real_c( numberOfLayers( nRad, level, polynomialOrderOfLagrangeDiscr ) ) *
                                             ( ( radius - rMin ) / ( rMax - rMin ) ) ) );
}

/// Computes the index of the nearest shell from a given radius.
inline uint_t nearestShellFromRadius( real_t radius, const std::vector< real_t >& shellRadii )
{
   uint_t nearestShell    = 0u;
   real_t nearestDistance = __DBL_MAX__;
   for ( uint_t iShell = 0u; iShell < shellRadii.size(); iShell++ )
   {
      if ( std::abs( radius - shellRadii[iShell] ) < nearestDistance )
      {
         nearestShell    = iShell;
         nearestDistance = std::abs( radius - shellRadii[iShell] );
      }
   }

   return nearestShell;
}

/// Computes the index of the nearest shell from a given radius.
inline uint_t nearestShellFromRadius( real_t                       radius,
                                      const std::vector< real_t >& layers,
                                      uint_t                       level,
                                      uint_t                       polynomialOrderOfLagrangeDiscr )
{
   auto shellRadii = computeShellRadii( layers, level, polynomialOrderOfLagrangeDiscr );

   return nearestShellFromRadius( radius, shellRadii );
}

/// Interpolates the radial shell ID at the nodes of the passed scalar P1 or P2 function.
template < typename ScalarFunctionType >
inline void interpolateRadialShellID( ScalarFunctionType& u, real_t rMin, real_t rMax, uint_t nRad, uint_t level )
{
   using ShellIDType = typename ScalarFunctionType::valueType;

   WALBERLA_CHECK(
       std::is_integral_v< ShellIDType >,
       "You should write the shell IDs to a function that is typed with some kind of integer. For instance P1Function< int16_t >." )

   const bool isP1Function = std::is_same_v< typename ScalarFunctionType::Tag, P1FunctionTag >;
   const bool isP2Function = std::is_same_v< typename ScalarFunctionType::Tag, P2FunctionTag >;

   WALBERLA_CHECK( isP1Function || isP2Function, "interpolateRadialShellID() only supported for scalar P1 and P2 functions." );

   std::function< ShellIDType( const Point3D& ) > radialShellID = [&]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      ShellIDType shellID = static_cast< ShellIDType >(
          nearestShellFromRadius( radius, rMin, rMax, nRad, level, polynomialDegreeOfBasisFunctions< ScalarFunctionType >() ) );

      // Returning the value to ensure that the values are not altered.
      return shellID;
   };

   u.interpolate( radialShellID, level );
}

/// Interpolates the radial shell ID at the nodes of the passed scalar P1 or P2 function.
template < typename ScalarFunctionType >
inline void interpolateRadialShellID( ScalarFunctionType& u, std::vector< real_t > layers, uint_t level )
{
   using ShellIDType = typename ScalarFunctionType::valueType;

   WALBERLA_CHECK(
       std::is_integral_v< ShellIDType >,
       "You should write the shell IDs to a function that is typed with some kind of integer. For instance P1Function< int16_t >." )

   const bool isP1Function = std::is_same_v< typename ScalarFunctionType::Tag, P1FunctionTag >;
   const bool isP2Function = std::is_same_v< typename ScalarFunctionType::Tag, P2FunctionTag >;

   WALBERLA_CHECK( isP1Function || isP2Function, "interpolateRadialShellID() only supported for scalar P1 and P2 functions." );

   auto shellRadii = computeShellRadii( layers, level, polynomialDegreeOfBasisFunctions< ScalarFunctionType >() );

   std::function< ShellIDType( const Point3D& ) > radialShellID = [&]( const Point3D& x ) {
      real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      ShellIDType shellID = static_cast< ShellIDType >( nearestShellFromRadius( radius, shellRadii ) );

      // Returning the value to ensure that the values are not altered.
      return shellID;
   };

   u.interpolate( radialShellID, level );
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
/// This avoids storing the points twice. The components must then be signalled through the key, all component indices will be
/// zero.)
///
/// No mesh information is stored - you get shell-wise point clouds.
///
/// Note that this stores only process-local data. Communication has to be performed explicitly via respective functions.
template < typename FunctionType >
class RadialShellData
{
 public:
   void addDataFromFunction( const FunctionType& u, real_t rMin, real_t rMax, std::vector< real_t > layers, uint_t level )
   {
      WALBERLA_CHECK_LESS( rMin, rMax, "The thick spherical seems to be degenerate :/" );

      WALBERLA_CHECK( u.getStorage()->hasGlobalCells(),
                      "The radial shell data can only be gathered in 3D on the spherical shell." )

      WALBERLA_CHECK_EQUAL( values_.count( u.getFunctionName() ), 0, "There already is data stored for that function name." )

      uint_t nRad = layers.size();

      WALBERLA_CHECK_GREATER( nRad, 0, "No layers?" );

      auto arePointsInitialized = points_.size() > 0;

      const auto numLayers = numberOfLayers( nRad, level, polynomialDegreeOfBasisFunctions< FunctionType >() );
      const auto numShells = numberOfShells( nRad, level, polynomialDegreeOfBasisFunctions< FunctionType >() );

      uint_t numComponents = 1;

      if constexpr ( std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > ||
                     std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag > )
      {
         numComponents = 3;
      }
      else if constexpr ( !std::is_same_v< typename FunctionType::Tag, P1FunctionTag > &&
                          !std::is_same_v< typename FunctionType::Tag, P2FunctionTag > )
      {
         WALBERLA_ABORT( "Currently only PxFunctions and PxVectorFunctions for x in {1, 2} are supported." );
      }

      // To allocate the correct amount of memory required to store all process-local points and values on each shell, we need to
      // count them first. Computing that number analytically may be possible, but is certainly very tricky.

      std::vector< uint_t > numLocalPointsPerShell( numShells, 0 );

      auto shellRadii = computeShellRadii( layers, level, polynomialDegreeOfBasisFunctions< FunctionType >() );

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > countNodes =
          [&]( const Point3D& x, const std::vector< real_t >& values ) {
             real_t radius      = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
             real_t scalarValue = values[0];

             uint_t shell = nearestShellFromRadius( radius, shellRadii );

             // Manual bounds checking.
             WALBERLA_ASSERT_LESS( shell, numShells );

             numLocalPointsPerShell[shell]++;

             // Returning the value to ensure that the values are not altered.
             return scalarValue;
          };

      if constexpr ( std::is_same_v< typename FunctionType::Tag, P1FunctionTag > ||
                     std::is_same_v< typename FunctionType::Tag, P2FunctionTag > )
      {
         u.interpolate( countNodes, { u }, level, All );
      }
      else if constexpr ( std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > ||
                          std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag > )
      {
         u[0].interpolate( countNodes, { u[0] }, level, All );
      }
      else
      {
         WALBERLA_ABORT( "Radial data cannot be collected for the selected function type." );
      }

      if ( !arePointsInitialized )
      {
         points_.resize( numShells );
         for ( uint_t shell = 0; shell < numShells; shell++ )
         {
            points_[shell].reserve( numLocalPointsPerShell[shell] );
         }
      }

      // Initialize/resize arrays
      values_[u.getFunctionName()].resize( numComponents );
      for ( uint_t component = 0; component < numComponents; component++ )
      {
         values_[u.getFunctionName()][component].resize( numShells );
         for ( uint_t shell = 0; shell < numShells; shell++ )
         {
            values_[u.getFunctionName()][component][shell].reserve( numLocalPointsPerShell[shell] );
         }
      }

      // Interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with total temperature and
      // number of DoFs.

      for ( uint_t component = 0; component < numComponents; component++ )
      {
         std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherValues =
             [&]( const Point3D& x, const std::vector< real_t >& values ) {
                real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

                real_t scalarValue = values[0];

                uint_t shell = nearestShellFromRadius( radius, shellRadii );

                // Manual bounds checking.
                WALBERLA_ASSERT_LESS( shell, numShells );

                // Set point
                if ( !arePointsInitialized )
                {
                   // No worries about push_back(), enough space has been reserved previously.
                   points_[shell].push_back( x );
                }

                // Add value to corresponding shell in data array.
                // No worries about push_back(), enough space has been reserved previously.
                values_[u.getFunctionName()][component][shell].push_back( scalarValue );

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

         arePointsInitialized = true;
      }
   }

   void addDataFromFunction( const FunctionType& u, real_t rMin, real_t rMax, uint_t nRad, uint_t level )
   {
      std::vector< real_t > layers( nRad, 0.0 );
      for ( uint_t layer = 0; layer < nRad; layer++ )
      {
         layers[layer] = radiusOfMacroMeshShell( layer, rMin, rMax, nRad );
      }

      addDataFromFunction( u, rMin, rMax, layers, level );
   }

   const std::vector< Point3D >& points( uint_t shellId ) const { return points_.at( shellId ); }

   const std::vector< real_t >& values( std::string functionName, uint_t component, uint_t shellId ) const
   {
      WALBERLA_CHECK_GREATER( values_.count( functionName ), 0, "Key not registered in RadialShellData instance." );
      return values_.at( functionName ).at( component ).at( shellId );
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
   std::vector< real_t > rms;
   std::vector< uint_t > numDoFsPerShell;
   std::vector< real_t > depthDim;

   /// Simple logging of all vectors to file.
   void logToFile( std::string fileName, std::string fieldName )
   {
      WALBERLA_CHECK_EQUAL( shellRadii.size(), depthDim.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), mean.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), max.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), min.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), rms.size() );
      WALBERLA_CHECK_EQUAL( shellRadii.size(), numDoFsPerShell.size() );

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream outFile( fileName );

         if ( outFile.fail() )
         {
            WALBERLA_ABORT( "Failed to open file \"" << fileName << "\" for logging radial profiles. " )
         }

         // only write out rms for velocity
         if ( fieldName == "velocity" )
         {
            outFile << std::string( "#Depth [km] " ) << fieldName << std::string( "_mean " ) << fieldName
                    << std::string( "_max " ) << fieldName << std::string( "_min " ) << fieldName << std::string( "_rms \n" );

            for ( uint_t shell = shellRadii.size(); shell-- > 0; )
            {
               outFile << walberla::format( "%7.2f  %9.4f  %9.4f  %9.4f  %9.4f \n",
                                            depthDim.at( shell ),
                                            mean.at( shell ),
                                            max.at( shell ),
                                            min.at( shell ),
                                            rms.at( shell ) );
            }
         }
         else
         {
            outFile << std::string( "#Depth [km] " ) << fieldName << std::string( "_mean " ) << fieldName
                    << std::string( "_max " ) << fieldName << std::string( "_min \n" );

            for ( uint_t shell = shellRadii.size(); shell-- > 0; )
            {
               outFile << walberla::format(
                   "%7.2f  %9.4f  %9.4f  %9.4f \n", depthDim.at( shell ), mean.at( shell ), max.at( shell ), min.at( shell ) );
            }
         }
         outFile.close();
      }
   };
};

/// Computes radial profiles from a scalar or vector-valued FE function.
///
/// Assumes that the underlying domain is the spherical shell.
///
/// It is also possible to use on a thin spherical shell with rMin == rMax == layers[0].
/// For example to calcualte these parameters on the surface of the Earth (from Plates):
/// rMin = rMax = nRad = 1.0.
///
/// Iterating over the coefficients of the passed FE function, this function computes and returns the min, max, mean, and rms of
/// the coefficients (for scalar functions) or of the magnitude (for vector-valued functions) in radial layers.
///
/// Difference of mean and rms (root mean square): let n be the number of coefficients per shell and x_j the jth data point
/// (either a scalar value or the magnitude of a vector coefficient (i.e., x_j = sqrt( z[0]^2 + z[1]^2 + z[2]^2 ), if z is the
/// vector at that point)).
///
/// We compute:
///
///   mean( x ) = (1/n) * sum_j x_j,
///
/// and
///
///   rms( x ) = sqrt( (1/n) * sum_j (x_j)^2 ).
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
RadialProfile computeRadialProfile(
    const FunctionType&                       u,
    real_t                                    rMin,
    real_t                                    rMax,
    std::vector< real_t >                     layers,
    uint_t                                    level,
    std::function< real_t( const Point3D& ) > radiusFunc = []( const Point3D& x ) { return x.norm(); } )
{
   WALBERLA_CHECK_LESS_EQUAL( rMin, rMax );

   RadialProfile profile;

   uint_t nRad = layers.size();

   const auto numLayers = numberOfLayers( nRad, level, polynomialDegreeOfBasisFunctions< FunctionType >() );
   const auto numShells = numberOfShells( nRad, level, polynomialDegreeOfBasisFunctions< FunctionType >() );

   profile.shellRadii.resize( numShells );
   profile.min.resize( numShells, std::numeric_limits< real_t >::max() );
   profile.max.resize( numShells, std::numeric_limits< real_t >::lowest() );
   profile.mean.resize( numShells );
   profile.rms.resize( numShells );
   profile.numDoFsPerShell.resize( numShells );
   profile.depthDim.resize( numShells );

   profile.shellRadii = computeShellRadii( layers, level, polynomialDegreeOfBasisFunctions< FunctionType >() );

   // WALBERLA_LOG_INFO_ON_ROOT("Shell Layer " << numShells - 1 << " = " << profile.shellRadii[numShells - 1]);

   // Interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with total temperature and
   // number of DoFs.

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherValues =
       [&]( const Point3D& x, const std::vector< real_t >& values ) {
          real_t radius = radiusFunc( x ); // std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

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

          uint_t shell = nearestShellFromRadius( radius, profile.shellRadii );

          // Manual bounds checking.
          WALBERLA_ASSERT_LESS( shell, numShells );

          // Add value to corresponding row in profile (using the shell number)
          profile.min[shell] = std::min( profile.min[shell], scalarValue );
          profile.max[shell] = std::max( profile.max[shell], scalarValue );
          profile.mean[shell] += scalarValue;
          profile.rms[shell] += scalarValue * scalarValue;
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
   walberla::mpi::allReduceInplace( profile.rms, walberla::mpi::SUM );
   walberla::mpi::allReduceInplace( profile.numDoFsPerShell, walberla::mpi::SUM );

   // Now compute mean with total / counter
   for ( uint_t shell = 0; shell < numShells; ++shell )
   {
      profile.mean[shell] /= real_c( profile.numDoFsPerShell[shell] );
      profile.rms[shell] = std::sqrt( profile.rms[shell] / real_c( profile.numDoFsPerShell[shell] ) );
   }

   return profile;
}

template < typename FunctionType >
RadialProfile computeRadialProfile(
    const FunctionType&                       u,
    real_t                                    rMin,
    real_t                                    rMax,
    uint_t                                    nRad,
    uint_t                                    level,
    std::function< real_t( const Point3D& ) > radiusFunc = []( const Point3D& x ) { return x.norm(); } )
{
   std::vector< real_t > layers( nRad, 0.0 );
   for ( uint_t layer = 0; layer < nRad; layer++ )
   {
      layers[layer] = radiusOfMacroMeshShell( layer, rMin, rMax, nRad );
   }

   return computeRadialProfile( u, rMin, rMax, layers, level, radiusFunc );
}

} //namespace terraneo
