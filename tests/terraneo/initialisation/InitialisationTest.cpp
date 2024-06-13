/*
 * Copyright (c) 2024 Eugenio D'Ascoli, Marcus Mohr.
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

#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "terraneo/helpers/RadialProfiles.hpp"
#include "terraneo/initialisation/TemperatureInitialisation.hpp"

using namespace hyteg;

/// Setup storage for a spherical shell
std::shared_ptr< PrimitiveStorage >
    setupSphericalShellStorage( const uint_t& nTan, const uint_t& nRad, const real_t& rMax, const real_t& rMin )
{
   std::shared_ptr< hyteg::MeshInfo > meshInfo;
   meshInfo = std::make_shared< hyteg::MeshInfo >( hyteg::MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax ) );

   hyteg::SetupPrimitiveStorage setupStorage( *meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   IcosahedralShellMap::setMap( setupStorage );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );
   return storage;
}

/// Setup storage for a spherical shell
std::shared_ptr< PrimitiveStorage >
    setupSphericalShellStorage( const uint_t& nTan, std::vector<real_t> layers )
{
   std::shared_ptr< hyteg::MeshInfo > meshInfo;
   meshInfo = std::make_shared< hyteg::MeshInfo >( hyteg::MeshInfo::meshSphericalShell( nTan, layers ) );

   hyteg::SetupPrimitiveStorage setupStorage( *meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   IcosahedralShellMap::setMap( setupStorage );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );
   return storage;
}

void checkProfile( const terraneo::RadialProfile& profile )
{
   WALBERLA_CHECK_EQUAL( profile.shellRadii.size(), profile.min.size() )
   WALBERLA_CHECK_EQUAL( profile.shellRadii.size(), profile.max.size() )
   WALBERLA_CHECK_EQUAL( profile.shellRadii.size(), profile.mean.size() )
   WALBERLA_CHECK_EQUAL( profile.shellRadii.size(), profile.numDoFsPerShell.size() )

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8s | %8s | %8s | %8s | %8s ", "radius", "min", "mean", "max", "count" ) );
   WALBERLA_LOG_INFO_ON_ROOT( " ---------+----------+----------+----------+---------- " );
   for ( uint_t i = 0; i < profile.shellRadii.size(); i++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8f | %8f | %8f | %8f | %8u ",
                                                   profile.shellRadii[i],
                                                   profile.min[i],
                                                   profile.mean[i],
                                                   profile.max[i],
                                                   profile.numDoFsPerShell[i] ) )

      const auto eps = std::abs( profile.mean[i] ) * real_c( 1e-8 );
      WALBERLA_CHECK_LESS_EQUAL( profile.min[i], profile.mean[i] + eps );
      WALBERLA_CHECK_LESS_EQUAL( profile.mean[i], profile.max[i] + eps );
   }
}

/// Tests multiple pieces of code from the terraneo module:
/// - Uses some of the temperature initialization functions.
/// - Computes radial profiles.
/// - Checks whether the radial profiles give reasonable numbers.
template < typename FunctionType, typename VectorFunctionType >
void runTest( const uint_t nTan,
              const uint_t nRad,
              const real_t rMax,
              const real_t rMin,
              const uint_t level,
              const real_t epsMeanTest )
{
   walberla::math::seedRandomGenerator( 42 );

   auto storage = setupSphericalShellStorage( nTan, nRad, rMax, rMin );

   FunctionType temperature( "temperature", storage, level, level );
   FunctionType temperatureDev( "temperatureDev", storage, level, level );

   VectorFunctionType velocity( "velocity", storage, level, level );

   real_t Tsurface = 300;
   real_t Tcmb     = 4200;
   real_t Tadb     = 1600;

   terraneo::TemperatureInitializationParameters tempInitParams( Tcmb, Tsurface, Tadb, real_c( 0.68 ), rMin, rMax );

   real_t noiseFactor                 = real_c( 0.05 );
   uint_t tempInit                    = 10;
   uint_t deg                         = 4;
   int    ord                         = 2;
   uint_t lmax                        = 25;
   uint_t lmin                        = 10;
   bool   superposition               = true;
   real_t buoyancyFactor              = real_c( 0.01 );
   real_t initialTemperatureSteepness = 10;

   auto temperatureReference = terraneo::temperatureReferenceExponential( tempInitParams );

   auto initTemperatureWhiteNoise = terraneo::temperatureWhiteNoise( tempInitParams, temperatureReference, noiseFactor );

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating ref + white noise" )
   temperature.interpolate( initTemperatureWhiteNoise, level );

   // This is pretty slow. Just here to be noticed by the compiler.
   if ( false )
   {
      auto initTemperatureSPH = terraneo::temperatureSPH( tempInitParams,
                                                          temperatureReference,
                                                          tempInit,
                                                          deg,
                                                          ord,
                                                          lmax,
                                                          lmin,
                                                          superposition,
                                                          buoyancyFactor,
                                                          initialTemperatureSteepness );

      WALBERLA_LOG_INFO_ON_ROOT( "Interpolating ref + SPH" )
      temperature.interpolate( initTemperatureSPH, level );
   }

   // Just interpolating something for testing.
   velocity.interpolate( initTemperatureWhiteNoise, level );

   std::string outputDirectory = "./output";
   std::string baseName        = "Test_T_field_Init";

   if ( !std::filesystem::exists( outputDirectory ) )
   {
      std::filesystem::create_directories( outputDirectory );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Computing profile - scalar function" )
   auto profile = terraneo::computeRadialProfile( temperature, rMin, rMax, nRad, level );
   checkProfile( profile );
   WALBERLA_LOG_INFO_ON_ROOT( "Computing profile - vector function" )
   auto profileVector = terraneo::computeRadialProfile( velocity, rMin, rMax, nRad, level );
   checkProfile( profileVector );

   // Just logging to check if that runs through.
   profile.logToFile( outputDirectory + "/" + baseName + ".txt", "temperature" );

   // The following code subtracts the mean from every shell.
   // This way we can test if the mean after that is zero on every shell.

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > temperatureDevFct =
       [&]( const Point3D& x, const std::vector< real_t >& T ) {
          auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
          auto shell  = terraneo::nearestShellFromRadius(
              radius, rMin, rMax, nRad, level, polynomialDegreeOfBasisFunctions< FunctionType >() );
          return T[0] - profile.mean[shell];
       };

   temperatureDev.interpolate( temperatureDevFct, { temperature }, level, All );

   auto profileTest = terraneo::computeRadialProfile( temperatureDev, rMin, rMax, nRad, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Test mean after subtraction" )
   for ( auto p : profileTest.mean )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" << std::showpos << std::scientific << p )
      WALBERLA_CHECK_LESS( std::abs( p ), epsMeanTest );
   }

   bool output = false;

   if ( output )
   {
      auto vtkOutput = hyteg::VTKOutput( outputDirectory, baseName, storage );
      vtkOutput.setVTKDataFormat( hyteg::vtk::DataFormat::BINARY );
      vtkOutput.add( temperature );
      vtkOutput.add( temperatureDev );
      vtkOutput.write( level );
   }
}

template < typename FunctionType, typename VectorFunctionType >
void runTest( const uint_t& nTan, std::vector<real_t> layers, const uint_t& level )
{
   walberla::math::seedRandomGenerator( 42 );

   real_t rMin = layers[0];
   real_t rMax = layers[layers.size() - 1];

   auto storage = setupSphericalShellStorage( nTan, layers );

   FunctionType temperature( "temperature", storage, level, level );
   FunctionType temperatureDev( "temperatureDev", storage, level, level );

   VectorFunctionType velocity( "velocity", storage, level, level );

   real_t Tsurface = 300;
   real_t Tcmb     = 4200;
   real_t Tadb     = 1600;

   terraneo::TemperatureInitializationParameters tempInitParams( Tcmb, Tsurface, Tadb, 0.68, rMin, rMax );

   real_t noiseFactor                 = 0.05;
   uint_t tempInit                    = 10;
   uint_t deg                         = 4;
   int    ord                         = 2;
   uint_t lmax                        = 25;
   uint_t lmin                        = 10;
   bool   superposition               = true;
   real_t buoyancyFactor              = 0.01;
   real_t initialTemperatureSteepness = 10;

   auto temperatureReference = terraneo::temperatureReferenceExponential( tempInitParams );

   auto initTemperatureWhiteNoise = terraneo::temperatureWhiteNoise( tempInitParams, temperatureReference, noiseFactor );

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating ref + white noise" )
   temperature.interpolate( initTemperatureWhiteNoise, level );

   // This is pretty slow. Just here to be noticed by the compiler.
   if ( false )
   {
      auto initTemperatureSPH = terraneo::temperatureSPH( tempInitParams,
                                                          temperatureReference,
                                                          tempInit,
                                                          deg,
                                                          ord,
                                                          lmax,
                                                          lmin,
                                                          superposition,
                                                          buoyancyFactor,
                                                          initialTemperatureSteepness );

      WALBERLA_LOG_INFO_ON_ROOT( "Interpolating ref + SPH" )
      temperature.interpolate( initTemperatureSPH, level );
   }

   // Just interpolating something for testing.
   velocity.interpolate( initTemperatureWhiteNoise, level );

   std::string outputDirectory = "./output";
   std::string baseName        = "Test_T_field_Init";

   if ( !std::filesystem::exists( outputDirectory ) )
   {
      std::filesystem::create_directories( outputDirectory );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Computing profile - scalar function" )
   auto profile = terraneo::computeRadialProfile( temperature, rMin, rMax, layers, level );
   checkProfile( profile );
   WALBERLA_LOG_INFO_ON_ROOT( "Computing profile - vector function" )
   auto profileVector = terraneo::computeRadialProfile( velocity, rMin, rMax, layers, level );
   checkProfile( profileVector );

   // Just logging to check if that runs through.
   profile.logToFile( outputDirectory + "/" + baseName + ".txt", "temperature" );

   // The following code subtracts the mean from every shell.
   // This way we can test if the mean after that is zero on every shell.

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > temperatureDevFct =
       [&]( const Point3D& x, const std::vector< real_t >& T ) {
          auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
          auto shell = terraneo::nearestShellFromRadius(
              radius, layers, level, polynomialDegreeOfBasisFunctions< FunctionType >() );
          return T[0] - profile.mean[shell];
       };

   temperatureDev.interpolate( temperatureDevFct, { temperature }, level, All );

   auto profileTest = terraneo::computeRadialProfile( temperatureDev, rMin, rMax, layers, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Test mean after subtraction" )
   for ( auto p : profileTest.mean )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" << std::showpos << std::scientific << p )
      WALBERLA_CHECK_LESS( std::abs( p ), epsMeanTest );
   }

   bool output = false;

   if ( output )
   {
      auto vtkOutput = hyteg::VTKOutput( outputDirectory, baseName, storage );
      vtkOutput.setVTKDataFormat( hyteg::vtk::DataFormat::BINARY );
      vtkOutput.add( temperature );
      vtkOutput.add( temperatureDev );
      vtkOutput.write( level );
   }
}

template < typename FunctionType, typename VectorFunctionType >
void runTest( const uint_t& nTan, std::vector<real_t> layers, const uint_t& level )
{
   walberla::math::seedRandomGenerator( 42 );

   real_t rMin = layers[0];
   real_t rMax = layers[layers.size() - 1];

   auto storage = setupSphericalShellStorage( nTan, layers );

   FunctionType temperature( "temperature", storage, level, level );
   FunctionType temperatureDev( "temperatureDev", storage, level, level );

   VectorFunctionType velocity( "velocity", storage, level, level );

   real_t Tsurface = 300;
   real_t Tcmb     = 4200;
   real_t Tadb     = 1600;

   terraneo::TemperatureInitializationParameters tempInitParams( Tcmb, Tsurface, Tadb, 0.68, rMin, rMax );

   real_t noiseFactor                 = 0.05;
   uint_t tempInit                    = 10;
   uint_t deg                         = 4;
   int    ord                         = 2;
   uint_t lmax                        = 25;
   uint_t lmin                        = 10;
   bool   superposition               = true;
   real_t buoyancyFactor              = 0.01;
   real_t initialTemperatureSteepness = 10;

   auto temperatureReference = terraneo::temperatureReferenceExponential( tempInitParams );

   auto initTemperatureWhiteNoise = terraneo::temperatureWhiteNoise( tempInitParams, temperatureReference, noiseFactor );

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating ref + white noise" )
   temperature.interpolate( initTemperatureWhiteNoise, level );

   // This is pretty slow. Just here to be noticed by the compiler.
   if ( false )
   {
      auto initTemperatureSPH = terraneo::temperatureSPH( tempInitParams,
                                                          temperatureReference,
                                                          tempInit,
                                                          deg,
                                                          ord,
                                                          lmax,
                                                          lmin,
                                                          superposition,
                                                          buoyancyFactor,
                                                          initialTemperatureSteepness );

      WALBERLA_LOG_INFO_ON_ROOT( "Interpolating ref + SPH" )
      temperature.interpolate( initTemperatureSPH, level );
   }

   // Just interpolating something for testing.
   velocity.interpolate( initTemperatureWhiteNoise, level );

   std::string outputDirectory = "./output";
   std::string baseName        = "Test_T_field_Init";

   if ( !std::filesystem::exists( outputDirectory ) )
   {
      std::filesystem::create_directories( outputDirectory );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Computing profile - scalar function" )
   auto profile = terraneo::computeRadialProfile( temperature, rMin, rMax, layers, level );
   checkProfile( profile );
   WALBERLA_LOG_INFO_ON_ROOT( "Computing profile - vector function" )
   auto profileVector = terraneo::computeRadialProfile( velocity, rMin, rMax, layers, level );
   checkProfile( profileVector );

   // Just logging to check if that runs through.
   profile.logToFile( outputDirectory + "/" + baseName + ".txt", "temperature" );

   // The following code subtracts the mean from every shell.
   // This way we can test if the mean after that is zero on every shell.

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > temperatureDevFct =
       [&]( const Point3D& x, const std::vector< real_t >& T ) {
          auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
          auto shell = terraneo::nearestShellFromRadius(
              radius, layers, level, polynomialDegreeOfBasisFunctions< FunctionType >() );
          return T[0] - profile.mean[shell];
       };

   temperatureDev.interpolate( temperatureDevFct, { temperature }, level, All );

   auto profileTest = terraneo::computeRadialProfile( temperatureDev, rMin, rMax, layers, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Test mean after subtraction" )
   for ( auto p : profileTest.mean )
   {
      WALBERLA_LOG_INFO_ON_ROOT( p )
      WALBERLA_CHECK_FLOAT_EQUAL( p, real_c( 0 ) );
   }

   bool output = true;

   if ( output )
   {
      auto vtkOutput = hyteg::VTKOutput( outputDirectory, baseName, storage );
      vtkOutput.setVTKDataFormat( hyteg::vtk::DataFormat::BINARY );
      vtkOutput.add( temperature );
      vtkOutput.add( temperatureDev );
      vtkOutput.write( level );
   }
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   real_t rMax  = real_c( 2.12 );
   real_t rMin  = real_c( 1.12 );
   uint_t level = 3;
   uint_t nTan  = 3;
   uint_t nRad  = 2;

   WALBERLA_LOG_INFO_ON_ROOT( "**************************************" );
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing with P1 type functions ***" );
   WALBERLA_LOG_INFO_ON_ROOT( "**************************************" );
   {
      real_t epsMeanTest = real_c( std::is_same_v< double, real_t > ? 1e-13 : 3e-5 );
      runTest< P1Function< real_t >, P1VectorFunction< real_t > >( nTan, nRad, rMax, rMin, level, epsMeanTest );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "**************************************" );
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing with P2 type functions ***" );
   WALBERLA_LOG_INFO_ON_ROOT( "**************************************" );
   {
      real_t epsMeanTest = real_c( std::is_same_v< double, real_t > ? 1e-13 : 2e-4 );
      runTest< P2Function< real_t >, P2VectorFunction< real_t > >( nTan, nRad, rMax, rMin, level, epsMeanTest );
   }

   runTest< P1Function< real_t >, P1VectorFunction< real_t > >( nTan, {rMin, rMin + (rMax - rMin) / 5, rMax}, level );
   runTest< P2Function< real_t >, P2VectorFunction< real_t > >( nTan, {rMin, rMin + (rMax - rMin) / 5, rMax}, level );
   
   return 0;
}
