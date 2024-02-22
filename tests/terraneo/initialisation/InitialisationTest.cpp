/*
 * Copyright (c) 2024 Eugenio D'Ascoli.
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

#include "terraneo/helpers/InitialisationTool.hpp"
#include "terraneo/helpers/RadialProfileTool.hpp"

using namespace hyteg;

// Setup storage for a spherical shell

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

// Test function to evaluate temperature initialisation using white noise (Gaussian White Noise (GWN)) or spherical
// harmonics superimposed on the background temperature.
// The noise factor will define how much temprature deviation [%] will be added as noise to the refernce temperature.
// For the spherical harmonics temperature intialisation a distinct degree and order of the spherical harmonics functions
// can be defined.

template < typename FunctionType >
void runTest( const uint_t& nTan,
              const uint_t& nRad,
              const real_t& rMax,
              const real_t& rMin,
              const uint_t& maxLevel,
              const uint_t& minLevel )
{
   auto storage = setupSphericalShellStorage( nTan, nRad, rMax, rMin );

   std::shared_ptr< FunctionType > temperature = std::make_shared< FunctionType >( "temperature", storage, minLevel, maxLevel );
   std::shared_ptr< FunctionType > temperatureDev =
       std::make_shared< FunctionType >( "temperatureDev", storage, minLevel, maxLevel );
   std::shared_ptr< FunctionType > tmp    = std::make_shared< FunctionType >( "tmp", storage, minLevel, maxLevel );
   std::shared_ptr< FunctionType > tmpDev = std::make_shared< FunctionType >( "tmpDev", storage, minLevel, maxLevel );

   real_t Tsurface = 300;
   real_t Tcmb     = 4200;
   real_t Tadb     = 1600;

   std::shared_ptr< terraneo::TemperaturefieldConv< FunctionType > > Temperaturefield =
       std::make_shared< terraneo::TemperaturefieldConv< FunctionType > >(
           temperature, Tcmb, Tsurface, Tadb, 0.68, rMax, rMin, maxLevel, minLevel );

   real_t noiseFactor                 = 0.05;
   uint_t tempInit                    = 10;
   uint_t deg                         = 4;
   int    ord                         = 2;
   uint_t lmax                        = 25;
   uint_t lmin                        = 10;
   bool   superposition               = true;
   real_t buoyancyFactor              = 0.01;
   real_t initialTemperatureSteepness = 10;
   bool   noiseInit                   = true;

   if ( noiseInit )
   {
      Temperaturefield->initialiseTemperatureWhiteNoise( noiseFactor );
   }
   else
   {
      Temperaturefield->initialiseTemperatureSPH(
          tempInit, deg, ord, lmax, lmin, superposition, buoyancyFactor, initialTemperatureSteepness );
   }

   std::string outputDirectory = "./output";
   std::string baseName        = "Test_T_field_Init";

   if ( !std::filesystem::exists( outputDirectory ) )
   {
      std::filesystem::create_directories( outputDirectory );
   }

   std::shared_ptr< terraneo::RadialProfileTool< FunctionType > > TemperatureProfileTool =
       std::make_shared< terraneo::RadialProfileTool< FunctionType > >( *temperature, *tmp, rMax, rMin, nRad, maxLevel );
   std::vector< real_t > TemperatureProfile = TemperatureProfileTool->getMeanProfile();
   auto                  numLayers          = 2 * ( nRad - 1 ) * ( levelinfo::num_microvertices_per_edge( maxLevel ) - 1 );

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > temperatureDevFct =
       [&]( const Point3D& x, const std::vector< real_t >& T ) {
          auto   radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
          real_t retVal;

          uint_t shell = static_cast< uint_t >( std::round( real_c( numLayers ) * ( ( radius - rMin ) / ( rMax - rMin ) ) ) );
          WALBERLA_ASSERT( shell < T.size() );

          retVal = ( T[0] - TemperatureProfile.at( shell ) );
          return retVal;
       };

   for ( uint_t l = minLevel; l <= maxLevel; ++l )
   {
      temperatureDev->interpolate( temperatureDevFct, { *temperature }, l, All );
   }

   std::shared_ptr< terraneo::RadialProfileTool< FunctionType > > TemperatureDevProfileTool =
       std::make_shared< terraneo::RadialProfileTool< FunctionType > >( *temperatureDev, *tmpDev, rMin, rMax, nRad, maxLevel );
   std::vector< real_t > TemperatureDevProfile = TemperatureDevProfileTool->getMeanProfile();

   // Evaluate that the mean of temperature deviation is zero.

   bool isZero = std::all_of( TemperatureDevProfile.begin(), TemperatureDevProfile.end(), []( uint_t i ) { return i == 0; } );
   WALBERLA_CHECK( isZero == true );

   bool output = true;

   if ( output )
   {
      auto vtkOutput = hyteg::VTKOutput( outputDirectory, baseName, storage );
      vtkOutput.setVTKDataFormat( hyteg::vtk::DataFormat::BINARY );
      vtkOutput.add( *temperature );
      vtkOutput.add( *temperatureDev );
      vtkOutput.write( maxLevel );
   }
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   real_t rMax     = 2.12;
   real_t rMin     = 1.12;
   uint_t minLevel = 0;
   uint_t maxLevel = 4;
   uint_t nTan     = 3;
   uint_t nRad     = 2;

   runTest< P2Function< real_t > >( nTan, nRad, rMax, rMin, maxLevel, minLevel );
   return 0;
}