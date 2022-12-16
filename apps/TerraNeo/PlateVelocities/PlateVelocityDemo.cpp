/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/ThinShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "terraneo/plates/PlateVelocityProvider.hpp"

using walberla::int_c;
using walberla::real_c;
using walberla::real_t;
using namespace hyteg;
using namespace terraneo;

typedef enum
{
   PLATE_IDS,
   VELOCITIES,
   VELOCITIES_AND_IDS
} job_t;

std::string genOutputFileName( real_t age )
{
   std::stringstream sstr;
   sstr << "PlateVelocitiesDemo_age=";
   sstr << std::setfill( '0' );
   sstr << std::setw( 5 );
   sstr << std::fixed;
   sstr << std::setprecision( 1 );
   sstr << age;
   return sstr.str();
}

std::tuple< uint_t, uint_t > decodeRangeString( std::string& rangeStr )
{
   size_t fPos = rangeStr.find( '-' );
   if ( fPos == std::string::npos )
   {
      WALBERLA_ABORT( "Cannot decode range string '" << rangeStr << "'" );
   }

   uint_t stageInit = stoul( rangeStr.substr( 0, fPos ) );
   uint_t stageStop = stoul( rangeStr.substr( fPos + 1 ) );

   return std::make_tuple( stageInit, stageStop );
}

// =================================================================================
//  Function to test computation of a velocity field from plate reconstruction data
//  for all DoFs on the surface of a sphere
// =================================================================================
template < typename feFuncType >
void performComputations( uint_t                                     level,
                          std::shared_ptr< hyteg::PrimitiveStorage > storage,
                          real_t                                     age,
                          job_t                                      jobType,
                          terraneo::plates::PlateVelocityProvider&   oracle,
                          std::string                                indent )
{
   // need that here, to capture it below ;-)
   uint_t                                           coordIdx = 0;
   terraneo::plates::StatisticsPlateNotFoundHandler handlerWithStatistics;

   // callback function for computing the velocity components
   std::function< real_t( const Point3D& ) > computeVelocityComponent =
       [&oracle, age, &coordIdx, &handlerWithStatistics]( const Point3D& point ) {
          vec3D coords{ point[0], point[1], point[2] };
          vec3D velocity =
              oracle.getPointVelocity( coords, age, terraneo::plates::LinearDistanceSmoother{ 0.015 }, handlerWithStatistics );
          return velocity[int_c( coordIdx )];
       };

   // callback function for determining plate IDs
   std::function< real_t( const Point3D& ) > findPlateID = [&oracle, age, &handlerWithStatistics]( const Point3D& point ) {
      vec3D  coords{ point[0], point[1], point[2] };
      uint_t id = oracle.findPlateID( coords, age );
      if ( id == oracle.idWhenNoPlateFound )
      {
         handlerWithStatistics( coords, age );
      }
      return id;
   };

   // use pointers so that we only request memory for function in the desired jobType
   std::shared_ptr< feFuncType >                               surfaceVelocity{ nullptr };
   std::shared_ptr< typename feFuncType::VectorComponentType > plateID{ nullptr };

   // for DoF locations on the surface determine their associate plate IDs
   if ( jobType == PLATE_IDS || jobType == VELOCITIES_AND_IDS )
   {
      plateID = std::make_shared< typename feFuncType::VectorComponentType >( "plateID", storage, level, level );
      plateID->interpolate( findPlateID, level, All );
   }

   // for DoF locations on the surface determine their associate velocity vectors
   if ( jobType == VELOCITIES || jobType == VELOCITIES_AND_IDS )
   {
      surfaceVelocity = std::make_shared< feFuncType >( "plateVelocities", storage, level, level, 3 );
      for ( coordIdx = 0; coordIdx < 3; ++coordIdx )
      {
         ( *surfaceVelocity )[coordIdx].interpolate( computeVelocityComponent, level, All );
      }
   }

   // export results
   std::string fName = genOutputFileName( age );
   WALBERLA_LOG_INFO_ON_ROOT( "" << indent << "Exporting simulation data to file with basename '" << fName << "'" );
   hyteg::VTKOutput vtkOutput( "./output", fName, storage );
   switch ( jobType )
   {
   case PLATE_IDS:
      vtkOutput.add( *plateID );
      break;
   case VELOCITIES:
      vtkOutput.add( *surfaceVelocity );
      break;
   case VELOCITIES_AND_IDS:
      vtkOutput.add( *plateID );
      vtkOutput.add( *surfaceVelocity );
      break;
   }
   vtkOutput.write( level, 0 );

   handlerWithStatistics.generateReport();
}

// ========
//  Driver
// ========
int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   // ------------
   //  Parameters
   // ------------
   WALBERLA_LOG_INFO_ON_ROOT( "*** STEP 1: Obtaining Steering Parameters" );

   // check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./PlateVelocityDemo.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle params = cfg->getBlock( "Parameters" );

   if ( walberla::MPIManager::instance()->worldRank() == 0 )
   {
      WALBERLA_LOG_INFO( "Running with the following steering parameters:" );
      params.listParameters();
   }

   // determine job type
   job_t             jobType;
   const std::string str = params.getParameter< std::string >( "jobType" );
   if ( str == "plateIDs" )
   {
      jobType = PLATE_IDS;
   }
   else if ( str == "velocities" )
   {
      jobType = VELOCITIES;
   }
   else if ( str == "both" )
   {
      jobType = VELOCITIES_AND_IDS;
   }
   else
   {
      WALBERLA_ABORT( "Invalid jobType in parameter file!" );
   }

   // determine age indices
   std::string rangeStr        = params.getParameter< std::string >( "ageRange" );
   auto [stageInit, stageStop] = decodeRangeString( rangeStr );

   // ---------
   //  Meshing
   // ---------
   WALBERLA_LOG_INFO_ON_ROOT( "*** STEP 2: Generating Mesh" );
   const uint_t level = params.getParameter< uint_t >( "level" );
   const uint_t nTan  = params.getParameter< uint_t >( "nTan" );
   const real_t radius = real_c( 1 );

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::meshThinSphericalShell( nTan, radius );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );
   ThinShellMap::setMap( setupStorage, radius );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );

   // --------
   //  Oracle
   // --------
   WALBERLA_LOG_INFO_ON_ROOT( "*** STEP 3: Generating an Oracle" );
   std::string                             dataDir{ "../../../data/terraneo/plates/" };
   std::string                             fnameTopologies      = dataDir + "topologies0-100Ma.geojson";
   std::string                             fnameReconstructions = dataDir + "Global_EarthByte_230-0Ma_GK07_AREPS.rot";
   terraneo::plates::PlateVelocityProvider oracle( fnameTopologies, fnameReconstructions );

   WALBERLA_LOG_INFO_ON_ROOT( "*** STEP 4: Checking plate stages to work with" );
   auto stages = oracle.getListOfPlateStages();
   if ( stageStop >= stages.size() )
   {
      WALBERLA_ABORT( "ageRange parameter does not work for available plate stages!" );
   }
   WALBERLA_LOG_INFO_ON_ROOT( " - Running from stage[" << stageInit << "] = " << stages[stageInit] << " Ma back to stage["
                                                       << stageStop << "] = " << stages[stageStop] << " Ma" );

   // ------------
   //  Delegation
   // ------------
   WALBERLA_LOG_INFO_ON_ROOT( "*** STEP 5: Running the actual computations" );
   std::string feSpace = params.getParameter< std::string >( "feSpace" );

   for ( uint_t currentStage = stageInit; currentStage <= stageStop; ++currentStage )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " - stage[" << std::setw( 3 ) << currentStage << "] = " << std::fixed << std::setprecision( 1 )
                                             << stages[currentStage] << " Ma" );
      if ( feSpace == "P1" )
      {
         performComputations< P1VectorFunction< real_t > >( level, storage, stages[currentStage], jobType, oracle, "   " );
      }
      else if ( feSpace == "P2" )
      {
         performComputations< P2VectorFunction< real_t > >( level, storage, stages[currentStage], jobType, oracle, "   " );
      }
   }

   return EXIT_SUCCESS;
}
