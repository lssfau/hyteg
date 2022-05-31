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
#include "hyteg/geometry/IcosahedralShellMap.hpp"
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

// =================================================================================
//  Function to test computation of a velocity field from plate reconstruction data
//  for all DoFs on the surface of a sphere
// =================================================================================
template < typename feFuncType >
void exportDemoFunc( uint_t level, std::shared_ptr< hyteg::PrimitiveStorage > storage, real_t age )
{
   // initialise an oracle
   std::string                             dataDir{ "../../../data/terraneo/plates/" };
   std::string                             fnameTopologies      = dataDir + "topologies0-100Ma.geojson";
   std::string                             fnameReconstructions = dataDir + "Global_EarthByte_230-0Ma_GK07_AREPS.rot";
   terraneo::plates::PlateVelocityProvider oracle( fnameTopologies, fnameReconstructions );

   // need that here, to capture it below ;-)
   uint_t coordIdx = 0;

   std::function< real_t( const Point3D& ) > computeVelocityComponent = [&oracle, age, &coordIdx]( const Point3D& point ) {
      vec3D coords{ point[0], point[1], point[2] };
      vec3D velocity = oracle.getPointVelocity( coords, age );
      return velocity[ int_c(coordIdx) ];
   };

   std::function< real_t( const Point3D& ) > findPlateID = [&oracle, age]( const Point3D& point ) {
      vec3D coords{ point[0], point[1], point[2] };
      return oracle.findPlateID( coords, age );
   };

   // set everything to zero in the interior
   feFuncType demoFunc( "demoFunc", storage, level, level );
   demoFunc.interpolate( { real_c( 0 ), real_c( 0 ), real_c( 0 ) }, level, Inner );

   typename feFuncType::VectorComponentType plates( "plateID", storage, level, level );
   plates.interpolate( real_c( 0 ), level, Inner );

   // now set plate velocities on boundary (both for the moment)
   for ( coordIdx = 0; coordIdx < 3; ++coordIdx )
   {
      // demoFunc[coordIdx].interpolate( computeVelocityComponent, level, Boundary );
   }
   plates.interpolate( findPlateID, level, Boundary );

   hyteg::VTKOutput vtkOutput( "./output", "PlateVelocities", storage );
   vtkOutput.add( plates );
   // vtkOutput.add( demoFunc );
   vtkOutput.write( level, 0 );
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

   // ---------
   //  Meshing
   // ---------
   const uint_t level = params.getParameter< uint_t >( "level" );
   const uint_t nRad  = params.getParameter< uint_t >( "nRad" );
   const uint_t nTan  = params.getParameter< uint_t >( "nTan" );

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, real_c(0.57), real_c(1.0) );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   IcosahedralShellMap::setMap( setupStorage );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   // ============
   //  Delegation
   // ============
   std::string feSpace = params.getParameter< std::string >( "feSpace" );
   const real_t age  = params.getParameter< real_t >( "age" );

   if ( feSpace == "P1" )
   {
      exportDemoFunc< P1VectorFunction< real_t > >( level, storage, age );
   }
   else if ( feSpace == "P2" )
   {
      exportDemoFunc< P2VectorFunction< real_t > >( level, storage, age );
   }

   return EXIT_SUCCESS;
}
