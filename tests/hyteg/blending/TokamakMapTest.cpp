/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/geometry/TokamakMap.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t minLevel = 2;
   const uint_t maxLevel = 3;
   const bool   writeVTK = true;

   // ITER configuration

   const uint_t                numToroidalSlices          = 8;
   const uint_t                numPoloidalSlices          = 6;
   const real_t                radiusOriginToCenterOfTube = 6.2;
   const std::vector< real_t > tubeLayerRadii             = { 3 };
   const real_t                torodialStartAngle         = 0.0;
   const real_t                polodialStartAngle         = 2.0 * pi / real_c( 2 * numPoloidalSlices );

   real_t delta = sin( 0.33 );
   real_t r1    = 2.0;
   real_t r2    = 3.7;

   // set to true to map an actual torus, not the tokamak
   const bool mapTorus = false;

   if ( mapTorus )
   {
      delta = 0;
      r1    = tubeLayerRadii.back();
      r2    = tubeLayerRadii.back();
   }

   const auto meshInfo = MeshInfo::meshTorus(
       numToroidalSlices, numPoloidalSlices, radiusOriginToCenterOfTube, tubeLayerRadii, torodialStartAngle, polodialStartAngle );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   TokamakMap::setMap( setupStorage,
                       numToroidalSlices,
                       numPoloidalSlices,
                       radiusOriginToCenterOfTube,
                       tubeLayerRadii,
                       torodialStartAngle,
                       polodialStartAngle,
                       delta,
                       r1,
                       r2 );

   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK( storage->hasGlobalCells() );

   writeDomainPartitioningVTK( storage, "../../output", "TokamakMapTest" );

   std::function< real_t( const Point3D& ) > exact = []( const Point3D& p ) -> real_t {
      return sin( p[0] ) * sinh( p[1] ) * p[2];
   };

   P1Function< real_t > u( "u", storage, minLevel, maxLevel );

   u.interpolate( exact, maxLevel, DoFType::All );

   VTKOutput vtkOutput( "../../output", "TokamakMapTest", storage );
   vtkOutput.add( u );

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   return 0;
}
