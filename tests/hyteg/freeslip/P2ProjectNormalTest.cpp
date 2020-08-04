/*
 * Copyright (c) 2020 Daniel Drzisga.
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
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_c;
using walberla::real_t;

using namespace hyteg;

static void testProjectNormal2D( )
{
   const bool   writeVTK   = true;
   const real_t errorLimit = 1e-13;
   const int level = 3;

   const auto  meshInfo = MeshInfo::meshAnnulus(0.5, 1.0, MeshInfo::CRISS, 6, 6);
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 3, 0, true );
   AnnulusMap::setMap( setupStorage );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   if ( writeVTK )
      writeDomainPartitioningVTK( storage, "../../output", "P2ProjectNormalTest2D_Domain" );

   auto normal_function = []( const Point3D& p, Point3D& n ) -> void {
     real_t norm = p.norm();
     real_t sign = (norm > 0.75) ? 1.0 : -1.0;

     n = sign/norm * p;
   };

   P2ProjectNormalOperator projectNormalOperator( storage, level, level, normal_function );

   P2P1TaylorHoodFunction< real_t > u( "u", storage, level, level );

   VTKOutput vtkOutput( "../../output", "P2ProjectNormalTest2D", storage );
   vtkOutput.add( u );

   u.interpolate( 1, level );
   projectNormalOperator.apply( u, level, FreeslipBoundary );

   if ( writeVTK )
      vtkOutput.write( level, 0 );
}

static void testProjectNormal3D( )
{
   const bool   writeVTK   = true;
   const real_t errorLimit = 1e-13;
   const int level = 2;

   auto meshInfo = MeshInfo::meshSphericalShell( 5, 2, 0.5, 1.0 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 3, 0, true );
   IcosahedralShellMap::setMap( setupStorage );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   if ( writeVTK )
      writeDomainPartitioningVTK( storage, "../../output", "P2ProjectNormalTest3D_Domain" );

   auto normal_function = []( const Point3D& p, Point3D& n ) -> void {
     real_t norm = p.norm();
     real_t sign = (norm > 0.75) ? 1.0 : -1.0;

     n = sign/norm * p;
   };

   P2ProjectNormalOperator projectNormalOperator( storage, level, level, normal_function );

   P2P1TaylorHoodFunction< real_t > u( "u", storage, level, level );

   VTKOutput vtkOutput( "../../output", "P2ProjectNormalTest3D", storage );
   vtkOutput.add( u );

   u.interpolate( 1, level );
   projectNormalOperator.apply( u, level, FreeslipBoundary );

   if ( writeVTK )
      vtkOutput.write( level, 0 );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   testProjectNormal2D( );
   testProjectNormal3D( );

   return 0;
}
