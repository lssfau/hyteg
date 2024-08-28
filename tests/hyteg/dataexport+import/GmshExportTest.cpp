/*
 * Copyright (c) 2024 Marcus Mohr.
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
#include <cstdlib>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/dataexport/Gmsh/ExportRefinedMesh.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/ThinShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// This test checks that exporting refined meshes to MSH4.1 files
// for Gmsh compiles and runs through without aborting.
//
// The test can optionally check that the meshfiles are valid by
// processing them with Gmsh (supply --verify_with_gmsh when invoking)

using namespace hyteg;

std::string outputDirectory{ "./GmshExportTest-Output" };

std::string decorateBaseFileName( const std::string& stringIn, uint_t level )
{
   std::stringstream stringStreamOut;
   stringStreamOut << outputDirectory << "/";
   stringStreamOut << stringIn;
   stringStreamOut << "-level=" << level;
   stringStreamOut << ".msh";
   return stringStreamOut.str();
}

void runTest( const SetupPrimitiveStorage& setupStorage, std::string baseFileName, uint_t level, bool verify_with_gmsh )
{
   std::string fileName = decorateBaseFileName( baseFileName, level );

   WALBERLA_LOG_INFO_ON_ROOT( "-> running with level = " << level );
   WALBERLA_LOG_INFO_ON_ROOT( "-> output will go to file '" << fileName << "'" );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   gmsh::exportRefinedMesh( storage, level, fileName );

   // feed meshfile through Gmsh to check its validity
   if( verify_with_gmsh ) {
     std::string command = "gmsh -parse_and_exit " + fileName + " > /dev/null";
     WALBERLA_LOG_INFO_ON_ROOT( "-> Checking validity with '" << command << "'" );

     int returnCode = std::system( command.c_str() );
     WALBERLA_CHECK_EQUAL( returnCode, 0 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "done" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   // a little hacky, but we do not want to work with *.prm file and overwritting settings
   int                   dummy{ 1 };
   walberla::Environment walberlaEnv( dummy, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   bool verify_with_gmsh = false;

   if ( argc != 1 && argc != 2 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" << argv[0] << " received " << argc << " command-line arguments:" );
      for ( uint_t i = 1; i <= argc; ++i )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "arg #" << i << " = '" << argv[i] << "'" );
         WALBERLA_ABORT( "Expecting either none or '--verify_with_gmsh'" );
      }
   }
   else if( argc == 2 )
   {
      if ( strcmp( argv[1], "--verify_with_gmsh" ) == 0 )
      {
         verify_with_gmsh = true;
      }
      else
      {
         WALBERLA_ABORT( "Do not understand command-line argument '" << argv[1] << "'" );
      }
   }

   // 2D mesh w/o blending
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing Mesh Export for LShape_6el.msh ***" );
   {
      MeshInfo              mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/LShape_6el.msh" ) );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      runTest( setupStorage, "LShape_6el", 4, verify_with_gmsh );
   }

   // 2D mesh w/ blending
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing Mesh Export for a blended annulus ***" );
   {
      MeshInfo              mesh = MeshInfo::meshAnnulus( real_c( 0.75 ), real_c( 1.5 ), MeshInfo::CROSS, 24, 2 );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      AnnulusMap::setMap( setupStorage );
      runTest( setupStorage, "BlendedAnnulus", 4, verify_with_gmsh );
   }

   // 3D mesh w/o blending
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing Mesh Export for cube_6el_offcenter ***" );
   {
      MeshInfo              mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/cube_6el_offcenter.msh" ) );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      runTest( setupStorage, "cube_6el_offcenter", 3, verify_with_gmsh );
   }

   // 3D mesh w/ blending
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing Mesh Export for Thick Spherical Shell ***" );
   {
      MeshInfo              mesh = MeshInfo::meshSphericalShell( 2, 2, real_c( 1 ), real_c( 2 ) );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      IcosahedralShellMap::setMap( setupStorage );
      runTest( setupStorage, "ThickSphericalShell", 3, verify_with_gmsh );
   }

   // 3D surface mesh
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing Mesh Export for Thin Spherical Shell ***" );
   {
      real_t                thinShellRadius{ real_c( 2.5 ) };
      MeshInfo              mesh = MeshInfo::meshThinSphericalShell( 2, thinShellRadius );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      ThinShellMap::setMap( setupStorage, thinShellRadius );
      runTest( setupStorage, "ThinSphericalShell", 5, verify_with_gmsh );
   }

   return EXIT_SUCCESS;
}
