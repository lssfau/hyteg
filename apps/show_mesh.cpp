/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include <sstream>

#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/ThinShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

void showUsage()
{
   std::cout
       << "\n --------\n  USAGE:\n --------\n\n"
       << " show_mesh demonstrates generation of a MeshInfo object by one of the\n following methods:\n\n"
       << "  1) by importing a file in Gmsh format\n"
       << "  2) by meshing a rectangle in a certain flavour\n"
       << "  3) by meshing a full or partial annulus\n"
       << "  4) by generating a strip of chained triangles\n"
       << "  5) by meshing a thick spherical shell (w/ or w/o blending)\n"
       << "  6) by meshing a thin spherical shell (w/ or w/o blending)\n"
       << "  7) by meshing a rectangular cuboid\n"
       << "  8) by meshing a symmetric rectangular cuboid\n"
       << "  9) by meshing a T-domain using the cubed domain generator\n"
       << " 10) by meshing a torus\n\n"
       << " This is steered by choosing one of the options below:\n\n"
       << "  --file <name of Gmsh file>\n"
       << "  --rect [criss|cross|crisscross|diamond]\n"
       << "  --annulus [full|partial]\n"
       << "  --face-chain [numFaces]\n"
       << "  --spherical-shell [ntan]\n"
       << "  --blended-spherical-shell [ntan]\n"
       << "  --thin-spherical-shell [ntan]\n"
       << "  --blended-thin-spherical-shell [ntan]\n"
       << "  --cuboid [nHint]\n"
       << "  --symm-cuboid [nSubCubes]\n"
       << "  --t-domain [nCubesInEachDirection]\n"
       << "  --torus\n\n"
       << " The generated base mesh will be tested be doing two levels of refinement.\n"
       << " Then it will be exported to a VTU file for visualisation.\n\n"
       << " Also visualization of the domain partitioning, mesh boundary flags and MPI rank assignment will be output.\n"
       << " The desired load balancing approach can be steered by:\n\n"
       << "  --load-balancing [allroot|roundrobin|roundrobinvolume|greedy|parmetis|diffusivecluster] (default: roundrobin)\n\n"
       << " Use  -v  to print info on mesh to console.\n\n"
       << " Use  --dofs <lvl>  to print number of DoFs on refinement level lvl.\n\n"
       << std::endl;
}

int main( int argc, char* argv[] )
{
   std::string meshFileName;
   std::string vtkFileName;
   typedef enum
   {
      FROM_FILE,
      RECTANGLE,
      ANNULUS,
      PARTIAL_ANNULUS,
      FACE_CHAIN,
      SPHERICAL_SHELL,
      THIN_SPHERICAL_SHELL,
      CUBOID,
      SYMM_CUBOID,
      T_DOMAIN,
      TORUS
   } meshDomainType;

   meshDomainType        meshDomain;
   MeshInfo::meshFlavour rectMeshType = MeshInfo::CROSS;
   MeshInfo*             meshInfo     = nullptr;
   bool                  beVerbose    = false;
   bool                  blending     = false;
   uint_t                ntan         = 5;
   uint_t                nHint        = 1;
   // std::vector< real_t > layers       = {0.5, 0.6, 0.7, 0.8};
   std::vector< real_t > layers          = { 1.0, 2.0 };
   uint_t                numFaces        = 2;
   real_t                thinShellRadius = real_t( 0.5 );

   typedef enum
   {
      ALL_ROOT,
      ROUND_ROBIN,
      ROUND_ROBIN_VOLUME,
      GREEDY,
      PARMETIS,
      DIFFUSIVE_CLUSTER
   } LoadBalancingType;
   LoadBalancingType loadBalancingType = ROUND_ROBIN;

   if ( argc < 3 || argc > 8 )
   {
      showUsage();
      WALBERLA_ABORT( "Please provide command-line parameters!" );
   }
   else if ( strcmp( argv[1], "--file" ) == 0 )
   {
      meshDomain   = FROM_FILE;
      meshFileName = std::string( argv[2] );
      auto pos     = meshFileName.find_last_of( '/' );
      pos == meshFileName.npos ? pos = 0 : ++pos;
      vtkFileName = meshFileName.substr( pos, meshFileName.length() - 4 );
   }
   else if ( strcmp( argv[1], "--rect" ) == 0 )
   {
      meshDomain = RECTANGLE;
      if ( strcmp( argv[2], "criss" ) == 0 )
      {
         rectMeshType = MeshInfo::CRISS;
         vtkFileName  = std::string( "rectMeshCriss" );
      }
      else if ( strcmp( argv[2], "cross" ) == 0 )
      {
         rectMeshType = MeshInfo::CROSS;
         vtkFileName  = std::string( "rectMeshCross" );
      }
      else if ( strcmp( argv[2], "crisscross" ) == 0 )
      {
         rectMeshType = MeshInfo::CRISSCROSS;
         vtkFileName  = std::string( "rectMeshCrissCross" );
      }
      else if ( strcmp( argv[2], "diamond" ) == 0 )
      {
         rectMeshType = MeshInfo::DIAMOND;
         vtkFileName  = std::string( "rectMeshDiamond" );
      }
      else
      {
         WALBERLA_ABORT( "Flavour for rect mesh not recognised!" );
      }
   }
   else if ( strcmp( argv[1], "--annulus" ) == 0 )
   {
      if ( strcmp( argv[2], "full" ) == 0 )
      {
         meshDomain  = ANNULUS;
         vtkFileName = std::string( "annulusMesh" );
      }
      else if ( strcmp( argv[2], "partial" ) == 0 )
      {
         meshDomain  = PARTIAL_ANNULUS;
         vtkFileName = std::string( "partialAnnulusMesh" );
      }
      else
      {
         WALBERLA_ABORT( "Subtype of --annulus not recognised!" );
      }
   }
   else if ( strcmp( argv[1], "--face-chain" ) == 0 )
   {
      numFaces    = uint_c( std::stoi( argv[2] ) );
      meshDomain  = FACE_CHAIN;
      vtkFileName = std::string( "faceChain" );
   }
   else if ( strcmp( argv[1], "--spherical-shell" ) == 0 )
   {
      ntan        = uint_c( std::stoi( argv[2] ) );
      meshDomain  = SPHERICAL_SHELL;
      blending    = false;
      vtkFileName = std::string( "sphericalShell" );
   }
   else if ( strcmp( argv[1], "--blended-spherical-shell" ) == 0 )
   {
      ntan        = uint_c( std::stoi( argv[2] ) );
      meshDomain  = SPHERICAL_SHELL;
      blending    = true;
      vtkFileName = std::string( "blendedSphericalShell" );
   }
   else if ( strcmp( argv[1], "--thin-spherical-shell" ) == 0 )
   {
      ntan        = uint_c( std::stoi( argv[2] ) );
      meshDomain  = THIN_SPHERICAL_SHELL;
      blending    = false;
      vtkFileName = std::string( "thinSphericalShell" );
   }
   else if ( strcmp( argv[1], "--blended-thin-spherical-shell" ) == 0 )
   {
      ntan        = uint_c( std::stoi( argv[2] ) );
      meshDomain  = THIN_SPHERICAL_SHELL;
      blending    = true;
      vtkFileName = std::string( "blendedThinSphericalShell" );
   }
   else if ( strcmp( argv[1], "--cuboid" ) == 0 )
   {
      nHint       = uint_c( std::stoi( argv[2] ) );
      meshDomain  = CUBOID;
      vtkFileName = std::string( "cuboidMesh" );
   }
   else if ( strcmp( argv[1], "--symm-cuboid" ) == 0 )
   {
      nHint       = uint_c( std::stoi( argv[2] ) );
      meshDomain  = SYMM_CUBOID;
      vtkFileName = std::string( "symmCuboidMesh" );
   }
   else if ( strcmp( argv[1], "--t-domain" ) == 0 )
   {
      nHint       = uint_c( std::stoi( argv[2] ) );
      meshDomain  = T_DOMAIN;
      vtkFileName = std::string( "tDomain" );
   }
   else if ( strcmp( argv[1], "--torus" ) == 0 )
   {
      meshDomain  = TORUS;
      vtkFileName = std::string( "torus" );
   }
   else
   {
      WALBERLA_ABORT( "Could not understand command-line args!" );
   }

   if ( argc > 4 && argc < 7 && strcmp( argv[3], "--load-balancing" ) == 0 )
   {
      if ( strcmp( argv[4], "allroot" ) == 0 )
      {
         loadBalancingType = ALL_ROOT;
      }
      else if ( strcmp( argv[4], "roundrobin" ) == 0 )
      {
         loadBalancingType = ROUND_ROBIN;
      }
      else if ( strcmp( argv[4], "roundrobinvolume" ) == 0 )
      {
         loadBalancingType = ROUND_ROBIN_VOLUME;
      }
      else if ( strcmp( argv[4], "greedy" ) == 0 )
      {
         loadBalancingType = GREEDY;
      }
      else if ( strcmp( argv[4], "diffusivecluster" ) == 0 )
      {
         loadBalancingType = DIFFUSIVE_CLUSTER;
      }
      else if ( strcmp( argv[4], "parmetis" ) == 0 )
      {
#ifdef WALBERLA_BUILD_WITH_PARMETIS
         loadBalancingType = PARMETIS;
#else
         WALBERLA_ABORT( "Framework was not built with ParMetis." );
#endif
      }
      else
      {
         WALBERLA_ABORT( "Could not understand command-line args. Possibly invalid load balancing approach." );
      }
   }

   if ( ( argc == 4 || argc == 6 ) && ( strcmp( argv[3], "-v" ) == 0 || strcmp( argv[5], "-v" ) == 0 ) )
   {
      beVerbose = true;
   }

   bool   reportDoFCount = false;
   uint_t dofLevel       = 0;
   for ( int k = 0; k < argc - 1; k++ )
   {
      if ( strcmp( argv[k], "--dofs" ) == 0 )
      {
         reportDoFCount = true;
         dofLevel       = (uint_t) atoi( argv[k + 1] );
      }
   }
   if ( strcmp( argv[argc - 1], "--dofs" ) == 0 )
   {
      WALBERLA_ABORT( "You need to give a level after '--dofs'" );
   }

   // ------------
   //  Let's rock
   // ------------
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "HyTeG Show Mesh Test\n" );

   switch ( meshDomain )
   {
   case FROM_FILE:
      meshInfo = new MeshInfo( MeshInfo::fromGmshFile( meshFileName ) );
      break;

   case RECTANGLE:
      meshInfo = new MeshInfo( MeshInfo::meshRectangle( Point2D( -2.0, 1.0 ), Point2D( 0.0, 3.0 ), rectMeshType, 16, 16 ) );
      break;

   case PARTIAL_ANNULUS:
      meshInfo = new MeshInfo( MeshInfo::meshAnnulus( 1.0, 2.0, 0.25 * pi, 0.75 * pi, MeshInfo::CRISSCROSS, 4, 2 ) );
      break;

   case ANNULUS:
      meshInfo = new MeshInfo( MeshInfo::meshAnnulus( 1.0, 2.0, MeshInfo::CROSS, 16, 5 ) );
      break;

   case FACE_CHAIN:
      meshInfo = new MeshInfo( MeshInfo::meshFaceChain( numFaces ) );
      break;

   case SPHERICAL_SHELL:
      meshInfo = new MeshInfo( MeshInfo::meshSphericalShell( ntan, layers ) );
      break;

   case THIN_SPHERICAL_SHELL:
      meshInfo = new MeshInfo( MeshInfo::meshThinSphericalShell( ntan, thinShellRadius ) );
      break;

   case CUBOID:
      meshInfo = new MeshInfo(
          MeshInfo::meshCuboid( Point3D( -1.0, -1.0, 0.0 ), Point3D( 2.0, 0.0, 2.0 ), nHint + 1, nHint + 1, nHint ) );
      break;

   case SYMM_CUBOID:
      meshInfo = new MeshInfo(
          MeshInfo::meshSymmetricCuboid( Point3D( -1.0, -1.0, -1.0 ), Point3D( 1.0, 1.0, 1.0 ), nHint, nHint, nHint ) );
      break;
   case T_DOMAIN: {
      std::set< std::array< int, 3 > > cubes;
      cubes.insert( { 0, 0, 0 } );
      for ( int i = 0; i <= walberla::int_c( nHint ); i++ )
      {
         cubes.insert( { -i, 0, 0 } );
         cubes.insert( { 0, i, 0 } );
         cubes.insert( { 0, -i, 0 } );
      }
      meshInfo = new MeshInfo( MeshInfo::meshCubedDomain( cubes, 1 ) );
      break;
   }
   case TORUS:

      const uint_t                toroidalResolution         = 8;
      const uint_t                poloidalResolution         = 6;
      const real_t                radiusOriginToCenterOfTube = 6.2;
      const std::vector< real_t > tubeLayerRadii             = { 3 };
      const real_t                torodialStartAngle         = 0.0;
      const real_t                polodialStartAngle         = 2.0 * pi / real_c( 2 * poloidalResolution );

      meshInfo = new MeshInfo( MeshInfo::meshTorus( toroidalResolution,
                                                    poloidalResolution,
                                                    radiusOriginToCenterOfTube,
                                                    tubeLayerRadii,
                                                    torodialStartAngle,
                                                    polodialStartAngle ) );
      break;
   }

   // ----------------
   //  Log mesh info
   // ----------------

   if ( beVerbose )
   {
      MeshInfo::VertexContainer verts = meshInfo->getVertices();
      std::ostringstream        msg;

      msg << "VERTEX INFO:\n";
      for ( const auto& it : verts )
      {
         msg << "node " << it.first << ": mesh boundary flag = " << it.second.getBoundaryFlag()
             << " | pos = " << it.second.getCoordinates() << "\n";
      }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str( "" );

      msg << "EDGE INFO: (vertex indices)\n";
      MeshInfo::EdgeContainer edges = meshInfo->getEdges();
      for ( const auto& it : edges )
      {
         std::array< MeshInfo::IDType, 2 > node = it.second.getVertices();
         msg << node[0] << " <--> " << node[1] << " : mesh boundary flag = " << it.second.getBoundaryFlag() << std::endl;
      }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str( "" );

      msg << "FACE INFO (vertex indices):\n";
      MeshInfo::FaceContainer faces = meshInfo->getFaces();
      for ( const auto& it : faces )
      {
         std::vector< MeshInfo::IDType > node = it.second.getVertices();
         msg << node[0] << " <--> " << node[1] << " <--> " << node[2] << " : mesh boundary flag = " << it.second.getBoundaryFlag()
             << std::endl;
      }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str( "" );
   }

   SetupPrimitiveStorage* setupStorage = nullptr;
   setupStorage = new SetupPrimitiveStorage( *meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // check, if blending for spherical shell needs to be done
   if ( blending )
   {
      if ( meshDomain == SPHERICAL_SHELL )
      {
         IcosahedralShellMap::setMap( *setupStorage );
         WALBERLA_LOG_INFO_ON_ROOT( "added geometry map for blending" );
      }
      else if ( meshDomain == THIN_SPHERICAL_SHELL )
      {
         ThinShellMap::setMap( *setupStorage, thinShellRadius );
         WALBERLA_LOG_INFO_ON_ROOT( "added geometry map for blending" );
      }
   }

   switch ( loadBalancingType )
   {
   case ALL_ROOT:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: all on root" );
      hyteg::loadbalancing::allPrimitivesOnRoot( *setupStorage );
      break;
   case ROUND_ROBIN:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: round robin" );
      hyteg::loadbalancing::roundRobin( *setupStorage );
      break;
   case ROUND_ROBIN_VOLUME:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: round robin volume" );
      hyteg::loadbalancing::roundRobinVolume( *setupStorage );
      break;
   case GREEDY:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: greedy" );
      hyteg::loadbalancing::greedy( *setupStorage );
      break;
   case DIFFUSIVE_CLUSTER:
      hyteg::loadbalancing::roundRobin( *setupStorage );
      break;
   default:
      break;
   }

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" << *setupStorage );
   }

   const size_t minLevel = 2;
   const size_t maxLevel = std::max( minLevel, (size_t) 2 );
   const size_t outLevel = minLevel;

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   switch ( loadBalancingType )
   {
   case DIFFUSIVE_CLUSTER:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: diffusive cluster" );
      hyteg::loadbalancing::distributed::diffusiveSmooth( *storage, 10, 1 );
      break;
   default:
      break;
   }

   if ( beVerbose )
   {
      std::string pInfo = storage->getGlobalInfo();
      WALBERLA_LOG_INFO_ON_ROOT( "" << pInfo );
   }

#ifdef WALBERLA_BUILD_WITH_PARMETIS
   if ( loadBalancingType == PARMETIS )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: parmetis" );
      hyteg::loadbalancing::distributed::parmetis( *storage );
   }
#endif

   hyteg::writeDomainPartitioningVTK( storage, "../output", vtkFileName + "_domain_partitioning" );
   WALBERLA_LOG_INFO_ON_ROOT(
       "Wrote domain partitioning (incl. rank assignment) and mesh boundary flags to files with base name: "
       << vtkFileName + "_domain_partitioning" );

   hyteg::VTKOutput                                 vtkOutput( "../output", vtkFileName, storage );
   hyteg::P1Function< real_t >                      someData( "test data", storage, minLevel, maxLevel );
   std::function< real_t( const hyteg::Point3D& ) > myFunc = []( const hyteg::Point3D& xx ) {
      return xx[0] * xx[0] - xx[1] * xx[1];
   };
   someData.interpolate( myFunc, maxLevel );
   vtkOutput.add( someData );
   WALBERLA_LOG_INFO_ON_ROOT( "Output goes to file with basename: " << vtkFileName );
   vtkOutput.write( outLevel );

   // -------------------
   //  Report DoF Counts
   // -------------------
   if ( reportDoFCount )
   {
      uint_t vDoFs = numberOfGlobalDoFs< VertexDoFFunctionTag >( *storage, dofLevel );
      uint_t eDoFs = numberOfGlobalDoFs< EdgeDoFFunctionTag >( *storage, dofLevel );
      WALBERLA_LOG_INFO_ON_ROOT( "\nDOF INFO:" );
      WALBERLA_LOG_INFO_ON_ROOT( "level ............ " << dofLevel );
      WALBERLA_LOG_INFO_ON_ROOT( "# vertex DoFs .... " << vDoFs );
      WALBERLA_LOG_INFO_ON_ROOT( "# edge DoFs ...... " << eDoFs );
   }

   delete meshInfo;
   delete setupStorage;

   return 0;
}
