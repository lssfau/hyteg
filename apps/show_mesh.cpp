/*
 * Copyright (c) 2017-2024 Christoph Schwarzmeier, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include <filesystem>
#include <sstream>
#include <utility>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/Gmsh/ExportRefinedMesh.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/geometry/AnnulusAlignedMap.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellAlignedMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/TorusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
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

/// \file show_mesh.cpp
///
/// \brief App demonstrating the use of mesh-reader and inline mesh-generators
///
/// The app demonstrates the use of various inline mesh-generators available in
/// HyTeG and can also be used to preview the resulting meshes or a GMSH file,
/// obtain info on the resulting primitive distribution for a given load balancing
/// scheme. Additionally one can query the number of degrees of freedom for a
/// mesh on a specific refinement level.

enum class LoadBalancingType
{
   ALL_ROOT,
   ROUND_ROBIN,
   ROUND_ROBIN_VOLUME,
   GREEDY,
   PARMETIS,
   DIFFUSIVE_CLUSTER
};

struct MeshBuilder
{
   MeshBuilder( std::string _flag, std::string _options, std::string _description )
   : flag( std::move( _flag ) )
   , options( std::move( _options ) )
   , description( std::move( _description ) )
   {}

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const = 0;

   virtual SetupPrimitiveStorage constructSetupStorage( const std::vector< std::string >& allArguments,
                                                        const MeshInfo&                   meshInfo ) const
   {
      WALBERLA_UNUSED( allArguments );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      return setupStorage;
   }

   virtual std::shared_ptr< PrimitiveStorage >
       constructPrimitiveStorage( const std::vector< std::string >& allArguments,
                                  const SetupPrimitiveStorage&      setupPrimitiveStorage ) const
   {
      WALBERLA_UNUSED( allArguments );
      auto primitiveStorage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage, 1 );
      return primitiveStorage;
   }

   const std::string flag;
   const std::string options;
   const std::string description;
};

struct MeshFromFile : public MeshBuilder
{
   MeshFromFile()
   : MeshBuilder( "--file", "<name of Gmsh file>", "importing a file in Gmsh format" )
   {}

   virtual ~MeshFromFile() = default;

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const
   {
      auto meshFileName = *std::next( std::find( allArguments.begin(), allArguments.end(), flag ) );
      return MeshInfo::fromGmshFile( meshFileName );
   }
};

struct MeshRectangle : public MeshBuilder
{
   MeshRectangle()
   : MeshBuilder( "--rect", "[criss|cross|crisscross|diamond]", "meshing rectangle in a certain flavour" )
   {}

   virtual ~MeshRectangle() = default;

   const std::map< std::string, MeshInfo::meshFlavour > rectMeshTypes = {
       { "criss", MeshInfo::CRISS },
       { "cross", MeshInfo::CROSS },
       { "crisscross", MeshInfo::CRISSCROSS },
       { "diamond", MeshInfo::DIAMOND },
   };

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const
   {
      auto flavour = *( std::next( std::find( allArguments.begin(), allArguments.end(), flag ) ) );
      WALBERLA_CHECK_EQUAL( rectMeshTypes.count( flavour ), 1, "Flavour for rect mesh not recognised!" );
      return MeshInfo::meshRectangle( Point2D( -2.0, 1.0 ), Point2D( 0.0, 3.0 ), rectMeshTypes.at( flavour ), 16, 16 );
   }
};

struct MeshCuboid : public MeshBuilder
{
   MeshCuboid()
   : MeshBuilder( "--cuboid", "[simple|symmetric] [number-of-cubes-in-each-direction]", "meshing cuboid in a certain flavour" )
   {}

   virtual ~MeshCuboid() = default;

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const
   {
      auto itFlag  = std::find( allArguments.begin(), allArguments.end(), flag );
      auto flavour = *( std::next( itFlag ) );
      WALBERLA_CHECK_UNEQUAL( std::next( itFlag, 2 ), allArguments.end(), "Parameter nCubes missing." );
      auto nCubes = uint_c( std::stoi( *std::next( itFlag, 2 ) ) );

      if ( flavour == "simple" )
      {
         return MeshInfo::meshCuboid( Point3D( -1.0, -1.0, 0.0 ), Point3D( 2.0, 0.0, 2.0 ), nCubes + 1, nCubes + 1, nCubes );
      }
      else if ( flavour == "symmetric" )
      {
         return MeshInfo::meshSymmetricCuboid( Point3D( -1.0, -1.0, -1.0 ), Point3D( 1.0, 1.0, 1.0 ), nCubes, nCubes, nCubes );
      }
      else
      {
         WALBERLA_ABORT( "Invalid parameters." )
      }
   }
};

struct MeshAnnulus : public MeshBuilder
{
   MeshAnnulus()
   : MeshBuilder( "--annulus",
                  "[full|partial] [affine|blended|blended-aligned]",
                  "meshing a full or partial annulus (potentially radially aligned)" )
   {}

   virtual ~MeshAnnulus() = default;

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const
   {
      auto itFlag  = std::find( allArguments.begin(), allArguments.end(), flag );
      auto flavour = *( std::next( itFlag ) );

      MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
      if ( flavour == "full" )
      {
         return meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, MeshInfo::CROSS, 16, 5 );
      }
      else if ( flavour == "partial" )
      {
         return MeshInfo::meshAnnulus( 1.0, 2.0, 0.25 * pi, 0.75 * pi, MeshInfo::CRISSCROSS, 4, 2 );
      }
      WALBERLA_ABORT( "Invalid annulus parameter." )
   }

   virtual SetupPrimitiveStorage constructSetupStorage( const std::vector< std::string >& allArguments,
                                                        const MeshInfo&                   meshInfo ) const
   {
      auto itFlag   = std::find( allArguments.begin(), allArguments.end(), flag );
      auto flavour  = *( std::next( itFlag ) );
      auto blending = *( std::next( itFlag, 2 ) );

      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      if ( blending == "affine" ) {}
      else if ( blending == "blended" )
      {
         WALBERLA_CHECK_EQUAL( flavour, "full", "Only the 'full' annulus can be blended." )
         AnnulusMap::setMap( setupStorage );
      }
      else if ( blending == "blended-aligned" )
      {
         WALBERLA_CHECK_EQUAL( flavour, "full", "Only the 'full' annulus can be blended." )
         AnnulusAlignedMap::setMap( setupStorage );
      }
      else
      {
         WALBERLA_ABORT( "Invalid annulus parameter." )
      }
      return setupStorage;
   }
};

struct MeshIcosahedralShell : public MeshBuilder
{
   MeshIcosahedralShell()
   : MeshBuilder( "--spherical-shell",
                  "[ntan] [nrad] [thin|thick] [affine|blended|blended-aligned]",
                  "meshing a thin or thick spherical shell (w/ or w/o blending, optionally radially aligned) - thick shell with "
                  "inner radius 0.5 and outer radius 1.0" )

   {}

   virtual ~MeshIcosahedralShell() = default;

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const
   {
      auto itFlag  = std::find( allArguments.begin(), allArguments.end(), flag );
      auto ntan    = uint_c( std::stoi( *( std::next( itFlag ) ) ) );
      auto nrad    = uint_c( std::stoi( *( std::next( itFlag, 2 ) ) ) );
      auto flavour = *( std::next( itFlag, 3 ) );

      const auto thinShellRadius = real_t( 0.5 );
      const auto rmin            = real_c( 0.5 );
      const auto rmax            = real_c( 1.0 );

      MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
      if ( flavour == "thin" )
      {
         return meshInfo = MeshInfo::meshThinSphericalShell( ntan, thinShellRadius );
      }
      else if ( flavour == "thick" )
      {
         return MeshInfo::meshSphericalShell( ntan, nrad, rmin, rmax );
      }
      WALBERLA_ABORT( "Invalid shell parameter." )
   }

   virtual SetupPrimitiveStorage constructSetupStorage( const std::vector< std::string >& allArguments,
                                                        const MeshInfo&                   meshInfo ) const
   {
      auto itFlag   = std::find( allArguments.begin(), allArguments.end(), flag );
      auto flavour  = *( std::next( itFlag, 3 ) );
      auto blending = *( std::next( itFlag, 4 ) );

      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      if ( blending == "affine" ) {}
      else if ( blending == "blended" )
      {
         WALBERLA_CHECK_EQUAL( flavour, "thick", "Only the 'thick' shell can be blended." )
         IcosahedralShellMap::setMap( setupStorage );
      }
      else if ( blending == "blended-aligned" )
      {
         WALBERLA_CHECK_EQUAL( flavour, "thick", "Only the 'thick' shell can be blended." )
         IcosahedralShellAlignedMap::setMap( setupStorage );
      }
      else
      {
         WALBERLA_ABORT( "Invalid shell parameter." )
      }
      return setupStorage;
   }
};

struct MeshTDomain : public MeshBuilder
{
   MeshTDomain()
   : MeshBuilder( "--t-domain", "[number-of-cubes-in-each-direction]", "meshing a T-domain using the cubed domain generator" )
   {}

   virtual ~MeshTDomain() = default;

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const
   {
      const auto nHint = uint_c( std::stoi( allArguments[2] ) );

      std::set< std::array< int, 3 > > cubes;
      cubes.insert( { 0, 0, 0 } );
      for ( int i = 0; i <= walberla::int_c( nHint ); i++ )
      {
         cubes.insert( { -i, 0, 0 } );
         cubes.insert( { 0, i, 0 } );
         cubes.insert( { 0, -i, 0 } );
      }
      return MeshInfo::meshCubedDomain( cubes, 1 );
   }
};

struct MeshTorus : public MeshBuilder
{
   MeshTorus()
   : MeshBuilder( "--torus", "[affine|blended]", "meshing a torus (with or without blending)" )
   {}

   const uint_t                toroidalResolution         = 8;
   const uint_t                poloidalResolution         = 6;
   const real_t                radiusOriginToCenterOfTube = 6.2;
   const std::vector< real_t > tubeLayerRadii             = { 3 };
   const real_t                torodialStartAngle         = 0.0;
   const real_t                polodialStartAngle         = 2.0 * pi / real_c( 2 * poloidalResolution );

   virtual ~MeshTorus() = default;

   virtual MeshInfo constructMeshInfo( const std::vector< std::string >& allArguments ) const
   {
      WALBERLA_UNUSED( allArguments );

      return MeshInfo::meshTorus( toroidalResolution,
                                  poloidalResolution,
                                  radiusOriginToCenterOfTube,
                                  tubeLayerRadii,
                                  torodialStartAngle,
                                  polodialStartAngle );
   }

   virtual SetupPrimitiveStorage constructSetupStorage( const std::vector< std::string >& allArguments,
                                                        const MeshInfo&                   meshInfo ) const
   {
      auto itFlag   = std::find( allArguments.begin(), allArguments.end(), flag );
      auto blending = *( std::next( itFlag, 1 ) );

      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      if ( blending == "affine" ) {}
      else if ( blending == "blended" )
      {
         TorusMap::setMap( setupStorage,
                           toroidalResolution,
                           poloidalResolution,
                           radiusOriginToCenterOfTube,
                           tubeLayerRadii,
                           torodialStartAngle,
                           polodialStartAngle );
      }
      else
      {
         WALBERLA_ABORT( "Invalid shell parameter." )
      }
      return setupStorage;
   }
};

std::vector< std::shared_ptr< MeshBuilder > > makeBuilders()
{
   std::vector< std::shared_ptr< MeshBuilder > > builders;
   builders.push_back( std::make_shared< MeshFromFile >() );
   builders.push_back( std::make_shared< MeshRectangle >() );
   builders.push_back( std::make_shared< MeshCuboid >() );
   builders.push_back( std::make_shared< MeshAnnulus >() );
   builders.push_back( std::make_shared< MeshIcosahedralShell >() );
   builders.push_back( std::make_shared< MeshTDomain >() );
   builders.push_back( std::make_shared< MeshTorus >() );
   return builders;
}

/// Print information on the app and available command-line parameters
void showUsage( const std::vector< std::shared_ptr< MeshBuilder > >& builders )
{
   std::stringstream usageStream;
   usageStream << "--------  USAGE: --------\n";
   usageStream << "show_mesh demonstrates generation of a MeshInfo object by one of the following methods:\n\n";
   for ( auto builder : builders )
   {
      usageStream << " " << builder->flag << " " << builder->options << "\n"
                  << "     " << builder->description << "\n\n";
   }
   usageStream
       << "The generated base mesh will be tested be doing --level levels of refinement.\n"
       << "Then it will be exported to a VTU file for visualisation.\n"
       << "Also visualization of the domain partitioning, mesh boundary flags and MPI rank assignment will be output.\n\n"
       << "The desired load balancing approach can be steered by:\n"
       << "  --load-balancing [allroot|roundrobin|roundrobinvolume|greedy|parmetis|diffusivecluster] (default: roundrobin)\n\n"
       << "Use -v                            to print info on mesh to console.\n\n"
       << "Use --level <level>               to specify the refinement level for all outputs (default = 2).\n\n"
       << "Use --export-fine-mesh            to store the mesh for refinement level lvl in a MESH4.1 file.\n\n"
       << "Use --parametric-map <map-degree> to write the node positions to the micro-mesh and plot via a linear or quadratic\n"
       << "                                  parametric map. Only degrees 1 and 2 are supported.";

   WALBERLA_LOG_INFO_ON_ROOT( "" << usageStream.str() );
}

/// The actual show_mesh app itself
int main( int argc, char* argv[] )
{
   walberla::mpi::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string       vtkFileName       = "show_mesh_output";
   LoadBalancingType loadBalancingType = LoadBalancingType::ROUND_ROBIN;
   bool              beVerbose         = false;
   uint_t            level             = 2;
   bool              exportFineMesh    = false;

   auto builders = makeBuilders();

   std::vector< std::string > allArguments( (size_t) argc );
   for ( size_t i = 0; i < (size_t) argc; i++ )
   {
      allArguments[i] = argv[i];
      if ( i > 0 )
      {
         vtkFileName += "_" + allArguments[i];
      }
   }

   std::replace( vtkFileName.begin(), vtkFileName.end(), '.', '_' );
   std::replace( vtkFileName.begin(), vtkFileName.end(), '/', '_' );

   if ( allArguments.size() == 1 || algorithms::contains( allArguments, "-h" ) || algorithms::contains( allArguments, "--help" ) )
   {
      showUsage( builders );
      return EXIT_SUCCESS;
   }

   std::shared_ptr< MeshBuilder > builder;
   for ( const auto& b : builders )
   {
      if ( algorithms::contains( allArguments, b->flag ) )
      {
         builder = b;
         WALBERLA_LOG_INFO_ON_ROOT( "Builder info: " << b->description )
         break;
      }
   }

   WALBERLA_CHECK_NOT_NULLPTR( builder, "Invalid mesh type." )

   const auto loadBalancingFlag = "--load-balancing";
   if ( algorithms::contains( allArguments, loadBalancingFlag ) )
   {
      if ( algorithms::contains( allArguments, "allroot" ) )
      {
         loadBalancingType = LoadBalancingType::ALL_ROOT;
      }
      else if ( algorithms::contains( allArguments, "roundrobin" ) )
      {
         loadBalancingType = LoadBalancingType::ROUND_ROBIN;
      }
      else if ( algorithms::contains( allArguments, "roundrobinvolume" ) )
      {
         loadBalancingType = LoadBalancingType::ROUND_ROBIN_VOLUME;
      }
      else if ( algorithms::contains( allArguments, "greedy" ) )
      {
         loadBalancingType = LoadBalancingType::GREEDY;
      }
      else if ( algorithms::contains( allArguments, "diffusivecluster" ) )
      {
         loadBalancingType = LoadBalancingType::DIFFUSIVE_CLUSTER;
      }
      else if ( algorithms::contains( allArguments, "parmetis" ) )
      {
#ifdef WALBERLA_BUILD_WITH_PARMETIS
         loadBalancingType = LoadBalancingType::PARMETIS;
#else
         WALBERLA_ABORT( "Framework was not built with ParMetis." );
#endif
      }
      else
      {
         WALBERLA_ABORT( "Could not understand command-line args. Possibly invalid load balancing approach." );
      }
   }

   if ( algorithms::contains( allArguments, "-v" ) )
   {
      beVerbose = true;
   }

   if ( algorithms::contains( allArguments, "--level" ) )
   {
      auto flagIt = std::find( allArguments.begin(), allArguments.end(), "--level" );
      WALBERLA_CHECK_UNEQUAL( std::next( flagIt ), allArguments.end(), "The option '--level' requires an integer parameter." )
      level = uint_c( std::stoi( *std::next( flagIt ) ) );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Output refinement level: " << level )

   exportFineMesh = algorithms::contains( allArguments, "--export-fine-mesh" );

   // ------------
   //  Let's rock
   // ------------

   auto meshInfo = builder->constructMeshInfo( allArguments );

   // ----------------
   //  Log mesh info
   // ----------------

   if ( beVerbose )
   {
      MeshInfo::VertexContainer verts = meshInfo.getVertices();
      std::ostringstream        msg;

      msg << "VERTEX INFO:\n";
      for ( const auto& it : verts )
      {
         msg << "node " << it.first << ":\n* mesh boundary flag = " << it.second.getBoundaryFlag() << "\n* coordinates =\n"
             << it.second.getCoordinates() << "\n";
      }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str( "" );

      msg << "EDGE INFO: (vertex indices)\n";
      MeshInfo::EdgeContainer edges = meshInfo.getEdges();
      for ( const auto& it : edges )
      {
         std::array< MeshInfo::IDType, 2 > node = it.second.getVertices();
         msg << node[0] << " <--> " << node[1] << " : mesh boundary flag = " << it.second.getBoundaryFlag() << std::endl;
      }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str( "" );

      msg << "FACE INFO (vertex indices):\n";
      MeshInfo::FaceContainer faces = meshInfo.getFaces();
      for ( const auto& it : faces )
      {
         std::vector< MeshInfo::IDType > node = it.second.getVertices();
         msg << node[0] << " <--> " << node[1] << " <--> " << node[2] << " : mesh boundary flag = " << it.second.getBoundaryFlag()
             << std::endl;
      }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str( "" );

      msg << "CELL INFO (vertex indices):\n";
      MeshInfo::CellContainer cells = meshInfo.getCells();
      if ( cells.size() > 0 )
      {
         for ( const auto& it : cells )
         {
            std::vector< MeshInfo::IDType > node = it.second.getVertices();
            WALBERLA_CHECK_EQUAL( node.size(), 4, "Cell seems to have " << node.size() << " vertices instead of 4!" );
            msg << node[0] << " <--> " << node[1] << " <--> " << node[2] << " <--> " << node[3]
                << " : mesh boundary flag = " << it.second.getBoundaryFlag() << std::endl;
         }
         WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      }
   }

   auto setupStorage = builder->constructSetupStorage( allArguments, meshInfo );

   switch ( loadBalancingType )
   {
   case LoadBalancingType::ALL_ROOT:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: all on root" );
      hyteg::loadbalancing::allPrimitivesOnRoot( setupStorage );
      break;
   case LoadBalancingType::ROUND_ROBIN:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: round robin" );
      hyteg::loadbalancing::roundRobin( setupStorage );
      break;
   case LoadBalancingType::ROUND_ROBIN_VOLUME:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: round robin volume" );
      hyteg::loadbalancing::roundRobinVolume( setupStorage );
      break;
   case LoadBalancingType::GREEDY:
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: greedy" );
      hyteg::loadbalancing::greedy( setupStorage );
      break;
   case LoadBalancingType::DIFFUSIVE_CLUSTER:
      hyteg::loadbalancing::roundRobin( setupStorage );
      break;
   default:
      break;
   }

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );
   }

   const size_t minLevel = level;
   const size_t maxLevel = level;
   const size_t outLevel = level;

   auto storage = builder->constructPrimitiveStorage( allArguments, setupStorage );

   if ( algorithms::contains( allArguments, "--parametric-map" ) )
   {
      auto flagIt              = std::find( allArguments.begin(), allArguments.end(), "--parametric-map" );
      auto parametricMapDegree = std::stoi( *std::next( flagIt, 1 ) );

      WALBERLA_CHECK_GREATER_EQUAL( parametricMapDegree, 1, "Parametric map degree has to be 1 or 2." );
      WALBERLA_CHECK_LESS_EQUAL( parametricMapDegree, 2, "Parametric map degree has to be 1 or 2." );

      WALBERLA_LOG_INFO_ON_ROOT(
          "Using parametric mapping.\n"
          "Essentially just copying over the node positions of the (possibly blended) elements to the micro-mesh.\n"
          "(Technically, we do not need the micro-mesh for the visualization here - but it showcases how to use it.)" )

      auto microMesh = std::make_shared< micromesh::MicroMesh >(
          storage, minLevel, maxLevel, parametricMapDegree, storage->hasGlobalCells() ? 3 : 2 );

      storage->setMicroMesh( microMesh );
   }

   switch ( loadBalancingType )
   {
   case LoadBalancingType::DIFFUSIVE_CLUSTER:
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
   if ( loadBalancingType == LoadBalancingType::PARMETIS )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Load balancing: parmetis" );
      hyteg::loadbalancing::distributed::parmetis( *storage );
   }
#endif

   const auto outputDirectory = "../output";

   WALBERLA_LOG_INFO_ON_ROOT( "Writing to directory: " << outputDirectory );

   hyteg::writeDomainPartitioningVTK( storage, outputDirectory, vtkFileName + "_domain_partitioning" );
   WALBERLA_LOG_INFO_ON_ROOT(
       "Wrote domain partitioning (incl. rank assignment) and mesh boundary flags to files with base name: "
       << vtkFileName + "_domain_partitioning" );

   hyteg::VTKOutput                                 vtkOutput( outputDirectory, vtkFileName, storage );
   hyteg::P1Function< real_t >                      someDataP1( "test data P1", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t >                      someDataP2( "test data P2", storage, minLevel, maxLevel );
   std::function< real_t( const hyteg::Point3D& ) > myFunc = []( const hyteg::Point3D& xx ) {
      return xx[0] * xx[0] - xx[1] * xx[1];
   };
   someDataP1.interpolate( myFunc, maxLevel );
   someDataP2.interpolate( myFunc, maxLevel );
   vtkOutput.add( someDataP1 );
   vtkOutput.add( someDataP2 );
   WALBERLA_LOG_INFO_ON_ROOT( "Output goes to file with basename: " << vtkFileName );
   vtkOutput.write( outLevel );

   // -------------------
   //  Report DoF Counts
   // -------------------
   uint_t vertexDoFs = numberOfGlobalDoFs< VertexDoFFunctionTag >( *storage, level );
   uint_t edgeDoFs   = numberOfGlobalDoFs< EdgeDoFFunctionTag >( *storage, level );
   uint_t volumeDoFs = numberOfGlobalDoFs< P0FunctionTag >( *storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of micro-vertices ... " << vertexDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of micro-edges ...... " << edgeDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of micro-volumes .... " << volumeDoFs );

   // ---------------------
   //  Export Refined Mesh
   // ---------------------
   if ( exportFineMesh )
   {
      std::stringstream streamExportFileName;
      streamExportFileName << vtkFileName << "-lvl=" << level << ".msh";
      gmsh::exportRefinedMesh( storage, level, streamExportFileName.str() );
   }

   return 0;
}
