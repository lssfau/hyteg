#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <sstream>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

void showUsage() {
  std::cout << "\n --------\n  USAGE:\n --------\n\n"
            << " show_mesh demonstrates generation of a MeshInfo object by one of three methods:\n\n"
            << " 1) by importing a file in Gmsh format\n"
            << " 2) by meshing a rectangle in a certain flavour\n"
            << " 3) by meshing a full or partial annulus\n\n"
            << " This is steered by choosing one of the options below:\n\n"
            << "  --file <name of Gmsh file>\n"
            << "  --rect [criss|cross|crisscross|diamond]\n"
            << "  --annulus [full|partial]\n\n"
            << " The generated base mesh will be tested be doing two levels of refinement.\n"
            << " Then it will be exported to a VTU file for visualisation.\n\n"
            << " Use  -v  to print info on mesh to console.\n\n"
            << std::endl;
}

int main( int argc, char* argv[] )
{

  std::string meshFileName;
  std::string vtkFileName;
  typedef enum { FROM_FILE, RECTANGLE, ANNULUS, PARTIAL_ANNULUS } meshDomainType;
  meshDomainType meshDomain;
  MeshInfo::meshFlavour rectMeshType;
  MeshInfo *meshInfo = nullptr;
  bool beVerbose = false;

  if( argc < 3 || argc > 4 )
    {
      showUsage();
      WALBERLA_ABORT( "Please provide command-line parameters!" );
    }
  else if( strcmp( argv[1], "--file" ) == 0 )
    {
      meshDomain = FROM_FILE;
      meshFileName = std::string( argv[2] );
      auto pos = meshFileName.find_last_of( '/' );
      pos == meshFileName.npos ? pos = 0 : ++pos;
      vtkFileName = meshFileName.substr( pos, meshFileName.length() - 4 );
    }
  else if( strcmp( argv[1], "--rect" ) == 0 )
    {
      meshDomain = RECTANGLE;
      if( strcmp( argv[2], "criss" ) == 0 )
        {
          rectMeshType = MeshInfo::CRISS;
          vtkFileName = std::string( "rectMeshCriss" );
        }
      else if( strcmp( argv[2], "cross" ) == 0 )
        {
          rectMeshType = MeshInfo::CROSS;
          vtkFileName = std::string( "rectMeshCross" );
        }
      else if( strcmp( argv[2], "crisscross" ) == 0 )
        {
          rectMeshType = MeshInfo::CRISSCROSS;
          vtkFileName = std::string( "rectMeshCrissCross" );
        }
      else if( strcmp( argv[2], "diamond" ) == 0 )
        {
          rectMeshType = MeshInfo::DIAMOND;
          vtkFileName = std::string( "rectMeshDiamond" );
        }
      else
        {
          WALBERLA_ABORT( "Flavour for rect mesh not recognised!" );
        }
    }
  else if( strcmp( argv[1], "--annulus" ) == 0 )
    {
      if( strcmp( argv[2], "full" ) == 0 )
        {
          meshDomain = ANNULUS;
          vtkFileName = std::string( "annulusMesh" );
        }
      else if( strcmp( argv[2], "partial" ) == 0 )
        {
          meshDomain = PARTIAL_ANNULUS;
          vtkFileName = std::string( "partialAnnulusMesh" );
        }
      else
        {
          WALBERLA_ABORT( "Subtype of --annulus not recognised!" );
        }
    }
  else
    {
      WALBERLA_ABORT( "Could not understand command-line args!" );
    }

  if( argc == 4 && strcmp( argv[3], "-v" ) == 0 )
    {
      beVerbose = true;
    }


  // ------------
  //  Let's rock
  // ------------
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  WALBERLA_LOG_INFO_ON_ROOT("TinyHHG Show Mesh Test\n");

  switch( meshDomain )
    {
    case FROM_FILE:
      meshInfo = new MeshInfo( MeshInfo::fromGmshFile( meshFileName ) );
      break;

    case RECTANGLE:
      meshInfo = new MeshInfo( MeshInfo::meshRectangle( Point2D( {-2.0, 1.0} ), Point2D( {0.0, 3.0} ),
                                                        rectMeshType, 3, 2 ) );
      break;

    case PARTIAL_ANNULUS:
      meshInfo = new MeshInfo( MeshInfo::meshAnnulus( 1.0, 2.0, 0.25*M_PI, 0.75*M_PI,
                                                      MeshInfo::CRISSCROSS, 4, 2 ) );
      break;

    case ANNULUS:
      meshInfo = new MeshInfo( MeshInfo::meshAnnulus( 1.0, 2.0, 15, 2 ) );
      break;
    }

  // ----------------
  //  Log mesh info
  // ----------------

  if( beVerbose )
    {
      MeshInfo::VertexContainer verts = meshInfo->getVertices();
      std::ostringstream msg;

      msg << "VERTEX INFO:\n";
      for ( const auto & it : verts )
        {
          msg << "node " << it.first << ": dofType = " << it.second.getDoFType() << "\n";
        }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str("");

      msg << "EDGE INFO: (vertex indices)\n";
      MeshInfo::EdgeContainer edges = meshInfo->getEdges();
      for ( const auto & it : edges )
        {
          std::array< MeshInfo::IDType, 2 > node = it.second.getVertices();
          msg << node[0] << " <--> " << node[1] << std::endl;
        }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str("");

      msg << "FACE INFO (vertex indices):\n";
      MeshInfo::FaceContainer faces = meshInfo->getFaces();
      for ( const auto & it : faces )
        {
          std::vector< MeshInfo::IDType > node = it.second.getVertices();
          msg << node[0] << " <--> " << node[1] << " <--> " << node[2] << std::endl;
        }
      WALBERLA_LOG_INFO_ON_ROOT( msg.str() );
      msg.str("");
    }

  SetupPrimitiveStorage *setupStorage;
  setupStorage = new SetupPrimitiveStorage( *meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( *setupStorage );

  const size_t minLevel = 0;
  const size_t maxLevel = std::max( minLevel, (size_t)2 );
  const size_t outLevel = minLevel;

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>( *setupStorage );

  hhg::VTKOutput vtkOutput(  "../output", vtkFileName );
  hhg::P1Function<real_t> someData( "test data", storage, minLevel, maxLevel );
  std::function<real_t( const hhg::Point3D& )> myFunc = []( const hhg::Point3D& xx )
    { return xx[0]*xx[0] - xx[1]*xx[1]; };
  someData.interpolate( myFunc, maxLevel );
  vtkOutput.add( &someData );
  WALBERLA_LOG_INFO_ON_ROOT( "Output goes to file with basename: '" << vtkFileName << "'\n" );
  vtkOutput.write( outLevel );

  delete meshInfo;

  return 0;
}
