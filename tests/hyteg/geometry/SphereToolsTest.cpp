
#include "hyteg/geometry/SphereTools.hpp"

#include <memory>

#include "core/Environment.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/types/pointnd.hpp"

namespace hyteg {

void writeSimpleVTKFile( const std::vector< Point3D >& points, std::string filename )
{
   const auto n = points.size();

   std::stringstream ss;
   ss << "<?xml version=\"1.0\"?>\n"
         "  <VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">\n"
         "    <UnstructuredGrid>\n";
   ss << "      <Piece NumberOfPoints=\"" << n << "\" NumberOfCells=\"" << n
      << "\">\n"
         "        <Points>\n"
         "          <DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
   for ( auto p : points )
   {
      ss << p[0] << " " << p[1] << " " << p[2] << "\n";
   }

   ss << "          </DataArray>\n"
         "        </Points>\n"
         "        <Cells>\n"
         "          <DataArray type=\"Int32\" Name=\"connectivity\">\n";

   for ( int i = 0; i < points.size(); i++ )
   {
      ss << i << "\n";
   }

   ss << "          </DataArray>\n"
         "          <DataArray type=\"Int32\" Name=\"offsets\">\n";

   for ( int i = 0; i < points.size(); i++ )
   {
      ss << i + 1 << "\n";
   }

   ss << "          </DataArray>\n"
         "          <DataArray type=\"UInt8\" Name=\"types\">\n";

   for ( int i = 0; i < points.size(); i++ )
   {
      ss << 1 << "\n";
   }

   ss << "          </DataArray>\n"
         "        </Cells>\n"
         "      </Piece>\n"
         "    </UnstructuredGrid>\n"
         "  </VTKFile>";

   std::ofstream fstream;
   fstream.open( filename );
   fstream << ss.str();
   fstream.close();
}

void referenceDomain( uint_t nTan, uint_t nRad, real_t rMin, real_t rMax )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   meshInfo = MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   loadbalancing::roundRobinVolume( *setupStorage );

   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   IcosahedralShellMap::setMap( *setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   writeDomainPartitioningVTK( *storage, "../../output", "SphereToolsTest_Domain" );

   const uint_t level = 2;

   P1Function< real_t > someFunction( "someFunction", storage, level, level );

   auto f = []( const Point3D& p ) { return p[0] * p[1] * p[2]; };
   someFunction.interpolate( f, level );

   VTKOutput vtk( "../../output", "SphereToolsTest_Function", storage );
   vtk.add( someFunction );
   vtk.write( level );
}

void testEvaluateUV()
{
   std::vector< Point3D > vertices;
   std::vector< int >     meridians;
   std::vector< int >     parallels;

   uvSphereSurfaceVertices( 1, 16, 16, vertices, meridians, parallels );
   writeSimpleVTKFile( vertices, "../../output/SphereToolsTest_UV_m16_p16.vtu" );
}

void testEvaluateIco()
{
   std::vector< Point3D >                 vertices;
   std::vector< std::array< uint_t, 3 > > triangles;

   icosahedralSurfaceTriangles( 1, 0, vertices, triangles );
   writeSimpleVTKFile( vertices, "../../output/SphereToolsTest_Ico_0_Refinements.vtu" );

   icosahedralSurfaceTriangles( 1, 1, vertices, triangles );
   writeSimpleVTKFile( vertices, "../../output/SphereToolsTest_Ico_1_Refinements.vtu" );

   icosahedralSurfaceTriangles( 1, 2, vertices, triangles );
   writeSimpleVTKFile( vertices, "../../output/SphereToolsTest_Ico_2_Refinements.vtu" );
}
} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::referenceDomain( 3, 3, 0.55, 1 );
   hyteg::testEvaluateUV();
   hyteg::testEvaluateIco();

   return 0;
}