#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = "../../data/meshes/tri_1el_neumann.msh";

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t level = 2;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P2Function< real_t >                      u( "u", storage, level, level );
   std::function< real_t( const hyteg::Point3D& ) > testExpression = []( const hyteg::Point3D& x ) { return x[1] + 1.0; };
   u.interpolate( testExpression, level, hyteg::All );

   // Sync interpolated function values
   u.getEdgeDoFFunction().communicate< Vertex, Edge >( level );
   u.getEdgeDoFFunction().communicate< Edge, Face >( level );
   u.getEdgeDoFFunction().communicate< Face, Edge >( level );
   u.getEdgeDoFFunction().communicate< Edge, Vertex >( level );

   // We assume that the bottom left vertex has following connectivity
   //
   // EdgeDoF(1)
   //    |
   //    |   EdgeDoF(2)
   //    |   /
   // Vertex(0) --- EdgeDoF(0)

   // Get bottom left vertex
   auto vertex = storage->getVertex( PrimitiveID( 0 ) );

   // Get vertex values
   auto vertexEdgeData = vertex->getData( u.getEdgeDoFFunction().getVertexDataID() )->getPointer( level );

   WALBERLA_CHECK_FLOAT_EQUAL( vertexEdgeData[0], 1.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( vertexEdgeData[1], 1.125 );
   WALBERLA_CHECK_FLOAT_EQUAL( vertexEdgeData[2], 1.125 );

   return EXIT_SUCCESS;
}
