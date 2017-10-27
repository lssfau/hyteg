
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

namespace hhg {

static void testEdgeDoFFunction()
{
  const uint_t minLevel = 2;
  const uint_t maxLevel = 3;

  MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto edgeDoFFunction = std::make_shared< EdgeDoFFunction< real_t > >( "x", storage, minLevel, maxLevel );

  auto vertexDataID = edgeDoFFunction->getVertexDataID();
  auto edgeDataID   = edgeDoFFunction->getEdgeDataID();
  auto faceDataID   = edgeDoFFunction->getFaceDataID();

  WALBERLA_UNUSED( vertexDataID );
  WALBERLA_UNUSED( edgeDataID );
  WALBERLA_UNUSED( faceDataID );

  std::function<real_t(const hhg::Point3D&)> expr = []( const Point3D & x ) -> real_t { return x[0]; };

  edgeDoFFunction->interpolate( expr, maxLevel, DoFType::All );

  std::vector< PrimitiveID > faces;
  storage->getFaceIDs( faces );
  Face * face = storage->getFace( faces[0] );
  real_t * faceData = face->getData( faceDataID )->getPointer( maxLevel );

  for ( const auto & it : indexing::edgedof::macroface::Iterator< maxLevel >() )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "col: " << it.col() << " row: " << it.row() << "\n"
                               "    E_HO = " << faceData[ indexing::edgedof::macroface::horizontalIndex< maxLevel >( it.col(), it.row() ) ] << "\n"
                               "    E_VE = " << faceData[ indexing::edgedof::macroface::verticalIndex< maxLevel >( it.col(), it.row() ) ]   << "\n"
                               "    E_DI = " << faceData[ indexing::edgedof::macroface::diagonalIndex< maxLevel >( it.col(), it.row() ) ]    );
  }


}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testEdgeDoFFunction();

   return EXIT_SUCCESS;
}
