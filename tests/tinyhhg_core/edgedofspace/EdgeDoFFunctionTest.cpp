
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

namespace hhg {

static void testEdgeDoFFunction()
{
  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto x = std::make_shared< EdgeDoFFunction< real_t > >( "x", storage, minLevel, maxLevel );
  auto y = std::make_shared< EdgeDoFFunction< real_t > >( "y", storage, minLevel, maxLevel );

  std::vector< PrimitiveID > faces;
  storage->getFaceIDs( faces );
  Face * face = storage->getFace( faces[0] );

  real_t * faceDataX = face->getData( x->getFaceDataID() )->getPointer( maxLevel );
  real_t * faceDataY = face->getData( y->getFaceDataID() )->getPointer( maxLevel );

  // Interpolate

  std::function<real_t(const hhg::Point3D&)> expr = []( const Point3D & ) -> real_t { return real_c( 2 ); };

  walberla::WcTimingPool timer;

  timer["Interpolate"].start();
  x->interpolate( expr, maxLevel, DoFType::All );
  timer["Interpolate"].end();

  for ( const auto & it : indexing::edgedof::macroface::Iterator< maxLevel, 1 >() )
  {
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataX[ indexing::edgedof::macroface::horizontalIndex< maxLevel >( it.col(), it.row() ) ], real_c( 2 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataX[ indexing::edgedof::macroface::diagonalIndex< maxLevel >( it.col(), it.row() ) ], real_c( 2 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataX[ indexing::edgedof::macroface::verticalIndex< maxLevel >( it.col(), it.row() ) ], real_c( 2 ) );
  }

  // Assign

  timer["Assign"].start();
  y->assign( {{ 3.0 }}, {{ x.get() }}, maxLevel, DoFType::All );
  timer["Assign"].end();

  for ( const auto & it : indexing::edgedof::macroface::Iterator< maxLevel, 1 >() )
  {
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[ indexing::edgedof::macroface::horizontalIndex< maxLevel >( it.col(), it.row() ) ], real_c( 6 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[ indexing::edgedof::macroface::diagonalIndex< maxLevel >( it.col(), it.row() ) ], real_c( 6 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[ indexing::edgedof::macroface::verticalIndex< maxLevel >( it.col(), it.row() ) ], real_c( 6 ) );
  }




  WALBERLA_LOG_INFO_ON_ROOT( timer );

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
