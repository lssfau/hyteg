
#include "hyteg/communication/Syncing.hpp"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"
#include "core/DataTypes.h"

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/VTKWriter.hpp"

using walberla::uint_t;
using walberla::real_c;

namespace hyteg {

static void testEdgeDoFFunction()
{
  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  MeshInfo mesh  = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
  MeshInfo mesh2 = MeshInfo::fromGmshFile( "../../data/meshes/annulus_coarse.msh" );

  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  SetupPrimitiveStorage setupStorage2( mesh2, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage2 = std::make_shared< PrimitiveStorage >( setupStorage2 );

  auto x = std::make_shared< EdgeDoFFunction< real_t > >( "x", storage, minLevel, maxLevel );
  auto y = std::make_shared< EdgeDoFFunction< real_t > >( "y", storage, minLevel, maxLevel );

  std::vector< PrimitiveID > faces;
  storage->getFaceIDs( faces );
  Face * face = storage->getFace( faces[0] );

  real_t * faceDataX = face->getData( x->getFaceDataID() )->getPointer( maxLevel );
  real_t * faceDataY = face->getData( y->getFaceDataID() )->getPointer( maxLevel );

  // Interpolate

  std::function<real_t(const hyteg::Point3D&)> expr = []( const Point3D & ) -> real_t { return real_c( 2 ); };

  walberla::WcTimingPool timer;

  timer["Interpolate"].start();
  x->interpolate( expr, maxLevel, DoFType::All );
  y->interpolate( expr, maxLevel, DoFType::All );
  timer["Interpolate"].end();

  hyteg::communication::syncFunctionBetweenPrimitives( *x, maxLevel );
  hyteg::communication::syncFunctionBetweenPrimitives( *y, maxLevel );

  for ( const auto & it : edgedof::macroface::Iterator( maxLevel ) )
  {
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataX[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row())], real_c( 2 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataX[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row())], real_c( 2 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataX[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row())], real_c( 2 ) );

    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row())], real_c( 2 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row())], real_c( 2 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row())], real_c( 2 ) );
  }

  // Assign

  timer["Assign"].start();
  y->assign( { 3.0, 2.0 }, { *x, *y }, maxLevel, DoFType::All );
  timer["Assign"].end();

  hyteg::communication::syncFunctionBetweenPrimitives( *x, maxLevel );
  hyteg::communication::syncFunctionBetweenPrimitives( *y, maxLevel );

  for ( const auto & it : edgedof::macroface::Iterator( maxLevel ) )
  {
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row())], real_c( 10 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row())], real_c( 10 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row())], real_c( 10 ) );
  }

  // Add

  timer["Add"].start();
  y->add( {{ 4.0, 3.0 }}, {{ *x, *y }}, maxLevel, DoFType::All );
  timer["Add"].end();

  hyteg::communication::syncFunctionBetweenPrimitives( *x, maxLevel );
  hyteg::communication::syncFunctionBetweenPrimitives( *y, maxLevel );

  for ( const auto & it : edgedof::macroface::Iterator( maxLevel ) )
  {
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row())], real_c( 48 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row())], real_c( 48 ) );
    WALBERLA_CHECK_FLOAT_EQUAL( faceDataY[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row())], real_c( 48 ) );
  }

  // Dot

  timer["Dot"].start();
  const real_t scalarProduct = y->dotLocal( *x, maxLevel, DoFType::All );
  timer["Dot"].end();

  WALBERLA_CHECK_FLOAT_EQUAL( scalarProduct, real_c( levelinfo::num_microedges_per_face( maxLevel ) * 48 * 2 ) );

  WALBERLA_LOG_INFO_ON_ROOT( timer );

  // Output interpolate VTK

  auto p1 = std::make_shared< P1Function< real_t > >( "p1", storage2, minLevel, maxLevel );
  std::function<real_t(const hyteg::Point3D&)> linearX = []( const Point3D & xx ) -> real_t { return xx[0] + xx[1]; };
  p1->interpolate( linearX, maxLevel, DoFType::All );

  auto z = std::make_shared< EdgeDoFFunction< real_t > >( "z", storage2, minLevel, maxLevel );
  z->interpolate( linearX, maxLevel, DoFType::All );

  auto dg = std::make_shared< DGFunction< real_t > >( "dg", storage2, minLevel, maxLevel );

  VTKOutput vtkOutput("../../output", "interpolate_test", storage);
  vtkOutput.add( *p1 );
  vtkOutput.add( *z );
  vtkOutput.add( *dg );
  vtkOutput.write( maxLevel );
}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testEdgeDoFFunction();

   return EXIT_SUCCESS;
}
