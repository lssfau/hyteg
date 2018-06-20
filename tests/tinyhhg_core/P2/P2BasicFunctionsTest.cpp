

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/communication/Syncing.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hhg {

static void testP2BasicFunctions()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto x = P2Function< real_t >( "x", storage, minLevel, maxLevel );
   auto y = P2Function< real_t >( "y", storage, minLevel, maxLevel );
   auto z = P2Function< real_t >( "y", storage, minLevel, maxLevel );

   std::vector< PrimitiveID > faces;
   storage->getFaceIDs( faces );
   Face* face = storage->getFace( faces[0] );

   real_t* faceEdgeDataX = face->getData( x.getEdgeDoFFunction()->getFaceDataID() )->getPointer( maxLevel );
   real_t* faceEdgeDataY = face->getData( y.getEdgeDoFFunction()->getFaceDataID() )->getPointer( maxLevel );
   real_t* faceEdgeDataZ = face->getData( z.getEdgeDoFFunction()->getFaceDataID() )->getPointer( maxLevel );

   real_t* faceVertexDataX = face->getData( x.getVertexDoFFunction()->getFaceDataID() )->getPointer( maxLevel );
   real_t* faceVertexDataY = face->getData( y.getVertexDoFFunction()->getFaceDataID() )->getPointer( maxLevel );
   real_t* faceVertexDataZ = face->getData( z.getVertexDoFFunction()->getFaceDataID() )->getPointer( maxLevel );

   // Interpolate

   std::function< real_t( const hhg::Point3D& ) > expr  = []( const Point3D& ) -> real_t { return real_c( 2 ); };
   std::function< real_t( const hhg::Point3D& ) > zeros = []( const Point3D& ) -> real_t { return real_c( 0 ); };
   std::function< real_t( const hhg::Point3D& ) > func  = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + xx[0] ) * ( 2.0 + xx[1] ) );
   };
   std::function< real_t( const hhg::Point3D& ) > func2 = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + ( xx[0] / 5.0 ) ) * ( 2.0 + xx[1] ) );
   };

   walberla::WcTimingPool timer;

   timer["Interpolate"].start();
   x.interpolate( expr, maxLevel, DoFType::All );
   y.interpolate( expr, maxLevel, DoFType::All );
   z.interpolate( func, maxLevel, DoFType::All );
   timer["Interpolate"].end();

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataX[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataX[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row() )], real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataX[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row() )], real_c( 2 ) );

      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row() )], real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row() )], real_c( 2 ) );
   }

   hhg::communication::syncFunctionBetweenPrimitives( x, maxLevel );
   hhg::communication::syncFunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataX[vertexdof::macroface::index( maxLevel, it.col(), it.row() )], real_c( 2 ) );

      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.col(), it.row() )], real_c( 2 ) );
   }

   // Assign

   timer["Assign"].start();
   y.assign( {3.0, 2.0}, {&x, &y}, maxLevel, DoFType::All );
   z.assign( {1.0, -1.0}, {&z, &z}, maxLevel, DoFType::All );
   timer["Assign"].end();

   hhg::communication::syncFunctionBetweenPrimitives( y, maxLevel );
   hhg::communication::syncFunctionBetweenPrimitives( z, maxLevel );

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 10 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 10 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 10 ) );
   }

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.col(), it.row() )], real_c( 10 ) );
   }

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 0 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row() )], real_c( 0 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row() )], real_c( 0 ) );
   }

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataZ[vertexdof::macroface::index( maxLevel, it.col(), it.row() )], real_c( 0 ) );
   }

   // Add

   timer["Add"].start();
   y.add( {{4.0, 3.0}}, {{&x, &y}}, maxLevel, DoFType::All );
   timer["Add"].end();

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 48 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 48 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row() )],
                                  real_c( 48 ) );
   }

   // Dot

   timer["Dot"].start();
   const real_t scalarProduct = y.dot( x, maxLevel, DoFType::All );
   timer["Dot"].end();

   uint_t totalPointsOnFace = levelinfo::num_microvertices_per_face( maxLevel ) + levelinfo::num_microedges_per_face( maxLevel );
   WALBERLA_CHECK_FLOAT_EQUAL( scalarProduct, real_c( totalPointsOnFace * 48 * 2 ) );

   WALBERLA_LOG_INFO_ON_ROOT( timer );

   ///Check assign against add

   x.interpolate( func, maxLevel, DoFType::All );
   y.interpolate( func2, maxLevel, DoFType::All );
   z.interpolate( zeros, maxLevel, DoFType::All );

   z.assign( {1.0, -1.0}, {&x, &y}, maxLevel );
   hhg::communication::syncFunctionBetweenPrimitives( z, maxLevel );
   x.add( {-1.0}, {&y}, maxLevel );

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row() )],
                                  faceEdgeDataX[edgedof::macroface::horizontalIndex( maxLevel, it.col(), it.row() )] );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row() )],
                                  faceEdgeDataX[edgedof::macroface::diagonalIndex( maxLevel, it.col(), it.row() )] );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row() )],
                                  faceEdgeDataX[edgedof::macroface::verticalIndex( maxLevel, it.col(), it.row() )] );
   }

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataZ[vertexdof::macroface::index( maxLevel, it.col(), it.row() )],
                                  faceVertexDataX[vertexdof::macroface::index( maxLevel, it.col(), it.row() )] );
   }
}

} // namespace hhg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testP2BasicFunctions();

   return EXIT_SUCCESS;
}
