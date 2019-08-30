#include <tinyhhg_core/communication/Syncing.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>

#include "core/Environment.h"
#include "core/math/Random.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

void test2D( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t minLevel = 2;
   const size_t maxLevel = 5;

   walberla::math::seedRandomGenerator( 12345678 );

   const uint_t numRandomEvaluations = 1000;

   auto testFunc = []( const Point3D& x ) {
      return 5.0 * x[0] * x[0] + 2.0 * x[0] * x[1] + 7.0 * x[1] * x[1] + 10.0 * x[0] + 3.0 * x[1] + 1.0;
   };
   auto testFuncDerivativeX = []( const Point3D& x ) { return 10.0 * x[0] + 2.0 * x[1] + 10.0; };
   auto testFuncDerivativeY = []( const Point3D& x ) { return 14.0 * x[1] + 2.0 * x[0] + 3.0; };

   P2Function< real_t > x( "x", storage, minLevel, maxLevel );
   x.interpolate( testFunc, maxLevel );
   // IMPORTANT
   communication::syncP2FunctionBetweenPrimitives( x, maxLevel );

   Point3D coordinates( {0.0, 0.5, 0.0} );
   real_t  eval = x.evaluate( coordinates, maxLevel );
   Point3D gradient;
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 2.0;
   coordinates[1] = 0.0;
   eval           = x.evaluate( coordinates, maxLevel );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 2.0;
   coordinates[1] = 1.0;
   eval           = x.evaluate( coordinates, maxLevel );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 0.0;
   coordinates[1] = 1.0;
   eval           = x.evaluate( coordinates, maxLevel );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   coordinates[0] = 0.5;
   coordinates[1] = 0.5;
   eval           = x.evaluate( coordinates, maxLevel );
   x.evaluateGradient( coordinates, maxLevel, gradient );
   WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
   WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );

   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      coordinates[0] = walberla::math::realRandom( 0.0, 1.0 );
      coordinates[1] = walberla::math::realRandom( 0.0, 1.0 );

      if ( coordinates[0] < 0.5 && coordinates[1] < 0.5 )
      {
         continue;
      }

      eval = x.evaluate( coordinates, maxLevel );
      x.evaluateGradient( coordinates, maxLevel, gradient );
      WALBERLA_CHECK_FLOAT_EQUAL( eval, testFunc( coordinates ) );
      WALBERLA_CHECK_FLOAT_EQUAL( gradient[0], testFuncDerivativeX( coordinates ) );
      WALBERLA_CHECK_FLOAT_EQUAL( gradient[1], testFuncDerivativeY( coordinates ) );
   }
}

int main( int argc, char** argv )
{
   test2D( argc, argv );
   return EXIT_SUCCESS;
}
