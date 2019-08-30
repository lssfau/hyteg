#include <vector>
#include <functional>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"


#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/PrimitiveID.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"

using walberla::real_t;

namespace hyteg {

static void testP2Swap()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;
   const uint_t testLevel = 3;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/quad_8el.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   auto timingTree = std::make_shared< walberla::WcTimingTree >();
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   P2Function< real_t > x( "x", storage, minLevel, maxLevel );
   P2Function< real_t > y( "y", storage, minLevel, maxLevel );
   P2Function< real_t > xCorrect( "xCorrect", storage, minLevel, maxLevel );
   P2Function< real_t > yCorrect( "yCorrect", storage, minLevel, maxLevel );
   P2Function< real_t > err( "err", storage, minLevel, maxLevel );
   

   // Interpolate

   std::function< real_t( const hyteg::Point3D& ) > funcX  = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + std::sin( xx[0] ) ) * ( 2.0 + xx[1] ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > funcY = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + ( xx[0] / 5.0 ) ) * ( 42.0 + xx[1] ) );
   };

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
     x.interpolate( funcX, level );
     y.interpolate( funcY, level );
     if ( level == testLevel )
     {
       xCorrect.interpolate( funcY, level );
       yCorrect.interpolate( funcX, level );
     }
     else
     {
       xCorrect.interpolate( funcX, level );
       yCorrect.interpolate( funcY, level );
     }
   }

   x.swap( y, testLevel );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
     err.assign( {1.0, -1.0}, {x, xCorrect}, level );
     const real_t errorX = err.dotGlobal( err, level );
     WALBERLA_LOG_DEVEL_ON_ROOT( errorX )
     WALBERLA_CHECK_LESS( errorX, 1e-16 );
     err.assign( {1.0, -1.0}, {y, yCorrect}, level );
     const real_t errorY = err.dotGlobal( err, level );
     WALBERLA_LOG_DEVEL_ON_ROOT( errorY )
     WALBERLA_CHECK_LESS( errorY, 1e-16 );
   }

}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2Swap();

   return EXIT_SUCCESS;
}
