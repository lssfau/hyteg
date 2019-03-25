#include <core/Environment.h>
#include <core/config/Config.h>

#include "core/timing/Timer.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/geometry/PolarCoordsMap.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1VariableOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::PI;

using namespace hhg;

int main( int argc, char* argv[] )
{
  // Setup enviroment
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  const size_t level = 2;

  // Generate annulus mesh in polar coordinates, so it's a rectangle
  real_t rmin = 1.0;
  real_t rmax = 2.0;

  Point2D cornerLL( { rmin, 0.0 } );
  Point2D cornerUR( { rmax, 2.0*PI } );

  MeshInfo meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CROSS, 1, 6 );
  WALBERLA_LOG_INFO_ON_ROOT( " *** Using Inline Mesher" );

  // Prepare storage and set geometry mapping
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  for( auto it : setupStorage.getFaces() )
    {
      setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
    }

  for( auto it : setupStorage.getEdges() )
    {
      setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
    }

  for( auto it : setupStorage.getVertices() )
    {
      setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
    }

  hhg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // Check surface area of annulus
  std::function< real_t( const hhg::Point3D& ) > one = []( const hhg::Point3D& ) { return 1.0; };

  uint_t minLevel = level;
  uint_t maxLevel = level;

  P1BlendingMassOperator massOp( storage, minLevel, maxLevel );

  P1Function< real_t > aux( "aux", storage, minLevel, maxLevel );
  P1Function< real_t > vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      vecOfOnes.interpolate( one, lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dotGlobal( aux, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "annulus area = " << std::scientific << measure );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, 3.0*PI );
    }

  return 0;
}
