// test that the product one^T*M*one with mass matrix M and vector of ones gives area of domain
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ElementwiseOperator.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"
#include "tinyhhg_core/communication/Syncing.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::PI;

using namespace hhg;

void checkArea( std::shared_ptr<PrimitiveStorage> storage, real_t area )
{

  const size_t minLevel = 2;
  const size_t maxLevel = 4;

  hhg::P1Function< real_t > microCoordX( "microCoordX", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > microCoordY( "microCoordY", storage, minLevel, maxLevel );

  std::function< real_t( const hhg::Point3D& ) > compX = []( const hhg::Point3D& pp )
    { return pp[0]; };
  std::function< real_t( const hhg::Point3D& ) > compY = []( const hhg::Point3D& pp )
    { return pp[1]; };

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      microCoordX.interpolate( compX, lvl );
      microCoordY.interpolate( compY, lvl );

      communication::syncFunctionBetweenPrimitives( microCoordX, lvl );
      communication::syncFunctionBetweenPrimitives( microCoordY, lvl );
    }

  P1ElementwisePolarMassOperator massOp( storage, {&microCoordX,&microCoordY}, minLevel, maxLevel );

  P1Function< real_t > aux( "aux", storage, minLevel, maxLevel );
  P1Function< real_t > vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );
  std::function< real_t( const Point3D& ) > ones = []( const Point3D& ) { return 1.0; };

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      vecOfOnes.interpolate( ones, lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dot( aux, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": measure = " << std::scientific << measure );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, area );
    }
}

int main(int argc, char **argv)
{
  walberla::debug::enterTestMode();

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  // Test with rectangle 1
  WALBERLA_LOG_INFO_ON_ROOT( "Testing with RECTANGLE 1 (full annulus)" );
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {1.0, 0.0} ), Point2D( {3.0, 2*PI} ),
                                               MeshInfo::CRISSCROSS, 1, 1 );
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);
  checkArea( storage, 8.0*PI );

  // Test with rectangle 2
  WALBERLA_LOG_INFO_ON_ROOT( "Testing with RECTANGLE 2 (partial annulus)" );
  meshInfo = MeshInfo::meshRectangle( Point2D( {1.0, 0.3*PI} ), Point2D( {3.0, 0.55*PI} ),
                                               MeshInfo::CRISSCROSS, 1, 1 );
  SetupPrimitiveStorage setupStorage2(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage2 = std::make_shared<PrimitiveStorage>(setupStorage2);
  checkArea( storage2, PI );

  return EXIT_SUCCESS;
}
