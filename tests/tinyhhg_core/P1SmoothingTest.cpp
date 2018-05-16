// -------------------------------------------------------------------------
//  Test that Jacobi, Gauss-Seidel and SOR smoothing give identical results
//  for the P1ConstantOperator and the P1ElementwiseOperator
// -------------------------------------------------------------------------

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ElementwiseOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"
#include "tinyhhg_core/communication/Syncing.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::PI;

using namespace hhg;

int main(int argc, char **argv)
{
  walberla::debug::enterTestMode();

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  WALBERLA_LOG_INFO_ON_ROOT( "Comparing Smoothing for P1Constant and P1Elementwise Operators" );

  // setup mesh and storage stuff
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {-1.0, -1.0} ), Point2D( {1.0, 1.0} ),
                                               MeshInfo::CRISSCROSS, 1, 1 );
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>( setupStorage );

  // how many levels
  const size_t minLevel = 2;
  const size_t maxLevel = 4;

  // prepare micro-coordinates for elementwise operator
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

      syncFunctionBetweenPrimitives( &microCoordX, lvl );
      syncFunctionBetweenPrimitives( &microCoordY, lvl );
    }

  // setup two Laplacians
  P1LaplaceOperator lapOpCO( storage, minLevel, maxLevel );
  P1ElementwiseLaplaceOperator lapOpEL( storage, {&microCoordX,&microCoordY}, minLevel, maxLevel );

  // setup auxilliary P1Functions
  P1Function< real_t > zeros( "zeros", storage, minLevel, maxLevel );
  P1Function< real_t > difference( "difference", storage, minLevel, maxLevel );
  P1Function< real_t > initialCO( "initCO", storage, minLevel, maxLevel );
  P1Function< real_t > initialEL( "initEL", storage, minLevel, maxLevel );
  P1Function< real_t > smoothCO( "smoothedCO", storage, minLevel, maxLevel );
  P1Function< real_t > smoothEL( "smoothedEL", storage, minLevel, maxLevel );

  std::function< real_t( const Point3D& ) > linear = []( const Point3D &pp ) { return pp[0] + 2.0*pp[1]; };
  std::function< real_t( const Point3D& ) > quadratic = []( const Point3D &pp ) { return pp[0]*pp[0] + pp[1]*pp[1]; };
  std::function< real_t( const Point3D& ) > one = []( const Point3D & ) { return real_t(1.0); };

  real_t value = real_t(0.0);

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {

      // Jacobi smoothing
      initialCO.interpolate( quadratic, lvl, All );
      initialEL.interpolate( quadratic, lvl, All );

      lapOpCO.smooth_jac( smoothCO, zeros, initialCO, lvl, All );
      lapOpEL.smooth_jac( smoothEL, zeros, initialEL, lvl, All );

      difference.assign( {1.0}, {&zeros}, lvl );
      difference.add( { 1.0, -1.0 }, { &smoothCO, &smoothEL }, lvl, All );
      value = sqrt( difference.dotGlobal( difference, lvl, All ) );

      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": JAC value = " << std::scientific << value );
      WALBERLA_CHECK_FLOAT_EQUAL( value, 0.0 );

      // Gauss-Seidel smoothing
      smoothCO.interpolate( quadratic, lvl, All );
      smoothEL.interpolate( quadratic, lvl, All );

      lapOpCO.smooth_gs( smoothCO, zeros, lvl, All );
      lapOpEL.smooth_gs( smoothEL, zeros, lvl, All );

      difference.assign( {1.0}, {&zeros}, lvl );
      difference.add( { 1.0, -1.0 }, { &smoothCO, &smoothEL }, lvl, All );
      value = sqrt( difference.dotGlobal( difference, lvl, All ) );

      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": GS value = " << std::scientific << value );
      WALBERLA_CHECK_FLOAT_EQUAL( value, 0.0 );

      // SOR smoothing
      smoothCO.interpolate( quadratic, lvl, All );
      smoothEL.interpolate( quadratic, lvl, All );

      lapOpCO.smooth_sor( smoothCO, zeros, real_t(1.25), lvl, All );
      lapOpEL.smooth_sor( smoothEL, zeros, real_t(1.25), lvl, All );

      difference.assign( {1.0}, {&zeros}, lvl );
      difference.add( { 1.0, -1.0 }, { &smoothCO, &smoothEL }, lvl, All );
      value = sqrt( difference.dotGlobal( difference, lvl, All ) );

      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": SOR value = " << std::scientific << value );
      WALBERLA_CHECK_FLOAT_EQUAL( value, 0.0 );
    }

  return EXIT_SUCCESS;
}
