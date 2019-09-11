/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
// -------------------------------------------------------------------------
//  Test that Jacobi, Gauss-Seidel and SOR smoothing give identical results
//  for the P1ConstantOperator and the P1ElementwiseOperator
// -------------------------------------------------------------------------

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ElementwiseOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

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
  hyteg::P1Function< real_t > microCoordX( "microCoordX", storage, minLevel, maxLevel );
  hyteg::P1Function< real_t > microCoordY( "microCoordY", storage, minLevel, maxLevel );

  std::function< real_t( const hyteg::Point3D& ) > compX = []( const hyteg::Point3D& pp )
    { return pp[0]; };
  std::function< real_t( const hyteg::Point3D& ) > compY = []( const hyteg::Point3D& pp )
    { return pp[1]; };

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      microCoordX.interpolate( compX, lvl );
      microCoordY.interpolate( compY, lvl );

      communication::syncFunctionBetweenPrimitives( microCoordX, lvl );
      communication::syncFunctionBetweenPrimitives( microCoordY, lvl );
    }

  // setup two Laplacians
  P1ConstantLaplaceOperator lapOpCO( storage, minLevel, maxLevel );
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

      difference.assign( {1.0}, {zeros}, lvl );
      difference.add( { 1.0, -1.0 }, { smoothCO, smoothEL }, lvl, All );
      value = sqrt( difference.dotGlobal( difference, lvl, All ) );

      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": JAC value = " << std::scientific << value );
      WALBERLA_CHECK_FLOAT_EQUAL( value, 0.0 );

      // Gauss-Seidel smoothing
      smoothCO.interpolate( quadratic, lvl, All );
      smoothEL.interpolate( quadratic, lvl, All );

      lapOpCO.smooth_gs( smoothCO, zeros, lvl, All );
      lapOpEL.smooth_gs( smoothEL, zeros, lvl, All );

      difference.assign( {1.0}, {zeros}, lvl );
      difference.add( { 1.0, -1.0 }, { smoothCO, smoothEL }, lvl, All );
      value = sqrt( difference.dotGlobal( difference, lvl, All ) );

      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": GS value = " << std::scientific << value );
      WALBERLA_CHECK_FLOAT_EQUAL( value, 0.0 );

      // SOR smoothing
      smoothCO.interpolate( quadratic, lvl, All );
      smoothEL.interpolate( quadratic, lvl, All );

      lapOpCO.smooth_sor( smoothCO, zeros, real_t(1.25), lvl, All );
      lapOpEL.smooth_sor( smoothEL, zeros, real_t(1.25), lvl, All );

      difference.assign( {1.0}, {zeros}, lvl );
      difference.add( { 1.0, -1.0 }, { smoothCO, smoothEL }, lvl, All );
      value = sqrt( difference.dotGlobal( difference, lvl, All ) );

      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": SOR value = " << std::scientific << value );
      WALBERLA_CHECK_FLOAT_EQUAL( value, 0.0 );
    }

  return EXIT_SUCCESS;
}
