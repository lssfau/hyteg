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
// test that the product one^T*M*one with mass matrix M and vector of ones gives area of domain
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ElementwiseOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "hyteg/communication/Syncing.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

void checkArea( std::shared_ptr<PrimitiveStorage> storage, real_t area )
{

  const size_t minLevel = 2;
  const size_t maxLevel = 4;

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

  P1ElementwisePolarMassOperator massOp( storage, {&microCoordX, &microCoordY}, minLevel, maxLevel );

  P1Function< real_t > aux( "aux", storage, minLevel, maxLevel );
  P1Function< real_t > vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );
  std::function< real_t( const Point3D& ) > ones = []( const Point3D& ) { return 1.0; };

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      vecOfOnes.interpolate( ones, lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dotGlobal( aux, lvl );
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
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {1.0, 0.0} ), Point2D( {3.0, 2*pi} ),
                                               MeshInfo::CRISSCROSS, 1, 1 );
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);
  checkArea( storage, 8.0*pi );

  // Test with rectangle 2
  WALBERLA_LOG_INFO_ON_ROOT( "Testing with RECTANGLE 2 (partial annulus)" );
  meshInfo = MeshInfo::meshRectangle( Point2D( {1.0, 0.3*pi} ), Point2D( {3.0, 0.55*pi} ),
                                               MeshInfo::CRISSCROSS, 1, 1 );
  SetupPrimitiveStorage setupStorage2(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage2 = std::make_shared<PrimitiveStorage>(setupStorage2);
  checkArea( storage2, pi );

  return EXIT_SUCCESS;
}
