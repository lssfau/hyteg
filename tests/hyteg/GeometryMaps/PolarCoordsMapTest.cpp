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
#include <core/Environment.h>
#include <core/config/Config.h>
#include <core/math/Constants.h>

#include "core/timing/Timer.h"

#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/VTKWriter.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

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
  Point2D cornerUR( { rmax, 2.0*pi } );

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

    hyteg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // Check surface area of annulus
  std::function< real_t( const hyteg::Point3D& ) > one = []( const hyteg::Point3D& ) { return 1.0; };

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
      WALBERLA_CHECK_FLOAT_EQUAL( measure, 3.0*pi );
    }

  return 0;
}
