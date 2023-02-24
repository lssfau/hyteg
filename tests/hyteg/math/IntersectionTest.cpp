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

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/TimingPool.h"
#include "core/logging/all.h"

#include "hyteg/geometry/Intersection.hpp"

namespace hyteg {

static void testIntersection()
{
  using walberla::uint_t;

  // triangle-circle intersection

  Point2D v1( 0, 0 );
  Point2D v2( 5, 0 );
  Point2D v3( 0, 5 );

  WALBERLA_CHECK(  circleTriangleIntersection( Point2D(   1, 1 ), 0.5, v1, v2, v3 ) );
  WALBERLA_CHECK(  circleTriangleIntersection( Point2D(   1, 1 ), 2.0, v1, v2, v3 ) );
  WALBERLA_CHECK(  circleTriangleIntersection( Point2D(   3, 3 ), 2.0, v1, v2, v3 ) );
  WALBERLA_CHECK( !circleTriangleIntersection( Point2D(  -1, 1 ), 0.5, v1, v2, v3 ) );

  // tet-sphere intersection

  Point3D tetV1( 0, 0, 0 );
  Point3D tetV2( 5, 0, 0 );
  Point3D tetV3( 0, 5, 0 );
  Point3D tetV4( 0, 0, 5 );

  WALBERLA_CHECK(  isSphereCompletelyInTetrahedron( Point3D( 1, 1, 1), 0.1, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK( !isSphereCompletelyInTetrahedron( Point3D( 1, 1, 1), 2.0, tetV1, tetV2, tetV3, tetV4 ) );

  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D( 1, 1, 1), 0.1, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D( 1, 1, 1), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D( 3, 0, 3), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D( 3, 0.1, 3), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK( !sphereTetrahedronIntersection( Point3D( 9, 1, 1), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK( !sphereTetrahedronIntersection( Point3D( -1, -1, -1), 1.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK( !sphereTetrahedronIntersection( Point3D(  0, 5, 5 ), 1.0, tetV1, tetV2, tetV3, tetV4 ) );

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hyteg::testIntersection();

  return EXIT_SUCCESS;
}
