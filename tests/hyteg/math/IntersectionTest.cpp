
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

  Point2D v1( {0, 0} );
  Point2D v2( {5, 0} );
  Point2D v3( {0, 5} );

  WALBERLA_CHECK(  circleTriangleIntersection( Point2D({  1, 1 }), 0.5, v1, v2, v3 ) );
  WALBERLA_CHECK(  circleTriangleIntersection( Point2D({  1, 1 }), 2.0, v1, v2, v3 ) );
  WALBERLA_CHECK(  circleTriangleIntersection( Point2D({  3, 3 }), 2.0, v1, v2, v3 ) );
  WALBERLA_CHECK( !circleTriangleIntersection( Point2D({ -1, 1 }), 0.5, v1, v2, v3 ) );

  // tet-sphere intersection

  Point3D tetV1( {0, 0, 0} );
  Point3D tetV2( {5, 0, 0} );
  Point3D tetV3( {0, 5, 0} );
  Point3D tetV4( {0, 0, 5} );

  WALBERLA_CHECK(  isSphereCompletelyInTetrahedron( Point3D({1, 1, 1}), 0.1, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK( !isSphereCompletelyInTetrahedron( Point3D({1, 1, 1}), 2.0, tetV1, tetV2, tetV3, tetV4 ) );

  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D({1, 1, 1}), 0.1, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D({1, 1, 1}), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D({3, 0, 3}), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK(  sphereTetrahedronIntersection( Point3D({3, 0.1, 3}), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK( !sphereTetrahedronIntersection( Point3D({9, 1, 1}), 2.0, tetV1, tetV2, tetV3, tetV4 ) );
  WALBERLA_CHECK( !sphereTetrahedronIntersection( Point3D({-1, -1, -1}), 1.0, tetV1, tetV2, tetV3, tetV4 ) );

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
