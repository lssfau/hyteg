/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

#include "hyteg/geometry/Polygons.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/debug/CheckFunctions.h"


using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::real_c;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   Point3D center( 1, 0, 42 );

   std::vector< Point3D > polygon = {
       Point3D(  2, 0, 42  ),
       Point3D(  2, 1, 42  ),
       Point3D(  1, 2, 42  ),
       Point3D(  0, 2, 42  ),
       Point3D(  0, -2, 42  ),
       Point3D(  1, -2, 42  ),
       Point3D(  2, -1, 42  ),
   };

   real_t r, angle;

   fractionalRadiusToPolygonBoundary( Point3D(  real_c(0.5), 0, 42  ), center, polygon, r, angle );
   WALBERLA_CHECK_FLOAT_EQUAL( r, real_c( real_c(0.5) ) );

   fractionalRadiusToPolygonBoundary( center, center, polygon, r, angle );
   WALBERLA_CHECK_FLOAT_EQUAL( r, real_c( 0 ) );

   fractionalRadiusToPolygonBoundary( Point3D(  real_c(1.5), 0, 42  ), center, polygon, r, angle );
   WALBERLA_CHECK_FLOAT_EQUAL( r, real_c( real_c(0.5) ) );

   for ( const auto & bp : polygon )
   {
      fractionalRadiusToPolygonBoundary( bp, center, polygon, r, angle );
      WALBERLA_CHECK_FLOAT_EQUAL( r, real_c( 1 ), "Point: " << bp );

      auto p = center + real_c(0.3) * (bp - center);
      fractionalRadiusToPolygonBoundary( p, center, polygon, r, angle );
      WALBERLA_CHECK_FLOAT_EQUAL( r, real_c( real_c(0.3) ), "Point: " << p );
   }

   std::vector< Point3D > polygon2 = {
       Point3D(  real_c(1.29493), real_c(0.114682), 0  ),
       Point3D(  real_c(1.29493), real_c(0.114682), real_c(0.175)  ),
       Point3D(  real_c(0.896491), real_c(0.0793954), real_c(0.525)  ),
       Point3D(  real_c(0.697271), real_c(0.061752), real_c(0.525)  ),
       Point3D(  real_c(0.697271), real_c(0.061752), -real_c(0.525)  ),
       Point3D(  real_c(0.896491), real_c(0.0793954), -real_c(0.525)  ),
       Point3D(  real_c(1.29493), real_c(0.114682), -real_c(0.175)  ),
   };

   Point3D center2(  real_c(0.796881), real_c(0.0705737), 0  );
   Point3D p2(  real_c(0.846686), real_c(0.0749846), -real_c(2.08167e-17)  );

   fractionalRadiusToPolygonBoundary( p2, center2, polygon2, r, angle );
   WALBERLA_LOG_DEVEL_ON_ROOT( "" << r << ", " << angle );

   return 0;
}
