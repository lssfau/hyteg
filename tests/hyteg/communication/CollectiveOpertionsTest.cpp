/*
* Copyright (c) 2017-2024 Dominik Thoennes.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

void runTest( const uint_t level )
{
   auto meshInfo = hyteg::MeshInfo::meshSymmetricCuboid(
       Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), uint_c( walberla::MPIManager::instance()->numProcesses() ), 1, 1 );

   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   hyteg::P1Function< double >                      func( "func", storage, level, level );
   std::function< real_t( const hyteg::Point3D& ) > random = []( const hyteg::Point3D& ) { return walberla::math::realRandom(); };
   func.interpolate( random, level );

   std::function< double( double, double ) > reduceMax = []( double reduced, double newValue ) {
      return std::max( reduced, newValue );
   };

   WALBERLA_CHECK_FLOAT_EQUAL( func.getMaxValue( level, All, false ),
                               func.reduceLocal( level, reduceMax, std::numeric_limits< double >::min(), All ) )
   WALBERLA_CHECK_FLOAT_EQUAL(
       func.getMaxValue( level, All, true ),
       func.reduceGlobal( level, reduceMax, std::numeric_limits< double >::min(), walberla::mpi::MAX, All ) )
   WALBERLA_CHECK_FLOAT_EQUAL( func.getMaxValue( level, All, true ),
                               func.reduceGlobal( level, reduceMax, std::numeric_limits< double >::min(), All ) )

   std::function< double( double, double ) > reduceMin = []( double reduced, double newValue ) {
      return std::min( reduced, newValue );
   };

   WALBERLA_CHECK_FLOAT_EQUAL( func.getMinValue( level, All, false ),
                               func.reduceLocal( level, reduceMin, std::numeric_limits< double >::max(), All ) )
   WALBERLA_CHECK_FLOAT_EQUAL(
       func.getMinValue( level, All, true ),
       func.reduceGlobal( level, reduceMin, std::numeric_limits< double >::max(), walberla::mpi::MIN, All ) )
   WALBERLA_CHECK_FLOAT_EQUAL( func.getMinValue( level, All, true ),
                               func.reduceGlobal( level, reduceMin, std::numeric_limits< double >::max(), All ) )

   std::function< double( double, double ) > maxMagnitued = []( double reduced, double newValue ) {
      return std::max( std::abs( reduced ), std::abs( newValue ) );
   };
   WALBERLA_CHECK_FLOAT_EQUAL( func.getMaxMagnitude( level, All, false ),
                               func.reduceLocal( level, maxMagnitued, std::numeric_limits< double >::min(), All ) )
   WALBERLA_CHECK_FLOAT_EQUAL(
       func.getMaxMagnitude( level, All, true ),
       func.reduceGlobal( level, maxMagnitued, std::numeric_limits< double >::min(), walberla::mpi::MAX, All ) )
   WALBERLA_CHECK_FLOAT_EQUAL( func.getMaxMagnitude( level, All, true ),
                               func.reduceGlobal( level, maxMagnitued, std::numeric_limits< double >::min(), All ) )
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   //walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::runTest( 2 );
   hyteg::runTest( 3 );
}