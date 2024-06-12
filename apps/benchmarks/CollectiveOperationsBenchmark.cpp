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

#include <hyteg/mesh/MeshInfo.hpp>
#include <hyteg/p1functionspace/P1Function.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/types/PointND.hpp"

namespace hyteg {

void runTest()
{
   auto meshInfo = hyteg::MeshInfo::meshSymmetricCuboid(
       Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), uint_c( walberla::MPIManager::instance()->numProcesses() ), 1, 1 );

   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   hyteg::P1Function< double > func( "func", storage, 4, 4 );

   timingTree->start( "getMaxValue" );
   double dummy = 0.0;
   for ( int i = 0; i < 100; ++i )
   {
      dummy += func.getMaxValue( 4, All, true );
   }
   timingTree->stop( "getMaxValue" );
   timingTree->start( "getMaxValueWoMPI" );
   for ( int i = 0; i < 100; ++i )
   {
      dummy += func.getMaxValue( 4, All, false );
   }
   timingTree->stop( "getMaxValueWoMPI" );

   std::function< double( double, double ) > reduceMax = []( double reduced, double newValue ) {
      return std::max( reduced, newValue );
   };
   timingTree->start( "getNewMaxValue" );
   for ( int i = 0; i < 100; ++i )
   {
      func.reduceGlobal( 4, reduceMax, std::numeric_limits< double >::lowest(), walberla::mpi::MAX, All );
   }
   timingTree->stop( "getNewMaxValue" );

   timingTree->start( "getNewMaxValueWoMPI" );
   for ( int i = 0; i < 100; ++i )
   {
      func.reduceLocal( 4, reduceMax, std::numeric_limits< double >::lowest(), All );
   }
   timingTree->stop( "getNewMaxValueWoMPI" );

   auto tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( dummy )
   WALBERLA_LOG_INFO_ON_ROOT( tt );
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   //walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::runTest();
}