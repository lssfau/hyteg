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

#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

static std::string getDateTimeID()
{
   std::vector< char > cTimeString( 64 );
   WALBERLA_ROOT_SECTION()
   {
      std::time_t t;
      std::time( &t );
      std::strftime( cTimeString.data(), 64, "%F_%H-%M-%S", std::localtime( &t ) );
   }

   walberla::mpi::broadcastObject( cTimeString );

   std::string timeString( cTimeString.data() );
   return timeString;
}

template < typename FunctionType, typename ConstantOperator >
void run()
{
   const uint_t minLevel        = 1;
   const uint_t maxLevel        = 8;
   const uint_t numEdgesPerSide = 1;
   const uint_t numIterations   = 100;

   const auto dateTimeID = getDateTimeID();

   std::string discretization = "P1";
   if ( std::is_same< FunctionType, P2Function< real_t > >::value )
   {
      discretization = "P2";
   }

   const auto basename = "SpMVvsGSBenchmark_" + discretization + "_" + dateTimeID;

   WALBERLA_LOG_INFO_ON_ROOT( "SpMV vs GS benchmark" );
   WALBERLA_LOG_INFO_ON_ROOT( " - ID:             " << basename );
   WALBERLA_LOG_INFO_ON_ROOT( " - discretization: " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( " - num processes:  " << walberla::mpi::MPIManager::instance()->numProcesses() );
   WALBERLA_LOG_INFO_ON_ROOT( " - min level:      " << minLevel );
   WALBERLA_LOG_INFO_ON_ROOT( " - max level:      " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( " - num iterations: " << numIterations );

   auto meshInfo = MeshInfo::meshSymmetricCuboid(
       Point3D(  0, 0, 0  ), Point3D(  1, 1, 1  ), numEdgesPerSide, numEdgesPerSide, numEdgesPerSide );

   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );

   FunctionType src( "src", storage, minLevel, maxLevel );
   FunctionType dst( "dst", storage, minLevel, maxLevel );

   ConstantOperator A( storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( " level | unknowns " );
   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      auto unknowns = numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, l );

      WALBERLA_LOG_INFO_ON_ROOT( "     " << l << ": " << unknowns );

      src.interpolate( 42, l );
      dst.interpolate( 42, l );
   }

   auto timingTree = storage->getTimingTree();

   WALBERLA_LOG_INFO_ON_ROOT( "SpMV ..." )
   timingTree->start( "SpMV" );

   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " - level " << l )
      timingTree->start( "level " + std::to_string( l ) );
      for ( uint_t i = 0; i < numIterations; i++ )
      {
         timingTree->start( "iteration " + std::to_string( i ) );
         A.apply( src, dst, l, Inner );
         timingTree->stop( "iteration " + std::to_string( i ) );
      }
      timingTree->stop( "level " + std::to_string( l ) );
   }

   timingTree->stop( "SpMV" );

   WALBERLA_LOG_INFO_ON_ROOT( "GS ..." )
   timingTree->start( "GS" );

   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " - level " << l )
      timingTree->start( "level " + std::to_string( l ) );
      for ( uint_t i = 0; i < numIterations; i++ )
      {
         timingTree->start( "iteration " + std::to_string( i ) );
         A.smooth_gs( src, dst, l, Inner );
         timingTree->stop( "iteration " + std::to_string( i ) );
      }
      timingTree->stop( "level " + std::to_string( l ) );
   }

   timingTree->stop( "GS" );

   writeTimingTreeJSON( *timingTree, basename + ".json" );
}

int main( int argc, char* argv[] )
{
   /// create enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   run< P1Function< real_t >, P1ConstantLaplaceOperator >();
   run< P2Function< real_t >, P2ConstantLaplaceOperator >();

   return 0;
}
