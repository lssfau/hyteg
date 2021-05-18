/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Marcus Mohr.
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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/composites/P2P1TaylorHoodBlockFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

// Perform a bacis compile and apply test for BlockOperator class

using namespace hyteg;

template < typename vfType, typename opType >
static void runTest( bool beVerbose, std::string tag, std::string opName )
{
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "checking " << opName );
   }

   // not using this currently
   WALBERLA_UNUSED( tag );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "=======================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing BlockOperator" );
   WALBERLA_LOG_INFO_ON_ROOT( "=======================" );

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t minLevel = 1;
   uint_t maxLevel = 1;

   typedef P2P1TaylorHoodBlockFunction< real_t > thType;
   BlockOperator< thType, thType > blockOper( storage, minLevel, maxLevel, 2, 3 );
   
   return EXIT_SUCCESS;
}
