/*
 * Copyright (c) 2021 Marcus Mohr
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
#include <vector>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/functions/FunctionMultiStore.hpp"

using walberla::uint_t;

using namespace hyteg;

// --------------------------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  uint_t theLevel = 3;

  // Generate mesh around origin
  MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
  meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CROSS, 1, 1 );

  // Generate primitives
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>( setupStorage );

  // =============
  //  P1Function 
  // =============
  hyteg::P1Function< double > dFunc1( "#1", storage, theLevel, theLevel );
  hyteg::P1Function< double > dFunc2( "#2", storage, theLevel, theLevel );
  // hyteg::P1Function< float > fFunc1( "#3", storage, theLevel, theLevel );
  hyteg::P1Function< int32_t > iFunc1( "#4", storage, theLevel, theLevel );
  hyteg::P1Function< int64_t > lFunc1( "#5", storage, theLevel, theLevel );

  FunctionMultiStore< P1Function > ms;

  ms.push_back( dFunc1 );
  ms.push_back( dFunc2 );
  ms.push_back( iFunc1 );
  ms.push_back( lFunc1 );

  WALBERLA_LOG_INFO_ON_ROOT( "FunctionMultiStore holds " << ms.size() << " P1Functions" ); 

  // ==================
  //  P2VectorFunction 
  // ==================
  hyteg::P2VectorFunction< double > dVecFunc1( "#1", storage, theLevel, theLevel );
  hyteg::P2VectorFunction< double > dVecFunc2( "#2", storage, theLevel, theLevel );
  // hyteg::P2VectorFunction< float > fVecFunc1( "#3", storage, theLevel, theLevel );
  hyteg::P2VectorFunction< int32_t > iVecFunc1( "#4", storage, theLevel, theLevel );
  hyteg::P2VectorFunction< int64_t > lVecFunc1( "#5", storage, theLevel, theLevel );
  hyteg::P2VectorFunction< int64_t > lVecFunc2( "#6", storage, theLevel, theLevel );

  FunctionMultiStore< P2VectorFunction > ms2;

  ms2.push_back( dVecFunc1 );
  ms2.push_back( dVecFunc2 );
  ms2.push_back( iVecFunc1 );
  ms2.push_back( lVecFunc1 );
  ms2.push_back( lVecFunc2 );

  WALBERLA_LOG_INFO_ON_ROOT( "FunctionMultiStore holds " << ms2.size() << " P2VectorFunctions" ); 

  // try obtaining a specific vector from the store
  auto dFuncs = ms2.getFunctions< double >();
  WALBERLA_LOG_INFO_ON_ROOT( " - of these " << dFuncs.size() << " have value type 'double'" ); 
  auto iFuncs = ms2.getFunctions< int32_t >();
  WALBERLA_LOG_INFO_ON_ROOT( " - of these " << iFuncs.size() << " have value type 'int32_t'" ); 

  return EXIT_SUCCESS;
}
