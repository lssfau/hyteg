/*
 * Copyright (c) 2017-2019 Nils Kohl.
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
#include "core/math/all.h"
#include "core/timing/all.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

namespace hyteg {

static void testParallelRoundRobin( const MeshInfo & mesh )
{
  const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
  const uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  SetupPrimitiveStorage setupStorage( mesh, numProcesses );
  const auto numGlobalPrimitives = setupStorage.getNumberOfPrimitives();

  // balance with round robin to test if both algorithms produce same results
  loadbalancing::roundRobin( setupStorage );  

  auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // get current distribution
  std::vector< PrimitiveID > distributionBefore, distributionAfter;
  storage->getPrimitiveIDs( distributionBefore );
  // balance
  loadbalancing::distributed::roundRobin( *storage );
  // compare distribution
  storage->getPrimitiveIDs( distributionAfter );
  WALBERLA_CHECK( distributionBefore == distributionAfter );
  // redistribute twice and double check
  if ( numProcesses > 1 )
  {
    loadbalancing::distributed::roundRobin( *storage, 1 );
    storage->getPrimitiveIDs( distributionAfter );
    WALBERLA_CHECK( distributionBefore != distributionAfter );
    if ( rank == 0 )
    {
      WALBERLA_CHECK( distributionAfter.size() == numGlobalPrimitives );
    }
    else
    {
      WALBERLA_CHECK( distributionAfter.empty() )
    }
  }
  // balance
  loadbalancing::distributed::roundRobin( *storage );
  // compare distribution
  storage->getPrimitiveIDs( distributionAfter );
  WALBERLA_CHECK( distributionBefore == distributionAfter );
}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testParallelRoundRobin( hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" ) );

   return EXIT_SUCCESS;
}

