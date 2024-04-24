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
#include <hyteg/primitivestorage/Visualization.hpp>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::uint_t;

namespace hyteg {


static void testFunctionMemorySerialization()
{
  const uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
  const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  const std::string meshFileName = "../../meshes/bfs_126el.msh";

  const uint_t minLevel = 2;
  const uint_t maxLevel = 3;

  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< PrimitiveStorage > storage( new PrimitiveStorage( setupStorage ) );

  P1Function< real_t > x("x", storage, minLevel, maxLevel);
  P1ConstantLaplaceOperator A(storage, minLevel, maxLevel);

  VTKOutput vtkOutputBefore("../../output/", "function_memory_serialization_test_data_before_migration", storage);
  vtkOutputBefore.add( x );

  VTKOutput vtkOutputAfter("../../output/", "function_memory_serialization_test_data_after_migration", storage);
  vtkOutputAfter.add( x );

  std::function<real_t(const hyteg::Point3D&)> gradient = [](const hyteg::Point3D& xx) { return xx[0]; };

  for ( uint_t level = minLevel; level <= maxLevel; level++ )
  {
    x.interpolate( gradient, level );
  }

  writeDomainPartitioningVTK( storage, "../../output/", "function_memory_serialization_test_domain_before_migration" );
  vtkOutputBefore.write( maxLevel );

  WALBERLA_LOG_INFO( "Number of local primitives (before migration): " << storage->getNumberOfLocalPrimitives() );

  std::map< hyteg::PrimitiveID, uint_t > primitivesToMigrate;
  std::vector< PrimitiveID > localPrimitiveIDs;
  storage->getPrimitiveIDs( localPrimitiveIDs );
  for ( const auto & id : localPrimitiveIDs )
  {
    primitivesToMigrate[ id ] = (rank + numProcesses / 2) % numProcesses;
  }

  const auto numReceivingPrimitives = getNumReceivingPrimitives( primitivesToMigrate );
  const MigrationInfo migrationInfo( primitivesToMigrate, numReceivingPrimitives );
  storage->migratePrimitives( migrationInfo );

  WALBERLA_LOG_INFO( "Number of local primitives (after migration): " << storage->getNumberOfLocalPrimitives() );

  writeDomainPartitioningVTK( storage, "../../output/", "function_memory_serialization_test_domain_after_migration" );
  vtkOutputAfter.write( maxLevel );

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();
  hyteg::testFunctionMemorySerialization();

  return EXIT_SUCCESS;
}
