/*
 * Copyright (c) 2017-2020 Nils Kohl.
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
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"
#include "core/math/Random.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/forms/form_fenics_generated/p1_tet_diffusion.h"
#include "hyteg/dataexport/VTKOutput.hpp"

using walberla::real_t;
using walberla::real_c;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

void TestFunctionMemoryAllocation()
{
   auto meshInfo = MeshInfo::fromGmshFile("../../data/meshes/3D/cube_24el.msh");
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P1Function< real_t > x( "x", storage, 3, 3 );

   auto globalMemoryBeforeDeletion =  FunctionMemory< real_t >::getGlobalAllocatedMemoryInBytes();
   WALBERLA_LOG_INFO_ON_ROOT( "Initial memory on level 3: " << globalMemoryBeforeDeletion );

   // now delete cell memory on level 3
   WALBERLA_CHECK_GREATER( storage->getNumberOfLocalCells(), 0 );
   auto cell = storage->getCell( storage->getCellIDs().front() );
   x.deleteMemory( 3, *cell );

   // calculate memory deleted
   auto memoryDeleted = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) * levelinfo::num_microvertices_per_cell( 3 ) * sizeof( real_t );
   auto globalMemoryAfterDeletion = FunctionMemory< real_t >::getGlobalAllocatedMemoryInBytes();
   WALBERLA_LOG_INFO_ON_ROOT( "Memory after deletion of 1 cell per process on level 3: " << globalMemoryAfterDeletion );
   WALBERLA_CHECK_EQUAL( globalMemoryAfterDeletion, globalMemoryBeforeDeletion - memoryDeleted );

   // now allocate memory on level 4
   x.allocateMemory( 4, *cell );

   // calculate memory allocated
   auto memoryAllocated = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) * levelinfo::num_microvertices_per_cell( 4 ) * sizeof( real_t );
   auto globalMemoryAfterReallocation = FunctionMemory< real_t >::getGlobalAllocatedMemoryInBytes();
   WALBERLA_LOG_INFO_ON_ROOT( "Memory after re-allocation of 1 cell per process on level 4: " << globalMemoryAfterReallocation );
   WALBERLA_CHECK_EQUAL( globalMemoryAfterReallocation, globalMemoryAfterDeletion + memoryAllocated );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   TestFunctionMemoryAllocation();
   return 0;
}
