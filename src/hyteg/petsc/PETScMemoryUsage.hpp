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

#pragma once

#ifdef HYTEG_BUILD_WITH_PETSC

#include "core/mpi/Reduce.h"

#include "hyteg/petsc/PETScWrapper.hpp"

namespace hyteg {

/// \brief Returns the current memory usage of the entire application on the current process in bytes.
///
/// This is not related to PETSc but uses the PETSc function PetscMemoryGetCurrentUsage() which
/// appears to not only include PETSc's memory usage but the total.
double getCurrentMemoryUsage()
{
   PetscLogDouble mem;
   PetscMemoryGetCurrentUsage( &mem );
   return static_cast< double >( mem );
}

/// \brief Prints the information on the current memory usage of the entire application.
///
/// Involves global reduction.
///
/// This is not related to PETSc but uses the PETSc function PetscMemoryGetCurrentUsage() which
/// appears to not only include PETSc's memory usage but the total.
void printCurrentMemoryUsage()
{
   double locallyAllocatedMemory = getCurrentMemoryUsage();

   const double globalActualAllocatedMemory =
       walberla::mpi::allReduce( locallyAllocatedMemory, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() ) /
       1e+09;
   const double minActualAllocatedMemory =
       walberla::mpi::allReduce( locallyAllocatedMemory, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() ) /
       1e+09;
   const double maxActualAllocatedMemory =
       walberla::mpi::allReduce( locallyAllocatedMemory, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() ) /
       1e+09;

   // print it all...
   WALBERLA_LOG_INFO_ON_ROOT( "========================= Memory Usage Info =========================" );
   WALBERLA_LOG_INFO_ON_ROOT( "                       +--------------+--------------+--------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "                       |          sum |          min |          max |" );
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------+--------------+--------------+--------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( " allocated memory (GB) | "
                              << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << globalActualAllocatedMemory << " | "
                              << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << minActualAllocatedMemory << " | "
                              << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << maxActualAllocatedMemory << " | " );
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------+--------------+--------------+--------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "=====================================================================" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

} // namespace hyteg

#endif