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

#include <cassert>
#include <iostream>

#if defined(_MSC_VER)
   #include <sys/timeb.h>
   #include <windows.h>
   #include <psapi.h>
   #pragma comment(lib, "psapi.lib")
#else
   #include <sys/time.h>
   #include <sys/resource.h>
#endif
#include <vector>

#include "core/Abort.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/mpi/Reduce.h"

#include "hyteg/petsc/PETScWrapper.hpp"

namespace hyteg {

enum class MemoryUsageDeterminationType
{
   /// uses the C function getrusage()
   C_RUSAGE,

   /// uses the PETSc function PetscMemoryGetCurrentUsage()
   /// which appears to not only include PETSc's memory usage but the total.
   PETSC
};

/// \brief Returns the current memory usage of the entire application on the current process in bytes.
///
/// This is not related to PETSc but
///
/// \param type How the memory usage is determined.
///
/// \return Memory usage in bytes.
inline double getCurrentMemoryUsage( MemoryUsageDeterminationType type = MemoryUsageDeterminationType::C_RUSAGE )
{
   if ( type == MemoryUsageDeterminationType::C_RUSAGE )
   {
      #if defined(_MSC_VER)
             PROCESS_MEMORY_COUNTERS info;
            GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
            return static_cast< double >(info.PeakWorkingSetSize);
      #else
         struct rusage usage;
         int           gru = getrusage( RUSAGE_SELF, &usage );
         WALBERLA_ASSERT( gru == 0, "getrusage() returned an error." );

         // printf("ru_maxrss: %ld (maximum resident set size -- MB)\n",usage.ru_maxrss / 1024);
         // assert(usage.ru_maxrss / 1024 > megabytes );

         return static_cast< double >( usage.ru_maxrss * 1024 );
      #endif
   }
   else if ( type == MemoryUsageDeterminationType::PETSC )
   {
#ifdef HYTEG_BUILD_WITH_PETSC
      PetscLogDouble mem;
      PetscMemoryGetCurrentUsage( &mem );
      return static_cast< double >( mem );
#else
      WALBERLA_ABORT( "Memory usage with PETSc requested, but HyTeG was not built with PETSc." )
#endif
   }
   else
   {
      WALBERLA_ABORT( "Invalid memory usage determination type." )
   }
}

/// \brief Prints and returns the information on the current memory usage of the entire application.
/// Order: sum, min, max
/// In gigabytes (GB, 1e9 bytes).
///
/// Involves global reduction.
inline std::tuple< double, double, double >
    printAndGetCurrentMemoryUsage( MemoryUsageDeterminationType type = MemoryUsageDeterminationType::C_RUSAGE )
{
   std::string method;
   switch ( type )
   {
   case MemoryUsageDeterminationType::C_RUSAGE:
      method = "rusage";
      break;
   case MemoryUsageDeterminationType::PETSC:
      method = "PETSc";
      break;
   }

   double locallyAllocatedMemory = getCurrentMemoryUsage( type );

   const double globalActualAllocatedMemory =
       walberla::mpi::allReduce( locallyAllocatedMemory, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() ) /
       1e+09;
   const double minActualAllocatedMemory =
       walberla::mpi::allReduce( locallyAllocatedMemory, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() ) /
       1e+09;
   const double maxActualAllocatedMemory =
       walberla::mpi::allReduce( locallyAllocatedMemory, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() ) /
       1e+09;

   // std::make_tuple( globalActualAllocatedMemory, minActualAllocatedMemory, maxActualAllocatedMemory );

   WALBERLA_LOG_INFO_ON_ROOT( "========================= Memory Usage Info =========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " method: " << method );
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

   return std::make_tuple( globalActualAllocatedMemory, minActualAllocatedMemory, maxActualAllocatedMemory );
}

/// \brief Prints the information on the current memory usage of the entire application.
///
/// Involves global reduction.
inline void printCurrentMemoryUsage( MemoryUsageDeterminationType type = MemoryUsageDeterminationType::C_RUSAGE )
{
   printAndGetCurrentMemoryUsage( type );
}

} // namespace hyteg