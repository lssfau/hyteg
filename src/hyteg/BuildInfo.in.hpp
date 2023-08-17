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

#ifdef WALBERLA_BUILD_WITH_MPI
#include <mpi.h>
#endif

#include <sstream>

#include "core/logging/Logging.h"

namespace hyteg {

inline std::string buildType()
{
   return "@HYTEG_BUILD_TYPE@";
}

inline std::string compilerFlags()
{
   return "@HYTEG_COMPILER_FLAGS@";
}

inline std::string compilerInfo()
{
   return "@HYTEG_COMPILER_INFO@";
}

inline std::string systemEndianess()
{
   return "@CMAKE_CXX_BYTE_ORDER@";
}

inline std::string mpiVersion()
{
#ifdef WALBERLA_BUILD_WITH_MPI

   // We need to do some C-style string handling here
   char mpiLibraryVersion[MPI_MAX_LIBRARY_VERSION_STRING];
   int  stringLength{ 0 };
   int  status{ 0 };

   // Run the query
   status = MPI_Get_library_version( mpiLibraryVersion, &stringLength );

   // Process results
   if ( status != MPI_SUCCESS )
   {
      WALBERLA_ABORT( "Problem with MPI_Get_library_version() detected!" );
   }

   return std::string( mpiLibraryVersion );

#else

   return std::string( "no MPI library info; sequential build!" );

#endif
}

inline void printBuildInfo()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Build info:" )
   WALBERLA_LOG_INFO_ON_ROOT( " - build type ....... " << buildType() );
   WALBERLA_LOG_INFO_ON_ROOT( " - compiler ......... " << compilerInfo() );
   WALBERLA_LOG_INFO_ON_ROOT( " - compiler flags ... " << compilerFlags() );
   WALBERLA_LOG_INFO_ON_ROOT( " - mpi version ...... " << mpiVersion() );
   WALBERLA_LOG_INFO_ON_ROOT( " - byte order ....... " << systemEndianess() );
}

} // namespace hyteg
