/*
 * Copyright (c) 2022 Marcus Mohr.
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

// Check that we get neither compile- nor run-time errors, when
// trying to print meta information on the build

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/BuildInfo.hpp"
#include "hyteg/Git.hpp"

using namespace hyteg;
using namespace hyteg::buildinfo;

// check for the unlikely event that the byte order specified during
// compilation is different to the one at executation time
void checkEndianess()
{
   // do some classical C magic
   uint16_t i = 1u;
   uint8_t* c = (uint8_t*) &i;

   std::string byteOrder = ( *c ) ? "LITTLE_ENDIAN" : "BIG_ENDIAN";
   WALBERLA_CHECK( byteOrder == systemEndianess(), "Byte order inconsistent between compile and execute architecture!" );
}

void checkFPType()
{
   std::string fpReported = fpType();

#ifdef WALBERLA_DOUBLE_ACCURACY
   std::string fpExpected = "double";
#else
   std::string fpExpected = "float";
#endif

   WALBERLA_CHECK( fpReported == fpExpected,
                   "Inconsistency in fpType detected! '" << fpReported << "' != '" << fpExpected << "'" );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::string separator{ "--------------------------------------------------" };

   printBuildInfo();
   WALBERLA_LOG_INFO_ON_ROOT( separator );
   WALBERLA_LOG_INFO_ON_ROOT( "buildType() returned ......... " << buildType() );
   WALBERLA_LOG_INFO_ON_ROOT( "compilerInfo() returned ...... " << compilerInfo() );
   WALBERLA_LOG_INFO_ON_ROOT( "compilerFlags() returned ..... " << compilerFlags() );
   WALBERLA_LOG_INFO_ON_ROOT( "mpiVersion() returned ........ " << mpiVersion() );
   WALBERLA_LOG_INFO_ON_ROOT( "systemEndianess() returned ... " << systemEndianess() );
   checkEndianess();
   WALBERLA_LOG_INFO_ON_ROOT( "fpType() returned ............ " << fpType() );
   checkFPType();

   WALBERLA_LOG_INFO_ON_ROOT( separator );
   printGitInfo();
   WALBERLA_LOG_INFO_ON_ROOT( separator );
   WALBERLA_LOG_INFO_ON_ROOT( "gitSHA1() returned ........... " << gitSHA1() );
   WALBERLA_LOG_INFO_ON_ROOT( "gitBranch() returned ......... " << gitBranch() );

   return 0;
}
