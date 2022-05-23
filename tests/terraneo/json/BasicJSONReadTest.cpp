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
#include "core/debug/CheckFunctions.h"
#include "core/extern/json.hpp"
#include "core/logging/Logging.h"

#include "hyteg/BuildInfo.hpp"
#include "hyteg/Git.hpp"

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::string separator{"--------------------------------------------------"};
   WALBERLA_LOG_INFO_ON_ROOT( separator );
   printBuildInfo();
   printGitInfo();
   WALBERLA_LOG_INFO_ON_ROOT( separator );

   std::string   jsonFileName{"../../data/terraneo/BasicJSONReadTest.json"};
   std::ifstream jsonFile;

   WALBERLA_LOG_INFO_ON_ROOT( "Opening '" << jsonFileName << "'" );
   jsonFile.open( jsonFileName.c_str() );

   WALBERLA_CHECK( !!jsonFile, "Error opening file: " << jsonFileName );

   WALBERLA_LOG_INFO_ON_ROOT( "Reading '" << jsonFileName << "'" );
   nlohmann::json jsonObject;
   jsonFile >> jsonObject;

   WALBERLA_LOG_INFO_ON_ROOT( "\n" << jsonObject.dump(4) );

   WALBERLA_CHECK_EQUAL( jsonObject[ "preconditionerType" ], "Multigrid" );
   WALBERLA_CHECK_EQUAL( jsonObject[ "preconditionerParams" ][ "preSmooth" ], 3 );

   return 0;
}
