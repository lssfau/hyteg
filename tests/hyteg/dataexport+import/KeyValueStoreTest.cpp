/*
 * Copyright (c) 2023 Marcus Mohr, Daniel Bauer.
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

#include "hyteg/dataexport/LaTeX/KeyValueStore.hpp"

#include <sstream>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"

using walberla::real_c;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::latex::KeyValueStore store;

   store.store( "/example/minLevel", 0 );
   store.store( "/example/maxLevel", 8 );
   store.store( "/example/string", "Hello, World!" );
   store.store( "/example/real_t", real_c( 42.0 ) );

   std::stringstream ss;
   ss << store;

   std::string expected =
R"(/example/maxLevel = 8
/example/minLevel = 0
/example/real_t   = 42
/example/string   = Hello, World!
)";

   WALBERLA_CHECK_EQUAL( ss.str(), expected );

   store.writePgfKeys( ".", "KeyValueStoreTest" );
}
