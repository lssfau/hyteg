/*
* Copyright (c) 2024 Marcus Mohr.
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

// check pipeline compilers for support of C++20 concepts

#include <concepts>
#include <core/Environment.h>

#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

template< typename T >
requires std::integral<T>
struct myClass {
  T attribute;
};

int main( int argc, char** argv )
{

   // environment stuff
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   myClass< uint > obj1;

   obj1.attribute = 2UL;

   WALBERLA_LOG_INFO_ON_ROOT( "If you see this message, compiling with concepts worked ;-)" );

   return EXIT_SUCCESS;
}
