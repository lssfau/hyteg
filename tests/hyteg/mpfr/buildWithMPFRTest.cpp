/*
* Copyright (c) 2023-2024 Michael Zikeli.
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

#include <mpfr.h>
#include <sstream>
#include <stdio.h>

#include "core/Environment.h"
#include "core/debug/Debug.h"

namespace hyteg::mpfrTest {

void callVersion()
{
   printf( "MPFR library: %-12s\nMPFR header:  %s (based on %d.%d.%d)\n",
           mpfr_get_version(),
           MPFR_VERSION_STRING,
           MPFR_VERSION_MAJOR,
           MPFR_VERSION_MINOR,
           MPFR_VERSION_PATCHLEVEL );
   WALBERLA_CHECK_NOT_NULLPTR( mpfr_get_version() );

#ifdef MPFR_VERSION_STRING
   WALBERLA_CHECK( true );
#else
   WALBERLA_CHECK( false );
#endif

#ifdef MPFR_VERSION_MAJOR
   WALBERLA_CHECK( true );
#else
   WALBERLA_CHECK( false );
#endif

#ifdef MPFR_VERSION_MINOR
   WALBERLA_CHECK( true );
#else
   WALBERLA_CHECK( false );
#endif

#ifdef MPFR_VERSION_PATCHLEVEL
   WALBERLA_CHECK( true );
#else
   WALBERLA_CHECK( false );
#endif

   return;
} //callVersion()
} // namespace hyteg::mpfrTest

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   hyteg::mpfrTest::callVersion();

   return EXIT_SUCCESS;
}
