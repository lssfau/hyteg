/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/types/Matrix.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   Eigen::Matrix< real_t, 3, 3 > m3;

   // clang-format off
   m3 <<  1.0,  2.0,  3.0,
         -4.0, -5.0, -6.0,
          0.1, -0.2, -1.3;
   // clang-format on

   Matrix3r matrix3{ m3 };

   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 0, 0 ), 1.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 0, 1 ), 2.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 0, 2 ), 3.0 )

   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 1, 0 ), -4.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 1, 1 ), -5.0 )
   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 1, 2 ), -6.0 )

   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 2, 0 ), 0.1 )
   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 2, 1 ), -0.2 )
   WALBERLA_CHECK_FLOAT_EQUAL( matrix3( 2, 2 ), -1.3 )

   return EXIT_SUCCESS;
}
