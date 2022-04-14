/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

// check for floating-point exceptions in invertElementwise()
// see issue #120
#include <cfenv>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

void logSectionHeader( const char* header )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separator( len + 2, '-' );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << "\n" << separator );
}

int main( int argc, char** argv )
{
#ifndef __APPLE__
   #ifndef _MSC_VER
   // should work with Intel and GCC compiler / not with windows
      feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
   #endif
#endif
   // environment stuff
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Set mesh and primitives
   // MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/unitsquare_with_circular_hole.msh" );
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/regular_octahedron_8el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t level = 2;

   // We test using P2 functions as these include P1 functions
   P2Function< real_t > tFunc1( "test vector #1", storage, level, level );
   P2Function< real_t > tFunc2( "test vector #2", storage, level, level );

   tFunc1.interpolate( real_c( 0.5 ), level, All );
   tFunc2.interpolate( real_c( 0.5 ), level, All );

   logSectionHeader( "Testing invertElementwise() after syncing halos:" );
   communication::syncP2FunctionBetweenPrimitives( tFunc1, level );
   bool workOnHalos = true;
   WALBERLA_LOG_INFO_ON_ROOT( "workOnHalos = " << std::boolalpha << workOnHalos );
   tFunc1.invertElementwise( level, All, workOnHalos );
   std::cout << std::flush;
   WALBERLA_LOG_INFO_ON_ROOT( "* Check! No floating point exceptions detected!" << std::endl );
   WALBERLA_MPI_BARRIER();

   logSectionHeader( "Testing invertElementwise() without syncing halos:" );
   workOnHalos = false;
   tFunc2.invertElementwise( level, All, workOnHalos );
   WALBERLA_LOG_INFO_ON_ROOT( "* Check! No floating point exceptions detected!" );

   return EXIT_SUCCESS;
}
