/*
 * Copyright (c) 2024 Andreas Burkhart.
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
#ifdef HYTEG_BUILD_WITH_PETSC
#ifdef PETSC_HAVE_HDF5
#define HYTEG_BUILD_WITH_PETSC_AND_HDF5
#endif
#endif

#include <core/timing/Timer.h>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
// clang-format off
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScHDF5FunctionSave.hpp"
// clang-format on
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   #ifdef HYTEG_BUILD_WITH_PETSC_AND_HDF5

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager;

   MeshInfo meshInfo =
       MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::meshFlavour::CRISSCROSS, 2, 2 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   // constants
   uint_t minLevel_ = 0;
   uint_t maxLevel_ = 2;

   // fem functions
   P1Function< real_t > f_( "f", storage, minLevel_, maxLevel_ );
   P1Function< real_t > f2_( "f2", storage, minLevel_, maxLevel_ );

   P1Function< real_t > f_comp_( "f_comp", storage, minLevel_, maxLevel_ );
   P1Function< real_t > f2_comp_( "f2_comp", storage, minLevel_, maxLevel_ );

   // interpolate
   std::function< real_t( const Point3D& ) > randFunc = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   f_.interpolate( randFunc, maxLevel_, All );
   f2_.interpolate( randFunc, maxLevel_, All );

   f_comp_.assign( { real_c( 1 ) }, { f_ }, maxLevel_, All );
   f2_comp_.assign( { real_c( 1 ) }, { f2_ }, maxLevel_, All );

   // create some random parameters to be saved
   real_t p1 = walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   real_t p2 = walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );

   real_t p1_comp = p1;
   real_t p2_comp = p2;

   // save functions
   saveMultipleFunctionsPETSc( maxLevel_, "HDF5SaveTest.dat", f_, f2_ );
   // save parameters

   std::vector< PetscScalar > parameters_save = { p1, p2 };
   saveParametersPETSc( parameters_save, "HDF5SaveTest.dat", FILE_MODE_APPEND );

   // clear values
   f_.interpolate( 0, maxLevel_, All );
   f2_.interpolate( 0, maxLevel_, All );
   p1 = real_c( 0.0 );
   p2 = real_c( 0.0 );

   // load functions
   loadMultipleFunctionsPETSc( maxLevel_, "HDF5SaveTest.dat", f_, f2_ );
   // load parameters
   std::vector< PetscScalar > parameters_load;
   loadParametersPETSc( parameters_load, "HDF5SaveTest.dat" );

   // compare
   f_comp_.assign( { real_c( 1 ), real_c( -1 ) }, { f_comp_, f_ }, maxLevel_, All );
   f2_comp_.assign( { real_c( 1 ), real_c( -1 ) }, { f2_comp_, f2_ }, maxLevel_, All );
   p1_comp = p1_comp - parameters_load[0];
   p2_comp = p2_comp - parameters_load[1];

   real_t f_sum  = f_comp_.dotGlobal( f_comp_, maxLevel_, All );
   real_t f2_sum = f2_comp_.dotGlobal( f2_comp_, maxLevel_, All );

   WALBERLA_CHECK_FLOAT_EQUAL( f_sum, real_c( 0.0 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( f2_sum, real_c( 0.0 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( p1_comp, real_c( 0.0 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( p2_comp, real_c( 0.0 ) );

   #endif

   return EXIT_SUCCESS;
}
