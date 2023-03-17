/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr.
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
#include <hyteg/communication/Syncing.hpp>

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include <cstdio>

using walberla::real_c;
using walberla::real_t;

using namespace hyteg;

void test_dotLocalGlobal_2D()
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_4el.msh" );
   // MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   const size_t minLevel = 0;
   const size_t maxLevel = 5;

   walberla::math::seedRandomGenerator( 12345678 );

   const uint_t numRandomEvaluations = 100;

   P0Function< real_t > x( "x", storage, minLevel, maxLevel );

   P0Function< real_t > y( "x", storage, minLevel, maxLevel );

   // IMPORTANT
   // communication::syncFunctionBetweenPrimitives( x, maxLevel );
   // communication::syncFunctionBetweenPrimitives( y, maxLevel );

   const uint_t numGlobalDoFs = x.getNumberOfGlobalDoFs(maxLevel);

   const uint_t numLocalDoFs = x.getNumberOfLocalDoFs(maxLevel);

   // const int rank = walberla::mpi::MPIManager::instance()->rank();

   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      real_t test_func_1_val = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
      real_t test_func_2_val = real_c( walberla::math::realRandom( 0.0, 1.0 ) );

      auto testFunc_1            = [&test_func_1_val]( const Point3D& ) { return real_c( test_func_1_val ); };
      auto testFunc_2            = [&test_func_2_val]( const Point3D& ) { return real_c( test_func_2_val ); };

      x.interpolate( testFunc_1, maxLevel );
      y.interpolate( testFunc_2, maxLevel );

      // Learnt that this is essential! 
      x.communicate(maxLevel);
      y.communicate(maxLevel);

      real_t dotvalue_local = x.dotLocal(y, maxLevel);
      real_t dotvalue_global = x.dotGlobal(y, maxLevel);

      /* For P0 function, the dot value must be direct multiplication of function values and the corresponding number of dofs */
      WALBERLA_CHECK_FLOAT_EQUAL( dotvalue_local, ((real_t)numLocalDoFs) * test_func_1_val * test_func_2_val );
      WALBERLA_CHECK_FLOAT_EQUAL( dotvalue_global, ((real_t)numGlobalDoFs) * test_func_1_val * test_func_2_val );

   }
   
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test_dotLocalGlobal_2D();

   return EXIT_SUCCESS;
}
