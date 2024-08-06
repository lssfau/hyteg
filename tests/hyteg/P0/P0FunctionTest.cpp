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
#include "hyteg/p0functionspace/P0Function.hpp"

#include <cstdio>
#include <hyteg/communication/Syncing.hpp>

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::real_t;

using namespace hyteg;

void test_dotLocalGlobal( std::string meshFileName )
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   const size_t minLevel = 0;
   const size_t maxLevel = 5;

   // walberla::math::seedRandomGenerator( 12345678 );

   const uint_t numRandomEvaluations = 1;

   P0Function< real_t > x( "x", storage, minLevel, maxLevel );

   P0Function< real_t > y( "x", storage, minLevel, maxLevel );

   const uint_t numGlobalDoFs = x.getNumberOfGlobalDoFs( maxLevel );

   const uint_t numLocalDoFs = x.getNumberOfLocalDoFs( maxLevel );

   real_t testFuncValX = 12.0;
   real_t testFuncValY = 27.0;

   // const int rank = walberla::mpi::MPIManager::instance()->rank();

   /* This tests using the constant interpolate function */
   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      x.interpolate( testFuncValX, maxLevel );
      y.interpolate( testFuncValY, maxLevel );

      real_t dotValueLocal  = x.dotLocal( y, maxLevel );
      real_t dotValueGlobal = x.dotGlobal( y, maxLevel );

      /* For P0 function, the dot value must be direct multiplication of function values and the corresponding number of dofs */
      WALBERLA_CHECK_FLOAT_EQUAL( dotValueLocal, ( (real_t) numLocalDoFs ) * testFuncValX * testFuncValY );
      WALBERLA_CHECK_FLOAT_EQUAL( dotValueGlobal, ( (real_t) numGlobalDoFs ) * testFuncValX * testFuncValY );
   }

   /* This tests using the expression interpolate function (in turn tests the 2D and 3D expr interpolate functions also) */
   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      auto testFuncX = [&testFuncValX]( const Point3D& ) { return real_c( testFuncValX ); };
      auto testFuncY = [&testFuncValY]( const Point3D& ) { return real_c( testFuncValY ); };

      x.interpolate( testFuncX, maxLevel );
      y.interpolate( testFuncY, maxLevel );

      real_t dotValueLocal  = x.dotLocal( y, maxLevel );
      real_t dotValueGlobal = x.dotGlobal( y, maxLevel );

      /* For P0 function, the dot value must be direct multiplication of function values and the corresponding number of dofs */
      WALBERLA_CHECK_FLOAT_EQUAL( dotValueLocal, ( (real_t) numLocalDoFs ) * testFuncValX * testFuncValY );
      WALBERLA_CHECK_FLOAT_EQUAL( dotValueGlobal, ( (real_t) numGlobalDoFs ) * testFuncValX * testFuncValY );
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test_dotLocalGlobal( prependHyTeGMeshDir( "tri_4el.msh" ) );        // 2D
   test_dotLocalGlobal( prependHyTeGMeshDir( "3D/pyramid_4el.msh" ) ); // 3D

   return EXIT_SUCCESS;
}
