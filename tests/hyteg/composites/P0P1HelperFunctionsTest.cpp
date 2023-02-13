/*
* Copyright (c) 2023 Marcus Mohr.
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

#include "hyteg/composites/P0P1HelperFunctions.hpp"

#include <cfenv>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/all.h"
#include "core/mpi/all.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

void runProjectionTestIn2D()
{
   // Set mesh and primitives
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // Setup functions
   uint_t               minLevel = 0;
   uint_t               maxLevel = 3;
   P1Function< real_t > p1Source( "source P1", storage, minLevel, maxLevel );
   P0Function< real_t > p0Destination( "destination P0", storage, minLevel, maxLevel );
   P0Function< real_t > p0Control( "control P0", storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > affineExpr = []( const Point3D& x ) {
      return real_c( 2.0 ) * x[0] - x[1] + real_c( 3.5 );
   };

   bool outputVTK = false;

#ifdef WALBERLA_DOUBLE_ACCURACY
   real_t tolerance = 1e-14;
#else
   real_t tolerance = 5e-6;
#endif

   // Test projection
   for ( uint_t currentLevel = minLevel; currentLevel <= maxLevel; ++currentLevel )
   {
      p1Source.interpolate( affineExpr, currentLevel, All );
      p0Control.interpolate( affineExpr, currentLevel, All );
      p0Destination.interpolate( real_c( -1 ), currentLevel, All );

      projectP1ToP0( p1Source, p0Destination, currentLevel, Replace );

      p0Control.assign( { real_c( 1.0 ), real_c( -1.0 ) }, { p0Control, p0Destination }, currentLevel, All );
      real_t error = std::sqrt( p0Control.dotGlobal( p0Control, currentLevel, All ) );
      WALBERLA_LOG_INFO_ON_ROOT( "Projection error in discrete L2-norm is " << error );
      WALBERLA_CHECK_LESS( error, tolerance );

      if ( outputVTK )
      {
         VTKOutput vtkOutput( "../../output", "P0P1HelperFunctionsTest2D", storage );
         vtkOutput.add( p1Source );
         vtkOutput.add( p0Control );
         vtkOutput.add( p0Destination );
         vtkOutput.write( currentLevel );
      }
   }
}

int main( int argc, char* argv[] )
{
#if defined (WALBERLA_DOUBLE_ACCURACY) || !defined(__INTEL_LLVM_COMPILER)
#ifndef __APPLE__
// should work with Intel, GCC, Clang and even MSVC compiler /nope not MSVC
#ifndef _MSC_VER
   feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#endif
#endif
#endif

   // environment stuff
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::debug::enterTestMode();

   // run test
   runProjectionTestIn2D();

   return 0;
}
