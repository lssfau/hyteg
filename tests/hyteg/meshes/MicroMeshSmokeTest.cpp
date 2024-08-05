/*
* Copyright (c) 2017-2024 Nils Kohl.
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

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

void test()
{
   const uint_t level = 4;

   const auto            meshInfo = MeshInfo::meshUnitSquare( 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto            storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const auto microMesh = std::make_shared< micromesh::MicroMesh >( storage, level, level, 1 );

   storage->setMicroMesh( microMesh );

   P1Function< real_t > u( "u", storage, level, level );
   u.interpolate( 1, level );

   writeDomainPartitioningVTK( storage, "../../output", "MicroMeshSmokeTestDomain" );

   VTKOutput vtkOutput( "../../output", "MicroMeshSmokeTest", storage );
   vtkOutput.add( u );
   vtkOutput.write( level, 0 );

   VTKOutput vtkOutputNewMesh( "../../output", "MicroMeshSmokeTest_NewMesh", storage );
   vtkOutputNewMesh.add( u );

   for ( uint_t i = 0; i < 100; i++ )
   {
      auto meshMoverY = [&]( const Point3D& x ) { return x[1] + 0.1 * sin( 2 * pi * ( x[0] + ( real_c( i ) / 100.0 ) ) ); };

      storage->getMicroMesh()->p1Mesh()->component( 1 ).interpolate( meshMoverY, level );
      vtkOutputNewMesh.write( level, i );
   }
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   test();

   return 0;
}
