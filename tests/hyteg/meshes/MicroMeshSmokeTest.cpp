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
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGrad.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

using Func_T = std::function< real_t( const hyteg::Point3D& ) >;

void improveMesh( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Improving mesh ..." )
   using Lap_T      = operatorgeneration::P2ElementwiseDivKGrad;
   using MeshFunc_T = P2VectorFunction< real_t >;

   // Can be used to steer blending of inner nodes.
   // Not sure yet if that's really useful.
   P2Function< real_t > k( "k", storage, level, level );
   auto                 kfunc = []( const Point3D& x ) { return 1; };
   k.interpolate( kfunc, level );

   Lap_T      B( storage, level, level, k );
   MeshFunc_T sol( "sol", storage, level, level );
   MeshFunc_T rhs( "rhs", storage, level, level );

   sol.assign( { 1.0 }, { *storage->getMicroMesh()->p2Mesh() }, level, DirichletBoundary );
   rhs.interpolate( 0, level );

   CGSolver< Lap_T > meshSolver( storage, level, level, 500 );

   const uint_t dim = storage->hasGlobalCells() ? 3 : 2;
   for ( uint_t i = 0; i < dim; i++ )
   {
      meshSolver.solve( B, sol.component( i ), rhs.component( i ), level );
   }
   storage->getMicroMesh()->p2Mesh()->assign( { 1.0 }, { sol }, level );
   WALBERLA_LOG_INFO_ON_ROOT( "Improving mesh ... done." )
}

void test( MeshInfo meshInfo, uint_t level, uint_t mappingDegree, std::vector< Func_T > blendingFunc )
{
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t dim = storage->hasGlobalCells() ? 3 : 2;

   writeDomainPartitioningVTK( storage,
                               "../../output",
                               "MicroMeshSmokeTest" + std::to_string( dim ) + "DPrimitiveStorage" +
                                   std::to_string( mappingDegree ) );

   const auto microMesh =
       std::make_shared< micromesh::MicroMesh >( storage, level, level, mappingDegree, storage->hasGlobalCells() ? 3 : 2 );

   micromesh::interpolateAndCommunicate( *microMesh, blendingFunc, level );
   storage->setMicroMesh( microMesh );

   P2Function< real_t > u( "u", storage, level, level );

   auto someGradient = []( const Point3D& x ) { return x.sum(); };
   u.interpolate( someGradient, level );

   VTKOutput vtkOutput(
       "../../output", "MicroMeshSmokeTest" + std::to_string( dim ) + "D" + "P" + std::to_string( mappingDegree ), storage );
   vtkOutput.add( u );
   vtkOutput.write( level, 0 );
   improveMesh( storage, level );
   vtkOutput.write( level, 1 );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::vector< Func_T > blendingFunc = {
       [&]( const Point3D& x ) { return x[0] + 0.1 * sin( 2 * pi * 1.5 * x[1] ); },
       [&]( const Point3D& x ) { return x[1] + 0.1 * sin( 2 * pi * x[0] ); },
       [&]( const Point3D& x ) { return x[2] + 0.1 * sin( 2 * pi * 1.1 * x[0] ); },
   };

   for ( uint_t level = 2; level <= 3; level++ )
   {
      test( MeshInfo::meshUnitSquare( 1 ), level, 2, blendingFunc );
      test( MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 ), level, 2, blendingFunc );
   }

   return 0;
}
