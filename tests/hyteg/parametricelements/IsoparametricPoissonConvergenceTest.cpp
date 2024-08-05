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
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/mesh/micro/DiffusionMicroMeshP1/P1ElementwiseDiffusionMicroMeshP1_const_fused_quadloops_float64.hpp"
#include "hyteg/mesh/micro/DiffusionMicroMeshP2/P2ElementwiseDiffusionMicroMeshP2_const_fused_quadloops_float64.hpp"
#include "hyteg/mesh/micro/MassMicroMeshP1/P1ElementwiseMassMicroMeshP1_const_fused_quadloops_float64.hpp"
#include "hyteg/mesh/micro/MassMicroMeshP2/P2ElementwiseMassMicroMeshP2_const_fused_quadloops_float64.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

using Func_T = std::function< real_t( const hyteg::Point3D& ) >;

template < typename Function_T, typename Operator_T, typename Mass_T >
real_t test( MeshInfo              meshInfo,
             uint_t                level,
             uint_t                mappingDegree,
             Func_T                solFunc,
             Func_T                rhsFunc,
             std::vector< Func_T > blendingFunc )
{
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage   = std::make_shared< PrimitiveStorage >( setupStorage );

   const auto microMesh = std::make_shared< micromesh::MicroMesh >( storage, level, level, mappingDegree );
   storage->setMicroMesh( microMesh );

   micromesh::interpolateAndCommunicate( storage, blendingFunc, level );

   Function_T u( "u", storage, level, level );
   Function_T f( "f", storage, level, level );
   Function_T error( "error", storage, level, level );
   Function_T exact( "exact", storage, level, level );

   auto mesh = *storage->getMicroMesh()->p2Mesh();

   Mass_T M( storage, level, level, mesh );
   exact.interpolate( rhsFunc, level );
   M.apply( exact, f, level, All );

   exact.interpolate( solFunc, level );
   u.interpolate( solFunc, level, DirichletBoundary );

#if 0
   Operator_T A( storage, level, level );
#else
   Operator_T A( storage, level, level, mesh );
#endif

   CGSolver< Operator_T > solver( storage, level, level );

   solver.solve( A, u, f, level );

   error.assign( { 1.0, -1.0 }, { u, exact }, level );

   auto errorL2 = std::sqrt( error.dotGlobal( error, level ) / real_c( numberOfGlobalDoFs( error, level ) ) );

   writeDomainPartitioningVTK(
       storage, "../../output", "IsoparametricPoissonConvergenceTest_mesh" + std::to_string( mappingDegree ) );

   VTKOutput vtkOutput( "../../output", "IsoparametricPoissonConvergenceTest_mesh" + std::to_string( mappingDegree ), storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( error );
   vtkOutput.add( exact );
   vtkOutput.write( level, 0 );

   return errorL2;
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::function< real_t( const hyteg::Point3D& ) > solFunc = []( const hyteg::Point3D& x ) {
      return std::exp( -x[0] - ( x[1] * x[1] ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhsFunc = []( const hyteg::Point3D& x ) {
      return -( 4 * x[1] * x[1] - 1 ) * std::exp( -x[0] - ( x[1] * x[1] ) );
   };

   std::vector< Func_T > blendingFunc = {
       [&]( const Point3D& x ) { return x[0] + 0.1 * sin( 2 * pi * 1.5 * x[1] ); },
       [&]( const Point3D& x ) { return x[1] + 0.1 * sin( 2 * pi * x[0] ); },
       [&]( const Point3D& x ) { return x[2]; },
   };
#if 0
   blendingFunc = {
       [&]( const Point3D& x ) { return x[0]; },
       [&]( const Point3D& x ) { return x[1]; },
       [&]( const Point3D& x ) { return x[2]; },
   };
#endif

#if 0
   WALBERLA_LOG_INFO_ON_ROOT( "P1 sol, P1 mesh" )

   for ( uint_t level = 2; level <= 5; level++ )
   {
      auto error =
          test< P1Function< real_t >, operatorgeneration::P1ElementwiseDiffusionMicroMeshP1_const_fused_quadloops_float64 >(
              MeshInfo::meshUnitSquare( 1 ), level, 1, solFunc, rhsFunc, blendingFunc );
      WALBERLA_LOG_DEVEL_ON_ROOT( error );
   }
#endif

#if 1
   WALBERLA_LOG_INFO_ON_ROOT( "P2 sol, P2 mesh" )

   for ( uint_t level = 2; level <= 5; level++ )
   {
      auto error = test< P2Function< real_t >,
                         operatorgeneration::P2ElementwiseDiffusionMicroMeshP2_const_fused_quadloops_float64,
                         operatorgeneration::P2ElementwiseMassMicroMeshP2_const_fused_quadloops_float64 >(
          MeshInfo::meshUnitSquare( 0 ), level, 2, solFunc, rhsFunc, blendingFunc );
      WALBERLA_LOG_DEVEL_ON_ROOT( error );
   }
#endif

   return 0;
}
