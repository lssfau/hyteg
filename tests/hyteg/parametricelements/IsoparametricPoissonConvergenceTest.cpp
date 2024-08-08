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
#include "hyteg/mesh/micro/Diffusion/P1ElementwiseDiffusion_ParametricMapP1_const_fused_quadloops_float64.hpp"
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
#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGrad.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGrad.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "mixed_operator/VectorLaplaceOperator.hpp"

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

   P2Function< real_t > k( "k", storage, level, level );
   // Can be used to steer blending of inner nodes.
   // Not sure yet if that's really useful.
   auto kfunc = []( const Point3D& x ) { return 1; };
   k.interpolate( kfunc, level );

   Lap_T      B( storage, level, level, k );
   MeshFunc_T sol( "sol", storage, level, level );
   MeshFunc_T rhs( "rhs", storage, level, level );

   sol.assign( { 1.0 }, { *storage->getMicroMesh()->p2Mesh() }, level, DirichletBoundary );
   rhs.interpolate( 0, level );

   CGSolver< Lap_T > meshSolver( storage, level, level, 500 );
   // meshSolver.setPrintInfo( true );

   const uint_t dim = storage->hasGlobalCells() ? 3 : 2;
   for ( uint_t i = 0; i < dim; i++ )
   {
      meshSolver.solve( B, sol.component( i ), rhs.component( i ), level );
   }
   storage->getMicroMesh()->p2Mesh()->assign( { 1.0 }, { sol }, level );
   WALBERLA_LOG_INFO_ON_ROOT( "Improving mesh ... done." )
}

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
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const auto microMesh =
       std::make_shared< micromesh::MicroMesh >( storage, level, level, mappingDegree, storage->hasGlobalCells() ? 3 : 2 );

   micromesh::interpolateAndCommunicate( *microMesh, blendingFunc, level );
   storage->setMicroMesh( microMesh );

   Function_T u( "u", storage, level, level );
   Function_T f( "f", storage, level, level );
   Function_T error( "error", storage, level, level );
   Function_T exact( "exact", storage, level, level );

   auto meshFunc = microMesh->p2Mesh();

   VTKOutput vtkOutput( "../../output", "IsoparametricPoissonConvergenceTest_mesh" + std::to_string( mappingDegree ), storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( error );
   vtkOutput.add( exact );
   vtkOutput.add( *meshFunc );
   vtkOutput.write( level, 0 );

   improveMesh( storage, level );

   vtkOutput.write( level, 1 );

   Mass_T M( storage, level, level, *meshFunc );
   exact.interpolate( rhsFunc, level );
   M.apply( exact, f, level, All );

   exact.interpolate( solFunc, level );
   u.interpolate( solFunc, level, DirichletBoundary );

   Operator_T A( storage, level, level, *meshFunc );

   CGSolver< Operator_T > solver( storage, level, level, 500 );
   // solver.setPrintInfo( true );

   solver.solve( A, u, f, level );

   error.assign( { 1.0, -1.0 }, { u, exact }, level );

   auto errorL2 = std::sqrt( error.dotGlobal( error, level ) / real_c( numberOfGlobalDoFs( error, level ) ) );

   writeDomainPartitioningVTK(
       storage, "../../output", "IsoparametricPoissonConvergenceTest_mesh" + std::to_string( mappingDegree ) );

   vtkOutput.write( level, 2 );

   return errorL2;
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

#if 0
   std::function< real_t( const hyteg::Point3D& ) > solFunc = []( const hyteg::Point3D& x ) {
      return std::exp( -x[0] - ( x[1] * x[1] ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhsFunc = []( const hyteg::Point3D& x ) {
      return -( 4 * x[1] * x[1] - 1 ) * std::exp( -x[0] - ( x[1] * x[1] ) );
   };
#endif

   std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
      return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] );
   };

   std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
      return 12 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) );
   };

   std::vector< Func_T > blendingFunc = {
       [&]( const Point3D& x ) { return x[0] + 0.1 * sin( 2 * pi * 1.5 * x[1] ); },
       [&]( const Point3D& x ) { return x[1] + 0.1 * sin( 2 * pi * x[2] ); },
       [&]( const Point3D& x ) { return x[2] + 0.1 * sin( 2 * pi * 1.1 * x[0] ); },
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

   for ( uint_t level = 2; level <= 4; level++ )
   {
      auto error = test< P1Function< real_t >,
                         operatorgeneration::P1ElementwiseDiffusion_ParametricMapP1_const_fused_quadloops_float64,
                         operatorgeneration::P1ElementwiseMassMicroMeshP1_const_fused_quadloops_float64 >(
          // MeshInfo::meshUnitSquare( 1 ), level, 1, solFunc, rhsFunc, blendingFunc );
          MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 ),
          level,
          1,
          solFunc,
          rhsFunc,
          blendingFunc );
      WALBERLA_LOG_DEVEL_ON_ROOT( error );
   }
#else
   WALBERLA_LOG_INFO_ON_ROOT( "P2 sol, P2 mesh" )

   for ( uint_t level = 2; level <= 4; level++ )
   {
      auto error = test< P2Function< real_t >,
                         operatorgeneration::P2ElementwiseDiffusionMicroMeshP2_const_fused_quadloops_float64,
                         operatorgeneration::P2ElementwiseMassMicroMeshP2_const_fused_quadloops_float64 >(
          // MeshInfo::meshUnitSquare( 0 ),
          MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 ),
          // MeshInfo::singleTetrahedron( { Point3D( 0, 0, 0 ), Point3D( 1, 0, 0 ), Point3D( 0, 1, 0 ), Point3D( 0, 0, 1 ) } ),
          level,
          2,
          solFunc,
          rhsFunc,
          blendingFunc );
      WALBERLA_LOG_DEVEL_ON_ROOT( error );
   }
#endif

   return 0;
}
