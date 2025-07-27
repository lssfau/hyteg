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

#include <algorithm>

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg_operators/operators/diffusion/P1ElementwiseDiffusionParametricP1Map.hpp"
#include "hyteg_operators/operators/diffusion/P2ElementwiseDiffusionParametricP2Map.hpp"
#include "hyteg_operators/operators/mass/P1ElementwiseMassParametricP1Map.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassParametricP2Map.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "mixed_operator/VectorLaplaceOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

using Func_T = std::function< real_t( const hyteg::Point3D& ) >;

template < typename ConstantLaplace_T, typename Mesh_T >
void improveMesh( const Mesh_T& mesh, uint_t level )
{
   auto storage = mesh.getStorage();

   ConstantLaplace_T B( storage, level, level );
   Mesh_T            rhs( "rhs", storage, level, level );

   CGSolver< ConstantLaplace_T > meshSolver( storage, level, level, 500 );

   const uint_t dim = storage->hasGlobalCells() ? 3 : 2;
   for ( uint_t i = 0; i < dim; i++ )
   {
      meshSolver.solve( B, mesh.component( i ), rhs.component( i ), level );
   }
}

template < typename Function_T,
           typename Operator_T,
           typename ConstantLaplace_T,
           typename Mass_T,
           typename Mesh_T,
           typename Restriction_T,
           typename Prolongation_T >
void test( const MeshInfo&              meshInfo,
           uint_t                       minLevel,
           uint_t                       maxLevel,
           uint_t                       mappingDegree,
           Func_T                       solFunc,
           Func_T                       rhsFunc,
           const std::vector< Func_T >& blendingFunc )
{
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Let's solve Laplace's equation for the mesh node displacement on the affine/computational mesh to improve the mesh quality.
   // That means we solve the vector Laplace equation
   //
   //   - Δ u = 0
   //
   // where u is a vector that represents the mesh displacement. The Dirichlet boundary conditions are imposed via the position at
   // the boundary. The domain is the computational (non-deformed) refined mesh. So we can use the constant Laplace operator.
   //
   // This is somewhat heuristic, and we need to make sure that this converges. But works quite nice here and just show how you
   // can work on the micro-mesh. There are more sophisticated ways to get a better mesh.

   auto mesh = std::make_shared< Mesh_T >( "mesh", storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Improving mesh ..." )
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      mesh->interpolate( blendingFunc, level );
      improveMesh< ConstantLaplace_T, Mesh_T >( *mesh, level );
      communication::syncVectorFunctionBetweenPrimitives( *mesh, level );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "Improving mesh ... done." )

   // From this point on we actually add the mesh to the storage.
   // After that we simply solve Poisson with a full multigrid solver on a domain with curved boundaries.

   auto microMesh = std::make_shared< micromesh::MicroMesh >( mesh );
   storage->setMicroMesh( microMesh );

   Function_T u( "u", storage, minLevel, maxLevel );
   Function_T f( "f", storage, minLevel, maxLevel );
   Function_T error( "error", storage, minLevel, maxLevel );
   Function_T exact( "exact", storage, minLevel, maxLevel );

   Operator_T A( storage, minLevel, maxLevel, *mesh );
   A.computeInverseDiagonalOperatorValues();

   auto prolongation = std::make_shared< Prolongation_T >();
   auto restriction  = std::make_shared< Restriction_T >();

   auto smoother = std::make_shared< ChebyshevSmoother< Operator_T > >( storage, minLevel, maxLevel );
   error.interpolate( []( const Point3D& ) { return walberla::math::realRandom(); }, maxLevel );
   smoother->setupCoefficients( 3, chebyshev::estimateRadius( A, maxLevel, 20, storage, error, exact ) );

   auto coarseGridSolver = std::make_shared< CGSolver< Operator_T > >( storage, minLevel, minLevel, 1000, real_c(0), 1e-10 );
   coarseGridSolver->setPrintInfo( false );

   auto gmgSolver = std::make_shared< GeometricMultigridSolver< Operator_T > >(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, 2, 2 );

   FullMultigridSolver< Operator_T > fmgSolver( storage, gmgSolver, prolongation, minLevel, maxLevel, 20 );

   Mass_T M( storage, minLevel, maxLevel, *mesh );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      exact.interpolate( rhsFunc, level );
      M.apply( exact, f, level, All );

      exact.interpolate( solFunc, level );

      u.interpolate( solFunc, level, DirichletBoundary );
   }

   real_t errorL2Prev = 0;
   real_t errorL2     = 0;
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      fmgSolver.solve( A, u, f, level );
      error.assign( { 1.0, -1.0 }, { u, exact }, level );
      errorL2Prev = errorL2;
      errorL2     = std::sqrt( error.dotGlobal( error, level ) / real_c( numberOfGlobalDoFs( error, level ) ) );
      auto rate   = errorL2Prev / errorL2;
      WALBERLA_LOG_INFO_ON_ROOT( "level " << level << " | error " << errorL2 << " | rate " << rate );
      if ( level > minLevel )
      {
         if ( ( std::is_same_v< real_t, float > && errorL2 > 3e-6 ) || std::is_same_v< real_t, double > )
         {
            WALBERLA_CHECK_GREATER( rate, 3.4 * real_c( polynomialDegreeOfBasisFunctions< Function_T >() ) );
         }
      }
   }

   writeDomainPartitioningVTK(
       storage, "../../output", "IsoparametricPoissonConvergenceTest_mesh" + std::to_string( mappingDegree ) );

   VTKOutput vtkOutput( "../../output", "IsoparametricPoissonConvergenceTest_mesh" + std::to_string( mappingDegree ), storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( error );
   vtkOutput.add( exact );
   vtkOutput.add( *mesh );
   vtkOutput.add( *A.getInverseDiagonalValues() );
   vtkOutput.write( maxLevel );
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
       [&]( const Point3D& ) { return 0; },
   };

   WALBERLA_LOG_INFO_ON_ROOT( "P1 sol, P1 mesh" )

   test< P1Function< real_t >,
         operatorgeneration::P1ElementwiseDiffusionParametricP1Map,
         P1ConstantLaplaceOperator,
         operatorgeneration::P1ElementwiseMassParametricP1Map,
         P1VectorFunction< real_t >,
         P1toP1LinearRestriction< real_t >,
         P1toP1LinearProlongation< real_t > >( MeshInfo::meshUnitSquare( 1 ), 2, 5, 1, solFunc, rhsFunc, blendingFunc );

   WALBERLA_LOG_INFO_ON_ROOT( "P2 sol, P2 mesh" )

   test< P2Function< real_t >,
         operatorgeneration::P2ElementwiseDiffusionParametricP2Map,
         P2ConstantLaplaceOperator,
         operatorgeneration::P2ElementwiseMassParametricP2Map,
         P2VectorFunction< real_t >,
         P2toP2QuadraticRestriction,
         P2toP2QuadraticProlongation >( MeshInfo::meshUnitSquare( 0 ), 2, 5, 2, solFunc, rhsFunc, blendingFunc );

   return 0;
}
