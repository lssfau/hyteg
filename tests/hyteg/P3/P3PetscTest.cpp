/*
 * Copyright (c) 2017-2026 Dominik Thoennes, Marcus Mohr.
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
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/eigen/EigenSparseDirectSolver.hpp"
#include "hyteg/elementwiseoperators/P3ElementwiseOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p3functionspace/P3Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( const uint_t& level, const std::string& meshFileName, const real_t& errEps, bool doVTKOutput = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "##### Mesh file: " << meshFileName << " / level: " << level << " #####" )

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P3Function< real_t > x( "x", storage, level, level );
   P3Function< real_t > x_exact( "x_exact", storage, level, level );
   P3Function< real_t > b( "b", storage, level, level );
   P3Function< real_t > err( "err", storage, level, level );
   P3Function< real_t > residuum( "residuum", storage, level, level );

   P3ElementwiseDiffusionOperator A( storage, level, level );

   std::function< real_t( const Point3D& ) > exact = []( const Point3D& xx ) { return std::sin( xx[0] ) * std::sinh( xx[1] ); };
   walberla::math::seedRandomGenerator( 0 );
   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) {
      return real_c( walberla::math::realRandom( 0.0, 1.0 ) );
   };

   x.interpolate( exact, level, DirichletBoundary );
   x.interpolate( rand, level, Inner );
   b.interpolate( exact, level, DirichletBoundary );
   x_exact.interpolate( exact, level );

   uint_t localDoFs  = numberOfLocalDoFs< P3FunctionTag >( *storage, level );
   uint_t globalDoFs = numberOfGlobalDoFs< P3FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO( "localDoFs: " << localDoFs << " globalDoFs: " << globalDoFs );

   // PETScLUSolver< P3ElementwiseDiffusionOperator > solver( storage, level );
   EigenSparseDirectSolver< P3ElementwiseDiffusionOperator > solver( storage, level );

   walberla::WcTimer timer;
   solver.solve( A, x, b, level );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, residuum, level, Inner );

   err.assign( { 1.0, -1.0 }, { x, x_exact }, level );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs );
   real_t residuum_l2  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error 1 = " << discr_l2_err );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2 );

   if ( doVTKOutput )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "****************************" );
      VTKOutput vtkOutput( ".", "P3PetscSolve", storage );
      vtkOutput.add( x );
      vtkOutput.add( x_exact );
      vtkOutput.add( err );
      vtkOutput.add( residuum );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_LESS( residuum_l2, real_c( 4e-15 ) );
   WALBERLA_CHECK_LESS( discr_l2_err, errEps );
}

void petscMassMatrixApplyTest( const uint_t& level )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Testing with Backward Facing Step" );

   MeshInfo meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/bfs_12el.msh" ) );
   real_t   area     = real_c( 1.75 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P3ElementwiseMassOperator massOp( storage, level, level );

   P3Function< real_t > aux( "aux", storage, level, level );
   P3Function< real_t > vecOfOnes( "vecOfOnes", storage, level, level );
   P3Function< idx_t >  enumerator( "enumerator", storage, level, level );

   vecOfOnes.interpolate( real_c( 1 ), level, All );
   enumerator.enumerate( level );

   WALBERLA_LOG_INFO_ON_ROOT( "* Assembling global mass matrix" );
   PETScManager                                   petscManager;
   PETScSparseMatrix< P3ElementwiseMassOperator > matrix( "GlobalMassMatrix" );

   matrix.createMatrixFromOperator( massOp, level, enumerator, All );

   WALBERLA_LOG_INFO_ON_ROOT( "* Assembling vector of ones" );
   PETScVector< real_t, P3Function > srcVector( vecOfOnes, enumerator, level, All, "src" );

   WALBERLA_LOG_INFO_ON_ROOT( "* Performing MatVec" );
   PETScVector< real_t, P3Function > dstVector( aux, enumerator, level, All, "dst" );
   auto                              petscErrorCode = MatMult( matrix.get(), srcVector.get(), dstVector.get() );
   WALBERLA_CHECK_EQUAL( petscErrorCode, PETSC_SUCCESS );

   WALBERLA_LOG_INFO_ON_ROOT( "* Computing dot product" );
   real_t measure = real_c( 0 );
   petscErrorCode = VecDot( srcVector.get(), dstVector.get(), &measure );
   WALBERLA_CHECK_EQUAL( petscErrorCode, PETSC_SUCCESS );

   WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure );
   WALBERLA_CHECK_FLOAT_EQUAL( measure, area );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   petscMassMatrixApplyTest( 2 );

   // petscSolveTest( 0, prependHyTeGMeshDir( "2D/quad_4el.msh" ), 3.0e-04, true );
   // petscSolveTest( 1, prependHyTeGMeshDir( "2D/quad_4el.msh" ), 2.0e-05 );
   // petscSolveTest( 3, prependHyTeGMeshDir( "2D/quad_4el.msh" ), 3.0e-07, true );

   return EXIT_SUCCESS;
}
