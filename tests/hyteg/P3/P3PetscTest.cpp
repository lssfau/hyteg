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
#include "core/math/Constants.h"
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
   WALBERLA_LOG_INFO_ON_ROOT( "\n*** PETSc Solve Test (cubic polynomial) ***" );

   WALBERLA_LOG_INFO_ON_ROOT( "-> loading mesh from: " << meshFileName );
   WALBERLA_LOG_INFO_ON_ROOT( "-> refinement level: " << level );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // setup FE functions
   WALBERLA_LOG_INFO_ON_ROOT( "-> setting up FE functions" );
   P3Function< real_t > uDiscr( "uDiscr", storage, level, level );
   P3Function< real_t > uExact( "uExact", storage, level, level );
   P3Function< real_t > b( "b", storage, level, level );
   P3Function< real_t > bInitial( "bInitial", storage, level, level );
   P3Function< real_t > bPostInner( "bPostInner", storage, level, level );
   P3Function< real_t > rhs( "rhs", storage, level, level );
   P3Function< real_t > err( "err", storage, level, level );

   // define and interpolate exact solution
   WALBERLA_LOG_INFO_ON_ROOT( "-> performing interpolation for exact solution and boundary values" );
   std::function< real_t( const Point3D& ) > exact = []( const Point3D& x ) { return x[0] * x[1] * x[1] + real_c( 2 ); };
   uExact.interpolate( exact, level, All );

   // this is helpful to check that nothing bad happens internally for facedofs in PETScLUSolver when assining stuff
   walberla::math::seedRandomGenerator( 0 );
   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) {
      return real_c( walberla::math::realRandom( 0.0, 1.0 ) );
   };
   uDiscr.interpolate( rand, level, Inner );

   // set Dirichlet BCs in discrete solution
   uDiscr.interpolate( exact, level, DirichletBoundary );

   // setup rhs
   WALBERLA_LOG_INFO_ON_ROOT( "-> computing RHS" );
   P3ElementwiseMassOperator massOpr( storage, level, level );

   std::function< real_t( const Point3D& ) > expressionRHS = []( const Point3D& x ) { return -real_c( 2 ) * x[0]; };

   rhs.interpolate( expressionRHS, level, All );
   massOpr.apply( rhs, b, level, All );
   bInitial.assign( { real_c( 1 ) }, { b }, level, All );

   // solve linear system
   WALBERLA_LOG_INFO_ON_ROOT( "-> solving linear system" );

   uint_t localDoFs  = numberOfLocalDoFs< P3FunctionTag >( *storage, level );
   uint_t globalDoFs = numberOfGlobalDoFs< P3FunctionTag >( *storage, level );
   WALBERLA_LOG_INFO( "   localDoFs: " << localDoFs << ", globalDoFs: " << globalDoFs );

   P3ElementwiseDiffusionOperator                  diffOpr( storage, level, level );
   PETScLUSolver< P3ElementwiseDiffusionOperator > solver( storage, level );

   solver.setVerbose( true );
   solver.assumeSymmetry( true );
   solver.solve( diffOpr, uDiscr, b, level );

   bPostInner.assign( { real_c( 1 ) }, { b }, level, Inner );

   WALBERLA_LOG_INFO_ON_ROOT( "-> checking error in discrete solution" );
   err.assign( { 1.0, -1.0 }, { uDiscr, uExact }, level );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs );

   WALBERLA_LOG_INFO_ON_ROOT( "-> discrete L2 error = " << discr_l2_err );

   if ( doVTKOutput )
   {
      std::string fname = "P3PetscSolve";
      WALBERLA_LOG_INFO_ON_ROOT( "-> Exporting functions to VTU file '" << fname << "'" );
      VTKOutput vtkOutput( ".", fname, storage );
      vtkOutput.add( uDiscr );
      vtkOutput.add( uExact );
      vtkOutput.add( b );
      vtkOutput.add( bInitial );
      vtkOutput.add( bPostInner );
      vtkOutput.add( rhs );
      vtkOutput.add( err );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_LESS( discr_l2_err, errEps );
}

void petscMassMatrixApplyTest( const uint_t& level )
{
   WALBERLA_LOG_INFO_ON_ROOT( "\n*** PETSc Mass Matrix Apply Test ***" );
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

   WALBERLA_LOG_INFO_ON_ROOT( "-> Assembling global mass matrix" );
   PETScSparseMatrix< P3ElementwiseMassOperator > matrix( "GlobalMassMatrix" );

   matrix.createMatrixFromOperator( massOp, level, enumerator, All );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Assembling vector of ones" );
   PETScVector< real_t, P3Function > srcVector( vecOfOnes, enumerator, level, All, "src" );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Performing MatVec" );
   PETScVector< real_t, P3Function > dstVector( aux, enumerator, level, All, "dst" );
   auto                              petscErrorCode = MatMult( matrix.get(), srcVector.get(), dstVector.get() );
   WALBERLA_CHECK_EQUAL( petscErrorCode, PETSC_SUCCESS );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Computing dot product" );
   real_t measure = real_c( 0 );
   petscErrorCode = VecDot( srcVector.get(), dstVector.get(), &measure );
   WALBERLA_CHECK_EQUAL( petscErrorCode, PETSC_SUCCESS );

   WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure );
   WALBERLA_CHECK_FLOAT_EQUAL( measure, area );
}

void petscDiffusionMatrixApplyTest( const uint_t& level, DoFType dofFlag, bool doVTKOutput )
{
   std::string tag;
   switch ( dofFlag )
   {
   case Inner:
      tag = "(Inner)";
      break;
   case All:
      tag = "(All)";
      break;
   default:
      WALBERLA_ABORT( "Can only handle dofFlag = Inner or All!" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "\n*** Diffusion Matrix Apply Test " << tag << " ***" );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing with unit square" );
   Point2D  lowerLeft( real_c( 0.0 ), real_c( 0.0 ) );
   Point2D  upperRight( real_c( 1.0 ), real_c( 1.0 ) );
   MeshInfo meshInfo = MeshInfo::meshRectangle( lowerLeft, upperRight, MeshInfo::meshFlavour::CRISS, 1, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up functions and operator" );
   P3ElementwiseDiffusionOperator diffusionOp( storage, level, level );

   P3Function< real_t > src( "Source function", storage, level, level );
   P3Function< real_t > dstWithHyTeG( "HyTeG result", storage, level, level );
   P3Function< real_t > dstWithPETSc( "PETSc result", storage, level, level );

   real_t freqX = real_c( 1 );
   real_t freqY = real_c( 2 );

   std::function< real_t( const hyteg::Point3D& ) > expression = [&freqX, &freqY]( const Point3D& x ) -> real_t {
      using walberla::math::pi;
      real_t retVal = std::sin( freqX * pi * x[0] ) * std::sin( freqY * pi * x[1] ) + real_c( 1 );
      return retVal;
   };

   src.interpolate( expression, level, All );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Assembling global diffusion matrix" );
   P3Function< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   PETScSparseMatrix< P3ElementwiseDiffusionOperator > matrix( "GlobalDiffusionMatrix" );

   matrix.createMatrixFromOperator( diffusionOp, level, enumerator, All );

   PETScVector< real_t, P3Function< real_t >::FunctionType > rhsVector(
       dstWithPETSc, enumerator, level, DirichletBoundary, "rhs" );

   if ( dofFlag == Inner )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "-> Eliminating Dirichlet Boundary Conditions" );
      matrix.applyDirichletBCSymmetrically( src, enumerator, rhsVector, level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "-> Applying HyTeG's elementwise operator" );
   dstWithHyTeG.assign( { real_c( 1 ) }, { src }, level, DirichletBoundary );
   diffusionOp.apply( src, dstWithHyTeG, level, dofFlag );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Converting src function to PETSc vector" );
   PETScVector< real_t, P3Function > srcVector( src, enumerator, level, All, "src" );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Performing MatVec" );
   PETScVector< real_t, P3Function > dstVector( dstWithPETSc, enumerator, level, All, "src" );
   auto                              petscErrorCode = MatMult( matrix.get(), srcVector.get(), dstVector.get() );
   WALBERLA_CHECK_EQUAL( petscErrorCode, PETSC_SUCCESS );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Performing VecAXPY to re-add Dirichlet Values for apply at neighbours" );
   petscErrorCode = VecAXPY( dstVector.get(), -1.0, rhsVector.get() );
   WALBERLA_CHECK_EQUAL( petscErrorCode, PETSC_SUCCESS );

   WALBERLA_LOG_INFO_ON_ROOT( "-> Converting PETSc result vector back" );
   dstVector.createFunctionFromVector( dstWithPETSc, enumerator, level, All );

   if ( dofFlag == Inner )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "-> Re-adding Dirichlet Values to PETSc result" );
      dstWithPETSc.add( { real_c( 1 ) }, { src }, level, DirichletBoundary );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "-> Computing difference between approaches" );
   P3Function< real_t > diff( "Differences", storage, level, level );
   diff.assign( { real_c( 1 ), real_c( -1 ) }, { dstWithHyTeG, dstWithPETSc }, level, All );

   if ( doVTKOutput )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "-> Exporting functions to VTU file" );
      VTKOutput vtkOutput( ".", "P3PetscTest-Diffusion", storage );
      vtkOutput.add( src );
      vtkOutput.add( dstWithHyTeG );
      vtkOutput.add( dstWithPETSc );
      vtkOutput.add( diff );
      vtkOutput.write( level );
   }

   real_t measure = diff.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure );
   WALBERLA_CHECK_LESS_EQUAL( measure, real_c( 1e-14 ) );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   petscMassMatrixApplyTest( 2 );
   petscDiffusionMatrixApplyTest( 5, All, false );
   petscDiffusionMatrixApplyTest( 5, Inner, false );
   petscSolveTest( 3, prependHyTeGMeshDir( "2D/quad_4el.msh" ), 1e-13 );

   return EXIT_SUCCESS;
}
