/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GKBSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "mixed_operator//P2P1TaylorHoodStokesOperator.hpp"
#include "mixed_operator/VectorLaplaceOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( const uint_t&   level,
                     const MeshInfo& meshInfo,
                     const uint_t&   solver,
                     const real_t&   nu,
                     const real_t&   resEps,
                     const real_t&   errEpsUSum,
                     const real_t&   errEpsP )
{
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2P1Stokes2DPetscSolve_Domain" );

   hyteg::P2P1TaylorHoodFunction< real_t > x( "x", storage, 0, level );
   hyteg::P2P1TaylorHoodFunction< real_t > x_exact( "x_exact", storage, 0, level );
   hyteg::P2P1TaylorHoodFunction< real_t > btmp( "btmp", storage, 0, level );
   hyteg::P2P1TaylorHoodFunction< real_t > b( "b", storage, 0, level );
   hyteg::P2P1TaylorHoodFunction< real_t > err( "err", storage, 0, level );
   hyteg::P2P1TaylorHoodFunction< real_t > residuum( "res", storage, 0, level );
   hyteg::P2P1TaylorHoodFunction< real_t > Ax( "Ax", storage, 0, level );

   hyteg::P2P1TaylorHoodStokesOperator A( storage, 0, level );

   std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& xx ) {
      return real_c( 20 ) * xx[0] * xx[1] * xx[1] * xx[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& xx ) {
      return real_c( 5 ) * xx[0] * xx[0] * xx[0] * xx[0] - real_c( 5 ) * xx[1] * xx[1] * xx[1] * xx[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& xx ) {
      return real_c( 60 ) * std::pow( xx[0], 2.0 ) * xx[1] - real_c( 20 ) * std::pow( xx[1], 3.0 );
   };
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return real_c( 0 ); };

   x.uvw().interpolate( { exactU, exactV }, level, hyteg::DirichletBoundary );
   x.p().interpolate( { exactP }, level, hyteg::DirichletBoundary );
   x_exact.uvw().interpolate( { exactU, exactV }, level );
   x_exact.p().interpolate( exactP, level );
   b.uvw().interpolate( { exactU, exactV }, level, DirichletBoundary );

   A.apply( x, residuum, level, hyteg::Inner | hyteg::NeumannBoundary );

   err.assign( { 1.0, -1.0 }, { x, x_exact }, level );
   uint_t localDoFs1  = hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );

   WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );

   GKBSolver_P2P1TH GKB_HOUSE_solver(
       storage, level, CGSolver< ALOP_P2P1TH >( storage, level, level, 1000, 1e-11 ), nu, 100, 1e-10 );

   PETScBlockPreconditionedStokesSolver< P2P1TaylorHoodStokesOperator > GKB_PETSC_solver( storage, level, 1e-10, 1000, 5, 0, 2 );

   ///////////////////////////////////////////////////////////////////////////// solve ///////////////////////////////////////////////
   real_t discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), level ) / (real_t) globalDoFs1 );
   real_t residuum_l2_1    = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

   walberla::WcTimer timer;
   if ( solver == 0 )
   {
      A.apply( x, residuum, level, hyteg::Inner | hyteg::NeumannBoundary );
      b.assign( { 1, -1 }, { b, residuum }, level, hyteg::Inner | hyteg::NeumannBoundary );
      GKB_HOUSE_solver.solve( A, x, b, level );
   }
   else if ( solver == 1 )
   {
      GKB_PETSC_solver.solve( A, x, b, level );
   }
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );

   hyteg::vertexdof::projectMean( x.p(), level );
   hyteg::vertexdof::projectMean( x_exact.p(), level );

   A.apply( x, residuum, level, hyteg::Inner | hyteg::NeumannBoundary );

   err.assign( { 1.0, -1.0 }, { x, x_exact }, level );

   discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level ) / (real_t) globalDoFs1 );
   discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level ) / (real_t) globalDoFs1 );
   discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), level ) / (real_t) globalDoFs1 );
   residuum_l2_1    = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

   //  vtkOutput.write( level, 1 );

   WALBERLA_CHECK_LESS( residuum_l2_1, resEps );
   WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v, errEpsUSum );
   WALBERLA_CHECK_LESS( discr_l2_err_1_p, errEpsP );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   petscSolveTest(
       4, hyteg::MeshInfo::fromGmshFile( "../../meshes/quad_center_at_origin_4el.msh" ), 0, 0, 2.2e-09, 0.00033, 0.0184 );
   petscSolveTest(
       4, hyteg::MeshInfo::fromGmshFile( "../../meshes/quad_center_at_origin_4el.msh" ), 1, 0, 2.2e-09, 0.00033, 0.0184 );
   //TODO AL with nu = 100 converges faster for velocity, but still has a problem for the pressure
   petscSolveTest(
       4, hyteg::MeshInfo::fromGmshFile( "../../meshes/quad_center_at_origin_4el.msh" ), 0, 100, 10, 0.00033, 10 );

   return EXIT_SUCCESS;
}