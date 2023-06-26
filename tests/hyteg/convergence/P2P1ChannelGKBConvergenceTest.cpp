/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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
#include "core/timing/Timer.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GKBSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

using walberla::math::pi;

void setRightBFSBoundaryNeumannPoiseuille( SetupPrimitiveStorage& setupStorage, const uint_t& channelLength )
{
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   const real_t eps = 0.001;

   for ( const auto& it : setupStorage.getVertices() )
   {
      if ( std::fabs( it.second->getCoordinates()[0] - 1.0 ) < eps && it.second->getCoordinates()[1] > -1.0 + eps &&
           it.second->getCoordinates()[1] < 1.0 - eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }

   for ( const auto& it : setupStorage.getEdges() )
   {
      const auto edgeCoordinates = it.second->getCoordinates();
      if ( std::fabs( edgeCoordinates[0][0] - static_cast< real_t >( channelLength ) / 2 ) < eps &&
           std::fabs( edgeCoordinates[1][0] - static_cast< real_t >( channelLength ) / 2 ) < eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }
}

void runBenchmark( const uint_t& level,
                   const uint_t& channelLength,
                   const uint_t& solver,
                   const real_t& resEps,
                   const real_t& errEpsU,
                   const real_t& errEpsP )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Poiseuille flow benchmark with channel length " << channelLength );

   /////////////////////////////////////////////////////////////////////////// Domain setup /////////////////////////////////////////////////////////
   //create a Rectangle as mesh with 4 triangles
   real_t halfLength = static_cast< real_t >( channelLength ) / 2;
   auto   meshInfo   = MeshInfo::meshRectangle(
       Point2D(  -halfLength, -1  ), Point2D(  halfLength, 1  ), MeshInfo::CRISSCROSS, channelLength, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setRightBFSBoundaryNeumannPoiseuille( setupStorage, channelLength );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   auto                                      storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   writeDomainPartitioningVTK( storage, "../../output", "P2P1ChannelTest" );

   P2P1TaylorHoodFunction< real_t > r( "r", storage, level, level );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, level, level );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, level, level );
   P2P1TaylorHoodFunction< real_t > Au( "Au", storage, level, level );
   P2P1TaylorHoodFunction< real_t > u_exact( "u_exact", storage, level, level );
   P2P1TaylorHoodFunction< real_t > err( "err", storage, level, level );

   hyteg::P2P1TaylorHoodStokesOperator A( storage, level, level );

   const auto setUVelocityBC = [halfLength]( const Point3D& x ) -> real_t {
      if ( x[0] < -halfLength + 1e-8 )
      {
         return real_c( ( 1 - x[1] * x[1] ) );
      }
      else
      {
         return real_c( 0 );
      }
   };

   const auto solutionU = []( const Point3D& x ) -> real_t { return real_c( 1 - x[1] * x[1] ); };

   const auto solutionP = [channelLength]( const Point3D& x ) -> real_t {
      return real_c( -2 * x[0] + static_cast< real_t >( channelLength ) );
   };

   u_exact.uvw()[0].interpolate( solutionU, level );
   u_exact.p().interpolate( solutionP, level );
   u.uvw()[0].interpolate( setUVelocityBC, level, DirichletBoundary );

   real_t localDoFs1 = static_cast< real_t >( hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level ) );
   ////////////////////////////////////////////////////////// VTK /////////////////////////////////////////////////////////////////////
   VTKOutput vtkOutput( "../../output", "P2P1ChannelTest", storage );
   bool      writeVTK = false;
   if ( writeVTK )
   {
      vtkOutput.add( u.uvw() );
      vtkOutput.add( u.p() );
      vtkOutput.add( u_exact.uvw() );
      vtkOutput.add( u_exact.p() );
      vtkOutput.write( level, 0 );
   }

   ////////////////////////////////////////////////////////// solver setup /////////////////////////////////////////////////////////////

   GKBSolver_P2P1TH GKB_HOUSE_solver(
       storage, level, CGSolver< ALOP_P2P1TH >( storage, level, level, 1000, 1e-11 ), 0, 100, 1e-10 );

   PETScBlockPreconditionedStokesSolver< P2P1TaylorHoodStokesOperator > GKB_PETSC_solver( storage, level, 1e-10, 1000, 5, 0, 2 );

   /////////////////////////////////////////////////////////////////// initial error residual /////////////////////////////////////////
   A.apply( u, Au, level, Inner | NeumannBoundary );

   err.assign( { 1.0, -1.0 }, { u, u_exact }, level, All );

   auto discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, Inner | NeumannBoundary ) / localDoFs1 );
   auto discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, Inner | NeumannBoundary ) / localDoFs1 );
   auto discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), level, Inner | NeumannBoundary ) / localDoFs1 );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_p );

   real_t residuum_l2_1;
   r.assign( { 1.0, -1.0 }, { f, Au }, level, hyteg::Inner | hyteg::NeumannBoundary );

   ///////////////////////////////////////////////////////////////////////////// solve ///////////////////////////////////////////////
   walberla::WcTimer timer;
   if ( solver == 0 )
   {
      A.apply( u, Au, level, hyteg::Inner | hyteg::NeumannBoundary );
      f.assign( { 1, -1 }, { f, Au }, level, hyteg::Inner | hyteg::NeumannBoundary );
      GKB_HOUSE_solver.solve( A, u, f, level );
      f.assign( { 0 }, { f }, level, hyteg::Inner | hyteg::NeumannBoundary );
   }
   else if ( solver == 1 )
   {
      GKB_PETSC_solver.solve( A, u, f, level );
   }
   timer.end();

   /////////////////////////////////////////////////////////////////////////////// final error and residual  //////////////////////////////////////////////

   err.assign( { 1.0, -1.0 }, { u, u_exact }, level, All );
   discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, Inner | NeumannBoundary ) / localDoFs1 );
   discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, Inner | NeumannBoundary ) / localDoFs1 );
   discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), level, Inner | NeumannBoundary ) / localDoFs1 );

   A.apply( u, Au, level, hyteg::Inner | hyteg::NeumannBoundary );
   r.assign( { 1.0, -1.0 }, { f, Au }, level, hyteg::Inner | hyteg::NeumannBoundary );
   residuum_l2_1 = std::sqrt( r.dotGlobal( r, level, hyteg::Inner | hyteg::NeumannBoundary ) );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );
   WALBERLA_CHECK_LESS( residuum_l2_1, resEps );
   WALBERLA_CHECK_LESS( discr_l2_err_u, errEpsU );
   WALBERLA_CHECK_LESS( discr_l2_err_p, errEpsP );

   /////////////////////////////////////////////////////////////////////////////// write solution //////////////////////////////////////////////
   if ( writeVTK )
   {
      vtkOutput.write( level, 1 );
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );
   runBenchmark( 2, 8, 0, 1e-8, 1e-9, 1e-9 );
   runBenchmark( 2, 8, 1, 1e-8, 1e-9, 1e-9 );

   return 0;
}
