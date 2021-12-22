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

#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
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

void setRightBFSBoundaryNeumannPoiseuille( SetupPrimitiveStorage& setupStorage, const real_t & channelLength )
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
      if ( std::fabs( edgeCoordinates[0][0] - channelLength/2 ) < eps && std::fabs( edgeCoordinates[1][0] - channelLength/2 ) < eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }
}



void runBenchmark( const uint_t & level, const uint_t & solverType, const uint_t & channelLength )
{
  
   WALBERLA_LOG_INFO_ON_ROOT( "Poiseuille flow benchmark with channel length " << channelLength);


   /////////////////////////////////////////////////////////////////////////// Domain setup /////////////////////////////////////////////////////////
   //create a Rectangle as mesh with 4 triangles
   real_t halfLength = static_cast<double>(channelLength)/2;
   auto meshInfo = MeshInfo::meshRectangle( Point2D( {-halfLength, -1} ), Point2D( {halfLength, 1} ), MeshInfo::CRISSCROSS, channelLength, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setRightBFSBoundaryNeumannPoiseuille( setupStorage, channelLength );
   
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   auto                                      storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   writeDomainPartitioningVTK( storage, "../../output", "P2P1ChannelTest" );

   P2P1TaylorHoodFunction< real_t > r( "r", storage, 0,level );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, 0, level );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, 0, level );
   P2P1TaylorHoodFunction< real_t > Au( "Au", storage, 0, level );
   P2P1TaylorHoodFunction< real_t > u_exact( "u_exact", storage, 0, level );
   P2P1TaylorHoodFunction< real_t > err( "err", storage, 0, level );

  hyteg::P2P1TaylorHoodStokesOperator A( storage, 0, level );


      const auto setUVelocityBC = [halfLength]( const Point3D& x ) -> real_t {
         if ( x[0] < -halfLength + 1e-8 )
         {
            return real_c( 1 - x[1] * x[1] );//*real_c(sin((fabs(x[1])/2 + 1/2)*pi));
         }
         else
         {
            return real_c( 0 );
         }
      };

      const auto solutionU = []( const Point3D& x ) -> real_t { return real_c( 1 - x[1] * x[1] );};//*real_c(sin((fabs(x[1])/2 + 1/2)*pi)); };

      const auto solutionP = [channelLength]( const Point3D& x ) -> real_t {  return real_c( -2 * x[0] + channelLength); }; // normalize x

      u.uvw[0].interpolate( setUVelocityBC, level, DirichletBoundary );
      u_exact.uvw[0].interpolate( solutionU, level );
      u_exact.p.interpolate( solutionP, level );


    real_t localDoFs1 =static_cast<real_t>(hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level ));
   ////////////////////////////////////////////////////////// VTK /////////////////////////////////////////////////////////////////////
   VTKOutput vtkOutput( "../../output", "P2P1ChannelTest", storage );
   vtkOutput.add( u.uvw );
   vtkOutput.add( u.p );
   vtkOutput.add( u_exact.uvw );
   vtkOutput.add( u_exact.p );
   bool writeVTK = true;
   if ( writeVTK )
   {
      vtkOutput.write( level, 0 );
   }

   ////////////////////////////////////////////////////////// solver setup /////////////////////////////////////////////////////////////
   PETScBlockPreconditionedStokesSolver< P2P1TaylorHoodStokesOperator > GKB_solver(
      storage, level, 1e-6, std::numeric_limits< PetscInt >::max(), 5, 0, 2 );
   PETScBlockPreconditionedStokesSolver< P2P1TaylorHoodStokesOperator > MINRES_solver(
      storage, level, 1e-6, std::numeric_limits< PetscInt >::max(), 3, 1, 0 );
   PETScBlockPreconditionedStokesSolver< P2P1TaylorHoodStokesOperator > FGMRES_solver(
      storage, level, 1e-6, std::numeric_limits< PetscInt >::max(), 5, 1, 1 );


   double uzawaSmootherOmega = 0.4;
   int uzawaPre = 3;
   int uzawaPost = 3;
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< P2P1TaylorHoodStokesOperator, P1LumpedInvMassOperator > >(
           storage, 0,0 );
   auto coarseGridSolver = std::make_shared< MinResSolver< P2P1TaylorHoodStokesOperator > >(
       storage, 0, 0, 1000, 1e-12, pressurePreconditioner );

   auto gaussSeidel        = std::make_shared< GaussSeidelSmoother< P2P1TaylorHoodStokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< P2P1TaylorHoodStokesOperator > >( storage, gaussSeidel );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< P2P1TaylorHoodStokesOperator > >(
       storage, uzawaVelocityPreconditioner, 0, level, uzawaSmootherOmega,  hyteg::Inner | hyteg::NeumannBoundary, 4 );

   auto restriction    = std::make_shared< P2P1StokesToP2P1StokesRestriction >();
   auto prolongation   = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   
   auto VCycle_solver = std::make_shared< GeometricMultigridSolver< P2P1TaylorHoodStokesOperator > >(
       storage, uzawaSmoother, coarseGridSolver, restriction, prolongation, 0, level, uzawaPre, uzawaPost, 4 );
   


   /////////////////////////////////////////////////////////////////// initial error residual /////////////////////////////////////////
   A.apply( u, Au, level, Inner | NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, Au}, level, Inner | NeumannBoundary );
   
   err.assign( {1.0, -1.0}, {u, u_exact}, level, All );

   auto discr_l2_err_u = std::sqrt( err.uvw[0].dotGlobal( err.uvw[0], level, Inner | NeumannBoundary )/  localDoFs1  );
   auto discr_l2_err_v = std::sqrt( err.uvw[1].dotGlobal( err.uvw[1], level, Inner | NeumannBoundary )  /  localDoFs1 );
   auto discr_l2_err_p = std::sqrt( err.p.dotGlobal( err.p, level, Inner | NeumannBoundary )/  localDoFs1   );

  real_t residuum_l2_1, residuum_l2_init, residuum_l2_reduction;
  residuum_l2_reduction = 1e-6;
  r.assign( {1.0, -1.0}, {f, Au}, level,  hyteg::Inner  | hyteg::NeumannBoundary);
  residuum_l2_init  = std::sqrt(r.dotGlobal( r, level, hyteg::Inner | hyteg::NeumannBoundary));

  WALBERLA_LOG_INFO_ON_ROOT( "initial discrete L2 error u = " << discr_l2_err_u );
  WALBERLA_LOG_INFO_ON_ROOT( "initial discrete L2 error v = " << discr_l2_err_v );
  WALBERLA_LOG_INFO_ON_ROOT( "initial discrete L2 error p = " << discr_l2_err_p );
  WALBERLA_LOG_INFO_ON_ROOT( "initial residuum 1 = " << residuum_l2_init );

  ///////////////////////////////////////////////////////////////////////////// solve /////////////////////////////////////////////// 
  residuum_l2_1 = residuum_l2_init;
   switch (solverType) {
     case 1:
      WALBERLA_LOG_INFO("GKB Solver");
        GKB_solver.solve( A, u,f , level );
        break;
     case 2:
      WALBERLA_LOG_INFO("MINRES Solver");
        MINRES_solver.solve( A,u, f, level );
        break;
     case 3:
      {
      WALBERLA_LOG_INFO("VCYCLES with Uzawa smoother");
        int c = 0;
        while( (c < 100 && (residuum_l2_1 / residuum_l2_init > residuum_l2_reduction) )) {
           VCycle_solver->solve( A, u, f, level );      
           
          A.apply( u, Au, level, hyteg::Inner  | hyteg::NeumannBoundary);
          r.assign( {1.0, -1.0}, {f, Au}, level,  hyteg::Inner  | hyteg::NeumannBoundary);
          residuum_l2_1  = std::sqrt( r.dotGlobal( r, level, hyteg::Inner | hyteg::NeumannBoundary )  );
           WALBERLA_LOG_INFO_ON_ROOT( "Cycle " << c << ", residuum 1 = " << residuum_l2_1  );
           ++c;
        }
        break;
      }
      case 4:
      WALBERLA_LOG_INFO("FGMRES Solver");
        FGMRES_solver.solve( A, u, f, level );
        break;
     default:
        WALBERLA_ABORT( "Invalid solver type." )
        break;
  }
   //vertexdof::projectMean(u.p, level);
   //vertexdof::projectMean(u_exact.p, level);

  /////////////////////////////////////////////////////////////////////////////// final error and residual  //////////////////////////////////////////////
   err.assign( {1.0, -1.0}, {u, u_exact}, level, All );
   discr_l2_err_u = std::sqrt( err.uvw[0].dotGlobal( err.uvw[0], level, Inner | NeumannBoundary ) /  localDoFs1 );
   discr_l2_err_v = std::sqrt( err.uvw[1].dotGlobal( err.uvw[1], level, Inner | NeumannBoundary ) /  localDoFs1 );
   discr_l2_err_p = std::sqrt( err.p.dotGlobal( err.p, level, Inner | NeumannBoundary )/  localDoFs1  );
   A.apply( u, Au, level, hyteg::Inner  | hyteg::NeumannBoundary);
   r.assign( {1.0, -1.0}, {f, Au}, level,  hyteg::Inner  | hyteg::NeumannBoundary);
   residuum_l2_1  = std::sqrt( r.dotGlobal( r, level, hyteg::Inner | hyteg::NeumannBoundary )  );
          
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

  /////////////////////////////////////////////////////////////////////////////// write solution //////////////////////////////////////////////
  if ( writeVTK )
   {
      vtkOutput.write( level, 1 );
   }
 
}


int main( int argc, char* argv[] ) {
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

  PETScManager petscManager( &argc, &argv );
   runBenchmark( 5,3,8);

   return 0;
}
