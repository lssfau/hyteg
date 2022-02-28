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

void setRightBFSBoundaryNeumannPoiseuille( SetupPrimitiveStorage& setupStorage )
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
      if ( std::fabs( edgeCoordinates[0][0] - 1.0 ) < eps && std::fabs( edgeCoordinates[1][0] - 1.0 ) < eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }
}

std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& x ) {
   return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) {
   return 2.5 * pi * std::cos( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) { return 0; };

std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& x ) {
   return -12.5 * pi * pi * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};

void runBenchmark( uint_t benchmark )
{
   Point2D leftBottom( {0, 0} );
   if ( benchmark == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Poiseuille flow benchmark" )
      leftBottom = Point2D( {-1, -1} );
   }
   else if ( benchmark == 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Convection cell benchmark" )
   }
   else
      WALBERLA_ABORT( "Invalid benchmark" )

   const uint_t minLevel       = 2;
   const uint_t maxLevel       = 4;
   const bool   writeVTK       = false;
   const uint_t uzawaPre       = 6;
   const uint_t uzawaPost      = 6;
   const uint_t innerJacSmooth = 2;
   const real_t uzawaOmega     = 0.37;
   const real_t jacobiOmega    = 0.66;
   const uint_t numIterations  = 8;

   WALBERLA_LOG_INFO_ON_ROOT( "Elementwise Stokes CC Uzawa GMG w/ Jacobi smoothing" )
   WALBERLA_LOG_INFO_ON_ROOT( "min level:           " << minLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "max level:           " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa pre smooth:    " << uzawaPre )
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa post smooth:   " << uzawaPost )
   WALBERLA_LOG_INFO_ON_ROOT( "inner Jacobi smooth: " << innerJacSmooth )
   WALBERLA_LOG_INFO_ON_ROOT( "Uzawa omega:         " << uzawaOmega )
   WALBERLA_LOG_INFO_ON_ROOT( "Jacobi omega:        " << jacobiOmega )

   //create a Rectangle as mesh with 4 triangles
   auto meshInfo = MeshInfo::meshRectangle( leftBottom, Point2D( {1, 1} ), MeshInfo::CRISSCROSS, 1, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( benchmark == 0 )
   {
      setRightBFSBoundaryNeumannPoiseuille( setupStorage );
   }
   else if ( benchmark == 1 )
   {
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   }

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   auto                                      storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   writeDomainPartitioningVTK( storage, "../../output", "UzawaConvergenceTestDomain" );

   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > Au( "Au", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > err( "err", storage, minLevel, maxLevel );

   typedef P2P1ElementwiseConstantCoefficientStokesOperator StokesOperator;
   typedef P2ElementwiseMassOperator                        MassOperator;

   StokesOperator L( storage, minLevel, maxLevel );
   MassOperator   M( storage, maxLevel, maxLevel );

   if ( benchmark == 0 )
   {
      const auto setUVelocityBC = []( const Point3D& x ) -> real_t {
         if ( x[0] < -1.0 + 1e-8 )
         {
            return real_c( 1 - x[1] * x[1] );
         }
         else
         {
            return real_c( 0 );
         }
      };

      const auto solutionU = []( const Point3D& x ) -> real_t { return real_c( 1 - x[1] * x[1] ); };

      const auto solutionP = []( const Point3D& x ) -> real_t { return real_c( -2.0 * x[0] + 2.0 ); };

      u.uvw()[0].interpolate( setUVelocityBC, maxLevel, DirichletBoundary );
      u_exact.uvw()[0].interpolate( solutionU, maxLevel );
      u_exact.p().interpolate( solutionP, maxLevel );
   }
   else if ( benchmark == 1 )
   {
      u.uvw().interpolate( { exactU, exactV }, maxLevel, DirichletBoundary );

      Au.uvw().interpolate( { rhsU, rhsV }, maxLevel, All );

      M.apply( Au.uvw()[0], f.uvw()[0], maxLevel, All );
      M.apply( Au.uvw()[1], f.uvw()[1], maxLevel, All );

      Au.uvw()[0].setToZero( maxLevel );
      Au.uvw()[1].setToZero( maxLevel );
      Au.p().setToZero( maxLevel );

      u_exact.uvw().interpolate( { exactU, exactV }, maxLevel, All );
      u_exact.p().interpolate( exactP, maxLevel, All );
   }

   // communication::syncP2FunctionBetweenPrimitives( u_exact.uvw()[0], maxLevel );
   // communication::syncP2FunctionBetweenPrimitives( u_exact.uvw()[1], maxLevel );
   communication::syncVectorFunctionBetweenPrimitives( u_exact.uvw(), maxLevel );
   communication::syncFunctionBetweenPrimitives( u_exact.p(), maxLevel );

   auto coarseGridSolver = solvertemplates::stokesMinResSolver< StokesOperator >( storage, minLevel, 1e-12, 1000 );

   auto fineGridSolver = solvertemplates::stokesMinResSolver< StokesOperator >( storage, maxLevel, 1e-12, 1000 );

   auto restriction    = std::make_shared< P2P1StokesToP2P1StokesRestriction >( benchmark == 1 );
   auto prolongation   = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto jacobiSmoother = std::make_shared< WeightedJacobiSmoother< StokesOperator::VelocityOperator_T > >(
       storage, minLevel, maxLevel, jacobiOmega );
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, jacobiSmoother );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< StokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaOmega, Inner | NeumannBoundary, innerJacSmooth );
   auto gmgSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >(
       storage, uzawaSmoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, uzawaPre, uzawaPost, 2 );

   const uint_t npoints = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   real_t       currRes = 0, oldRes = 0;

   L.apply( u, Au, maxLevel, Inner | NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, Au}, maxLevel, Inner | NeumannBoundary );
   oldRes = std::sqrt( r.dotGlobal( r, maxLevel, All ) ) / real_c( npoints );

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel, All );

   auto discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel, Inner | NeumannBoundary ) /
                                    real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   auto discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], maxLevel, Inner | NeumannBoundary ) /
                                    real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   auto discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel, Inner | NeumannBoundary ) /
                                    real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "%15s|%15s|%15s|%15s|%15s|%15s", "iteration", "residual", "residual redct", "L2 error u", "L2 error v", "L2 error p" ) )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "%15s|%15e|%15s|%15e|%15e|%15e", "initial", oldRes, "-", discr_l2_err_u, discr_l2_err_v, discr_l2_err_p ) )

   VTKOutput vtkOutput( "../../output", "P2P1ElementwiseUzawaConvergence", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( err );

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   for ( uint_t j = 0; j < numIterations; ++j )
   {
      gmgSolver->solve( L, u, f, maxLevel );

      L.apply( u, Au, maxLevel, Inner | NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, Au}, maxLevel, Inner | NeumannBoundary );
      currRes = std::sqrt( r.dotGlobal( r, maxLevel, All ) ) / real_c( npoints );

      auto residualReduction = currRes / oldRes;
      WALBERLA_CHECK_LESS( residualReduction, 0.2 );
      oldRes = currRes;

      err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
      discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel, Inner | NeumannBoundary ) /
                                  real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
      discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], maxLevel, Inner | NeumannBoundary ) /
                                  real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
      discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel, Inner | NeumannBoundary ) /
                                  real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          "%15d|%15e|%15e|%15e|%15e|%15e", j + 1, currRes, residualReduction, discr_l2_err_u, discr_l2_err_v, discr_l2_err_p ) )

      if ( writeVTK )
      {
         vtkOutput.write( maxLevel, j + 1 );
      }
   }

   if ( benchmark == 0 )
   {
      WALBERLA_CHECK_LESS( currRes, 1.0e-10 );
      WALBERLA_CHECK_LESS( discr_l2_err_u, 5e-09 );
      WALBERLA_CHECK_LESS( discr_l2_err_v, 5e-09 );
      WALBERLA_CHECK_LESS( discr_l2_err_p, 3e-07 );
   }
   else if ( benchmark == 1 )
   {
      WALBERLA_CHECK_LESS( currRes, 1.0e-10 );
      WALBERLA_CHECK_LESS( discr_l2_err_u, 3.0e-05 );
      WALBERLA_CHECK_LESS( discr_l2_err_v, 5.0e-05 );
      WALBERLA_CHECK_LESS( discr_l2_err_p, 4.0e-02 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Solving down with MINRES ..." )

   // No idea why this is necessary to make MINRES converge properly here...
   u.p().interpolate( 0, maxLevel, All );

   fineGridSolver->solve( L, u, f, maxLevel );

   L.apply( u, Au, maxLevel, Inner | NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, Au}, maxLevel, Inner | NeumannBoundary );
   currRes = std::sqrt( r.dotGlobal( r, maxLevel, All ) ) / real_c( npoints );

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
   discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel, Inner | NeumannBoundary ) /
                               real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], maxLevel, Inner | NeumannBoundary ) /
                               real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel, Inner | NeumannBoundary ) /
                               real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "%15s|%15e|%15s|%15e|%15e|%15e", "MINRES", currRes, "-", discr_l2_err_u, discr_l2_err_v, discr_l2_err_p ) )

   if ( benchmark == 0 )
   {
      WALBERLA_CHECK_LESS( currRes, 1.0e-14 );
      WALBERLA_CHECK_LESS( discr_l2_err_u, 1e-12 );
      WALBERLA_CHECK_LESS( discr_l2_err_v, 1e-12 );
      WALBERLA_CHECK_LESS( discr_l2_err_p, 2e-12 );
   }
   else if ( benchmark == 1 )
   {
      WALBERLA_CHECK_LESS( currRes, 1.0e-14 );
      WALBERLA_CHECK_LESS( discr_l2_err_u, 3.0e-05 );
      WALBERLA_CHECK_LESS( discr_l2_err_v, 5.0e-05 );
      WALBERLA_CHECK_LESS( discr_l2_err_p, 4.0e-02 );
   }

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, numIterations + 1 );
   }
}
int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   runBenchmark( 0 );
   runBenchmark( 1 );

   return 0;
}
