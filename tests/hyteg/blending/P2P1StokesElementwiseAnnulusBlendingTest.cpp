/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
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

std::function< real_t( const hyteg::Point3D& ) > radius = []( const hyteg::Point3D& x ) {
   return std::sqrt( x[0] * x[0] + x[1] * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > angle = []( const hyteg::Point3D& x ) { return std::atan2( x[1], x[0] ); };

/// \brief Benchmark taken from Master's thesis of Leonard Schlag.
void runTest( bool preCompute )
{
   WALBERLA_LOG_INFO_ON_ROOT( "P2-P1 Stokes elementwise blending operator on annulus with Uzawa GMG solver" )

   const uint_t minLevel       = 2;
   const uint_t maxLevel       = 3;
   const bool   writeVTK       = false;
   const uint_t uzawaPre       = 10;
   const uint_t uzawaPost      = 10;
   const uint_t innerJacSmooth = 4;
   const real_t uzawaOmega     = 0.37;
   const real_t jacobiOmega    = 0.66;
   const uint_t numIterations  = 2;

   const uint_t nTan = 12;
   const uint_t nRad = 2;

   const real_t rMin = pi;
   const real_t rMax = 2.0 * pi;
   const real_t k    = 8.0;

   WALBERLA_LOG_INFO_ON_ROOT( "Domain parameters:" )

   WALBERLA_LOG_INFO_ON_ROOT( "- min level:           " << minLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "- max level:           " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "- nTan:                " << nTan )
   WALBERLA_LOG_INFO_ON_ROOT( "- nRad:                " << nRad )
   WALBERLA_LOG_INFO_ON_ROOT( "- rMin:                " << rMin )
   WALBERLA_LOG_INFO_ON_ROOT( "- rMax:                " << rMax )

   WALBERLA_LOG_INFO_ON_ROOT( "Solver parameters:" )

   WALBERLA_LOG_INFO_ON_ROOT( "- Uzawa pre smooth:    " << uzawaPre )
   WALBERLA_LOG_INFO_ON_ROOT( "- Uzawa post smooth:   " << uzawaPost )
   WALBERLA_LOG_INFO_ON_ROOT( "- inner Jacobi smooth: " << innerJacSmooth )
   WALBERLA_LOG_INFO_ON_ROOT( "- Uzawa omega:         " << uzawaOmega )
   WALBERLA_LOG_INFO_ON_ROOT( "- Jacobi omega:        " << jacobiOmega )

   auto                  meshInfo = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CRISS, nTan, nRad );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   AnnulusMap::setMap( setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "P2P1StokesElementwiseAnnulusBlendingTest_Domain" );

   std::function< real_t( const hyteg::Point3D& ) > uSolution = [&]( const hyteg::Point3D& x ) {
      const auto r       = radius( x );
      const auto phi     = angle( x );
      const auto sinPhi  = std::sin( phi );
      const auto cosPhi  = std::cos( phi );
      const auto sinKPhi = std::sin( k * phi );
      const auto cosKPhi = std::cos( k * phi );
      const auto sinR    = std::sin( r );
      const auto cosR    = std::cos( r );
      return cosPhi * ( k / r ) * sinR * cosKPhi + sinPhi * cosR * sinKPhi;
   };

   std::function< real_t( const hyteg::Point3D& ) > vSolution = [&]( const hyteg::Point3D& x ) {
      const auto r       = radius( x );
      const auto phi     = angle( x );
      const auto sinPhi  = std::sin( phi );
      const auto cosPhi  = std::cos( phi );
      const auto sinKPhi = std::sin( k * phi );
      const auto cosKPhi = std::cos( k * phi );
      const auto sinR    = std::sin( r );
      const auto cosR    = std::cos( r );
      return sinPhi * ( k / r ) * sinR * cosKPhi - cosPhi * cosR * sinKPhi;
   };

   std::function< real_t( const hyteg::Point3D& ) > pSolution = [&]( const hyteg::Point3D& x ) {
      const auto r       = radius( x );
      const auto phi     = angle( x );
      const auto cosKPhi = std::cos( k * phi );
      const auto sinR    = std::sin( r );
      const auto cosR    = std::cos( r );
      return ( 2.0 / ( r * r ) ) * k * sinR * cosKPhi - ( 1.0 / k ) * sinR * cosKPhi -
             ( ( k * k + r * r + 1.0 ) / ( k * r ) ) * cosR * cosKPhi;
   };

   std::function< real_t( const hyteg::Point3D& ) > TSolution = [&]( const hyteg::Point3D& x ) {
      const auto r       = radius( x );
      const auto phi     = angle( x );
      const auto cosKPhi = std::cos( k * phi );
      const auto sinR    = std::sin( r );
      const auto cosR    = std::cos( r );
      const auto velU    = ( k / r ) * sinR * cosKPhi;
      return ( k * k + r * r - 4.0 + ( ( r * r ) / ( k * k ) ) * ( k * k + r * r + 1.0 ) ) * ( velU / ( r * r ) ) +
             ( k + ( ( k * k - 2.0 * r * r + 1.0 ) / ( k ) ) ) * ( ( cosR * cosKPhi ) / ( r * r ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > uRhs = [&]( const hyteg::Point3D& x ) {
      const auto phi    = angle( x );
      const auto cosPhi = std::cos( phi );
      return cosPhi * TSolution( x );
   };

   std::function< real_t( const hyteg::Point3D& ) > vRhs = [&]( const hyteg::Point3D& x ) {
      const auto phi    = angle( x );
      const auto sinPhi = std::sin( phi );
      return sinPhi * TSolution( x );
   };

   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f_strong( "f_strong", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > Au( "Au", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > err( "err", storage, minLevel, maxLevel );
   P2Function< real_t >             T( "T", storage, minLevel, maxLevel );

   typedef P2P1ElementwiseBlendingStokesOperator StokesOperator;
   typedef P2ElementwiseBlendingMassOperator     MassOperator;

   StokesOperator L( storage, minLevel, maxLevel );
   if ( preCompute )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " Precomputing local element Matrices." )
      L.computeAndStoreLocalElementMatrices();
   }
   MassOperator M( storage, maxLevel, maxLevel );

   u.uvw().interpolate( { uSolution, vSolution }, maxLevel, DirichletBoundary );
   f_strong.uvw().interpolate( { uRhs, vRhs }, maxLevel, All );
   T.interpolate( TSolution, maxLevel, All );

   M.apply( f_strong.uvw()[0], f.uvw()[0], maxLevel, All );
   M.apply( f_strong.uvw()[1], f.uvw()[1], maxLevel, All );

   Au.uvw()[0].setToZero( maxLevel );
   Au.uvw()[1].setToZero( maxLevel );
   Au.p().setToZero( maxLevel );

   u_exact.uvw()[0].interpolate( uSolution, maxLevel, All );
   u_exact.uvw()[1].interpolate( vSolution, maxLevel, All );
   u_exact.p().interpolate( pSolution, maxLevel, All );

   communication::syncP2FunctionBetweenPrimitives( u_exact.uvw()[0], maxLevel );
   communication::syncP2FunctionBetweenPrimitives( u_exact.uvw()[1], maxLevel );
   communication::syncFunctionBetweenPrimitives( u_exact.p(), maxLevel );

   auto coarseGridSolver = solvertemplates::stokesMinResSolver< StokesOperator >( storage, minLevel, 1e-08, 500 );

   auto fineGridSolver = solvertemplates::stokesMinResSolver< StokesOperator >( storage, maxLevel, 1e-08, 100 );

   auto restriction    = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );
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

   VTKOutput vtkOutput( "../../output", "P2P1StokesElementwiseAnnulusBlendingTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( err );
   vtkOutput.add( f.uvw() );
   vtkOutput.add( f_strong.uvw() );
   vtkOutput.add( T );

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

      WALBERLA_CHECK_LESS( residualReduction, 2.0e-01 );

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

   WALBERLA_CHECK_LESS( currRes, 1.0e-06 );
   WALBERLA_CHECK_LESS( discr_l2_err_u, 4.0e-04 );
   WALBERLA_CHECK_LESS( discr_l2_err_v, 4.0e-04 );
   WALBERLA_CHECK_LESS( discr_l2_err_p, 4.0e-02 );

   WALBERLA_LOG_INFO_ON_ROOT( "Solving with MINRES ..." )

   fineGridSolver->solve( L, u, f, maxLevel );

   L.apply( u, Au, maxLevel, Inner | NeumannBoundary );
   r.assign( { 1.0, -1.0 }, { f, Au }, maxLevel, Inner | NeumannBoundary );
   currRes = std::sqrt( r.dotGlobal( r, maxLevel, All ) ) / real_c( npoints );

   err.assign( { 1.0, -1.0 }, { u, u_exact }, maxLevel );
   discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel, Inner | NeumannBoundary ) /
                               real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], maxLevel, Inner | NeumannBoundary ) /
                               real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel, Inner | NeumannBoundary ) /
                               real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "%15s|%15e|%15s|%15e|%15e|%15e", "MINRES", currRes, "-", discr_l2_err_u, discr_l2_err_v, discr_l2_err_p ) )

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, numIterations + 1 );
   }

   WALBERLA_CHECK_LESS( currRes, 5.0e-9 );
   WALBERLA_CHECK_LESS( discr_l2_err_u, 4.0e-04 );
   WALBERLA_CHECK_LESS( discr_l2_err_v, 4.0e-04 );
   WALBERLA_CHECK_LESS( discr_l2_err_p, 4.0e-02 );
}
int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   runTest( false );
   runTest( true );

   return 0;
}
