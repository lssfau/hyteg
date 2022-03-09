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
#include "core/timing/Timer.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

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

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t minLevel  = 2;
   const uint_t maxLevel  = 4;
   const bool   writeVTK  = true;

   //create a Rectangle as mesh with 4 triangles
   auto meshInfo =
       hyteg::MeshInfo::meshRectangle( hyteg::Point2D( {-1, -1} ), hyteg::Point2D( {1, 1} ), hyteg::MeshInfo::CRISSCROSS, 1, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setRightBFSBoundaryNeumannPoiseuille( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   auto                                      storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   hyteg::writeDomainPartitioningVTK( storage, "../../output", "UzawaConvergenceTestDomain" );

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > Au( "Au", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > err( "err", storage, minLevel, maxLevel );

   const auto setUVelocityBC = []( const hyteg::Point3D& x ) -> real_t {
      if ( x[0] < -1.0 + 1e-8 )
      {
         return real_c( 1 - x[1] * x[1] );
      }
      else
      {
         return real_c( 0 );
      }
   };

   const auto solutionU = []( const hyteg::Point3D& x ) -> real_t { return real_c( 1 - x[1] * x[1] ); };

   const auto solutionP = []( const hyteg::Point3D& x ) -> real_t { return real_c( -2.0 * x[0] + 2.0 ); };

   hyteg::P1P1StokesOperator L( storage, minLevel, maxLevel );

   u.uvw()[0].interpolate( setUVelocityBC, maxLevel, hyteg::DirichletBoundary );
   u_exact.uvw()[0].interpolate( solutionU, maxLevel );
   u_exact.p().interpolate( solutionP, maxLevel );

   hyteg::communication::syncFunctionBetweenPrimitives( u_exact.uvw()[0], maxLevel );
   hyteg::communication::syncFunctionBetweenPrimitives( u_exact.p(), maxLevel );

   auto gmgSolver = solvertemplates::stokesGMGUzawaSolver< P1P1StokesOperator >( storage, minLevel, maxLevel, 3, 3, 0.37 );

   const uint_t npoints = hyteg::numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, maxLevel );
   real_t       discr_l2_err_u, discr_l2_err_p, currRes, oldRes = 0;

   L.apply( u, Au, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
   oldRes = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::All ) / real_c( npoints ) );

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel, hyteg::All );

   discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel, hyteg::Inner ) /
                               real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );
   discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel, hyteg::Inner ) /
                               real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Totalpoints         = " << npoints );
   WALBERLA_LOG_INFO_ON_ROOT( "initial Residual    = " << oldRes );
   WALBERLA_LOG_INFO_ON_ROOT( "initial L2 error U  = " << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "initial L2 error P  = " << discr_l2_err_p );

   hyteg::VTKOutput vtkOutput( "../../output", "P1P1UzawaConvergence", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( err );

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   for ( int j = 0; j < 8; ++j )
   {
      gmgSolver->solve( L, u, f, maxLevel );

      L.apply( u, Au, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      currRes = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::All ) /
                           real_c( numberOfGlobalDoFs< P1StokesFunctionTag >( *storage, maxLevel ) ) );

      WALBERLA_LOG_INFO_ON_ROOT( "current Residual = " << currRes );
      // WALBERLA_CHECK_LESS( currRes / oldRes, 0.6 );
      oldRes = currRes;

      err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
      discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel, hyteg::Inner ) /
                                  real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );
      discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel, hyteg::Inner ) /
                                  real_c( numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel ) ) );

      WALBERLA_LOG_INFO_ON_ROOT( "current Err u = " << discr_l2_err_u )
      WALBERLA_LOG_INFO_ON_ROOT( "current Err p = " << discr_l2_err_p )
   }

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel, 1 );
   }

   WALBERLA_CHECK_LESS( currRes, 1.0e-10 );
   WALBERLA_CHECK_LESS( discr_l2_err_u, 4e-3 );
   WALBERLA_CHECK_LESS( discr_l2_err_p, 4e-2 );

   return 0;
}
