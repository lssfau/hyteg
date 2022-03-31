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
#include <cmath>

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager manager;

   const std::string meshFileName  = "../../data/meshes/3D/cube_24el.msh";
   const uint_t      minLevel      = 2;
   const uint_t      maxLevel      = 5;
   const uint_t      maxIterations = 3;
   const bool        writeVTK      = false;

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::writeDomainPartitioningVTK( storage, "../../output", "P1P1_Stokes_3D_Uzawa_convergence_partitioning" );

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > uExact( "uExact", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );

   hyteg::VTKOutput vtkOutput( "../../output", "P1P1_Stokes_3D_Uzawa_convergence", storage );
   vtkOutput.add( u );
   vtkOutput.add( uExact );
   vtkOutput.add( err );

   hyteg::P1P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > inflowPoiseuille = []( const hyteg::Point3D& x ) {
      if ( x[2] < 1e-8 )
      {
         return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
      }
      else
      {
         return 0.0;
      }
   };

   std::function< real_t( const hyteg::Point3D& ) > solutionPoiseuille = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > collidingFlow_x = []( const hyteg::Point3D& x ) {
      return real_c( 20 ) * x[0] * x[1] * x[1] * x[1];
   };

   std::function< real_t( const hyteg::Point3D& ) > collidingFlow_y = []( const hyteg::Point3D& x ) {
      return real_c( 5 ) * x[0] * x[0] * x[0] * x[0] - real_c( 5 ) * x[1] * x[1] * x[1] * x[1];
   };

   std::function< real_t( const hyteg::Point3D& ) > collidingFlow_p = []( const hyteg::Point3D& xx ) {
      return real_c( 60 ) * std::pow( xx[0], 2.0 ) * xx[1] - real_c( 20 ) * std::pow( xx[1], 3.0 );
   };

   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.uvw().interpolate( {collidingFlow_x, collidingFlow_y, zero}, maxLevel, hyteg::DirichletBoundary );
   uExact.uvw().interpolate( {collidingFlow_x, collidingFlow_y, zero}, maxLevel );
   uExact.p().interpolate( collidingFlow_p, maxLevel );

   if ( writeVTK )
      vtkOutput.write( maxLevel, 0 );

   typedef hyteg::StokesPressureBlockPreconditioner< hyteg::P1P1StokesOperator, hyteg::P1LumpedInvMassOperator >
        PressurePreconditioner_T;
   auto pressurePrec = std::make_shared< PressurePreconditioner_T >( storage, minLevel, maxLevel );

   auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1P1StokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner = std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< hyteg::P1P1StokesOperator > >( storage, gaussSeidel );
   auto smoother =
       std::make_shared< hyteg::UzawaSmoother< hyteg::P1P1StokesOperator > >( storage, uzawaVelocityPreconditioner, minLevel, maxLevel, 0.3 );
   auto restriction      = std::make_shared< hyteg::P1P1StokesToP1P1StokesRestriction >( true );
   auto prolongation     = std::make_shared< hyteg::P1P1StokesToP1P1StokesProlongation >();
   auto coarseGridSolver = std::make_shared< hyteg::PETScLUSolver< hyteg::P1P1StokesOperator > >( storage, minLevel );
   hyteg::GeometricMultigridSolver< hyteg::P1P1StokesOperator > solver(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, 3, 3, 2 );

   const uint_t globalDoFsVelocity = hyteg::numberOfGlobalDoFs< hyteg::P1FunctionTag >( *storage, maxLevel );
   const uint_t globalDoFsPressure = hyteg::numberOfGlobalDoFs< hyteg::P1FunctionTag >( *storage, maxLevel );

   L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
   err.assign( {1.0, -1.0}, {u, uExact}, maxLevel );
   real_t lastResidual =
       std::sqrt( r.dotGlobal( r, maxLevel ) / ( 3 * (real_t) globalDoFsVelocity + real_c( globalDoFsPressure ) ) );

   real_t discr_l2_err_1_u;
   real_t discr_l2_err_1_v;
   real_t discr_l2_err_1_w;
   real_t discr_l2_err_1_p;
   real_t residuum_l2_1;

   for ( uint_t i = 1; i <= maxIterations; i++ )
   {
      solver.solve( L, u, f, maxLevel );

      hyteg::vertexdof::projectMean( u.p(), maxLevel );
      hyteg::vertexdof::projectMean( uExact.p(), maxLevel );

      L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

      err.assign( {1.0, -1.0}, {u, uExact}, maxLevel );

      if ( writeVTK )
      {
         vtkOutput.write( maxLevel, i );
      }

      discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel ) / (real_t) globalDoFsVelocity );
      discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], maxLevel ) / (real_t) globalDoFsVelocity );
      discr_l2_err_1_w = std::sqrt( err.uvw()[2].dotGlobal( err.uvw()[2], maxLevel ) / (real_t) globalDoFsVelocity );
      discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel ) / (real_t) globalDoFsPressure );
      residuum_l2_1 =
          std::sqrt( r.dotGlobal( r, maxLevel ) / ( 3 * (real_t) globalDoFsVelocity + real_c( globalDoFsPressure ) ) );

      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error w = " << discr_l2_err_1_w );
      WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
      WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

      const real_t discrResConvRate = residuum_l2_1 / lastResidual;
      WALBERLA_CHECK_LESS( discrResConvRate, 1.4e-01 )

      lastResidual = residuum_l2_1;

      WALBERLA_LOG_INFO_ON_ROOT( "After " << std::setw( 3 ) << i << " VCycles: Residual: " << std::scientific << residuum_l2_1
                                          << " | convRate: " << discrResConvRate << " | Error L2 u: " << discr_l2_err_1_u
                                          << " | Error L2 p: " << discr_l2_err_1_p );
   }

   auto tt        = storage->getTimingTree();
   auto ttreduced = tt->getReduced().getCopyWithRemainder();
   WALBERLA_LOG_INFO_ON_ROOT( ttreduced );

   nlohmann::json ttjson = nlohmann::json( ttreduced );
   std::ofstream  o( "/tmp/uzawa.json" );
   o << ttjson;
   o.close();

   WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v + discr_l2_err_1_w, 2.8e-03 );
   WALBERLA_CHECK_LESS( discr_l2_err_1_p, 0.13 );
   WALBERLA_CHECK_LESS( residuum_l2_1, 4.0e-06 );

   return EXIT_SUCCESS;
}
