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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

template < typename P2P1P1StokesOperator >
void stokesMinResConvergenceTest()
{
   const std::string meshFileName  = "../../data/meshes/3D/cube_24el.msh";
   const uint_t      minLevel      = 2;
   const uint_t      maxLevel      = 2;
   const uint_t      maxIterations = 20;
   const real_t      tolerance     = 1e-13;

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

#if 0
   // Get all primitive IDs of primitives at the outflow boundary (z == 1)
  const real_t eps = 1e-8;
  const real_t zBoundary = 1.0;
  for ( const auto & it : setupStorage.getVertices() )
  {
    if ( std::fabs( it.second->getCoordinates()[2] - zBoundary ) < eps &&
         std::fabs( it.second->getCoordinates()[0] ) > eps &&
         std::fabs( it.second->getCoordinates()[1] ) > eps &&
         std::fabs( it.second->getCoordinates()[0] - 1.0 ) > eps &&
         std::fabs( it.second->getCoordinates()[1] - 1.0 ) > eps )
      setupStorage.setMeshBoundaryFlag( it.first, 2 );
  }
  for ( const auto & it : setupStorage.getEdges() )
  {
    if ( std::fabs( it.second->getCoordinates()[0][2] - zBoundary ) < eps &&
         std::fabs( it.second->getCoordinates()[1][2] - zBoundary ) < eps )

      setupStorage.setMeshBoundaryFlag( it.first, 2 );

    if ( std::fabs( it.second->getCoordinates()[0][0] ) < eps &&
         std::fabs( it.second->getCoordinates()[1][0] ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );

    if ( std::fabs( it.second->getCoordinates()[0][0] - 1.0 ) < eps &&
         std::fabs( it.second->getCoordinates()[1][0] - 1.0 ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );

    if ( std::fabs( it.second->getCoordinates()[0][1] ) < eps &&
         std::fabs( it.second->getCoordinates()[1][1] ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );

    if ( std::fabs( it.second->getCoordinates()[0][1] -1.0 ) < eps &&
         std::fabs( it.second->getCoordinates()[1][1] -1.0 ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 1 );
  }
  for ( const auto & it : setupStorage.getFaces() )
  {
    if ( std::fabs( it.second->getCoordinates()[0][2] - zBoundary ) < eps &&
         std::fabs( it.second->getCoordinates()[1][2] - zBoundary ) < eps &&
         std::fabs( it.second->getCoordinates()[2][2] - zBoundary ) < eps )
      setupStorage.setMeshBoundaryFlag( it.first, 2 );
  }
#endif

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::writeDomainPartitioningVTK( storage, "../../output", "P2P1_Stokes_3D_MinRes_convergence_partitioning" );

   hyteg::P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P2P1TaylorHoodFunction< real_t > uExact( "uExact", storage, minLevel, maxLevel );
   hyteg::P2P1TaylorHoodFunction< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P2P1TaylorHoodFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );

//   hyteg::VTKOutput vtkOutput( "../../output", "P2P1_Stokes_3D_MinRes_convergence", storage );
//   vtkOutput.add( u.u );
//   vtkOutput.add( u.v );
//   // vtkOutput.add( u.w );
//   vtkOutput.add( u.p() );
//   vtkOutput.add( uExact.u );
//   vtkOutput.add( uExact.v );
//   // vtkOutput.add( uExact.w );
//   vtkOutput.add( uExact.p() );
//   vtkOutput.add( err.u );
//   vtkOutput.add( err.v );
//   // vtkOutput.add( err.w );
//   vtkOutput.add( err.p() );

   P2P1P1StokesOperator L( storage, minLevel, maxLevel );

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

#if 0
   u.w.interpolate( inflowPoiseuille, maxLevel, hyteg::DirichletBoundary );
#else
   u.uvw().interpolate( {collidingFlow_x, collidingFlow_y, zero}, maxLevel, hyteg::DirichletBoundary );
   uExact.uvw().interpolate( {collidingFlow_x, collidingFlow_y, zero}, maxLevel );
   uExact.p().interpolate( collidingFlow_p, maxLevel );
#endif

   //   vtkOutput.write( maxLevel, 0 );

   typedef hyteg::StokesPressureBlockPreconditioner< P2P1P1StokesOperator, hyteg::P1LumpedInvMassOperator >
        PressurePreconditioner_T;
   auto pressurePrec = std::make_shared< PressurePreconditioner_T >( storage, minLevel, maxLevel );

   auto solver = hyteg::MinResSolver< P2P1P1StokesOperator >( storage, minLevel, maxLevel, maxIterations, tolerance, pressurePrec );

   solver.solve( L, u, f, maxLevel );

   hyteg::vertexdof::projectMean( u.p(), maxLevel );
   hyteg::vertexdof::projectMean( uExact.p(), maxLevel );

   L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

   err.assign( {1.0, -1.0}, {u, uExact}, maxLevel );

   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< hyteg::P2P1TaylorHoodFunctionTag >( *storage, maxLevel );

   real_t discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], maxLevel ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_w = std::sqrt( err.uvw()[2].dotGlobal( err.uvw()[2], maxLevel ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel ) / (real_t) globalDoFs1 );
   real_t residuum_l2_1    = std::sqrt( r.dotGlobal( r, maxLevel ) / (real_t) globalDoFs1 );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error w = " << discr_l2_err_1_w );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

//   vtkOutput.write( maxLevel, 1 );

   auto tt        = storage->getTimingTree();
   auto ttreduced = tt->getReduced().getCopyWithRemainder();
   WALBERLA_LOG_INFO_ON_ROOT( ttreduced );

   WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v + discr_l2_err_1_w, 0.6417 );
   WALBERLA_CHECK_LESS( discr_l2_err_1_p, 4.003 );
   WALBERLA_CHECK_LESS( residuum_l2_1, 0.01934 );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "Stokes CC stencil-based" )
   stokesMinResConvergenceTest< hyteg::P2P1TaylorHoodStokesOperator >();
   WALBERLA_LOG_INFO_ON_ROOT( "Stokes CC elementwise" )
   stokesMinResConvergenceTest< hyteg::P2P1ElementwiseConstantCoefficientStokesOperator >();
   return EXIT_SUCCESS;
}
