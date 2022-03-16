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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;


int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   const std::string meshFileName = "../../data/meshes/3D/cube_24el.msh";
   const uint_t minLevel         =  2;
   const uint_t maxLevel         =  3;
   const uint_t maxIterations    =  5;
   const real_t tolerance = 1e-16;

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

#if 1
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

   hyteg::writeDomainPartitioningVTK( storage, "../../output", "P1_Stokes_3D_MinRes_convergence_partitioning" );

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > uExact( "uExact", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );

//   hyteg::VTKOutput vtkOutput( "../../output", "P1_Stokes_3D_MinRes_convergence", storage );
//   vtkOutput.add( u.u );
//   vtkOutput.add( u.v );
//   vtkOutput.add( u.w );
//   vtkOutput.add( u.p );
//   vtkOutput.add( uExact.u );
//   vtkOutput.add( uExact.v );
//   vtkOutput.add( uExact.w );
//   vtkOutput.add( uExact.p );

   hyteg::P1P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > inflowPoiseuille = []( const hyteg::Point3D& x )
   {
      if( x[2] < 1e-8 )
      {
        return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
      }
      else
      {
         return 0.0;
      }
   };


  std::function< real_t( const hyteg::Point3D& ) > solutionPoiseuille = []( const hyteg::Point3D& x )
  {
      return ( 1.0 / 16.0 ) * x[0] * ( 1 - x[0] ) * x[1] * ( 1.0 - x[1] );
  };

  std::function< real_t( const hyteg::Point3D& ) > collidingFlow_x = []( const hyteg::Point3D& x )
  {
    return real_c(20) * x[0] * x[1] * x[1] * x[1];
  };

  std::function< real_t( const hyteg::Point3D& ) > collidingFlow_y = []( const hyteg::Point3D& x )
  {
      return real_c(5) * x[0] * x[0] * x[0] * x[0] - real_c(5) * x[1] * x[1] * x[1] * x[1];
  };

   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

#if 1
   u.uvw()[2].interpolate( inflowPoiseuille, maxLevel, hyteg::DirichletBoundary );
#else
   u.uvw().u.interpolate( { collidingFlow_x, collidingFlow_y}, maxLevel, hyteg::DirichletBoundary );
   uExact.uvw().interpolate( { collidingFlow_x, collidingFlow_y }, maxLevel );
#endif

//   vtkOutput.write( maxLevel, 0 );
#if 1

   typedef hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator > CoarseGridSolver_T;
   typedef hyteg::GeometricMultigridSolver< hyteg::P1ConstantLaplaceOperator  > GMGSolver_T;
   typedef hyteg::StokesBlockDiagonalPreconditioner< hyteg::P1P1StokesOperator, hyteg::P1LumpedInvMassOperator > Preconditioner_T;

   auto coarseGridSolver = std::make_shared< CoarseGridSolver_T  >( storage, minLevel, maxLevel );
   auto smoother = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1ConstantLaplaceOperator>  >();
   auto prolongationOperator = std::make_shared< hyteg::P1toP1LinearProlongation >();
   auto restrictionOperator = std::make_shared< hyteg::P1toP1LinearRestriction >();
   auto gmgSolver = std::make_shared< GMGSolver_T >( storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2 );
   //hyteg::P1LumpedInvMassOperator massOperator( storage, minLevel, maxLevel );
   //Preconditioner_T prec( storage, minLevel, maxLevel, 2, gmgSolver );
   auto prec = std::make_shared< Preconditioner_T >( storage, minLevel, maxLevel, 2, gmgSolver );

   auto solver = hyteg::MinResSolver< hyteg::P1P1StokesOperator >( storage, minLevel, maxLevel, maxIterations, tolerance, prec );
   // auto solver = hyteg::MinResSolver< hyteg::P1StokesFunction< real_t >, hyteg::P1P1StokesOperator, PressurePreconditioner_T >( storage, minLevel, maxLevel, pressurePrec );
   // auto solver = hyteg::MinResSolver< hyteg::P1StokesFunction< real_t >, hyteg::P1P1StokesOperator >( storage, minLevel, maxLevel );

   solver.solve( L, u, f, maxLevel );
#else
   auto         numerator  = std::make_shared< hyteg::P1StokesFunction< idx_t > >( "numerator", storage, level, level );
   uint_t globalSize = 0;
   const uint_t localSize = numerator->enumerate(level, globalSize);
   PETScManager petscManager( &argc, &argv );
   PETScLUSolver< real_t, hyteg::P1StokesFunction, hyteg::P1P1StokesOperator > petScLUSolver( numerator, localSize, globalSize );
   f.u.assign( {1.0}, {u.u}, level, DirichletBoundary );
   f.v.assign( {1.0}, {u.v}, level, DirichletBoundary );
   f.w.assign( {1.0}, {u.w}, level, DirichletBoundary );
   petScLUSolver.solve( L, u, f, r, level, tolerance, maxIterations, Inner | NeumannBoundary );
#endif
//   vtkOutput.write( maxLevel, 1 );

   L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
   real_t final_residual = r.dotGlobal( r, maxLevel, hyteg::Inner ) / real_c( hyteg::numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, maxLevel ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residual: " << final_residual )
   WALBERLA_CHECK_LESS( final_residual, 3.1e-12 );

   return EXIT_SUCCESS;
}
