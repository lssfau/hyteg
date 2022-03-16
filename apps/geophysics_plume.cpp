/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include <core/Environment.h>
#include <core/timing/Timer.h>

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGUpwindOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1HelperFunctions.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::DETAIL );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   timingTree->start( "Global" );

   std::string meshFileName = "../data/meshes/annulus_fine.msh";

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   const uint_t minLevel      = 2;
   const uint_t maxLevel      = 6;
   const uint_t solverMaxiter = 100;

   std::function< real_t( const hyteg::Point3D& ) > initialConcentration = []( const hyteg::Point3D& x ) {
      if ( sqrt( x[0] * x[0] + x[1] * x[1] ) < 1.1 )
      {
         return 1.0;
      }
      else
      {
         return 0.0;
      }
   };

   std::function< real_t( const hyteg::Point3D&, const std::vector< real_t >& ) > expr_f =
       []( const hyteg::Point3D&, const std::vector< real_t >& val ) { return 25.0 * val[0]; };

   std::function< real_t( const hyteg::Point3D& ) > expr_n_x = []( const hyteg::Point3D& x ) {
      return std::cos( std::atan2( x[1], x[0] ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > expr_n_y = []( const hyteg::Point3D& x ) {
      return std::sin( std::atan2( x[1], x[0] ) );
   };

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
   loadbalancing::distributed::parmetis( *storage );
#endif

   // Setting up Functions
   auto c_old = std::make_shared< hyteg::DGFunction< real_t > >( "c", storage, minLevel, maxLevel );
   auto c     = std::make_shared< hyteg::DGFunction< real_t > >( "c", storage, minLevel, maxLevel );

   auto f_dg = std::make_shared< hyteg::DGFunction< real_t > >( "f_dg", storage, minLevel, maxLevel );

   auto r = std::make_shared< hyteg::P1StokesFunction< real_t > >( "r", storage, minLevel, maxLevel );
   auto f = std::make_shared< hyteg::P1StokesFunction< real_t > >( "f", storage, minLevel, maxLevel );
   auto u = std::make_shared< hyteg::P1StokesFunction< real_t > >( "u", storage, minLevel, maxLevel );

   auto n_x = std::make_shared< hyteg::P1Function< real_t > >( "n_x", storage, maxLevel, maxLevel );
   auto n_y = std::make_shared< hyteg::P1Function< real_t > >( "n_y", storage, maxLevel, maxLevel );

   auto tmp = std::make_shared< hyteg::P1Function< real_t > >( "tmp", storage, minLevel, maxLevel );

   // Setting up Operators
   std::array< hyteg::P1Function< real_t >, 2 >           velocity{ u->uvw()[0], u->uvw()[1] };
   hyteg::DGUpwindOperator< hyteg::P1Function< real_t > > N( storage, velocity, minLevel, maxLevel );
   hyteg::P1P1StokesOperator                                L( storage, minLevel, maxLevel );
   hyteg::P1ConstantMassOperator                          M( storage, minLevel, maxLevel );

   real_t       estimatedMaxVelocity = P1::getApproximateEuclideanNorm< 2 >( { { &u->uvw()[0], &u->uvw()[1] } }, maxLevel );
   const real_t minimalEdgeLength    = hyteg::MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "minimalEdgeLength: " << minimalEdgeLength )
   real_t dt = std::min( 1.0, 0.25 * minimalEdgeLength / estimatedMaxVelocity );
   WALBERLA_LOG_INFO_ON_ROOT( "dt: " << dt )
   const real_t finalTime = 100000.0;
   //  const real_t plotEach = 2.0;
   const auto timesteps = (uint_t) std::ceil( finalTime / dt );
   //  const uint_t plotModulo = (uint_t) std::ceil(plotEach/dt);
   const uint_t plotModulo = 10;
   real_t       time       = 0.0;

   // Interpolate normal components
   n_x->interpolate( expr_n_x, maxLevel );
   n_y->interpolate( expr_n_y, maxLevel );

   // Interpolate initial functions
   c_old->interpolate( initialConcentration, maxLevel );
   c->assign( { 1.0 }, { *c_old }, maxLevel );

   auto pressurePreconditioner =
       std::make_shared< hyteg::StokesPressureBlockPreconditioner< hyteg::P1P1StokesOperator, hyteg::P1LumpedInvMassOperator > >(
           storage, minLevel, maxLevel );
   auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1P1StokesOperator::VelocityOperator_T > >();
   auto uzawaVelocitySmoother =
       std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< hyteg::P1P1StokesOperator > >( storage,
                                                                                                             gaussSeidel );
   auto smoother = std::make_shared< hyteg::UzawaSmoother< hyteg::P1P1StokesOperator > >(
       storage, uzawaVelocitySmoother, minLevel, maxLevel, 0.37 );
   auto coarseGridSolver = std::make_shared< hyteg::MinResSolver< hyteg::P1P1StokesOperator > >(
       storage, minLevel, minLevel, solverMaxiter, 1e-16, pressurePreconditioner );
   auto restrictionOperator  = std::make_shared< hyteg::P1P1StokesToP1P1StokesRestriction >();
   auto prolongationOperator = std::make_shared< hyteg::P1P1StokesToP1P1StokesProlongation >();

   auto solver = hyteg::GeometricMultigridSolver< hyteg::P1P1StokesOperator >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

   WALBERLA_LOG_DETAIL_ON_ROOT( "Total number of faces: " << storage->getNumberOfLocalFaces() );

   uint_t totalGlobalDofsStokes = 0;
   uint_t totalGlobalDofsTemp   = 0;
   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      uint_t tmpDofStokes = numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, lvl );
      uint_t tmpDofTemp   = numberOfGlobalDoFs< hyteg::DGFunctionTag >( *storage, lvl );
      WALBERLA_LOG_DETAIL_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
      WALBERLA_LOG_DETAIL_ON_ROOT( "Temperature DoFs on level " << lvl << " : " << tmpDofTemp );
      totalGlobalDofsStokes += tmpDofStokes;
      totalGlobalDofsTemp += tmpDofTemp;
   }
   WALBERLA_LOG_INFO_ON_ROOT( "Total Stokes DoFs on all level :" << totalGlobalDofsStokes );
   WALBERLA_LOG_INFO_ON_ROOT( "Total Temperature DoFs on all level :" << totalGlobalDofsTemp );
   WALBERLA_LOG_INFO_ON_ROOT( "Total DoFs on all level :" << ( totalGlobalDofsTemp + totalGlobalDofsStokes ) );

   hyteg::VTKOutput vtkOutput( "../output", "plume", storage, plotModulo );
   vtkOutput.add( *u );
   vtkOutput.add( f->uvw() );
   vtkOutput.add( *c );

   uint_t plotIter = 0;
   for ( uint_t t = 0; t <= timesteps; ++t )
   {
      WALBERLA_LOG_PROGRESS_ON_ROOT( "Current timestep: " << time )

      if ( t % 3 == 0 )
      {
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Solving Stokes system..." )

         f_dg->interpolate( expr_f, { *c_old }, maxLevel );

         f->uvw()[0].integrateDG( *f_dg, *n_x, maxLevel, hyteg::All );
         f->uvw()[1].integrateDG( *f_dg, *n_y, maxLevel, hyteg::All );

         for ( uint_t outer = 0; outer < 2; ++outer )
         {
            solver.solve( L, *u, *f, maxLevel );
            hyteg::vertexdof::projectMean( u->p(), maxLevel );

            L.apply( *u, *r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
            hyteg::vertexdof::projectMean( u->p(), maxLevel );

            r->assign( { 1.0, -1.0 }, { *f, *r }, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
            real_t residuum = std::sqrt( r->dotGlobal( *r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary ) );
            WALBERLA_LOG_PROGRESS_ON_ROOT( "[Uzawa] residuum: " << std::scientific << residuum )
         }
      }

      if ( t % plotModulo == 0 )
      {
         timingTree->start( "VTK" );
         vtkOutput.write( maxLevel, plotIter );
         ++plotIter;
         timingTree->stop( "VTK" );
      }

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Advecting temperature..." )
      N.apply( *c_old, *c, maxLevel, hyteg::Inner, Replace );
      c->assign( { 1.0, -dt }, { *c_old, *c }, maxLevel, hyteg::Inner );

      c_old.swap( c );
      time += dt;

      // compute new dt by CFL condition
      estimatedMaxVelocity = P1::getApproximateEuclideanNorm< 2 >( { { &u->uvw()[0], &u->uvw()[1] } }, maxLevel );
      dt                   = std::min( 1.0, 0.25 * minimalEdgeLength / estimatedMaxVelocity );
   }

   timingTree->stop( "Global" );
   auto reduced_tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( reduced_tt )

   return EXIT_SUCCESS;
}
