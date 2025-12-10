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
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using walberla::math::e;
using walberla::math::pi;

namespace hyteg {

class Solution
{
 public:
   Solution( real_t diffusivity, uint_t setting )
   : diffusivity_( diffusivity )
   , setting_( setting )
   , t_( 0 )
   {}

   real_t operator()( const Point3D& p )
   {
      switch ( setting_ )
      {
      case 0:
         return diffusivity_ * ( 1 - std::pow( e, -t_ ) ) * std::sin( p[0] ) * std::cos( p[1] );
      case 1:
         return std::pow( e, -t_ ) * std::sin( pi * p[0] ) * std::sin( pi * p[1] );
      default:
         WALBERLA_ABORT( "Invalid setting" )
      }
   }

   void inc( const real_t& dt ) { t_ += dt; }

 private:
   real_t diffusivity_;
   uint_t setting_;
   real_t t_;
};

class Rhs
{
 public:
   Rhs( real_t diffusivity, uint_t setting )
   : diffusivity_( diffusivity )
   , setting_( setting )
   , t_( 0 )
   {}

   real_t operator()( const Point3D& p )
   {
      switch ( setting_ )
      {
      case 0:
         return diffusivity_ * ( 2 - std::pow( e, -t_ ) ) * std::sin( p[0] ) * std::cos( p[1] );
      case 1:
         return ( 2 * pi * pi - 1 ) * std::pow( e, -t_ ) * std::sin( pi * p[0] ) * std::sin( pi * p[1] );
      default:
         WALBERLA_ABORT( "Invalid setting" )
      }
   }

   void inc( const real_t& dt ) { t_ += dt; }

 private:
   real_t diffusivity_;
   uint_t setting_;
   real_t t_;
};

void P2UnsteadyDiffusionBlendingTest( const uint_t minLevel,
                                      const uint_t maxLevel,
                                      const uint_t testSolution,
                                      const uint_t steps,
                                      const real_t tMax,
                                      const uint_t timeSteppingScheme,
                                      const real_t discrL2Eps )
{
   const auto            meshInfo = MeshInfo::meshAnnulus( real_c( 0.35 ) * pi, pi, MeshInfo::CRISS, 8, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   AnnulusMap::setMap( setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2UnsteadyDiffusionBlendingTest_domain" );

   const real_t dt          = tMax / real_c( steps );
   const bool   vtk         = false;
   const real_t diffusivity = 1.0;

   hyteg::P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > uOld( "uOld", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > fOld( "fOld", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > uExact( "uExact", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > error( "error", storage, minLevel, maxLevel );

   DiffusionTimeIntegrator timeIntegrator =
       timeSteppingScheme == 0 ? DiffusionTimeIntegrator::ImplicitEuler : DiffusionTimeIntegrator::CrankNicolson;

   P2ElementwiseUnsteadyDiffusionOperator diffusionOperator( storage, minLevel, maxLevel, dt, diffusivity, timeIntegrator );
   P2ElementwiseBlendingLaplaceOperator   L( storage, minLevel, maxLevel );
   P2ElementwiseBlendingMassOperator      M( storage, minLevel, maxLevel );

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager manager;
   auto         solver = std::make_shared< PETScLUSolver< P2ElementwiseUnsteadyDiffusionOperator > >( storage, maxLevel );
#else
   auto coarseGridSolver =
       std::make_shared< CGSolver< P2ElementwiseUnsteadyDiffusionOperator > >( storage, minLevel, minLevel, 1000 );
   auto smoother =
       std::make_shared< WeightedJacobiSmoother< P2ElementwiseUnsteadyDiffusionOperator > >( storage, minLevel, maxLevel, 0.66 );
   auto restriction  = std::make_shared< P2toP2QuadraticRestriction >();
   auto prolongation = std::make_shared< P2toP2QuadraticProlongation >();
   auto gmgSolver    = std::make_shared< GeometricMultigridSolver< P2ElementwiseUnsteadyDiffusionOperator > >(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, 6, 6 );
   auto solver = std::make_shared< SolverLoop< P2ElementwiseUnsteadyDiffusionOperator > >( gmgSolver, 4 );
#endif
   UnsteadyDiffusion< P2Function< real_t >,
                      P2ElementwiseUnsteadyDiffusionOperator,
                      P2ElementwiseBlendingLaplaceOperator,
                      P2ElementwiseBlendingMassOperator >
       diffusionSolver( storage, minLevel, maxLevel, solver );

   hyteg::VTKOutput vtkOutput( "../../output", "P2UnsteadyDiffusionBlendingTest", storage );
   vtkOutput.add( uExact );
   vtkOutput.add( u );
   vtkOutput.add( error );

   Solution solution( diffusivity, testSolution );
   Rhs      rhs( diffusivity, testSolution );

   uExact.interpolate( solution, maxLevel, All );
   f.interpolate( rhs, maxLevel, All );
   u.interpolate( solution, maxLevel, DirichletBoundary );
   error.assign( { 1.0, -1.0 }, { u, uExact }, maxLevel, All );
   real_t l2Error    = std::sqrt( error.dotGlobal( error, maxLevel, Inner ) /
                               real_c( numberOfGlobalInnerDoFs< P2FunctionTag >( *storage, maxLevel ) ) );
   real_t l2Residual = 0.0;
   WALBERLA_LOG_INFO_ON_ROOT( " timestep | time total | discr. L2 error | discr. L2 residual" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------+------------+-----------------+--------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8d | %10.2e | %15.2e | %15.2e", 0, 0, l2Error, l2Residual ) );

   if ( vtk )
      vtkOutput.write( maxLevel, 0 );

   real_t timeTotal = 0;
   for ( uint_t step = 1; step <= steps; step++ )
   {
      solution.inc( dt );
      rhs.inc( dt );
      timeTotal += dt;

      uOld.assign( { 1.0 }, { u }, maxLevel, All );
      fOld.assign( { 1.0 }, { f }, maxLevel, All );

      uExact.interpolate( solution, maxLevel, All );
      f.interpolate( rhs, maxLevel, All );
      u.interpolate( solution, maxLevel, DirichletBoundary );

      diffusionSolver.step( diffusionOperator, L, M, u, uOld, f, fOld, maxLevel, Inner );

      error.assign( { 1.0, -1.0 }, { u, uExact }, maxLevel, All );
      l2Error = std::sqrt( error.dotGlobal( error, maxLevel, Inner ) /
                           real_c( numberOfGlobalInnerDoFs< P2FunctionTag >( *storage, maxLevel ) ) );

      diffusionSolver.calculateResidual( diffusionOperator, L, M, u, uOld, f, fOld, r, maxLevel, Inner );
      l2Residual = std::sqrt( r.dotGlobal( r, maxLevel, Inner ) /
                              real_c( numberOfGlobalInnerDoFs< P2FunctionTag >( *storage, maxLevel ) ) );

      if ( vtk )
         vtkOutput.write( maxLevel, step );
   }
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8d | %10.2e | %15.2e | %15.2e", steps, timeTotal, l2Error, l2Residual ) );

   WALBERLA_CHECK_LESS( l2Error, discrL2Eps );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   uint_t testSolution   = 0;
   real_t tMax           = 1;
   uint_t maxLevel       = 4;
   uint_t timeIntegrator = 0;
   hyteg::P2UnsteadyDiffusionBlendingTest( 2, maxLevel, testSolution, 10, tMax, timeIntegrator, real_c( 3.0e-03 ) );
   hyteg::P2UnsteadyDiffusionBlendingTest( 2, maxLevel, testSolution, 20, tMax, timeIntegrator, real_c( 1.5e-03 ) );
   hyteg::P2UnsteadyDiffusionBlendingTest( 2, maxLevel, testSolution, 40, tMax, timeIntegrator, real_c( 7.5e-04 ) );

   timeIntegrator = 1;
   hyteg::P2UnsteadyDiffusionBlendingTest( 2, maxLevel, testSolution, 10, tMax, timeIntegrator, real_c( 5.0e-05 ) );
   hyteg::P2UnsteadyDiffusionBlendingTest( 2, maxLevel, testSolution, 20, tMax, timeIntegrator, real_c( 1.5e-05 ) );
   hyteg::P2UnsteadyDiffusionBlendingTest( 2, maxLevel, testSolution, 40, tMax, timeIntegrator, real_c( 4.1e-06 ) );
}
