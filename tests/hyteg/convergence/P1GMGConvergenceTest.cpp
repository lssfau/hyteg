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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t      minLevel                  = 0;
   const uint_t      maxLevel                  = 5;
   const std::string meshFile                  = "../../data/meshes/quad_8el.msh";
   const real_t      coarseGridSolverTolerance = real_c( 1e-16 );
   const uint_t      maxCoarseGridSolverIter   = 10000;
   const uint_t      numVCycles                = 10;
   const bool        writeVTK                  = false;

   auto meshInfo = MeshInfo::fromGmshFile( meshFile );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   hyteg::P1Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > Au( "Au", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   hyteg::P1ConstantMassOperator    M( storage, minLevel, maxLevel );
   hyteg::P1ConstantLaplaceOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, maxLevel, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );
   npoints_helper.interpolate( rhs, maxLevel );
   M.apply( npoints_helper, f, maxLevel, hyteg::All );

   auto smoother         = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1ConstantLaplaceOperator > >();
   auto coarseGridSolver = std::make_shared< hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator > >(
       storage, minLevel, minLevel, maxCoarseGridSolverIter, coarseGridSolverTolerance );
   auto restrictionOperator  = std::make_shared< hyteg::P1toP1LinearRestriction<> >();
   auto prolongationOperator = std::make_shared< hyteg::P1toP1LinearProlongation<> >();

   auto gmgSolver = hyteg::GeometricMultigridSolver< hyteg::P1ConstantLaplaceOperator >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

   npoints_helper.interpolate( ones, maxLevel );
   const real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

   // init residual once
   L.apply( u, Au, maxLevel, hyteg::Inner );
   r.assign( { 1.0, -1.0 }, { f, Au }, maxLevel, hyteg::Inner );

   real_t discr_l2_err;
   real_t discr_l2_res           = std::sqrt( r.dotGlobal( r, maxLevel, DoFType::Inner ) / npoints );
   real_t discr_l2_res_last_step = discr_l2_res;

   for ( uint_t vCycleCount = 0; vCycleCount < numVCycles; vCycleCount++ )
   {
      gmgSolver.solve( L, u, f, maxLevel );

      err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );

      L.apply( u, Au, maxLevel, hyteg::Inner );
      r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hyteg::Inner );

      discr_l2_err           = std::sqrt( err.dotGlobal( err, maxLevel, DoFType::All ) / npoints );
      discr_l2_res_last_step = discr_l2_res;
      discr_l2_res           = std::sqrt( r.dotGlobal( r, maxLevel, DoFType::Inner ) / npoints );

      const real_t convRate = discr_l2_res / discr_l2_res_last_step;

      WALBERLA_LOG_INFO_ON_ROOT( "After " << vCycleCount
                                          << " V-cycles: "
                                             "discrete L2 error = "
                                          << discr_l2_err << ", discrete L2 residual = " << discr_l2_res
                                          << ", conv rate = " << convRate )
   }

   if ( writeVTK )
   {
      VTKOutput vtkOutput( "../../output", "P1GMGConvergenceTest", storage );
      vtkOutput.add( u );
      vtkOutput.add( u_exact );
      vtkOutput.add( err );
      vtkOutput.write( maxLevel );
   }

   bool dp = std::is_same< real_t, double >();
   WALBERLA_CHECK_LESS( discr_l2_res, dp ? 3.0e-14 : 5e-8 );
   WALBERLA_CHECK_LESS( discr_l2_err, 2.9e-06 );

   return 0;
}
