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
#include "core/Format.hpp"
#include "core/config/Config.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"
#include "constantStencilOperator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "../../data/param/gmg_P2.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t minLevel         = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel         = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t max_outer_iter   = parameters.getParameter< uint_t >( "max_outer_iter" );
   const uint_t max_cg_iter      = parameters.getParameter< uint_t >( "max_cg_iter" );
   const real_t mg_tolerance     = parameters.getParameter< real_t >( "mg_tolerance" );
   const real_t coarse_tolerance = parameters.getParameter< real_t >( "coarse_tolerance" );
   const uint_t nuPre            = 3;
   const uint_t nuPost           = 3;

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( parameters.getParameter< std::string >( "mesh" ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   hyteg::P1Function< real_t >  r_p1( "r_p1", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t >  f_p1( "f_p1", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t >  u_p1( "u_p1", storage, minLevel, maxLevel );

   hyteg::P2Function< real_t > r_p2( "r_p2", storage, maxLevel, maxLevel );
   hyteg::P2Function< real_t > f_p2( "f_p2", storage, maxLevel, maxLevel );
   hyteg::P2Function< real_t > u_p2( "u_p2", storage, maxLevel, maxLevel );
   hyteg::P2Function< real_t > Lu_p2( "Lu_p2", storage, maxLevel, maxLevel );
   hyteg::P2Function< real_t > tmp_p2( "tmp_p2", storage, maxLevel, maxLevel );
   hyteg::P2Function< real_t > u_exact_p2( "u_exact_p2", storage, maxLevel, maxLevel );
   hyteg::P2Function< real_t > err_p2( "err_p2", storage, maxLevel, maxLevel );
   hyteg::P2Function< real_t > npoints_helper_p2( "npoints_helper_p2", storage, maxLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > zero  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating u" );
   u_p2.interpolate( exact, maxLevel, hyteg::DirichletBoundary );

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating exact function" );
   u_exact_p2.interpolate( exact, maxLevel );
   //  WALBERLA_LOG_INFO_ON_ROOT("Interpolating and integrating rhs");
   //  npoints_helper.interpolate(rhs, maxLevel);
   //  M.apply(npoints_helper, f, maxLevel, hyteg::All);

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up stiffness operator" );
   auto                             start = walberla::timing::getWcTime();
   hyteg::P1ConstantLaplaceOperator L_p1( storage, minLevel, maxLevel );
   hyteg::P2ConstantLaplaceOperator L_p2( storage, maxLevel, maxLevel );
   auto                             end       = walberla::timing::getWcTime();
   real_t                           setupTime = end - start;

   npoints_helper_p2.interpolate( ones, maxLevel );
   real_t npoints = npoints_helper_p2.dotGlobal( npoints_helper_p2, maxLevel );

   auto smoother         = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1ConstantLaplaceOperator > >();
   auto coarseGridSolver = std::make_shared< hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator > >(
       storage, minLevel, minLevel, max_cg_iter, coarse_tolerance );
   auto restrictionOperator  = std::make_shared< hyteg::P1toP1LinearRestriction<> >();
   auto prolongationOperator = std::make_shared< hyteg::P1toP1LinearProlongation<> >();

   auto gmgSolver = hyteg::GeometricMultigridSolver< hyteg::P1ConstantLaplaceOperator >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

   WALBERLA_LOG_INFO_ON_ROOT( "Starting V cycles" );
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "%6s|%10s|%10s|%10s|%10s|%10s", "iter", "abs_res", "rel_res", "conv", "L2-error", "Time" ) );

   real_t rel_res = 1.0;

   L_p2.apply( u_p2, Lu_p2, maxLevel, hyteg::Inner );
   r_p2.assign( {1.0, -1.0}, {f_p2, Lu_p2}, maxLevel, hyteg::Inner );

   real_t begin_res   = std::sqrt( r_p2.dotGlobal( r_p2, maxLevel, hyteg::Inner ) );
   real_t abs_res_old = begin_res;

   err_p2.assign( {1.0, -1.0}, {u_p2, u_exact_p2}, maxLevel );
   real_t discr_l2_err = std::sqrt( err_p2.dotGlobal( err_p2, maxLevel ) / npoints );

   //WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  -", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err));
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res / abs_res_old, discr_l2_err, 0 ) )

   real_t       solveTime              = real_c( 0.0 );
   real_t       averageConvergenceRate = real_c( 0.0 );
   const uint_t convergenceStartIter   = 3;

   uint_t i = 0;
   for( ; i < max_outer_iter; ++i )
   {
      start = walberla::timing::getWcTime();

      // pre-smooth
      for( size_t nu = 0; nu < nuPre; ++nu )
      {
         L_p2.smooth_gs( u_p2, f_p2, maxLevel, hyteg::Inner );
      }

      // compute residuum
      L_p2.apply( u_p2, Lu_p2, maxLevel, hyteg::Inner );
      r_p2.assign( {1.0, -1.0}, {f_p2, Lu_p2}, maxLevel, hyteg::Inner );

      // restrict
      r_p2.restrictP2ToP1( f_p1, maxLevel, hyteg::Inner );

      u_p1.interpolate( zero, maxLevel );

      // Apply P1 geometric multigrid solver
      gmgSolver.solve(L_p1, u_p1, f_p1, maxLevel );

      // prolongate
      tmp_p2.assign( {1.0}, {u_p2}, maxLevel, hyteg::Inner );
      u_p2.prolongateP1ToP2( u_p1, maxLevel, hyteg::Inner );
      u_p2.add( {1.0}, {tmp_p2}, maxLevel, hyteg::Inner );

      // post-smooth
      for( size_t nu = 0; nu < nuPost; ++nu )
      {
         L_p2.smooth_gs( u_p2, f_p2, maxLevel, hyteg::Inner );
      }

      end = walberla::timing::getWcTime();

      L_p2.apply( u_p2, Lu_p2, maxLevel, hyteg::Inner );
      r_p2.assign( {1.0, -1.0}, {f_p2, Lu_p2}, maxLevel, hyteg::Inner );
      real_t abs_res = std::sqrt( r_p2.dotGlobal( r_p2, maxLevel, hyteg::Inner ) );
      rel_res        = abs_res / begin_res;
      err_p2.assign( {1.0, -1.0}, {u_p2, u_exact_p2}, maxLevel );
      discr_l2_err = std::sqrt( err_p2.dotGlobal( err_p2, maxLevel ) / npoints );

      //WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err, end-start));
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", i + 1, abs_res, rel_res, abs_res / abs_res_old, discr_l2_err, end - start ) )
      solveTime += end - start;

      if( i >= convergenceStartIter )
      {
         averageConvergenceRate += abs_res / abs_res_old;
      }

      abs_res_old = abs_res;

      if( rel_res < mg_tolerance )
      {
         break;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Setup time: " << std::defaultfloat << setupTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Solve time " << std::defaultfloat << solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Time to solution: " << std::defaultfloat << setupTime + solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Avg. convergence rate: " << std::scientific
                                                        << averageConvergenceRate / real_c( i - convergenceStartIter ) );
   WALBERLA_LOG_INFO_ON_ROOT( "L^2 error: " << std::scientific << discr_l2_err );
   WALBERLA_LOG_INFO_ON_ROOT( "DoFs: " << (uint_t) npoints );

   WALBERLA_CHECK_LESS( discr_l2_err, 9e-08 );

   if( parameters.getParameter< bool >( "vtkOutput" ) )
   {
      VTKOutput vtkOutput("../output", "gmg_P2", storage);
      vtkOutput.add( u_p2 );
      vtkOutput.add( u_exact_p2 );
      vtkOutput.add( f_p2 );
      vtkOutput.add( r_p2 );
      vtkOutput.add( err_p2 );
      vtkOutput.add( npoints_helper_p2 );
      vtkOutput.write( maxLevel );
   }

   if( parameters.getParameter< bool >( "printTiming" ) )
   {
      walberla::WcTimingTree tt = timingTree->getReduced();
      WALBERLA_LOG_INFO_ON_ROOT( tt );
   }

   return 0;
}
