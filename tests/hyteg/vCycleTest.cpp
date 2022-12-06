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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/Format.hpp"
#include "core/logging/Logging.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/types/PointND.hpp"
#include "hyteg/types/types.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::real_c;

enum class CycleType
{
   VCYCLE,
   WCYCLE
};

template < class O, class F, class CSolver, class RestrictionOperator, class ProlongationOperator >
void cscycle( size_t    level,
              size_t    minLevel,
              CSolver&  csolver,
              RestrictionOperator restrictionOperator,
              ProlongationOperator prolongationOperator,
              O&        A,
              F&        x,
              F&        ax,
              F&        b,
              F&        r,
              F&        tmp,
              real_t    coarse_tolerance,
              size_t    coarse_maxiter,
              size_t    nu_pre,
              size_t    nu_post,
              CycleType cycleType = CycleType::VCYCLE )
{
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };

   if( level == minLevel )
   {
      csolver.solve( A, x, b, minLevel );
   } else
   {
      // pre-smooth
      for( size_t i = 0; i < nu_pre; ++i )
      {
         A.smooth_gs( x, b, level, hyteg::Inner );
      }

      A.apply( x, ax, level, hyteg::Inner );
      r.assign( {1.0, -1.0}, {b, ax}, level, hyteg::Inner );

      // restrict
      restrictionOperator.restrict( r, level, hyteg::Inner | hyteg::NeumannBoundary );
      b.assign( {1.0}, {r}, level - 1, hyteg::Inner );

      x.interpolate( zero, level - 1 );

      cscycle( level - 1, minLevel, csolver, restrictionOperator, prolongationOperator,
               A, x, ax, b, r, tmp, coarse_tolerance, coarse_maxiter, nu_pre, nu_post, cycleType );

      if( cycleType == CycleType::WCYCLE )
      {
         cscycle(
             level - 1, minLevel, csolver, restrictionOperator, prolongationOperator,
             A, x, ax, b, r, tmp, coarse_tolerance, coarse_maxiter, nu_pre, nu_post, cycleType );
      }

      // prolongate
      tmp.assign( {1.0}, {x}, level, hyteg::Inner );
      prolongationOperator.prolongate( x, level - 1, hyteg::Inner | hyteg::NeumannBoundary );
      x.add( {1.0}, {tmp}, level, hyteg::Inner );

      // post-smooth
      for( size_t i = 0; i < nu_post; ++i )
      {
         A.smooth_gs( x, b, level, hyteg::Inner );
      }
   }
}

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   if( walberlaEnv.config() == nullptr )
   {
      auto defaultFile = "../../data/param/vCycleTest.prm";
      cfg->readParameterFile( defaultFile );
      if( !*cfg )
      {
         WALBERLA_ABORT( "could not open default file: " << defaultFile );
      }
   } else
   {
      cfg = walberlaEnv.config();
   }

   auto parameters = cfg->getOneBlock( "Parameters" );

   WALBERLA_LOG_INFO_ON_ROOT( "HyTeG FMG Test" );

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( parameters.getParameter< std::string >( "mesh" ) );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   size_t minLevel = parameters.getParameter< size_t >( "minlevel" );
   size_t maxLevel = parameters.getParameter< size_t >( "maxlevel" );
   size_t nu_pre   = parameters.getParameter< size_t >( "nu_pre" );
   size_t nu_post  = parameters.getParameter< size_t >( "nu_post" );
   size_t outer    = parameters.getParameter< size_t >( "outer_iter" );

   size_t coarse_maxiter   = 100;
   real_t coarse_tolerance = 1e-6;
   real_t mg_tolerance     = 1e-8;

   hyteg::P1Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > b( "b", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > x( "x", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > x_exact( "x_exact", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > ax( "ax", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );

   hyteg::P1ConstantLaplaceOperator A( storage, minLevel, maxLevel );
   hyteg::P1toP1LinearRestriction restrictionOperator;
   hyteg::P1toP1LinearProlongation prolongationOperator;

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) { return xx[0] * xx[0] - xx[1] * xx[1]; };
   std::function< real_t( const hyteg::Point3D& ) > zeros = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };
   std::function< real_t( const hyteg::Point3D& ) > rand  = []( const hyteg::Point3D& ) {
      return static_cast< real_t >( std::rand() ) / static_cast< real_t >( RAND_MAX );
   };

   x.interpolate( zeros, maxLevel, hyteg::DirichletBoundary );
   x.interpolate( rand, maxLevel, hyteg::Inner );
   x_exact.interpolate( zeros, maxLevel );

   tmp.interpolate( ones, maxLevel );
   real_t npoints = tmp.dotGlobal( tmp, maxLevel );

   auto csolver = hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator >( storage, minLevel, minLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Num dofs = {}" << uint_c( npoints ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Starting V cycles" );
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "%6s|%10s|%10s|%10s|%10s|%10s", "iter", "abs_res", "rel_res", "conv", "L2-error", "Time" ) )

   real_t rel_res = 1.0;

   A.apply( x, ax, maxLevel, hyteg::Inner );
   r.assign( {1.0, -1.0}, {b, ax}, maxLevel, hyteg::Inner );

   real_t begin_res   = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) );
   real_t abs_res_old = begin_res;

   err.assign( {1.0, -1.0}, {x, x_exact}, maxLevel );
   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res / abs_res_old, discr_l2_err, 0 ) )

   real_t       totalTime              = real_c( 0.0 );
   real_t       averageConvergenceRate = real_c( 0.0 );
   const uint_t convergenceStartIter   = 3;

   LIKWID_MARKER_START( "Compute" );
   uint_t i = 0;
   for( ; i < outer; ++i )
   {
      auto start = walberla::timing::getWcTime();
      cscycle( maxLevel,
               minLevel,
               csolver,
               restrictionOperator,
               prolongationOperator,
               A,
               x,
               ax,
               b,
               r,
               tmp,
               coarse_tolerance,
               coarse_maxiter,
               nu_pre,
               nu_post,
               CycleType::VCYCLE );
      auto end = walberla::timing::getWcTime();
      A.apply( x, ax, maxLevel, hyteg::Inner );
      r.assign( {1.0, -1.0}, {b, ax}, maxLevel, hyteg::Inner );
      real_t abs_res = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) );
      rel_res        = abs_res / begin_res;
      err.assign( {1.0, -1.0}, {x, x_exact}, maxLevel );
      discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e",
                                              i + 1,
                                              begin_res,
                                              rel_res,
                                              begin_res / abs_res_old,
                                              discr_l2_err,
                                              end - start ) )
      totalTime += end - start;

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
   LIKWID_MARKER_STOP( "Compute" );

   WALBERLA_LOG_INFO_ON_ROOT( "Time to solution: " << std::scientific << totalTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Avg. convergence rate: " << std::scientific
                                                        << averageConvergenceRate / real_c( i + 1 - convergenceStartIter ) );

   WALBERLA_CHECK_LESS( i, outer );

   //  hyteg::VTKWriter< hyteg::P1Function< real_t >, hyteg::DGFunction<real_t> >({ x, &b, x_exact }, {}, maxLevel, "../../output", "test");
   LIKWID_MARKER_CLOSE;
   return EXIT_SUCCESS;
}
