/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "core/math/Constants.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::math::pi;

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   if( walberlaEnv.config() == nullptr )
   {
      auto defaultFile = "../data/param/fullMGTest.prm";
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

   hyteg::P1Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > b( "b", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > x( "x", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > x_exact( "x_exact", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > ax( "ax", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );

   hyteg::P1ConstantLaplaceOperator A( storage, minLevel, maxLevel );
   hyteg::P1ConstantMassOperator    M( storage, minLevel, maxLevel );

   hyteg::P1toP1LinearRestriction restrictionOperator;
   hyteg::P1toP1LinearProlongation prolongationOperator;
   hyteg::P1toP1QuadraticProlongation quadraticProlongationOperator;

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   r.enableTiming( timingTree );
   b.enableTiming( timingTree );
   x.enableTiming( timingTree );
   x_exact.enableTiming( timingTree );
   ax.enableTiming( timingTree );
   tmp.enableTiming( timingTree );
   err.enableTiming( timingTree );

   A.enableTiming( timingTree );
   M.enableTiming( timingTree );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return sin( pi * xx[0] ) * sin( pi * xx[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& xx ) {
      return 2 * pi * pi * sin( pi * xx[0] ) * sin( pi * xx[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   for( size_t ll = minLevel; ll <= maxLevel; ++ll )
   {
      x.interpolate( zero, ll, hyteg::Inner );
      x.interpolate( exact, ll, hyteg::DirichletBoundary );
      x_exact.interpolate( exact, ll );
      tmp.interpolate( rhs, ll, hyteg::Inner );
      M.apply( tmp, b, ll, hyteg::Inner );
   }

   auto solver =
       hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator >( storage, minLevel, minLevel, coarse_maxiter, real_c(0), coarse_tolerance );

   std::function< void( size_t ) > cscycle;

   cscycle = [&]( size_t level ) {
      if( level == minLevel )
      {
         solver.solve( A, x, b, minLevel );
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
         restrictionOperator.restrict( r, level, hyteg::Inner );
         b.assign( {1.0}, {r}, level - 1, hyteg::Inner );

         x.interpolate( zero, level - 1 );

         cscycle( level - 1 );

         // prolongate
         tmp.assign( {1.0}, {x}, level, hyteg::Inner );
         prolongationOperator.prolongate( x, level - 1, hyteg::Inner );
         x.add( {1.0}, {tmp}, level, hyteg::Inner );

         // post-smooth
         for( size_t i = 0; i < nu_post; ++i )
         {
            A.smooth_gs( x, b, level, hyteg::Inner );
         }
      }
   };
   //hyteg::VTKWriter< hyteg::P1Function< real_t > >({ &x, &x_exact, &b }, minLevel, "../output", "fullmg");

   LIKWID_MARKER_START( "Compute" );
   for( size_t ll = minLevel; ll <= maxLevel; ++ll )
   {
      tmp.interpolate( ones, ll );
      real_t npoints = tmp.dotGlobal( tmp, ll );
      WALBERLA_LOG_INFO_ON_ROOT( "Level = " << ll );
      WALBERLA_LOG_INFO_ON_ROOT( "Num dofs = " << npoints );
      WALBERLA_LOG_INFO_ON_ROOT( "Starting V cycles" );
      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%6s|%10s|%10s|%10s|%10s|%10s", "iter", "abs_res", "rel_res", "conv", "L2-error", "H1-semi" ) );
      real_t rel_res = 1.0;

      A.apply( x, ax, ll, hyteg::Inner );
      r.assign( {1.0, -1.0}, {b, ax}, ll, hyteg::Inner );
      real_t abs_res_old = std::sqrt( r.dotGlobal( r, ll, hyteg::Inner ) );
      real_t begin_res   = abs_res_old;
      err.assign( {1.0, -1.0}, {x, x_exact}, ll );
      real_t discr_l2_err = std::sqrt( err.dotGlobal( err, ll ) / npoints );
      A.apply( err, tmp, ll, hyteg::Inner );
      real_t discr_h1_err = std::sqrt( err.dotGlobal( tmp, ll ) );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
          "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res / abs_res_old, discr_l2_err, discr_h1_err ) )

      for( size_t i = 0; i < outer; ++i )
      {
         cscycle( ll );
         A.apply( x, ax, ll, hyteg::Inner );
         r.assign( {1.0, -1.0}, {b, ax}, ll, hyteg::Inner );
         real_t abs_res = std::sqrt( r.dotGlobal( r, ll, hyteg::Inner ) );
         rel_res        = abs_res / begin_res;
         err.assign( {1.0, -1.0}, {x, x_exact}, ll );
         discr_l2_err = std::sqrt( err.dotGlobal( err, ll ) / npoints );
         A.apply( err, tmp, ll, hyteg::Inner );
         discr_h1_err = std::sqrt( err.dotGlobal( tmp, ll ) );

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e",
                                                 i + 1,
                                                 begin_res,
                                                 rel_res,
                                                 begin_res / abs_res_old,
                                                 discr_l2_err,
                                                 discr_h1_err ) )
         abs_res_old = abs_res;
      }
      if( ll < maxLevel )
      {
         quadraticProlongationOperator.prolongate( x, ll, hyteg::Inner );
      }
   }
   LIKWID_MARKER_STOP( "Compute" );

   walberla::WcTimingTree tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt );

   LIKWID_MARKER_CLOSE;
   return EXIT_SUCCESS;
}
