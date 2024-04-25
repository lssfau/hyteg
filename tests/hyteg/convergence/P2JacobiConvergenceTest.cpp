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
#include "core/Format.hpp"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

template< typename P2LaplaceOperator_T >
void jacobiTest()
{
   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "./P2JacobiConvergenceTest.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
   size_t                        level      = parameters.getParameter< size_t >( "level" );
   size_t                        maxiter    = parameters.getParameter< size_t >( "maxiter" );
   MeshInfo                      meshInfo   = MeshInfo::fromGmshFile( parameters.getParameter< std::string >( "mesh" ) );
   SetupPrimitiveStorage         setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   P2LaplaceOperator_T L( storage, level, level );

   hyteg::P2Function< real_t > residuum( "residuum", storage, level, level );
   hyteg::P2Function< real_t > rhs( "rhs", storage, level, level );
   hyteg::P2Function< real_t > p2function( "p2Function", storage, level, level );
   hyteg::P2Function< real_t > Lu( "Lu", storage, level, level );
   hyteg::P2Function< real_t > p2Exact( "p2Exact", storage, level, level );
   hyteg::P2Function< real_t > error( "error", storage, level, level );
   hyteg::P2Function< real_t > helperFun( "helperFun", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exactFunction = []( const hyteg::Point3D& x ) {
     return sin( x[0] ) * sinh( x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > ones   = []( const hyteg::Point3D& ) { return 1.0; };

   p2function.interpolate( exactFunction, level, hyteg::DirichletBoundary );
   helperFun.interpolate( exactFunction, level, hyteg::DirichletBoundary );
   p2Exact.interpolate( exactFunction, level );

   real_t begin_res, abs_res_old, rel_res = 0, abs_res = 0;

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%10s|%10s|%10s", "iter", "abs_res", "rel_res", "conv" ) );

   L.apply( p2function, Lu, level, hyteg::Inner );
   residuum.assign( {1.0, -1.0}, {rhs, Lu}, level, hyteg::Inner );
   begin_res   = std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner ) );
   abs_res_old = begin_res;

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res / abs_res_old ) )
   walberla::WcTimer timer;
   for ( uint_t i = 0; i < maxiter; ++i )
   {
      helperFun.assign( {1.0}, {p2function}, level, hyteg::Inner );
      L.smooth_jac( p2function, rhs, helperFun, 2.0 / 3.0, level, hyteg::Inner );
      L.apply( p2function, Lu, level, hyteg::Inner );
      residuum.assign( {1.0, -1.0}, {rhs, Lu}, level, hyteg::Inner );
      abs_res = std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner ) );
      rel_res = abs_res / begin_res;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%10.3e|%10.3e|%10.3e", i + 1, abs_res, rel_res, abs_res / abs_res_old ) )
      WALBERLA_CHECK_LESS( abs_res, abs_res_old );
      abs_res_old = abs_res;
   }
   timer.end();

   if ( parameters.getParameter< bool >( "vtkOutput" ) )
   {
      VTKOutput vtkOutput( "../../output", "gs_P2", storage );
      vtkOutput.add( p2function );
      vtkOutput.add( p2Exact );
      vtkOutput.add( rhs );
      vtkOutput.add( residuum );
      vtkOutput.add( error );
      vtkOutput.add( helperFun );
      vtkOutput.write( level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   error.assign( {1.0, -1.0}, {p2function, p2Exact}, level );

   helperFun.interpolate( ones, level );
   real_t npoints = helperFun.dotGlobal( helperFun, level );

   real_t discr_l2_err = std::sqrt( error.dotGlobal( error, level ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   if ( parameters.getParameter< bool >( "printTiming" ) )
   {
      walberla::WcTimingTree tt = timingTree->getReduced();
      WALBERLA_LOG_INFO_ON_ROOT( tt );
   }

   WALBERLA_CHECK_LESS( discr_l2_err, 2.0e-06 );
   WALBERLA_CHECK_LESS( abs_res,      8.0e-06 );

}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   jacobiTest< hyteg::P2ConstantLaplaceOperator >();
   jacobiTest< hyteg::P2ElementwiseLaplaceOperator >();

   return 0;
}
