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
#include "core/timing/Timer.h"
#include "core/math/Constants.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1ElementwiseOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = "../../data/meshes/quad_2el.msh";

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 4;

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
   loadbalancing::distributed::parmetis( *storage );
#endif

   hyteg::P1Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   auto coordX = std::make_shared< hyteg::P1Function< real_t > >( "x", storage, minLevel, maxLevel );
   auto coordY = std::make_shared< hyteg::P1Function< real_t > >( "y", storage, minLevel, maxLevel );

   auto coords    = std::make_shared< std::array< const hyteg::P1Function< real_t >*, 2 > >();
   ( *coords )[0] = coordX.get();
   ( *coords )[1] = coordY.get();

   //  hyteg::P1MassOperator M(storage, minLevel, maxLevel);
   hyteg::P1ConstantLaplaceOperator    Ltest( storage, minLevel, maxLevel );
   hyteg::P1ElementwiseLaplaceOperator L( storage, *coords, minLevel, maxLevel );
   hyteg::P1ElementwiseMassOperator    M( storage, *coords, minLevel, maxLevel );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   r.enableTiming( timingTree );
   f.enableTiming( timingTree );
   u.enableTiming( timingTree );
   u_exact.enableTiming( timingTree );
   err.enableTiming( timingTree );
   npoints_helper.enableTiming( timingTree );
   coordX->enableTiming( timingTree );
   coordY->enableTiming( timingTree );

   L.enableTiming( timingTree );

   std::function< real_t( const hyteg::Point3D& ) > map_x = []( const hyteg::Point3D& x ) {
      return x[0] + 0.1 * x[0] * std::sin( 3.0 * walberla::math::pi * x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > map_y = []( const hyteg::Point3D& x ) { return x[1]; };

   std::function< real_t( const hyteg::Point3D& ) > exact = [&map_x, &map_y]( const hyteg::Point3D& x ) {
      Point3D xt{{{map_x( x ), map_y( x ), 0.0}}};
      return ( 1.0 / 2.0 ) * sin( 2 * xt[0] ) * sinh( xt[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs = [&map_x, &map_y]( const hyteg::Point3D& x ) {
      Point3D xt{{{map_x( x ), map_y( x ), 0.0}}};
      return ( 3.0 / 2.0 ) * sin( 2 * xt[0] ) * sinh( xt[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   for( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      coordX->interpolate( map_x, level );
      coordY->interpolate( map_y, level );
   }

   // make sure that all coordinates are synchronized over all primitives and levels
   for( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      coordX->communicate< Vertex, Edge >( level );
      coordX->communicate< Edge, Face >( level );
      coordX->communicate< Face, Edge >( level );
      coordX->communicate< Edge, Vertex >( level );

      coordY->communicate< Vertex, Edge >( level );
      coordY->communicate< Edge, Face >( level );
      coordY->communicate< Face, Edge >( level );
      coordY->communicate< Edge, Vertex >( level );
   }

   u.interpolate( exact, maxLevel, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );
   npoints_helper.interpolate( rhs, maxLevel );
   M.apply( npoints_helper, f, maxLevel, hyteg::All );

   auto solver = hyteg::CGSolver< hyteg::P1ElementwiseLaplaceOperator >( storage, minLevel, maxLevel );
   walberla::WcTimer timer;
   solver.solve( L, u, f, maxLevel);
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );

   npoints_helper.interpolate( ones, maxLevel );
   real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   //  hyteg::VTKWriter<hyteg::P1Function<real_t>, hyteg::DGFunction<real_t >>({ u, u_exact, &f, &r, &err }, {}, maxLevel,
   //                                                                    "../output", "varcoords", coords);

   walberla::WcTimingTree tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt );

   /// a higher level would significantly decrase the error but also the runtime
   WALBERLA_CHECK_LESS( discr_l2_err, 3e-04 );

   return EXIT_SUCCESS;
}
