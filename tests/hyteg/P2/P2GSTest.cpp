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
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   /// create enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const size_t level = 3;

   /// create timingTree
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   /// read mesh file and create storage
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   /// create operator and functions
   hyteg::P2ConstantLaplaceOperator L( storage, level, level );

   hyteg::P2Function< real_t > residuumFunction( "residuumFunction", storage, level, level );
   hyteg::P2Function< real_t > rightHandSide( "rightHandSide", storage, level, level );
   hyteg::P2Function< real_t > u( "u", storage, level, level );
   hyteg::P2Function< real_t > u_exact( "u_exact", storage, level, level );
   hyteg::P2Function< real_t > error( "error", storage, level, level );
   hyteg::P2Function< real_t > npoints_helper( "npoints_helper", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   real_t residuum = 100;

   uint_t iterations;
   /// use different checks for debug and release since the runtime in debug would be to high otherwise
   WALBERLA_DEBUG_SECTION() { iterations = 100; }
   else { iterations = 5000; }

   /// apply Gauss Seidl smoother i times
   for( uint_t i = 0; i < iterations; ++i )
   {
      L.smooth_gs( u, rightHandSide, level, hyteg::Inner );
   }

   /// calculate residuum and check
   L.apply( u, residuumFunction, level, hyteg::Inner );
   residuumFunction.add( {-1}, {rightHandSide}, level, hyteg::Inner );
   residuum = std::sqrt( residuumFunction.dotGlobal( residuumFunction, level, hyteg::Inner ) );
   WALBERLA_LOG_INFO_ON_ROOT( "residual: = " << residuum );
   WALBERLA_DEBUG_SECTION() { WALBERLA_CHECK_GREATER( 0.2, residuum ); }
   else { WALBERLA_CHECK_GREATER( 1e-14, residuum ); }

   /// calculate and print error
   error.assign( {1.0, -1.0}, {u, u_exact}, level );
   npoints_helper.interpolate( ones, level );
   real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   real_t discr_l2_err = std::sqrt( error.dotGlobal( error, level ) / npoints );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   walberla::WcTimingTree tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt );

   return 0;
}
