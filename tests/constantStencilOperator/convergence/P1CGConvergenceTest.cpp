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
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

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

   const uint_t      level    = 4;
   const std::string meshFile = "../../data/meshes/quad_8el.msh";

   auto meshInfo = MeshInfo::fromGmshFile( meshFile );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   hyteg::P1Function< real_t > r( "r", storage, level, level );
   hyteg::P1Function< real_t > f( "f", storage, level, level );
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > u_exact( "u_exact", storage, level, level );
   hyteg::P1Function< real_t > err( "err", storage, level, level );
   hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, level, level );

   hyteg::P1ConstantMassOperator    M( storage, level, level );
   hyteg::P1ConstantLaplaceOperator L( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );
   npoints_helper.interpolate( rhs, level );
   M.apply( npoints_helper, f, level, hyteg::All );

   auto solver = hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator >( storage, level, level );

   solver.solve( L, u, f, level );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );
   npoints_helper.interpolate( ones, level );

   const real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_CHECK_LESS( discr_l2_err, 1.2e-5 );

//   hyteg::VTKOutput vtkOutput( "../../output", "P2CGConvergenceTest", storage );
//   vtkOutput.add( u );
//   vtkOutput.write( level );

   walberla::WcTimingTree tt = storage->getTimingTree()->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt );

   return 0;
}
