/*
 * Copyright (c) 2017-2012 Dominik Thoennes, Marcus Mohr.
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
#include "core/mpi/MPIManager.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   WALBERLA_LOG_INFO_ON_ROOT( "HyTeG CG Test\n" );

   std::string meshFileName = "../../data/meshes/quad_4el.msh";

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 2;
   //size_t maxiter  = 10000;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t >       npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   hyteg::P1ConstantVectorLaplaceOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return xx[0] * xx[0] - xx[1] * xx[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.uvw().interpolate( { exact, exact }, maxLevel, hyteg::DirichletBoundary );
   u_exact.uvw().interpolate( { exact, exact }, maxLevel );

   auto              solver = hyteg::CGSolver< hyteg::P1ConstantVectorLaplaceOperator >( storage, minLevel, maxLevel );
   walberla::WcTimer timer;
   solver.solve( L, u.uvw(), f.uvw(), maxLevel );
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   err.assign( { 1.0, -1.0 }, { u, u_exact }, maxLevel );

   npoints_helper.interpolate( ones, maxLevel );

   real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   WALBERLA_CHECK_LESS( discr_l2_err, 1.0e-16 )

   //  hyteg::VTKWriter({ u, u_exact, &f, &r, &err }, maxLevel, "../output", "test");
   return 0;
}
