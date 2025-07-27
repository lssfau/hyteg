/*
 * Copyright (c) 2024 Andreas Burkhart.
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
#include <string>

#include "core/mpi/MPIManager.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = hyteg::prependHyTeGMeshDir( "2D/quad_4el.msh" );

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 3;
   //size_t maxiter = 1000;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P2Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P2Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   hyteg::P2ConstantLaplaceOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D x ) -> real_t {
      return x[0] * x[0] - x[1] * x[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, maxLevel, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );

   auto solver = hyteg::GMRESSolver< hyteg::P2ConstantLaplaceOperator >(
       storage,
       minLevel,
       maxLevel,
       1000,
       real_c( 0 ),
       1e-16,
       std::make_shared< hyteg::IdentityPreconditioner< hyteg::P2ConstantLaplaceOperator > >(),
       1000,
       1e-16,
       0 );
   //solver.setPrintInfo( true );

   solver.solve( L, u, f, maxLevel );

   err.assign( { 1.0, -1.0 }, { u, u_exact }, maxLevel );

   npoints_helper.interpolate( ones, maxLevel );
   real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << std::scientific << discr_l2_err );
   bool dp = std::is_same< real_t, double >();
   WALBERLA_CHECK_LESS( discr_l2_err, dp ? 1e-13 : 1e-6 );
   return EXIT_SUCCESS;
}
