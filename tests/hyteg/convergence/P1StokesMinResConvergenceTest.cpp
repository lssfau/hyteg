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
#include <cmath>

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

using walberla::real_c;
using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = "../../data/meshes/quad_4el_neumann.msh";

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 2;
   //size_t maxiter  = 10000;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );

   hyteg::P1P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > bc_x = []( const hyteg::Point3D& x ) {
      if( x[0] < 1e-8 )
      {
         return 16.0 * ( x[1] - 0.5 ) * ( 1.0 - x[1] );
      } else
      {
         return 0.0;
      }
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.uvw().interpolate( {bc_x, zero}, maxLevel, hyteg::DirichletBoundary );

   auto solver = hyteg::MinResSolver< hyteg::P1P1StokesOperator >( storage, minLevel, maxLevel );
   solver.solve( L, u, f, maxLevel );

   L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
   real_t final_residuum = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) ) /
                           real_c( hyteg::numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, maxLevel ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residuum: " << final_residuum )

   WALBERLA_CHECK_LESS( final_residuum, 9.1e-07 );
   return EXIT_SUCCESS;
}
