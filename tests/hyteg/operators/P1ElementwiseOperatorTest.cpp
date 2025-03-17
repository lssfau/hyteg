/*
 * Copyright (c) 2025 Benjamin Mann.
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

#include <core/DataTypes.h>
#include <core/math/Random.h>
#include <core/mpi/MPIManager.h>
#include <hyteg/elementwiseoperators/P1ElementwiseOperator.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp>

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using namespace hyteg;

/* compare elementwise with constant stencil operator
 */
template < uint8_t DIM >
void P1ElementwiseOperatorTest( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t level )
{
   double epsilon, errorMax;

   // operators
   P1ElementwiseLaplaceOperator A( storage, level, level );
   P1ConstantLaplaceOperator    B( storage, level, level );

   // functions
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > Au( "Au", storage, level, level );
   hyteg::P1Function< real_t > Bu( "(Bu", storage, level, level );
   hyteg::P1Function< real_t > err( "(A-B)u", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > initialU = []( const hyteg::Point3D& x ) {
      return cos( 2 * M_PI * x[0] ) * cos( 2 * M_PI * x[1] ) * cos( 2 * M_PI * x[2] );
   };
   u.interpolate( initialU, level );

   epsilon = std::is_same< real_t, double >() ? 2e-12 : 1e-5;

   // apply operators
   A.apply( u, Au, level, All, Replace );
   B.apply( u, Bu, level, All, Replace );
   err.assign( { 1.0, -1.0 }, { Au, Bu }, level, All );
   errorMax = err.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%37s = %e", "||(A - B)u||_inf", errorMax ) )
   WALBERLA_CHECK_LESS( errorMax, epsilon, "||(A - B)u||_inf" );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // ----------------------------
   //  Prepare setup for 3D tests
   // ----------------------------
   MeshInfo              meshInfo3d = MeshInfo::meshCuboid( Point3D( 0.0, 0.0, 0.0 ), Point3D( 1.0, 1.0, 1.0 ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage3d( meshInfo3d, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage3d );
   std::shared_ptr< PrimitiveStorage > storage3d = std::make_shared< PrimitiveStorage >( setupStorage3d );

   // -------------------
   //  Run tests
   // -------------------
   for ( uint_t lvl = 3; lvl <= 5; ++lvl )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "level=%d, 2d", lvl ) );
      P1ElementwiseOperatorTest< 2 >( storage, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "level=%d, 3d", lvl ) );
      P1ElementwiseOperatorTest< 3 >( storage3d, lvl );
   }
   return 0;
}
