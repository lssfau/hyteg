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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.1
 */

#include <core/DataTypes.h>
#include <core/math/Random.h>
#include <core/mpi/MPIManager.h>
#include <hyteg/elementwiseoperators/P1ElementwiseSurrogateOperator.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp>

#include "operators/P1ElementwiseDiffusion.hpp"
#include "operators/P1ElementwiseDiffusion_cubes_const_vect_polycse.hpp"

using walberla::real_t;
using namespace hyteg;

void benchmark( const std::shared_ptr< PrimitiveStorage >& storage, const uint8_t q_max, const uint_t level )
{
   double epsilon, errorMax;

   // setup pde coefficient k ∈ P_q
   uint8_t                                            dim = storage->hasGlobalCells() ? 3 : 2;
   hyteg::surrogate::polynomial::Polynomial< real_t > k_poly( dim, q_max + 1 );
   for ( auto& c : k_poly )
   {
      c = walberla::math::realRandom();
   }
   auto k = [&]( const hyteg::Point3D& x ) { return k_poly.eval_naive( x ); };

   forms::p1_div_k_grad_affine_q3 form( k, k );

   // constant operators
   operatorgeneration::P1ElementwiseDiffusion                          A( storage, level, level );
   operatorgeneration::P1ElementwiseDiffusion_cubes_const_vect_polycse A_opt( storage, level, level );

   // functions
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > Au( "Au", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > initialU = []( const hyteg::Point3D& x ) {
      return cos( 2 * M_PI * x[0] ) * cos( 2 * M_PI * x[1] ) * cos( 2 * M_PI * x[2] );
   };
   u.interpolate( initialU, level );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s | t in sec for apply()", "operator type" ) );
   walberla::WcTimer timer;

   // apply constant operator
   timer.start();
   A.apply( u, Au, level, All, Replace );
   timer.end();
   auto t = timer.total();
   timer.reset();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s | %e", "generated", t ) );

   // apply optimized operator
   timer.start();
   A_opt.apply( u, Au, level, All, Replace );
   timer.end();
   t = timer.total();
   timer.reset();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s | %e", "generated-optimized", t ) );

   // apply surrogate operators
   for ( uint8_t q = 0; q <= q_max; ++q )
   {
      P1ElementwiseSurrogateAffineDivKGradOperator A_q( storage, level, level, form );
      A_q.init( q, 0, "", false );
      timer.start();
      A_q.apply( u, Au, level, All, Replace );
      timer.end();
      t = timer.total();
      timer.reset();
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%18s %d | %e", "surrogate q =", q, t ) );
   }
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
   //  Run benchmarks
   // -------------------
   uint_t  lvl2   = 10;
   uint_t  lvl3   = 7;
   uint8_t q_max = 5;
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "2d, level=%d", lvl2 ) );
   benchmark( storage, q_max, lvl2 );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "3d, level=%d", lvl3 ) );
   benchmark( storage3d, q_max, lvl3 );

   return 0;
}
