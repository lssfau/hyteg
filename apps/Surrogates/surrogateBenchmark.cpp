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
#include "operators/P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab.hpp"
#include "operators/P1ElementwiseDiffusion_cubes_const_vect_vect512_fused_quadloops_tab.hpp"

using walberla::real_t;
using namespace hyteg;

// we want roughly one element per process
template < uint_t DIM >
std::array< uint_t, DIM > compute_domain_size( uint_t n_procs )
{
   uint_t N = n_procs / ( ( DIM == 2 ) ? 2 : 6 ); // number of cubes/quads

   std::vector< uint_t > primeFactors;
   for ( uint_t i = 2; i * i <= N; ++i )
   {
      while ( N % i == 0 )
      {
         primeFactors.push_back( i );
         N /= i;
      }
   }
   if ( N > 1 )
   {
      primeFactors.push_back( N );
   }

   std::array< uint_t, DIM > n;
   n.fill( 1 );

   for ( auto factor : primeFactors )
   {
      auto min_index = std::distance( n.begin(), std::min_element( n.begin(), n.end() ) );
      n[min_index] *= factor;
   }

   return n;
}

template < typename Operator >
double apply( const Operator&                    A,
              const hyteg::P1Function< real_t >& u,
              hyteg::P1Function< real_t >&       Au,
              const uint_t                       level,
              const uint_t                       iter )
{
   walberla::WcTimer timer;
   timer.start();
   for ( uint_t i = 0; i < iter; ++i )
   {
      A.apply( u, Au, level, All, Replace );
   }
   timer.end();
   return timer.total() / iter;
}

template < uint8_t DEGREE >
double apply_surrogate( const std::shared_ptr< PrimitiveStorage >& storage,
                        const hyteg::P1Function< real_t >&         u,
                        hyteg::P1Function< real_t >&               Au,
                        const uint_t                               level,
                        const forms::p1_div_k_grad_affine_q3&      form,
                        const uint_t                               iter )
{
   P1ElementwiseSurrogateAffineDivKGradOperator< DEGREE > A_q( storage, level, level, form );
   A_q.init( 0, "", false );
   return apply( A_q, u, Au, level, iter );
}

template < uint8_t MAX_SURROGATE_DEGREE >
void benchmark( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t level, const uint_t iter )
{
   uint8_t dim  = storage->hasGlobalCells() ? 3 : 2;
   uint_t  n_el = ( dim == 2 ) ? storage->getNumberOfGlobalFaces() : storage->getNumberOfGlobalCells();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "number of elements: %d", n_el ) );

   // setup pde coefficient k ∈ P_q
   hyteg::surrogate::polynomial::Polynomial< real_t, 3, MAX_SURROGATE_DEGREE + 1 > k_poly;
   for ( auto& c : k_poly )
   {
      c = walberla::math::realRandom();
   }
   auto k = [&]( const hyteg::Point3D& x ) { return k_poly.eval_naive( x ); };

   forms::p1_div_k_grad_affine_q3 form( k, k );

   // functions
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > Au( "Au", storage, level, level );

   auto   n_dof = u.getNumberOfGlobalDoFs( level );
   real_t t;
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "number of global DoF: %1.1e", n_dof ) );

   std::function< real_t( const hyteg::Point3D& ) > initialU = []( const hyteg::Point3D& x ) {
      return cos( 2 * M_PI * x[0] ) * cos( 2 * M_PI * x[1] ) * cos( 2 * M_PI * x[2] );
   };
   u.interpolate( initialU, level );

   WALBERLA_LOG_INFO_ON_ROOT( "apply() benchmark" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s | %10s | %10s |", "operator type", "t in sec", "MDoF/s" ) );
   walberla::WcTimer timer;

   // apply constant operator
   operatorgeneration::P1ElementwiseDiffusion A( storage, level, level );
   t = apply( A, u, Au, level, iter );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s | %10.1e | %10d |", "generated", t, int( n_dof / t * 1e-6 ) ) );

   // apply optimized operator
   operatorgeneration::P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab A_opt( storage, level, level );
   t = apply( A_opt, u, Au, level, iter );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s | %10.1e | %10d |", "gen. optimized", t, int( n_dof / t * 1e-6 ) ) );

   // // apply optimized operator with avx512
   // operatorgeneration::P1ElementwiseDiffusion_cubes_const_vect_vect512_fused_quadloops_tab A_512( storage, level, level );
   // t = apply(A_512, u, Au, level, iter);
   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s | %10.1e | %10d |", "gen. opt. AVX512", t, int( n_dof / t * 1e-6 ) ) );

   // apply surrogate operators
   for ( uint8_t q = 0; q <= MAX_SURROGATE_DEGREE; ++q )
   {
      switch ( q )
      {
      case 0:
         t = apply_surrogate< 0 >( storage, u, Au, level, form, iter );
         break;
      case 1:
         t = apply_surrogate< 1 >( storage, u, Au, level, form, iter );
         break;
      case 2:
         t = apply_surrogate< 2 >( storage, u, Au, level, form, iter );
         break;
      case 3:
         t = apply_surrogate< 3 >( storage, u, Au, level, form, iter );
         break;
      case 4:
         t = apply_surrogate< 4 >( storage, u, Au, level, form, iter );
         break;
      case 5:
         t = apply_surrogate< 5 >( storage, u, Au, level, form, iter );
         break;
      case 6:
         t = apply_surrogate< 6 >( storage, u, Au, level, form, iter );
         break;
      default:
         WALBERLA_LOG_WARNING_ON_ROOT( "Unsupported surrogate degree" );
         break;
      }
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%19s%d | %10.1e | %10d |", "surrogate q=", q, t, int( n_dof / t * 1e-6 ) ) );
   }
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto n_procs = walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   auto     n2       = compute_domain_size< 2 >( n_procs );
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, n2[0], n2[1] );
   SetupPrimitiveStorage setupStorage( meshInfo, n_procs );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // ----------------------------
   //  Prepare setup for 3D tests
   // ----------------------------
   auto     n3         = compute_domain_size< 3 >( n_procs );
   MeshInfo meshInfo3d = MeshInfo::meshCuboid( Point3D( 0.0, 0.0, 0.0 ), Point3D( 1.0, 1.0, 1.0 ), n3[0], n3[1], n3[2] );
   SetupPrimitiveStorage setupStorage3d( meshInfo3d, n_procs );
   loadbalancing::roundRobin( setupStorage3d );
   std::shared_ptr< PrimitiveStorage > storage3d = std::make_shared< PrimitiveStorage >( setupStorage3d );

   // -------------------
   //  Run benchmarks
   // -------------------
   uint_t            lvl2  = 10;
   uint_t            lvl3  = 7;
   constexpr uint8_t q_max = 3;
   uint_t            iter  = 10;
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "2d, level=%d", lvl2 ) );
   benchmark< q_max >( storage, lvl2, iter );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "3d, level=%d", lvl3 ) );
   benchmark< q_max >( storage3d, lvl3, iter );

   return 0;
}
