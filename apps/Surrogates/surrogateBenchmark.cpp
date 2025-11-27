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
#include <hyteg/geometry/IcosahedralShellMap.hpp>
#include <hyteg/p1functionspace/P1SurrogateOperator_optimized.hpp>
#include <hyteg/primitivestorage/PrimitiveStorage.hpp>
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include <hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp>

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "operators/P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab.hpp"
#include "operators/P1ElementwiseDivKGrad_cubes_const_vect_fused_quadloops_tab.hpp"

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

template < class Operator, bool OperatorNeedsForm, bool OperatorNeedsCoeff, bool SurrogateOperator, bool smooth_gs, class Form >
void benchmark( const uint_t                                     dim,
                const uint_t                                     level,
                const uint_t                                     iter,
                const Form&                                      form,
                const std::function< real_t( const Point3D& ) >& coeff,
                const char*                                      opName,
                const bool                                       verbose,
                const bool                                       printTimingTree )
{
   auto n_procs = walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   // create domain and PrimitiveStorage
   std::shared_ptr< PrimitiveStorage > storage;
   if ( dim == 2 )
   {
      auto     n        = compute_domain_size< 2 >( n_procs );
      MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, n[0], n[1] );
      SetupPrimitiveStorage setupStorage( meshInfo, n_procs );
      loadbalancing::roundRobin( setupStorage );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }
   else // dim == 3
   {
      auto     n        = compute_domain_size< 3 >( n_procs );
      MeshInfo meshInfo = MeshInfo::meshCuboid( Point3D( 0.0, 0.0, 0.0 ), Point3D( 1.0, 1.0, 1.0 ), n[0], n[1], n[2] );
      SetupPrimitiveStorage setupStorage( meshInfo, n_procs );
      loadbalancing::roundRobin( setupStorage );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }

   // functions
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > Au( "Au", storage, level, level );

   // domain info
   uint_t n_el  = ( dim == 2 ) ? storage->getNumberOfGlobalFaces() : storage->getNumberOfGlobalCells();
   auto   n_dof = u.getNumberOfGlobalDoFs( level );
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%dd, level=%d", dim, level ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "number elements: %d, global DoF: %1.1f M", n_el, n_dof * 1e-6 ) );
   }

   std::function< real_t( const hyteg::Point3D& ) > initialU = []( const hyteg::Point3D& x ) {
      return cos( 2 * M_PI * x[0] ) * cos( 2 * M_PI * x[1] ) * cos( 2 * M_PI * x[2] );
   };
   u.interpolate( initialU, level );

   // operator
   std::shared_ptr< Operator > A;
   if constexpr ( SurrogateOperator )
   {
      A = std::make_shared< Operator >( storage, level, level, form, 0, "svd" );
   }
   else if constexpr ( OperatorNeedsForm )
   {
      A = std::make_shared< Operator >( storage, level, level, form );
   }
   else if constexpr ( OperatorNeedsCoeff )
   {
      hyteg::P1Function< real_t > k( "k", storage, level, level );
      k.interpolate( coeff, level );
      A = std::make_shared< Operator >( storage, level, level, k );
   }
   else
   {
      A = std::make_shared< Operator >( storage, level, level );
   }

   // benchmark
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%30s | %10s | %10s |", "operator type", "t in sec", "MDoF/s" ) );
   }
   walberla::WcTimer timer;
   timer.start();
   for ( uint_t i = 0; i < iter; ++i )
   {
      if constexpr ( smooth_gs )
      {
         A->smooth_gs( u, Au, level, All );
      }
      else
      {
         A->apply( u, Au, level, All, Add );
      }
   }
   timer.end();
   auto t = timer.total() / iter;
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%30s | %10.1e | %10d |", opName, t, int( n_dof / t * 1e-6 ) ) );

   if ( printTimingTree )
   {
      auto& timingTree = *( storage->getTimingTree() );
      if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
      {
         std::cout << "\n" << timingTree << "\n";
      }
   }
}

template < uint_t q >
void precompute_SVDs( uint_t lvl, uint_t dim )
{
   std::shared_ptr< PrimitiveStorage > storage;
   if ( dim == 2 )
   {
      MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, 1, 1 );
      SetupPrimitiveStorage setupStorage( meshInfo, 1 );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }
   else // dim == 3
   {
      MeshInfo              meshInfo = MeshInfo::meshCuboid( Point3D( 0.0, 0.0, 0.0 ), Point3D( 1.0, 1.0, 1.0 ), 1, 1, 1 );
      SetupPrimitiveStorage setupStorage( meshInfo, 1 );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }

   P1SurrogateBlendingLaplaceOperator< q > op( storage, lvl, lvl, 0 );
   op.store_svd( "svd" );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   //  parameters
   uint_t lvl2 = 10;
   uint_t lvl3 = 7;
   uint_t iter = 100;

   // precompute SVDs
   WALBERLA_LOG_INFO_ON_ROOT( "precompute SVDs" );
   WALBERLA_ROOT_SECTION()
   {
      std::filesystem::create_directory( "svd" );
   }
   for ( uint_t d : { 2, 3 } )
   {
      const auto lvl = ( d == 2 ) ? lvl2 : lvl3;
      precompute_SVDs< 0 >( lvl, d );
      precompute_SVDs< 1 >( lvl, d );
      precompute_SVDs< 2 >( lvl, d );
      precompute_SVDs< 3 >( lvl, d );
      precompute_SVDs< 4 >( lvl, d );
      precompute_SVDs< 5 >( lvl, d );
      precompute_SVDs< 6 >( lvl, d );
      precompute_SVDs< 7 >( lvl, d );
      precompute_SVDs< 8 >( lvl, d );
      precompute_SVDs< 9 >( lvl, d );
      precompute_SVDs< 10 >( lvl, d );
   }

   auto k = [&]( const hyteg::Point3D& x ) { return sin( 30 * x[0] ) * cos( 40 * x[1] ) * cos( 20 * x[2] ); };

   // run benchmarks
   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark: A.apply()" );
   for ( int d : { 2, 3 } )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      forms::p1_div_k_grad_affine_q3 form( k, k );
      auto                           lvl = ( d == 2 ) ? lvl2 : lvl3;
      benchmark< P1ConstantLaplaceOperator, false, false, false, false >( d, lvl, iter, form, k, "const. stencil", true, true );
      benchmark< operatorgeneration::P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab, false, false, false, false >(
          d, lvl, iter, form, k, "const. elwise.", false, true );
      benchmark< operatorgeneration::P1ElementwiseDivKGrad_cubes_const_vect_fused_quadloops_tab, false, true, false, false >(
          d, lvl, iter, form, k, "DivKGrad elwise.", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 0 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=0", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 1 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=1", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 2 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=2", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 3 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=3", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 4 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=4", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 5 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=5", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 6 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=6", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 7 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=7", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 8 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=8", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 9 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=9", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 10 >, true, false, true, false >(
          d, lvl, iter, form, k, "surrogate q=10", false, true );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark: A.smooth_gs()" );
   for ( int d : { 2, 3 } )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      forms::p1_div_k_grad_affine_q3 form( k, k );
      auto                           lvl = ( d == 2 ) ? lvl2 : lvl3;
      benchmark< P1ConstantLaplaceOperator, false, false, false, true >( d, lvl, iter, form, k, "const. stencil", true, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 0 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=0", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 1 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=1", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 2 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=2", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 3 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=3", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 4 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=4", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 5 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=5", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 6 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=6", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 7 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=7", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 8 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=8", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 9 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=9", false, true );
      benchmark< P1SurrogateAffineDivKGradOperator< 10 >, true, false, true, true >(
          d, lvl, iter, form, k, "surrogate q=10", false, true );
   }

   WALBERLA_ROOT_SECTION()
   {
      std::filesystem::remove_all( "svd" );
   }

   return 0;
}
