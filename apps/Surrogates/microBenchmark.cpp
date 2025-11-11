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
#include <core/logging/Logging.h>
#include <core/timing/Timer.h>

#define RESTRICT WALBERLA_RESTRICT

using walberla::real_t;
using walberla::uint_t;

constexpr uint_t stencil_size = 15;
template <typename T>
using Stencil = std::array<T,stencil_size>;
template <int q>
using poly = std::array<std::array<real_t, stencil_size>, q+1>; // stencil-valued polynomial
template <int q>
using poly_transposed = std::array<std::array<real_t, q+1>, stencil_size>; // stencil of polynomials

template <int q>
real_t apply_standard(const uint_t iter, const uint_t n_dof, const uint_t boundary, const poly<q> surrogate, const Stencil<int> offsets, const real_t* RESTRICT const src, real_t* RESTRICT dst)
{
   walberla::WcTimer timer;
   timer.start();
   for ( uint_t i = 0; i < iter; ++i )
   {
      Stencil<real_t> stencil{};
      for (uint_t j = boundary; j < n_dof - boundary; ++j)
      {
         const auto x = real_t(j)/real_t(n_dof);

         // eval poly
         for ( uint_t d = 0; d < stencil_size; ++d )
         {
            stencil[d] = surrogate[0][d];
         }
         auto xpow = x;
         for ( uint_t k = 1; k <= q; ++k )
         {
            for ( uint_t d = 0; d < stencil_size; ++d )
            {
               stencil[d] += surrogate[k][d] * xpow;
            }
            xpow *= x;
         }

         // apply stencil
         for (uint_t dir = 0; dir < stencil_size; ++dir)
         {
            dst[j] += stencil[dir] * src[j + offsets[dir]];
         }
      }
   }
   timer.end();
   if (q > 17) std::cout << dst[n_dof/2]; // prevent the loop from being optimized out;
   return  timer.total() / iter;
}

template <int q>
real_t apply_unrolled(const uint_t iter, const uint_t n_dof, const uint_t boundary, const poly<q> surrogate, const Stencil<int> offsets, const real_t* RESTRICT const src, real_t* RESTRICT dst)
{
   walberla::WcTimer timer;
   timer.start();
   for ( uint_t i = 0; i < iter; ++i )
   {
      Stencil<real_t> stencil{};
      for (uint_t j = boundary; j < n_dof - boundary; ++j)
      {
         const auto x = real_t(j)/real_t(n_dof);

         // eval poly
         for ( uint_t d = 0; d < stencil_size; ++d )
         {
            stencil[d] = surrogate[0][d];
         }
         auto xpow = x;
         for ( uint_t k = 1; k <= q; ++k )
         {
            stencil[0]  += surrogate[k][0]  * xpow;
            stencil[1]  += surrogate[k][1]  * xpow;
            stencil[2]  += surrogate[k][2]  * xpow;
            stencil[3]  += surrogate[k][3]  * xpow;
            stencil[4]  += surrogate[k][4]  * xpow;
            stencil[5]  += surrogate[k][5]  * xpow;
            stencil[6]  += surrogate[k][6]  * xpow;
            stencil[7]  += surrogate[k][7]  * xpow;
            stencil[8]  += surrogate[k][8]  * xpow;
            stencil[9]  += surrogate[k][9]  * xpow;
            stencil[10] += surrogate[k][10] * xpow;
            stencil[11] += surrogate[k][11] * xpow;
            stencil[12] += surrogate[k][12] * xpow;
            stencil[13] += surrogate[k][13] * xpow;
            stencil[14] += surrogate[k][14] * xpow;
            xpow *= x;
         }

         // apply stencil
         for (uint_t dir = 0; dir < stencil_size; ++dir)
         {
            dst[j] += stencil[dir] * src[j + offsets[dir]];
         }
      }
   }
   timer.end();
   if (q > 17) std::cout << dst[n_dof/2]; // prevent the loop from being optimized out;
   return  timer.total() / iter;
}

template <int q>
real_t apply_transposed(const uint_t iter, const uint_t n_dof, const uint_t boundary, const poly_transposed<q> surrogate, const Stencil<int> offsets, const real_t* RESTRICT const src, real_t* RESTRICT dst)
{
   walberla::WcTimer timer;
   timer.start();
   for ( uint_t i = 0; i < iter; ++i )
   {
      Stencil<real_t> stencil{};
      for (uint_t j = boundary; j < n_dof - boundary; ++j)
      {
         const auto x = real_t(j)/real_t(n_dof);

         // eval poly
         for ( uint_t d = 0; d < stencil_size; ++d )
         {
            stencil[d] = surrogate[d][0];
         }
         auto xpow = x;
         for ( uint_t k = 1; k <= q; ++k )
         {
            for ( uint_t d = 0; d < stencil_size; ++d )
            {
               stencil[d] += surrogate[d][k] * xpow;
            }
            xpow *= x;
         }

         // apply stencil
         for (uint_t dir = 0; dir < stencil_size; ++dir)
         {
            dst[j] += stencil[dir] * src[j + offsets[dir]];
         }
      }
   }
   timer.end();
   if (q > 17) std::cout << dst[n_dof/2]; // prevent the loop from being optimized out;
   return  timer.total() / iter;
}

template <int q>
real_t apply_fd(const uint_t iter, const uint_t n_dof, const uint_t boundary, const poly<q> surrogate, const Stencil<int> offsets, const real_t* RESTRICT const src, real_t* RESTRICT dst)
{
   walberla::WcTimer timer;
   timer.start();
   for ( uint_t i = 0; i < iter; ++i )
   {
      Stencil<real_t> stencil{};
      for (uint_t j = boundary; j < n_dof - boundary; ++j)
      {
         const auto x = real_t(j)/real_t(n_dof);

         // eval poly
         //todob

         // apply stencil
         for (uint_t dir = 0; dir < stencil_size; ++dir)
         {
            dst[j] += stencil[dir] * src[j + offsets[dir]];
         }
      }
   }
   timer.end();
   if (q > 17) std::cout << dst[n_dof/2]; // prevent the loop from being optimized out;
   return  timer.total() / iter;
}

template <int q>
void benchmark(const uint_t n_edge, const uint_t iter)
{

   const auto n_dof = n_edge*n_edge*n_edge;
   const auto n_dof_global = walberla::MPIManager::instance()->numProcesses() * n_dof;

   WALBERLA_LOG_INFO_ON_ROOT( "microbenchmark, q=" << q << ", global DoF: " << n_dof_global);

   const auto dy = int(n_edge);
   const auto dz = int(n_edge*n_edge);

   const Stencil<int> offsets{0, -1, 1, -dy, -dy+1, dy-1, dy, dz, dz-1, dz-dy, dz-dy+1, -dz, -dz+1, -dz+dy-1, -dz+dy};
   poly<q> surrogate{};
   poly_transposed<q> surrogate_transposed{};
   for (uint_t dir = 0; dir < 15; ++dir)
   {
      for (uint_t k = 0; k <= q; ++k)
      {
         surrogate[k][dir] = walberla::math::realRandom();
         surrogate_transposed[dir][k] = surrogate[k][dir];
      }
   }

   const auto boundary = uint_t(1 + dy + dz);

   std::vector<real_t> src(2*boundary + n_dof);
   std::vector<real_t> dst(2*boundary + n_dof);
   for (auto& val : src) {
      val = walberla::math::realRandom();
   }


   auto t = apply_standard<q>(iter, n_dof, boundary, surrogate, offsets, src.data(), dst.data());
   WALBERLA_LOG_INFO_ON_ROOT( "standard:    " << t << " s, " <<  int( n_dof_global / t * 1e-6 ) << " MDoF/s");

   t = apply_unrolled<q>(iter, n_dof, boundary, surrogate, offsets, src.data(), dst.data());
   WALBERLA_LOG_INFO_ON_ROOT( "unrolled:    " << t << " s, " <<  int( n_dof_global / t * 1e-6 ) << " MDoF/s");

   t = apply_transposed<q>(iter, n_dof, boundary, surrogate_transposed, offsets, src.data(), dst.data());
   WALBERLA_LOG_INFO_ON_ROOT( "transposed:  " << t << " s, " <<  int( n_dof_global / t * 1e-6 ) << " MDoF/s");

   t = apply_fd<q>(iter, n_dof, boundary, surrogate, offsets, src.data(), dst.data());
   WALBERLA_LOG_INFO_ON_ROOT( "forw. diff.: " << t << " s, " <<  int( n_dof_global / t * 1e-6 ) << " MDoF/s");
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   uint_t edge = 67;
   if (argc > 1)
   {
      edge = atoi(argv[1]);
   }


   benchmark<0>(edge, 100);
   benchmark<1>(edge, 100);
   benchmark<2>(edge, 100);
   benchmark<3>(edge, 100);
   benchmark<4>(edge, 100);
   benchmark<5>(edge, 100);
   benchmark<6>(edge, 100);

   return 0;
}
