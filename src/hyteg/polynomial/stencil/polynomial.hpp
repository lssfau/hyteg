/*
 * Copyright (c) 2024-2025 Benjamin Mann.
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

#pragma once
#include <hyteg/p1functionspace/globalIndices.hpp>
#include <hyteg/polynomial/new/polynomial.hpp>
#include <simd/SIMD.h>
#include <waLBerlaDefinitions.h>

namespace hyteg {
namespace p1 {
namespace stencil {
namespace surrogate {

template < uint8_t DIM_domain, uint8_t DIM_primitive, uint8_t DEGREE, typename FLOAT = real_t >
class Polynomial
: public std::array< StencilData< DIM_domain, FLOAT >, hyteg::surrogate::polynomial::dimP( DIM_primitive, DEGREE ) >
{
 public:
   static constexpr uint_t n_coeff   = hyteg::surrogate::polynomial::dimP( DIM_primitive, DEGREE );
   static constexpr uint_t n_stencil = stencilSize( DIM_domain );
   using Stencil                     = StencilData< DIM_domain, FLOAT >;

   inline Polynomial()
   : std::array< Stencil, n_coeff >{}
   {}

   inline void set_coefficients( p1::stencil::Dir d, auto& coeffs )
   {
      if ( uint_t( coeffs.size() ) != n_coeff )
      {
         WALBERLA_ABORT( "Number of coefficients must match dimension of polynomial space" );
      }
      using coeff_idx_t = decltype( coeffs.size() );

      for ( uint_t i = 0; i < n_coeff; ++i ) // polynomial coefficients
      {
         ( *this )[i][d] = static_cast< FLOAT >( coeffs[static_cast< coeff_idx_t >( i )] );
      }
   }

   // fix z coordinate s.th. only 2d polynomial must be evaluated
   inline Polynomial< DIM_domain, DIM_primitive - 1, DEGREE, FLOAT > fix_z( const FLOAT z ) const
   {
      static_assert( DIM_primitive == 3, "fix_z(z) only available in 3D!" );
      return fix_coord( z );
   }

   // fix y coordinate s.th. only 1d polynomial must be evaluated
   inline Polynomial< DIM_domain, DIM_primitive - 1, DEGREE, FLOAT > fix_y( const FLOAT y ) const
   {
      static_assert( DIM_primitive == 2, "fix_y(y) only available in 2D!" );
      return fix_coord( y );
   }

   /** @brief Evaluate the 1d polynomial at x
    * @param x The x-coordinate.
    */
   inline Stencil eval( const FLOAT x ) const
   {
      static_assert( DIM_primitive == 1, "eval(x) only available in 1D!" );

      Stencil result = ( *this )[0];

      auto xpow = x;
      for ( uint_t i = 1; i <= DEGREE; ++i )
      {
         for ( uint_t d = 0; d < n_stencil; ++d )
         {
            result[d] += ( *this )[i][d] * xpow;
         }
         xpow *= x;
      }

      return result;
   }

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

   /** @brief Evaluate the 1d polynomial at x
    * @param x The x-coordinates.
    */
   inline StencilData< DIM_domain, walberla::simd::double4_t > eval_vec( const std::array< FLOAT, 4 >& x ) const
   {
      static_assert( DIM_primitive == 1, "eval(x) only available in 1D!" );
      static_assert( std::is_same< FLOAT, double >::value, "eval_vec only available for double precision!" );

      StencilData< DIM_domain, walberla::simd::double4_t > result;
      for ( uint_t d = 0; d < n_stencil; ++d )
      {
         result[d] = walberla::simd::make_double4( ( *this )[0][d] );
      }

      const auto vx   = walberla::simd::load_unaligned( x.data() );
      auto       xpow = vx;
      for ( uint_t i = 1; i <= DEGREE; ++i )
      {
         for ( uint_t d = 0; d < n_stencil; ++d )
         {
            const auto ci = walberla::simd::make_double4( ( *this )[i][d] );
            result[d]     = result[d] + ci * xpow;
         }
         xpow = xpow * vx;
      }

      return result;
   }
#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

   // evaluate polynomial by summing up basis functions, only use for debugging or testing
   Stencil eval_naive( const PointND< FLOAT, 3 >& x ) const
   {
      // monomial basis
      constexpr hyteg::surrogate::polynomial::Basis< DEGREE > phi;

      Stencil px{};
      for ( uint_t i = 0; i < n_coeff; ++i ) // polynomial coefficients
      {
         auto phi_i_x = phi[i].eval( x );
         for ( uint_t d = 0; d < n_stencil; ++d ) // stencil direction
         {
            px[d] += ( *this )[i][d] * phi_i_x;
         }
      }
      return px;
   }

   std::string to_string( const stencil::Dir& dir ) const
   {
      // monomial basis
      constexpr hyteg::surrogate::polynomial::Basis< DEGREE > phi;

      std::ostringstream oss;
      for ( uint_t n = 0; n < n_coeff; ++n )
      {
         const auto [i, j, k] = phi[n].expand();
         const auto c         = ( *this )[n][dir];

         if ( n > 0 && c >= 0 )
         {
            oss << " + ";
         }
         else if ( c < 0 )
         {
            oss << " - ";
         }
         oss << std::scientific << std::abs( c );

         if ( i > 0 )
         {
            oss << "*x";
            if ( i > 1 )
            {
               oss << "^" << i;
            }
         }

         if ( j > 0 )
         {
            oss << "*y";
            if ( j > 1 )
            {
               oss << "^" << j;
            }
         }

         if ( k > 0 )
         {
            oss << "*z";
            if ( k > 1 )
            {
               oss << "^" << k;
            }
         }
      }
      return oss.str();
   }

 private:
   // fix coordinate s.th. only lower dimensional polynomial must be evaluated
   inline Polynomial< DIM_domain, DIM_primitive - 1, DEGREE, FLOAT > fix_coord( const FLOAT z ) const
   {
      // monomial basis
      constexpr hyteg::surrogate::polynomial::Basis< DEGREE > phi;

      // z^k for k=0,...,q
      std::array< FLOAT, DEGREE + 1 > z_pow;
      z_pow[0] = 1.0;
      for ( uint_t k = 1; k <= DEGREE; ++k )
      {
         z_pow[k] = z_pow[k - 1] * z;
      }

      // lower dimensional polynomial
      Polynomial< DIM_domain, DIM_primitive - 1, DEGREE, FLOAT > result;

      // first index where the 3d extension starts
      auto ijk = result.n_coeff;
      // iterate over coefficients of 2d polynomial
      for ( uint_t ij = 0; ij < result.n_coeff; ++ij )
      {
         const auto max_k = DEGREE - phi[ij].degree();

         for ( uint_t d = 0; d < n_stencil; ++d )
         {
            result[ij][d] = ( *this )[ij][d]; // k=0
         }
         for ( int k = 1; k <= max_k; ++k )
         {
            for ( uint_t d = 0; d < n_stencil; ++d )
            {
               result[ij][d] += ( *this )[ijk][d] * z_pow[k];
            }
            ++ijk;
         }
      }

      return result;
   }
};

} // namespace surrogate
} // namespace stencil
} // namespace p1
} // namespace hyteg
