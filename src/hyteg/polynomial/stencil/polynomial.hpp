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
#include <hyteg/polynomial/elementwise/polynomial.hpp>

namespace hyteg {
namespace p1 {
namespace stencil {
namespace surrogate {

template < uint8_t DIM_domain, uint8_t DIM_primitive, uint8_t DEGREE >
class Polynomial : public std::array< StencilData< DIM_domain >, hyteg::surrogate::polynomial::dimP( DIM_primitive, DEGREE ) >
{
 public:
   static constexpr uint_t n_coeff   = hyteg::surrogate::polynomial::dimP( DIM_primitive, DEGREE );
   static constexpr uint_t n_stencil = stencilSize( DIM_domain );
   template < typename T = real_t >
   using Stencil = StencilData< DIM_domain, T >;

   inline Polynomial()
   : std::array< Stencil<>, n_coeff >{}
   {}

   /**
    * @brief Sets the coefficients of the polynomial.
    *
    * @param all_coeffs A stencil of vectors containing the coefficients.
    *          Must provide functions `size()` and `operator[]`
    */
   template < typename CoeffVector >
   inline void set_coefficients( const Stencil< CoeffVector >& all_coeffs )
   {
      for ( auto& coeffs : all_coeffs )
      {
         if ( uint_t( coeffs.size() ) != n_coeff )
         {
            WALBERLA_ABORT( "Number of coefficients must match dimension of polynomial space" );
         }
      }
      using coeff_idx_t = decltype( all_coeffs[0].size() );

      for ( uint_t i = 0; i < n_coeff; ++i ) // polynomial coefficients
      {
         for ( uint_t d = 0; d < n_stencil; ++d ) // stencil direction
         {
            ( *this )[i][d] = static_cast< real_t >( all_coeffs[d][static_cast< coeff_idx_t >( i )] );
         }
      }
   }

   template < typename CoeffVector >
   inline void set_coefficients( p1::stencil::Dir d, CoeffVector& coeffs )
   {
      if ( uint_t( coeffs.size() ) != n_coeff )
      {
         WALBERLA_ABORT( "Number of coefficients must match dimension of polynomial space" );
      }
      using coeff_idx_t = decltype( coeffs.size() );

      for ( uint_t i = 0; i < n_coeff; ++i ) // polynomial coefficients
      {
         ( *this )[i][d] = static_cast< real_t >( coeffs[static_cast< coeff_idx_t >( i )] );
      }
   }

   // fix z coordinate s.th. only 2d polynomial must be evaluated
   inline Polynomial< DIM_domain, DIM_primitive - 1, DEGREE > fix_z( const real_t z ) const
   {
      static_assert( DIM_primitive == 3, "fix_z(z) only available in 3D!" );
      return fix_coord( z );
   }

   // fix y coordinate s.th. only 1d polynomial must be evaluated
   inline Polynomial< DIM_domain, DIM_primitive - 1, DEGREE > fix_y( const real_t y ) const
   {
      static_assert( DIM_primitive == 2, "fix_y(y) only available in 2D!" );
      return fix_coord( y );
   }

   /** @brief Evaluate the 1d polynomial at x
    * @param x The x-coordinate.
    */
   inline void eval( const real_t x ) const
   {
      static_assert( DIM_primitive == 1, "eval(x) only available in 1D!" );
      for ( uint_t d = 0; d < n_stencil; ++d )
      {
         _result[d] = ( *this )[0][d];
      }
      auto xpow = x;
      for ( uint_t i = 1; i < n_coeff; ++i )
      {
         for ( uint_t d = 0; d < n_stencil; ++d )
         {
            _result[d] += ( *this )[i][d] * xpow;
         }
         xpow *= x;
      }
   }

   /** @brief Vectorized version of eval(). Evaluate the 1d polynomial at x1,x2,x3,x4
    * @param x The x-coordinates.
    */
   inline void eval_vec( const std::array< real_t, 4 >& x ) const
   {
      static_assert( DIM_primitive == 1, "eval(x) only available in 1D!" );
      for ( uint_t d = 0; d < n_stencil; ++d )
      {
         _result_vec[d][0] = ( *this )[0][d];
         _result_vec[d][1] = ( *this )[0][d];
         _result_vec[d][2] = ( *this )[0][d];
         _result_vec[d][3] = ( *this )[0][d];
      }
      auto xpow = x;
      for ( uint_t i = 1; i < n_coeff; ++i )
      {
         for ( uint_t d = 0; d < n_stencil; ++d )
         {
            _result_vec[d][0] += ( *this )[i][d] * xpow[0];
            _result_vec[d][1] += ( *this )[i][d] * xpow[1];
            _result_vec[d][2] += ( *this )[i][d] * xpow[2];
            _result_vec[d][3] += ( *this )[i][d] * xpow[3];
         }
         xpow[0] *= x[0];
         xpow[1] *= x[1];
         xpow[2] *= x[2];
         xpow[3] *= x[3];
      }
   }

   // return the result of last call of eval(x)
   inline const Stencil<>& px() const
   {
      static_assert( DIM_primitive == 1, "Final result can only be obtained from 1D polynomial!" );
      return _result;
   }

   // return the result of last call of eval_vec(x)
   inline const StencilData< DIM_domain, std::array< real_t, 4 > >& px_vec() const
   {
      static_assert( DIM_primitive == 1, "Final result can only be obtained from 1D polynomial!" );
      return _result_vec;
   }

   // evaluate polynomial by summing up basis functions, only use for debugging or testing
   Stencil<> eval_naive( const Point3D& x ) const
   {
      // monomial basis
      constexpr hyteg::surrogate::polynomial::Basis< DEGREE > phi;

      Stencil<> px{};
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
   inline Polynomial< DIM_domain, DIM_primitive - 1, DEGREE > fix_coord( const real_t z ) const
   {
      // monomial basis
      constexpr hyteg::surrogate::polynomial::Basis< DEGREE > phi;

      // z^k for k=0,...,q
      std::array< real_t, DEGREE + 1 > z_pow;
      z_pow[0] = 1.0;
      for ( uint_t k = 1; k <= DEGREE; ++k )
      {
         z_pow[k] = z_pow[k - 1] * z;
      }

      // lower dimensional polynomial
      Polynomial< DIM_domain, DIM_primitive - 1, DEGREE > result;

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

   // location to write result
   mutable Stencil<>                                          _result;
   mutable StencilData< DIM_domain, std::array< real_t, 4 > > _result_vec;
};

// Base case for _restriction
template < uint8_t DIM_domain, uint8_t DEGREE >
class Polynomial< DIM_domain, 0, DEGREE >
{};

} // namespace surrogate
} // namespace stencil
} // namespace p1
} // namespace hyteg
