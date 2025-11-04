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
   inline void fix_z( const real_t z ) const
   {
      static_assert( DIM_primitive == 3, "fix_z can only be used in 3d!" );
      fix_coord( z );
   }

   // fix y coordinate s.th. only 1d polynomial must be evaluated
   inline void fix_y( const real_t y ) const
   {
      if constexpr ( DIM_primitive == 2 )
      {
         fix_coord( y );
      }
      else // DIM == 3
      {
         _restriction.fix_y( y );
      }
   }

   /** @brief Evaluate the 1d polynomial at x
    *
    * Usage: `p.fix_z(z);` `p.fix_y(y);` `p.eval(x);`
    *
    * @param x The x-coordinate.
    */
   inline void eval( const real_t x ) const
   {
      if constexpr ( DIM_primitive == 1 )
      {
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
      else
      {
         _restriction.eval( x );
      }
   }

   // return the result of last call of eval(x)
   inline const Stencil<>& px() const
   {
      if constexpr ( DIM_primitive == 1 )
      {
         return _result;
      }
      else
      {
         return _restriction.px();
      }
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
   inline void fix_coord( const real_t z ) const
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

      // first index where the 3d extension starts
      auto ijk = _restriction.size();
      // iterate over coefficients of 2d polynomial
      for ( uint_t ij = 0; ij < _restriction.size(); ++ij )
      {
         const auto max_k = DEGREE - phi[ij].degree();

         for ( uint_t d = 0; d < n_stencil; ++d )
         {
            _restriction[ij][d] = ( *this )[ij][d]; // k=0
         }
         for ( int k = 1; k <= max_k; ++k )
         {
            for ( uint_t d = 0; d < n_stencil; ++d )
            {
               _restriction[ij][d] += ( *this )[ijk][d] * z_pow[k];
            }
            ++ijk;
         }
      }
   }

   // restriction to lower spatial dimension
   mutable Polynomial< DIM_domain, DIM_primitive - 1, DEGREE > _restriction;
   // location to write result
   mutable Stencil<> _result;
};

// Base case for _restriction
template < uint8_t DIM_domain, uint8_t DEGREE >
class Polynomial< DIM_domain, 0, DEGREE >
{};

} // namespace surrogate
} // namespace stencil
} // namespace p1
} // namespace hyteg
