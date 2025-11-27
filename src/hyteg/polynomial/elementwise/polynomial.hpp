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
#include <core/Abort.h>
#include <core/DataTypes.h>
#include <hyteg/indexing/Common.hpp>
#include <hyteg/types/PointND.hpp>

namespace hyteg {
namespace surrogate {

using walberla::uint_t;

namespace polynomial {

/**
 * @brief Compute the dimension of P_q(R^d)
 */
static constexpr uint_t dimP( uint8_t d, uint8_t q )
{
   switch ( d )
   {
   case 1:
      return q + 1;
   case 2:
      return ( q + 1 ) * ( q + 2 ) / 2;
   case 3:
      return ( q + 1 ) * ( q + 2 ) * ( q + 3 ) / 6;
   default:
      return 0;
   }
}

// store x^i*y^j*z^k by exponents i,j,k compressed into 16 bit
struct Monomial
{
 private:
   using compression_t = int16_t;

   // bit masks for exponents
   static constexpr int SHIFT_J = 5;
   static constexpr int SHIFT_K = 10;
   static constexpr int I       = 1;
   static constexpr int J       = 1 << SHIFT_J;
   static constexpr int K       = 1 << SHIFT_K;

   // internal data
   compression_t ijk_compressed;

 public:
   constexpr inline Monomial()
   : ijk_compressed( 0 )
   {}

   constexpr inline Monomial( const uint_t i, const uint_t j, const uint_t k )
   : ijk_compressed( compression_t( I * i | J * j | K * k ) )
   {}

   // extract exponent i of x^i*y^j*z^k
   constexpr inline int i() const { return ijk_compressed & ( J - 1 ); }

   // extract exponent j of x^i*y^j*z^k
   constexpr inline int j() const { return ( ijk_compressed & ( K - 1 ) ) >> SHIFT_J; }

   // extract exponent k of x^i*y^j*z^k
   constexpr inline int k() const { return ijk_compressed >> SHIFT_K; }

   // extract exponents i,j,k from compressed basis function
   constexpr inline std::array< int, 3 > expand() const { return { i(), j(), k() }; }

   // total polynomial degree of this basis function
   constexpr inline int degree() const { return i() + j() + k(); }

   // evaluate basis function at x
   template < typename FLOAT >
   constexpr inline FLOAT eval( const PointND< FLOAT, 3 >& x ) const
   {
      const auto [i, j, k] = expand();

      FLOAT val = static_cast< FLOAT >( 1.0 );

      for ( int p = 0; p < i; ++p )
         val *= x[0]; // compute x^i
      for ( int p = 0; p < j; ++p )
         val *= x[1]; // compute y^j
      for ( int p = 0; p < k; ++p )
         val *= x[2]; // compute z^k

      return val;
   }
};

// monomial basis of polynomial space, given by the exponents of the monomials
template < uint8_t DEGREE >
struct Basis : public std::array< Monomial, dimP( 3, DEGREE ) >
{
   constexpr Basis()
   {
      // 1d basis
      for ( uint8_t i = 0; i <= DEGREE; ++i )
      {
         ( *this )[i] = Monomial( i, 0, 0 );
      }
      // extend to 2d
      auto ij = dimP( 1, DEGREE );
      for ( uint8_t i = 0; i < DEGREE; ++i )
      {
         for ( uint8_t j = 1; i + j <= DEGREE; ++j )
         {
            ( *this )[ij++] = Monomial( i, j, 0 );
         }
      }
      // extend to 3d
      auto ijk = dimP( 2, DEGREE );
      for ( uint8_t i = 0; i < DEGREE; ++i )
      {
         for ( uint8_t k = 1; i + k <= DEGREE; ++k )
         {
            ( *this )[ijk++] = Monomial( i, 0, k );
         }
      }
      for ( uint8_t i = 0; i < DEGREE; ++i )
      {
         for ( uint8_t j = 1; i + j < DEGREE; ++j )
         {
            for ( uint8_t k = 1; i + j + k <= DEGREE; ++k )
            {
               ( *this )[ijk++] = Monomial( i, j, k );
            }
         }
      }
   }
};

/* Domain of polynomial space = [-1,3]
 * @brief Affine mapping from {v}_v ∈ macro-edge to [-1,3]
 */
template < typename FLOAT >
struct Domain
{
   constexpr Domain( uint_t lvl )
   : scaling( FLOAT( 4.0 ) / FLOAT( ( 1 << lvl ) - 1 ) )
   {}
   // convert index i to coordinate x ∈ [-1,3]
   inline constexpr FLOAT operator[]( idx_t i ) const { return scaling * FLOAT( i ) - shift; }

   // convert index [i,j,k] to coordinate x ∈ [-1,3]^3
   inline PointND< FLOAT, 3 > operator()( const indexing::Index& idx ) const
   {
      return { ( *this )[idx.x()], ( *this )[idx.y()], ( *this )[idx.z()] };
   }

   // scaling factor to convert index i ∈ {0,..., 2^lvl - 1} to coordinate x ∈ [-1,3]
   const FLOAT            scaling;
   static constexpr FLOAT shift = FLOAT( 1.0 );
};

template < typename FLOAT, uint8_t DIM, uint8_t DEGREE >
class Polynomial : public std::array< FLOAT, dimP( DIM, DEGREE ) >
{
 public:
   inline Polynomial()
   : std::array< FLOAT, dimP( DIM, DEGREE ) >{}
   {}

   template < typename CoeffVector >
   inline Polynomial( const CoeffVector& coeffs )
   {
      set_coefficients( coeffs );
   }

   /**
    * @brief Sets the coefficients of the polynomial.
    *
    * @param coeffs A vector containing the coefficients.
    *          Must provide functions `size()` and `operator[]`
    */
   template < typename CoeffVector >
   inline void set_coefficients( const CoeffVector& coeffs )
   {
      if ( uint_t( coeffs.size() ) != this->size() )
      {
         WALBERLA_ABORT( "Number of coefficients must match dimension of polynomial space" );
      }

      using coeff_idx_t = decltype( coeffs.size() );

      for ( uint_t i = 0; i < this->size(); ++i )
      {
         ( *this )[i] = static_cast< FLOAT >( coeffs[static_cast< coeff_idx_t >( i )] );
      }
   }

   // fix z coordinate s.th. only 2d polynomial must be evaluated
   inline void fix_z( const FLOAT z ) const
   {
      static_assert( DIM == 3, "fix_z can only be used in 3d!" );
      fix_coord( z );
   }

   // fix y coordinate s.th. only 1d polynomial must be evaluated
   inline void fix_y( const FLOAT y ) const
   {
      if constexpr ( DIM == 2 )
      {
         fix_coord( y );
      }
      else // DIM == 3
      {
         _restriction.fix_y( y );
      }
   }

   /** @brief Evaluate the 1d polynomial at x using Horner's method.
    *
    * Usage: `p.fix_z(z);` `p.fix_y(y);` `p.eval(x);`
    *
    * @param x The x-coordinate.
    * @return p|_zy(x)
    */
   inline FLOAT eval( const FLOAT x ) const { return _restriction.eval( x ); }

   // evaluate polynomial by summing up basis functions, only use for debugging or testing
   FLOAT eval_naive( const PointND< FLOAT, 3 >& x ) const
   {
      // monomial basis
      constexpr Basis< DEGREE > phi;

      FLOAT p_xyz = 0.0;
      for ( uint_t i = 0; i < this->size(); ++i )
      {
         p_xyz += ( *this )[i] * phi[i].eval( x );
      }
      return p_xyz;
   }

   std::string to_string() const
   {
      // monomial basis
      constexpr Basis< DEGREE > phi;

      std::ostringstream oss;
      for ( uint_t n = 0; n < this->size(); ++n )
      {
         const auto [i,j,k] = phi[n].expand();
         const auto c = ( *this )[n];

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

   inline const Polynomial< FLOAT, 2, DEGREE >& get_2d_restriction() const
   {
      if constexpr ( DIM == 2 )
      {
         return ( *this );
      }
      else // DIM == 3
      {
         return _restriction;
      }
   }

   inline const Polynomial< FLOAT, 1, DEGREE >& get_1d_restriction() const
   {
      if constexpr ( DIM == 2 )
      {
         return _restriction;
      }
      else // DIM == 3
      {
         return _restriction.get_1d_restriction();
      }
   }

 private:
   // fix coordinate s.th. only lower dimensional polynomial must be evaluated
   inline void fix_coord( const FLOAT z ) const
   {
      // monomial basis
      constexpr Basis< DEGREE > phi;

      // z^k for k=0,...,q
      std::array< FLOAT, DEGREE + 1 > z_pow;
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

         _restriction[ij] = ( *this )[ij]; // k=0
         for ( int k = 1; k <= max_k; ++k )
         {
            _restriction[ij] += ( *this )[ijk++] * z_pow[k];
         }
      }
   }

   // restriction to lower spatial dimension
   mutable Polynomial< FLOAT, DIM - 1, DEGREE > _restriction;
};

// partial specialization for 1D polynomial
template < typename FLOAT, uint8_t DEGREE >
class Polynomial< FLOAT, 1, DEGREE > : public std::array< FLOAT, dimP( 1, DEGREE ) >
{
 public:
   /** @brief Evaluate the 1d polynomial at x using Horner's method.
    * @param x The x-coordinate.
    * @return p(x)
   */
   inline FLOAT eval( const FLOAT x ) const
   {
      if constexpr ( DEGREE == 0 )
      {
         return ( *this )[0];
      }
      else
      {
         auto px = ( *this )[this->size() - 1];
         for ( int i = this->size() - 2; i >= 0; --i )
         {
            px = px * x + ( *this )[static_cast< uint_t >( i )];
         }
         return px;
      }
   }
};

} // namespace polynomial
} // namespace surrogate
} // namespace hyteg
