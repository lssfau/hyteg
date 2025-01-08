/*
 * Copyright (c) 2024 Benjamin Mann.
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
#include <core/DataTypes.h>

namespace hyteg {
namespace surrogate {

using walberla::uint_t;

namespace polynomial {

/**
 * @brief Compute the dimension of a polynomial space of degree q over a d-dimensional domain
 */
static constexpr uint_t dimP( uint8_t d, uint8_t q )
{
   return ( d == 0 ) ? 0 : ( d == 1 ) ? q + 1 : ( d == 2 ) ? ( q + 1 ) * ( q + 2 ) / 2 : ( q + 1 ) * ( q + 2 ) * ( q + 3 ) / 6;
}

// store x^i*y^j*z^k by exponents i,j,k compressed into 16 bit
struct Monomial
{
   using compression_t = uint16_t;

   constexpr inline Monomial()
   : ijk_compressed( 0 )
   {}

   constexpr inline Monomial( const compression_t i, const compression_t j, const compression_t k )
   : ijk_compressed( I * i | J * j | K * k )
   {}

   // extract exponents i,j,k from compressed basis function
   constexpr inline std::array< compression_t, 3 > expand() const { return { i(), j(), k() }; }
   // extract exponent i of x^i*y^j*z^k
   constexpr inline compression_t i() const { return ijk_compressed & compression_t( J - 1 ); }
   // extract exponent j of x^i*y^j*z^k
   constexpr inline compression_t j() const { return ( ijk_compressed & compression_t( K - 1 ) ) >> 5; }
   // extract exponent k of x^i*y^j*z^k
   constexpr inline compression_t k() const { return ijk_compressed >> 10; }

   // polynomial degree of this basis function
   constexpr inline compression_t degree() const { return i() + j() + k(); }

   // evaluate basis function at x
   constexpr inline double eval( const std::array< double, 3 >& x ) const
   {
      const auto [i, j, k] = expand();

      double val = 1.0;

      for ( uint8_t p = 0; p < i; ++p )
         val *= x[0]; // compute x^i
      for ( uint8_t p = 0; p < j; ++p )
         val *= x[1]; // compute y^j
      for ( uint8_t p = 0; p < k; ++p )
         val *= x[2]; // compute z^k

      return val;
   }

 private:
   compression_t ijk_compressed;
   // bit masks for exponents
   static constexpr compression_t I = 1;
   static constexpr compression_t J = 1 << 5;
   static constexpr compression_t K = 1 << 10;
};

// monomial basis of polynomial space, given by the exponents of the monomials
struct Basis : public std::vector< Monomial >
{
   const uint8_t _q; // polynomial degree

   Basis( uint8_t q = 0 )
   : std::vector< Monomial >( dimP( 3, q ) )
   , _q( q )
   {
      // 1d basis
      for ( uint8_t i = 0; i <= _q; ++i )
      {
         ( *this )[i] = Monomial( i, 0, 0 );
      }
      // extend to 2d
      auto ij = dimP( 1, _q );
      for ( uint8_t i = 0; i < _q; ++i )
      {
         for ( uint8_t j = 1; i + j <= _q; ++j )
         {
            ( *this )[ij++] = Monomial( i, j, 0 );
         }
      }
      // extend to 3d
      auto ijk = dimP( 2, _q );
      for ( uint8_t i = 0; i < _q; ++i )
      {
         for ( uint8_t k = 1; i + k <= _q; ++k )
         {
            ( *this )[ijk++] = Monomial( i, 0, k );
         }
      }
      for ( uint8_t i = 0; i < _q; ++i )
      {
         for ( uint8_t j = 1; i + j < _q; ++j )
         {
            for ( uint8_t k = 1; i + j + k <= _q; ++k )
            {
               ( *this )[ijk++] = Monomial( i, j, k );
            }
         }
      }
   }
};

/* Coordinate system for polynomial evaluation.
 ! Different from the coordinates of PDE domain !
*/
struct Coordinates
{
   constexpr Coordinates( uint_t lvl )
   : scaling( 4.0 / double( ( 1 << lvl ) - 1 ) )
   {}
   // convert index i to coordinate x ∈ [-1,3]
   constexpr double x( uint_t i ) const { return scaling * i - 1.0; }
   // convert indices i,j,k to coordinate x ∈ [-1,3]^3
   constexpr std::array< double, 3 > x( uint_t i, uint_t j, uint_t k ) const { return { x( i ), x( j ), x( k ) }; }

   // convert index i ∈ {0,..., 2^lvl - 1} to coordinate x ∈ [-1,3]
   const double scaling;
};

class Polynomial : public std::vector< double >
{
 public:
   static constexpr uint8_t X = 0;
   static constexpr uint8_t Y = 1;
   static constexpr uint8_t Z = 2;

   inline Polynomial( uint8_t d, uint8_t q )
   : std::vector< double >( dimP( d, q ) )
   , _d( d )
   , _q( q )
   , _restriction( ( d > 1 ) ? std::make_unique< Polynomial >( d - 1, q ) : nullptr )
   , _basis( basis_q[q] )
   {}

   // set coefficients
   inline Polynomial& operator=( const Eigen::Matrix< double, -1, 1, Eigen::ColMajor >& coeffs )
   {
      WALBERLA_ASSERT( coeffs.size() == size(), "Number of coefficients must match dimension of polynomial space" );

      for ( uint_t i = 0; i < size(); ++i )
      {
         ( *this )[i] = coeffs( i );
      }

      return *this;
   }

   // fix z coordinate s.th. only 2d polynomial must be evaluated
   void fix_z( const double z ) const
   {
      WALBERLA_ASSERT( _d == 3, "fix_z can only be used in 3d" );

      return fix_coord( z );
   }

   // fix y coordinate s.th. only 1d polynomial must be evaluated
   void fix_y( const double y ) const
   {
      WALBERLA_ASSERT( _d == 2 || _d == 3, "fix_z can only be used in 2d and 3d" );
      if ( _d == 2 )
      {
         return fix_coord( y );
      }
      else
      {
         return _restriction->fix_y( y );
      }
   }

   /** @brief Evaluate the 1d polynomial at x using Horner's method.
    *
    * Usage: `p.fix_z(z);` `p.fix_y(y);` `p.eval(x);`
    *
    * @param x The x-coordinate.
    * @return p|_zy(x)
    */
   inline double eval( const double x ) const
   {
      if ( _d == 1 )
      {
         auto px = ( *this )[size() - 1];
         for ( int i = int( size() - 2 ); i >= 0; --i )
         {
            px = px * x + ( *this )[i];
         }
         return px;
      }
      else if ( _d == 0 )
      {
         WALBERLA_ABORT( "0-d polynomial can't be evaluated!" )
      }
      else
      {
         return _restriction->eval( x );
      }
   }

   // evaluate polynomial by summing up basis functions, only use for debugging or testing
   double eval_naive( const std::array< double, 3 >& x ) const
   {
      double p_xyz = 0.0;
      for ( uint_t i = 0; i < size(); ++i )
      {
         p_xyz += ( *this )[i] * _basis[i].eval( x );
      }
      return p_xyz;
   }

   const Basis& basis() const { return _basis; }

 private:
   // fix coordinate s.th. only lower dimensional polynomial must be evaluated
   inline void fix_coord( const double z ) const
   {
      // z^k for k=0,...,q
      std::vector< double > z_pow( _q + 1, 1.0 );
      for ( uint_t k = 0; k < _q; ++k )
      {
         z_pow[k + 1] = z_pow[k] * z;
      }

      // first index where the 3d extension starts
      uint_t ijk = _restriction->size();
      // iterate over coefficients of 2d polynomial
      for ( uint_t ij = 0; ij < _restriction->size(); ++ij )
      {
         const auto d = _basis[ij].degree();

         ( *_restriction )[ij] = ( *this )[ij]; // k=0
         for ( uint_t k = 1; k <= _q - d; ++k )
         {
            ( *_restriction )[ij] += ( *this )[ijk++] * z_pow[k];
         }
      }
   }

   // dimension of domain
   const uint8_t _d;
   // degree of polynomial
   const uint8_t _q;
   // restriction to lower dimension
   std::unique_ptr< Polynomial > _restriction;
   // basis functions
   const Basis& _basis;

   // basis functions for q = 0,1,...
   static const std::vector< Basis > basis_q;
};

/** Initialization of bases
 *  todo: add bases for higher degree spaces if needed
 */
const std::vector< Basis > Polynomial::basis_q = { Basis( 0 ),
                                                   Basis( 1 ),
                                                   Basis( 2 ),
                                                   Basis( 3 ),
                                                   Basis( 4 ),
                                                   Basis( 5 ),
                                                   Basis( 6 ),
                                                   Basis( 7 ),
                                                   Basis( 8 ),
                                                   Basis( 9 ),
                                                   Basis( 10 ),
                                                   Basis( 11 ),
                                                   Basis( 12 ) };


// RxC matrix of polynomials
template <uint_t R, uint_t C=R>
using Matrix = std::array<std::array<Polynomial, C>, R>;

} // namespace polynomial
} // namespace surrogate
} // namespace hyteg
