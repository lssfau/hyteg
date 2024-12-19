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

// helper function to compute binomial coefficients
template < uint_t n, uint_t k >
static constexpr uint_t binom = binom< n - 1, k - 1 > * n / k;
template < uint_t n >
static constexpr uint_t binom< n, 0 > = 1;
template < uint_t n >
static constexpr uint_t binom< n, n > = 1;

// dimension of polynomial space of degree Q over R^D
template < uint_t D, uint_t Q >
static constexpr uint_t dimP = binom< D + Q, D >;

// store x^i*y^j*z^k by exponents i,j,k compressed into 16 bit
struct Compression
{
   using basis_t = uint16_t;

   // return compressed basis function phi = x^i*y^j*z^k
   static constexpr inline basis_t compress( const basis_t i, const basis_t j, const basis_t k ) { return I * i | J * j | K * k; }

   // extract exponents i,j,k from compressed basis function phi
   static constexpr inline std::array< basis_t, 3 > expand( const basis_t phi ) { return { i( phi ), j( phi ), k( phi ) }; }
   // extract exponent i of x^i*y^j*z^k
   static constexpr inline basis_t i( const basis_t phi ) { return phi & basis_t( J - 1 ); }
   // extract exponent j of x^i*y^j*z^k
   static constexpr inline basis_t j( const basis_t phi ) { return ( phi & basis_t( K - 1 ) ) >> 5; }
   // extract exponent k of x^i*y^j*z^k
   static constexpr inline basis_t k( const basis_t phi ) { return phi >> 10; }

   // polynomial degree of this basis function
   static constexpr inline basis_t degree( const basis_t phi ) { return i( phi ) + j( phi ) + k( phi ); }

   // evaluate phi at x
   static constexpr inline double eval( const basis_t phi, const std::array< double, 3 >& x )
   {
      const auto [i, j, k] = expand( phi );

      double val = 1.0;

      for ( uint8_t p = 0; p < i; ++p )
         val *= x[0]; // compute x^i
      for ( uint8_t p = 0; p < j; ++p )
         val *= x[1]; // compute y^j
      for ( uint8_t p = 0; p < k; ++p )
         val *= x[2]; // compute z^k

      return val;
   }

   // bit masks for exponents
   static constexpr basis_t I = 1;
   static constexpr basis_t J = 1 << 5;
   static constexpr basis_t K = 1 << 10;
};

// monomial basis of polynomial space, given by the exponents of the monomials
template < uint8_t D, uint8_t Q >
struct Basis : public std::array< Compression::basis_t, dimP< D, Q > >
{
   using basis_t = Compression::basis_t;

   constexpr Basis()
   : std::array< basis_t, dimP< D, Q > >{}
   {
      // 1d basis
      for ( basis_t i = 0; i <= Q; ++i )
      {
         ( *this )[i] = Compression::compress( i, 0, 0 );
      }
      if constexpr ( D == 1 )
      {
         return;
      }
      // extend to 2d
      auto ij = dimP< 1, Q >;
      for ( basis_t i = 0; i < Q; ++i )
      {
         for ( basis_t j = 1; i + j <= Q; ++j )
         {
            ( *this )[ij++] = Compression::compress( i, j, 0 );
         }
      }
      if constexpr ( D == 2 )
      {
         return;
      }
      // extend to 3d
      auto ijk = dimP< 2, Q >;
      for ( basis_t i = 0; i < Q; ++i )
      {
         for ( basis_t k = 1; i + k <= Q; ++k )
         {
            ( *this )[ijk++] = Compression::compress( i, 0, k );
         }
      }
      for ( basis_t i = 0; i < Q; ++i )
      {
         for ( basis_t j = 1; i + j < Q; ++j )
         {
            for ( basis_t k = 1; i + j + k <= Q; ++k )
            {
               ( *this )[ijk++] = Compression::compress( i, j, k );
            }
         }
      }
   }

   // evaluate the n-th basis function at x
   constexpr inline double eval( uint_t n, const std::array< double, 3 >& x ) const
   {
      return Compression::eval( ( *this )[n], x );
   }

   // return the n-th basis function as array of exponents [i,j,k]
   constexpr inline auto operator()( uint_t n ) const { return Compression::expand( ( *this )[n] ); }

   // compute the degree of the n-th basis function
   constexpr inline basis_t degree( uint_t n ) const { return Compression::degree( ( *this )[n] ); }
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

// polynomials of degree Q over R^D
template < uint8_t D, uint8_t Q >
class Polynomial
{
 public:
   static constexpr uint8_t X = 0;
   static constexpr uint8_t Y = 1;
   static constexpr uint8_t Z = 2;

   static constexpr uint_t nc = dimP< D, Q >; // number of coefficients

   inline constexpr Polynomial()
   : c{}
   , restriction()
   {}

   inline constexpr Polynomial( const Eigen::Matrix< double, nc, 1, Eigen::ColMajor >& coeffs )
   : c{ coeffs }
   , restriction()
   {}

   // fix z coordinate s.th. only 2d polynomial must be evaluated
   void fix_z( const double z, bool use_for_y = false ) const
   {
      if constexpr ( D == 2 )
      {
         WALBERLA_ASSERT( use_for_y, "fix_z can only be used in 3d" );
      }
      else if constexpr ( D != 3 )
      {
         WALBERLA_ABORT( "fix_z can only be used in 3d" );
      }

      // z^k for k=0,...,Q
      std::array< double, Q + 1 > z_pow{ 1.0 };
      for ( uint_t k = 0; k < Q; ++k )
      {
         z_pow[k + 1] = z_pow[k] * z;
      }

      // first index where the 3d extension starts
      uint_t ijk = restriction.nc;
      // iterate over coefficients of 2d polynomial
      for ( uint_t ij = 0; ij < restriction.nc; ++ij )
      {
         const auto d = basis.degree( ij );

         restriction.c( ij ) = c( ij ); // k=0
         for ( uint_t k = 1; k <= Q - d; ++k )
         {
            restriction.c( ij ) += c[ijk++] * z_pow[k];
         }
      }
   }

   // fix y coordinate s.th. only 1d polynomial must be evaluated
   void fix_y( const double y ) const
   {
      if constexpr ( D == 2 )
      {
         return fix_z( y, true );
      }
      else if constexpr ( D == 3 )
      {
         return restriction.fix_y( y );
      }
      else
      {
         WALBERLA_ABORT( "fix_y can only be used in 2d and 3d" );
      }
   }

   // evaluate 1d polynomial using Horner's method

   /** @brief Evaluate the 1d polynomial at x using Horner's method.
    *
    * Usage: p.fix_z(z); p.fix_y(y); p.eval(x);
    *
    * @param x The x-coordinate.
    * @return p(x)
    */
   inline double eval( const double x ) const
   {
      if constexpr ( D == 1 )
      {
         auto val = c[nc - 1];
         for ( int i = int( nc - 2 ); i >= 0; --i )
         {
            val = val * x + c( i );
         }
         return val;
      }
      else
      {
         return restriction.eval( x );
      }
   }

   // evaluate polynomial by summing up basis functions, only use for debugging or testing
   double eval_naive( const std::array< double, 3 >& x ) const
   {
      double val = 0.0;
      for ( uint_t i = 0; i < nc; ++i )
      {
         val += c( i ) * basis.eval( i, x );
      }
      return val;
   }

   static constexpr auto basis = Basis< D, Q >();

 private:
   // coefficients
   Eigen::Matrix< double, nc, 1, Eigen::ColMajor > c;
   // restriction to lower dimension
   mutable Polynomial< D - 1, Q > restriction;
   // make all polynomials friends
   template < uint8_t d, uint8_t q >
   friend class Polynomial;
};

template < uint8_t Q >
class Polynomial< 0, Q >
{};

} // namespace polynomial
} // namespace surrogate
} // namespace hyteg
