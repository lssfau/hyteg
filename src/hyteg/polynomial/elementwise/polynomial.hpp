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
#include <vector>

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
 private:
   using compression_t = int16_t;

   // bit masks for exponents
   static constexpr int           SHIFT_J = 5;
   static constexpr int           SHIFT_K = 10;
   static constexpr compression_t I       = 1;
   static constexpr compression_t J       = 1 << SHIFT_J;
   static constexpr compression_t K       = 1 << SHIFT_K;

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
   constexpr inline uint_t i() const { return ijk_compressed & compression_t( J - 1 ); }
   // extract exponent j of x^i*y^j*z^k
   constexpr inline uint_t j() const { return ( ijk_compressed & compression_t( K - 1 ) ) >> SHIFT_J; }
   // extract exponent k of x^i*y^j*z^k
   constexpr inline uint_t k() const { return ijk_compressed >> SHIFT_K; }
   // extract exponents i,j,k from compressed basis function
   constexpr inline std::array< uint_t, 3 > expand() const { return { i(), j(), k() }; }
   // total polynomial degree of this basis function
   constexpr inline uint_t degree() const { return i() + j() + k(); }

   // evaluate basis function at x
   constexpr inline double eval( const PointND< double, 3 >& x ) const
   {
      const auto [i, j, k] = expand();

      double val = 1.0;

      for ( uint_t p = 0; p < i; ++p )
         val *= x[0]; // compute x^i
      for ( uint_t p = 0; p < j; ++p )
         val *= x[1]; // compute y^j
      for ( uint_t p = 0; p < k; ++p )
         val *= x[2]; // compute z^k

      return val;
   }
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
   inline constexpr double operator()( idx_t i ) const { return scaling * double( i ) - 1.0; }
   // convert index [i,j,k] to coordinate x ∈ [-1,3]^3
   inline PointND< double, 3 > operator()( const indexing::Index& idx ) const
   {
      return { ( *this )( idx.x() ), ( *this )( idx.y() ), ( *this )( idx.z() ) };
   }

   // scaling factor to convert index i ∈ {0,..., 2^lvl - 1} to coordinate x ∈ [-1,3]
   const double scaling;
};

class Polynomial : public std::vector< double >
{
 public:
   static constexpr uint8_t X = 0;
   static constexpr uint8_t Y = 1;
   static constexpr uint8_t Z = 2;

   inline Polynomial( uint8_t d = 0, uint8_t q = 0 )
   : std::vector< double >( dimP( d, q ) )
   , _d( d )
   , _q( q )
   , _restriction( ( d > 1 ) ? std::make_unique< Polynomial >( d - 1, q ) : nullptr )
   , _basis( &basis_q[q] )
   {}

   inline Polynomial( uint8_t d, uint8_t q, const Eigen::Matrix< double, -1, 1, Eigen::ColMajor >& coeffs )
   : Polynomial( d, q )
   {
      set_coefficients( coeffs );
   }

   // copy constructor
   inline Polynomial( const Polynomial& other )
   : std::vector< double >( other )
   , _d( other._d )
   , _q( other._q )
   , _restriction( other._restriction ? std::make_unique< Polynomial >( *other._restriction ) : nullptr )
   , _basis( other._basis )
   {}
   inline Polynomial( Polynomial&& other ) = default;

   // copy assignment
   inline Polynomial& operator=( const Polynomial& other )
   {
      std::vector< double >::operator=( other );
      _d           = other._d;
      _q           = other._q;
      _restriction = other._restriction ? std::make_unique< Polynomial >( *other._restriction ) : nullptr;
      _basis       = other._basis;
      return *this;
   }
   inline Polynomial& operator=( Polynomial&& other ) = default;

   // get i-th coefficient w.r.t. monomial basis
   inline double& c( idx_t i ) { return ( *this )[uint_t( i )]; }
   // get i-th coefficient w.r.t. monomial basis
   inline const double& c( idx_t i ) const { return ( *this )[uint_t( i )]; }
   // get number of coefficients ( = dimension of polynomial space )
   inline idx_t n_coefficients() const { return idx_t( size() ); }

   // set coefficients
   inline void set_coefficients( const Eigen::Matrix< double, -1, 1, Eigen::ColMajor >& coeffs )
   {
      if ( coeffs.size() != n_coefficients() )
      {
         WALBERLA_ABORT( "Number of coefficients must match dimension of polynomial space" );
      }

      for ( idx_t i = 0; i < n_coefficients(); ++i )
      {
         c( i ) = coeffs( i );
      }
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
         auto px = c( n_coefficients() - 1 );
         for ( idx_t i = n_coefficients() - 2; i >= 0; --i )
         {
            px = px * x + c( i );
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
   double eval_naive( const PointND< double, 3 >& x ) const
   {
      double p_xyz = 0.0;
      for ( idx_t i = 0; i < n_coefficients(); ++i )
      {
         p_xyz += c( i ) * phi( i ).eval( x );
      }
      return p_xyz;
   }

   // get basis of polynomial space
   const Basis& basis() const { return *_basis; }
   // get i-th basis function
   const Monomial& phi( uint_t i ) const { return basis()[i]; }

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
      auto ijk = _restriction->n_coefficients();
      // iterate over coefficients of 2d polynomial
      for ( idx_t ij = 0; ij < _restriction->n_coefficients(); ++ij )
      {
         const auto max_k = _q - phi( ij ).degree();

         _restriction->c( ij ) = c( ij ); // k=0
         for ( uint_t k = 1; k <= max_k; ++k )
         {
            _restriction->c( ij ) += c( ijk++ ) * z_pow[k];
         }
      }
   }

   // dimension of domain
   uint8_t _d;
   // degree of polynomial
   uint8_t _q;
   // restriction to lower dimension
   std::unique_ptr< Polynomial > _restriction;
   // basis of the this polynomial's space
   const Basis* _basis;

   // basis functions for q = 0,1,...
   static const std::vector< Basis > basis_q;
};

// RxC matrix of polynomials
template < uint_t R, uint_t C = R >
using Matrix = Eigen::Matrix< Polynomial, R, C, Eigen::RowMajor >;
// using Matrix = std::array< std::array< Polynomial, C >, R >;

} // namespace polynomial
} // namespace surrogate
} // namespace hyteg
