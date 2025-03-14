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

   // get i-th basis function
   const inline Monomial& phi( const idx_t i ) const { return ( *this )[uint_t( i )]; }
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
   inline constexpr FLOAT operator[]( idx_t i ) const { return scaling * FLOAT( i ) - FLOAT( 1.0 ); }

   // convert index [i,j,k] to coordinate x ∈ [-1,3]^3
   inline PointND< FLOAT, 3 > operator()( const indexing::Index& idx ) const
   {
      return { ( *this )[idx.x()], ( *this )[idx.y()], ( *this )[idx.z()] };
   }

   // scaling factor to convert index i ∈ {0,..., 2^lvl - 1} to coordinate x ∈ [-1,3]
   const FLOAT scaling;
};

template < typename FLOAT >
class Polynomial : public std::vector< FLOAT >
{
 public:
   static constexpr uint8_t X = 0;
   static constexpr uint8_t Y = 1;
   static constexpr uint8_t Z = 2;

   inline Polynomial( uint8_t d = 0, uint8_t q = 0 )
   : std::vector< FLOAT >( dimP( d, q ) )
   , _d( d )
   , _q( q )
   , _restriction( ( d > 1 ) ? std::make_unique< Polynomial >( d - 1, q ) : nullptr )
   {
      // initialize basis of P_q
      if ( q > 0 && _basis.count( q ) == 0 )
      {
         _basis.emplace( std::pair< uint8_t, Basis >{ q, Basis( q ) } );
      }
   }

   template < typename CoeffVector >
   inline Polynomial( uint8_t d, uint8_t q, const CoeffVector& coeffs )
   : Polynomial( d, q )
   {
      set_coefficients( coeffs );
   }

   // copy constructor
   inline Polynomial( const Polynomial& other )
   : std::vector< FLOAT >( other )
   , _d( other._d )
   , _q( other._q )
   , _restriction( other._restriction ? std::make_unique< Polynomial >( *other._restriction ) : nullptr )
   {}
   inline Polynomial( Polynomial&& other ) = default;

   // copy assignment
   inline Polynomial& operator=( const Polynomial& other )
   {
      std::vector< FLOAT >::operator=( other );
      _d           = other._d;
      _q           = other._q;
      _restriction = other._restriction ? std::make_unique< Polynomial >( *other._restriction ) : nullptr;
      return *this;
   }
   inline Polynomial& operator=( Polynomial&& other ) = default;

   // get i-th coefficient w.r.t. monomial basis
   inline FLOAT& c( idx_t i ) { return ( *this )[uint_t( i )]; }

   // get i-th coefficient w.r.t. monomial basis
   inline const FLOAT& c( idx_t i ) const { return ( *this )[uint_t( i )]; }

   // get number of coefficients ( = dimension of polynomial space )
   inline idx_t n_coefficients() const { return idx_t( this->size() ); }

   /**
    * @brief Sets the coefficients of the polynomial.
    *
    * @param coeffs A vector containing the coefficients.
    *          Must provide functions `size()` and `operator[]`
    */
   template < typename CoeffVector >
   inline void set_coefficients( const CoeffVector& coeffs )
   {
      if ( idx_t( coeffs.size() ) != n_coefficients() )
      {
         WALBERLA_ABORT( "Number of coefficients must match dimension of polynomial space" );
      }

      using coeff_idx_t = decltype( coeffs.size() );

      for ( idx_t i = 0; i < n_coefficients(); ++i )
      {
         c( i ) = static_cast< FLOAT >( coeffs[static_cast< coeff_idx_t >( i )] );
      }
   }

   // fix z coordinate s.th. only 2d polynomial must be evaluated
   void fix_z( const FLOAT z ) const
   {
      WALBERLA_ASSERT( _d == 3, "fix_z can only be used in 3d" );

      return fix_coord( z );
   }

   // fix y coordinate s.th. only 1d polynomial must be evaluated
   const Polynomial< FLOAT >& fix_y( const FLOAT y ) const
   {
      WALBERLA_ASSERT( _d == 2 || _d == 3, "fix_z can only be used in 2d and 3d" );
      if ( _d == 2 )
      {
         fix_coord( y );
         return *_restriction;
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
   inline FLOAT eval( const FLOAT x ) const
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
   FLOAT eval_naive( const PointND< FLOAT, 3 >& x ) const
   {
      FLOAT p_xyz = 0.0;
      for ( idx_t i = 0; i < n_coefficients(); ++i )
      {
         p_xyz += c( i ) * phi( i ).eval( x );
      }
      return p_xyz;
   }

   // get i-th basis function of P_q
   const inline Monomial& phi( const idx_t i ) const { return _basis[_q].phi( i ); }

 private:
   // fix coordinate s.th. only lower dimensional polynomial must be evaluated
   inline void fix_coord( const FLOAT z ) const
   {
      // z^k for k=0,...,q
      Eigen::Vector< FLOAT, -1 > z_pow( _q + 1 );
      z_pow[0] = 1.0;
      for ( int k = 1; k <= _q; ++k )
      {
         z_pow[k] = z_pow[k - 1] * z;
      }

      // first index where the 3d extension starts
      auto ijk = _restriction->n_coefficients();
      // iterate over coefficients of 2d polynomial
      for ( idx_t ij = 0; ij < _restriction->n_coefficients(); ++ij )
      {
         const auto max_k = _q - phi( ij ).degree();

         _restriction->c( ij ) = c( ij ); // k=0
         for ( int k = 1; k <= max_k; ++k )
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
   // basis functions of P_q for all required q
   static std::map< uint8_t, Basis > _basis;
};

template < typename FLOAT >
std::map< uint8_t, Basis > Polynomial< FLOAT >::_basis;

} // namespace polynomial
} // namespace surrogate
} // namespace hyteg
