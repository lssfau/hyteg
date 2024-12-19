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

#include <core/logging/Logging.h>
#include <hyteg/eigen/EigenWrapper.hpp>

#include "polynomial.hpp"

namespace hyteg {
namespace surrogate {

using walberla::uint_t;

namespace interpolation {

/**
 * @brief Computes the number of vertices in a triangle with edge length n
 *
 * @param n number of vertices along an edge
 * @return The triangular number of n.
 */
static constexpr inline uint_t tri( uint_t n )
{
   return n * ( n + 1 ) / 2;
}

/**
 * @brief Computes the number of vertices in a tetrahedron with edge length n
 *
 * @param n number of vertices along an edge
 * @return The tetrahedral number of n.
 */
static constexpr inline uint_t tet( uint_t n )
{
   return n * ( n + 1 ) * ( n + 2 ) / 6;
}

/**
 * @brief Computes the downsampled number of vertices along an edge on a given level
 *
 * @param lvl The level of refinement.
 * @param downsampling The downsampling factor.
 * @return The number of downsampled vertices along the edge.
 */
static constexpr inline uint_t n_edge( uint_t lvl, uint_t downsampling )
{
   return ( ( ( 1 << lvl ) - 1 ) / downsampling ) + 1; // ceil(n/d) with n = 2^lvl, d = downsampling
}

/**
 * @brief Computes the downsampled number of vertices in a volume element
 *
 * @tparam D The dimension (2 for 2D, 3 for 3D).
 * @param lvl The level of refinement.
 * @param downsampling The downsampling factor.
 * @return The number of vertices in the volume.
 */
template < uint8_t D >
static constexpr inline uint_t n_volume( uint_t n_edge )
{
   return ( D == 2 ) ? tri( n_edge ) : tet( n_edge );
}

/**
 * @brief Computes the number of vertices in a volume element with given edge length
 *
 * @tparam D The dimension (2 for 2D, 3 for 3D).
 * @param n_edge The number of vertices along an edge.
 * @return The number of vertices in the volume.
 */
template < uint8_t D >
static constexpr inline uint_t n_volume( uint_t lvl, uint_t downsampling )
{
   return n_volume< D >( n_edge( lvl, downsampling ) );
}

} // namespace interpolation

template < uint8_t D, uint8_t Q >
class LeastSquares
{
 private:
   /**
    * @class Iterator
    * @brief Iterator class to iterate over sample points
    */
   class Iterator
   {
    public:
      /**
       * @brief Constructs an Iterator object.
       *
       * @param lvl The level of refinement.
       * @param downsampling The downsampling factor.
       */
      inline Iterator( uint_t lvl, uint_t downsampling )
      : _n( 0 )
      , _i( 0 )
      , _j( 0 )
      , _k( 0 )
      , _ijk_max( interpolation::n_edge( lvl, 1 ) )
      , _n_max( interpolation::n_volume< D >( lvl, downsampling ) )
      , _stride( downsampling )
      {}

      /**
       * @brief Prefix increment operator.
       *
       * @return Iterator& Reference to the incremented iterator.
       */
      inline Iterator& operator++()
      {
         _i += _stride;
         if ( _i >= _ijk_max - _j - _k )
         {
            _i = 0;
            _j += _stride;
            if ( _j >= _ijk_max - _k )
            {
               _j = 0;
               _k += _stride;
            }
         }
         ++_n;
         return *this;
      }

      /**
       * @brief Getter for vertex index in x direction.
       *
       * @return uint_t Vertex index in x direction.
       */
      inline constexpr uint_t i() const { return _i; }
      /**
       * @brief Getter for vertex index in y direction.
       *
       * @return uint_t Vertex index in y direction.
       */
      inline constexpr uint_t j() const { return _j; }
      /**
       * @brief Getter for vertex index in z direction.
       *
       * @return uint_t Vertex index in z direction.
       */
      inline constexpr uint_t k() const { return _k; }
      /**
       * @brief Getter for index of sample point.
       *
       * @return uint_t Index of sample point = row of Vandermonde matrix.
       */
      inline constexpr uint_t operator()() const { return _n; }
      /**
       * @brief end of this iterator. Useage:
       * while ( it != it.end() ) { ... }
       *
       * @return uint_t end of iterator.
       */
      inline constexpr uint_t end() const { return _n_max; }
      inline constexpr bool   operator!=( const uint_t other ) const { return _n != other; }
      inline constexpr bool   operator==( const uint_t other ) const { return _n == other; }

    private:
      uint_t _n; // index of sample point (row in A)
      uint_t _i; // index of x-coordinate
      uint_t _j; // index of y-coordinate
      uint_t _k; // index of z-coordinate

      uint_t _ijk_max; // number of microedges along an edge (only coincides with sample points along edge if downsampling=1)
      uint_t _n_max;   // number of sample points
      uint_t _stride;  // downsampling factor
   };

   // compute downsampling factor s.th. number of sample points exceeds number of coefficients
   static constexpr uint_t compute_automatic_downsampling( uint_t lvl )
   {
      uint_t d = 3;
      for ( ; d > 1; --d )
      {
         if ( interpolation::n_volume< D >( lvl, d ) > cols )
         {
            break;
         }
      }
      return d;
   }

   // load matrix from file
   template < typename MatrixType >
   static bool load_matrix( const std::string& filename, MatrixType& M )
   {
      std::ofstream out( filename, std::ios::binary );
      if ( out.is_open() )
      {
         out.write( reinterpret_cast< const char* >( M.data() ), M.size() * sizeof( double ) );
         out.close();
         return true;
      }
      else
      {
         return false;
      }
   }

   // store matrix in file
   template < typename MatrixType >
   static bool store_matrix( const std::string& filename, const MatrixType& M )
   {
      std::ofstream out( filename, std::ios::binary );
      if ( out.is_open() )
      {
         out.write( reinterpret_cast< const char* >( M.data() ), M.size() * sizeof( double ) );
         out.close();
         return true;
      }
      else
      {
         return false;
      }
   }

   void compute_svd()
   {
      // Setup Vandermonde matrix
      constexpr polynomial::Basis< D, Q > phi;
      const polynomial::Coordinates       coords( _lvl );

      auto it = samplingIterator();
      while ( it != it.end() )
      {
         auto x = coords.x( it.i(), it.j(), it.k() );

         for ( uint_t col = 0; col < cols; ++col )
         {
            A( it(), col ) = phi.eval( col, x );
         }
         ++it;
      }

      // compute SVD
      Eigen::JacobiSVD< Matrix > svd( A, Eigen::ComputeThinU | Eigen::ComputeThinV );
      Uh = svd.matrixU().adjoint();
      Si = svd.singularValues().cwiseInverse();
      V  = svd.matrixV();
   }

   bool load_svd( const std::string& path )
   {
      std::string ext = walberla::format( "_%d_%d_%d_%d.dat", D, Q, _lvl, _downsampling );
      std::string d   = "/";
      return load_matrix( path + d + "A" + ext, A ) && load_matrix( path + d + "Uh" + ext, Uh ) &&
             load_matrix( path + d + "Si" + ext, Si ) && load_matrix( path + d + "V" + ext, V );
   }

   bool store_svd( const std::string& path ) const
   {
      std::string ext = walberla::format( "_%d_%d_%d_%d.dat", D, Q, _lvl, _downsampling );
      std::string d   = "/";
      return store_matrix( path + d + "A" + ext, A ) && store_matrix( path + d + "Uh" + ext, Uh ) &&
             store_matrix( path + d + "Si" + ext, Si ) && store_matrix( path + d + "V" + ext, V );
   }

   LeastSquares( uint_t lvl, uint_t downsampling, bool use_precomputed, const std::string& path_to_svd = "" )
   : _lvl( lvl )
   , _downsampling( ( downsampling == 0 ) ? compute_automatic_downsampling( _lvl ) : downsampling )
   , rows( interpolation::n_volume< D >( _lvl, _downsampling ) )
   , A( rows, cols )
   , Uh( cols, rows )
   , Si( cols )
   , V( cols, cols )
   , b( rows )
   , c( cols )
   {
      if ( use_precomputed )
      {
         if ( load_svd( path_to_svd ) )
         {
            return;
         }
         else
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "Could not load SVD. Compute SVD instead." );
         }
      }
      compute_svd();
   }

 public:
   /**
    * @brief Constructs a LeastSquares object.
    *
    * @param lvl The level of refinement.
    * @param downsampling The downsampling factor.
    * (0: choose automatically, 1: no downsampling, n: only use every n-th vertex in each direction)
    */
   LeastSquares( uint_t lvl, uint_t downsampling = 1 )
   : LeastSquares( lvl, downsampling, false, "" )
   {}

   /**
    * @brief Constructs a LeastSquares object using a precomputed SVD.
    *
    * @param path_to_svd Path to a the files containing the SVD of the Vandermonde matrix.
    * @param lvl The level of refinement.
    * @param downsampling The downsampling factor.
    * (0: choose automatically, 1: no downsampling, n: only use every n-th vertex in each direction)
    */
   LeastSquares( const std::string& path_to_svd, uint_t lvl, uint_t downsampling = 1 )
   : LeastSquares( lvl, downsampling, true, path_to_svd )
   {}

   void write_to_file( const std::string& path ) const
   {
      if ( !store_svd( path ) )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Could not open files for writing. SVD not stored" );
      }
   }

   /**
    * @brief Getter for sampling iterator.
    *
    * @return An iterator over all sampling points. Use to set right-hand side.
    */
   Iterator samplingIterator() const { return Iterator( _lvl, _downsampling ); }

   /**
    * @brief Set n-th coefficient of right hand side vector. Should be used together with it = samplingIterator().
    *
    * @param n Index of the coefficient, ie, n = it()
    * @param f_xyz f( it.i(), it.j(), it.k() ), where f is the function to be approximated.
    */
   void setRHS( const uint_t n, const double f_xyz ) { b( n ) = f_xyz; }

   /**
    * @brief Solve the least squares problem. Use after setting all right-hand side coefficients.
    *
    * @return The polynomial object with the computed coefficients.
    */
   const polynomial::Polynomial< D, Q > solve()
   {
      c = V * ( Si.cwiseProduct( Uh * b ) );
      return polynomial::Polynomial< D, Q >( c );
   }

   /**
    * @brief Compute the residual of the least squares problem.
    *
    * @return weighted residual norm ||r|| = sqrt(∑_{n=1}^N r_n^2 / N)
    */
   const double residual() const { return ( A * c - b ).norm() / std::sqrt( rows ); }

   using Matrix = Eigen::Matrix< double, -1, -1, Eigen::RowMajor >;
   using Vector = Eigen::Matrix< double, -1, 1, Eigen::ColMajor >;

 private:
   // level of refinement
   const uint_t _lvl;
   // downsampling factor
   const uint_t _downsampling;

 public:
   // number of interpolation points
   const int rows;
   // dimension of polynomial space
   static constexpr int cols = polynomial::dimP< D, Q >;

 private:
   // Vandermonde matrix
   Matrix A;
   // SVD of Vandermonde matrix
   Matrix Uh; // left singular vectors
   Vector Si; // inverse of singular values
   Matrix V;  // right singular vectors
   // rhs of least squares problem
   Vector b;
   // coefficients of the polynomial
   Vector c;
};

} // namespace surrogate

} // namespace hyteg