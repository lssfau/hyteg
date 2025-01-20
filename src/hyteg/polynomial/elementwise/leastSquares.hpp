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

#include <core/Format.hpp>
#include <core/logging/Logging.h>
#include <hyteg/eigen/EigenWrapper.hpp>
#include <hyteg/indexing/Common.hpp>

#include "polynomial.hpp"

// todo change from double to template <typename Float>

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
 * @brief Computes the number of vertices in a volume element with given edge length
 *
 * @param d The dimension (2 for 2D, 3 for 3D).
 * @param n The number of vertices along an edge.
 * @return The number of vertices in the volume.
 */
static constexpr inline uint_t n_volume( uint_t d, uint_t n )
{
   return ( d == 2 ) ? tri( n ) : tet( n );
}

/**
 * @brief Computes the downsampled number of vertices in a volume element
 *
 * @param d The dimension (2 for 2D, 3 for 3D).
 * @param lvl The level of refinement.
 * @param downsampling The downsampling factor.
 * @return The number of vertices in the volume.
 */
static constexpr inline uint_t n_volume( uint_t d, uint_t lvl, uint_t downsampling )
{
   return n_volume( d, n_edge( lvl, downsampling ) );
}

} // namespace interpolation

/**
 * @class LeastSquares
 * @brief Class to compute polynomial approximations p ∈ P_q(T) of functions f:T→ℝ,
 *          where T is the reference macro element of dimension d (triangle or tetrahedron).
 *          The approximation is computed by solving a least squares problem using the Vandermonde matrix.
 *          The Sample points are a subset of the micro vertices corresponding to a given level.
 */
template < typename FLOAT >
class LeastSquares
{
 public:
   using Matrix = Eigen::Matrix< FLOAT, -1, -1, Eigen::RowMajor >;
   using Vector = Eigen::Matrix< FLOAT, -1, 1, Eigen::ColMajor >;

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
      inline Iterator( uint_t d, uint_t lvl, uint_t downsampling )
      : _n( idx_t( 0 ) )
      , _i( idx_t( 0 ) )
      , _j( idx_t( 0 ) )
      , _k( idx_t( 0 ) )
      , _ijk_max( idx_t( interpolation::n_edge( lvl, 1 ) ) )
      , _n_max( interpolation::n_volume( d, lvl, downsampling ) )
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
       * @return idx_t Vertex index in x direction.
       */
      inline constexpr idx_t i() const { return _i; }
      /**
       * @brief Getter for vertex index in y direction.
       *
       * @return idx_t Vertex index in y direction.
       */
      inline constexpr idx_t j() const { return _j; }
      /**
       * @brief Getter for vertex index in z direction.
       *
       * @return idx_t Vertex index in z direction.
       */
      inline constexpr idx_t k() const { return _k; }
      /**
       * @brief Getter for vertex index in x,y,z direction.
       *
       * @return Index
       */
      inline indexing::Index ijk() const { return { _i, _j, _k }; }
      /**
       * @brief Getter for index of sample point.
       *
       * @return uint_t Index of sample point = row of Vandermonde matrix.
       */
      inline constexpr idx_t operator()() const { return _n; }
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
      idx_t _n; // index of sample point (row in A)
      idx_t _i; // index of x-coordinate
      idx_t _j; // index of y-coordinate
      idx_t _k; // index of z-coordinate

      idx_t  _ijk_max; // number of microedges along an edge (only coincides with sample points along edge if downsampling=1)
      uint_t _n_max;   // number of sample points
      uint_t _stride;  // downsampling factor
   };

   // compute downsampling factor s.th. number of sample points exceeds number of coefficients
   constexpr uint_t compute_automatic_downsampling() const
   {
      uint_t downsampling = 3;
      for ( ; downsampling > 1; --downsampling )
      {
         if ( interpolation::n_volume( _dim, _lvl, downsampling ) > polynomial::dimP( _dim, _q ) )
         {
            break;
         }
      }
      return downsampling;
   }

   // load matrix from file
   template < typename MatrixType >
   static bool load_matrix( const std::string& filename, MatrixType& M )
   {
      std::ifstream file( filename, std::ios::binary );
      if ( file.is_open() )
      {
         file.read( reinterpret_cast< char* >( M.data() ), M.size() * sizeof( FLOAT ) );
         file.close();
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
      std::ofstream file( filename, std::ios::binary );
      if ( file.is_open() )
      {
         file.write( reinterpret_cast< const char* >( M.data() ), M.size() * sizeof( FLOAT ) );
         file.close();
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
      const polynomial::Basis           phi( _q );
      const polynomial::Domain< FLOAT > X( _lvl );

      auto it = samplingIterator();
      while ( it != it.end() )
      {
         auto x = X( it.ijk() );

         for ( uint_t col = 0; col < cols; ++col )
         {
            A( it(), col ) = phi[col].eval( x );
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
      std::string ext = walberla::format( "_d%d_q%d_l%d_s%d_f%d.dat", _dim, _q, _lvl, _downsampling, 8*sizeof( FLOAT ) );
      std::string d   = "/";
      return load_matrix( path + d + "A" + ext, A ) && load_matrix( path + d + "Uh" + ext, Uh ) &&
             load_matrix( path + d + "Si" + ext, Si ) && load_matrix( path + d + "V" + ext, V );
   }

   bool store_svd( const std::string& path ) const
   {
      std::string ext = walberla::format( "_d%d_q%d_l%d_s%d_f%d.dat", _dim, _q, _lvl, _downsampling, 8*sizeof( FLOAT ) );
      std::string d   = "/";
      return store_matrix( path + d + "A" + ext, A ) && store_matrix( path + d + "Uh" + ext, Uh ) &&
             store_matrix( path + d + "Si" + ext, Si ) && store_matrix( path + d + "V" + ext, V );
   }

   LeastSquares( uint_t             dim,
                 uint_t             degree,
                 uint_t             lvl,
                 uint_t             downsampling,
                 bool               use_precomputed,
                 const std::string& path_to_svd )
   : _dim( dim )
   , _q( degree )
   , _lvl( lvl )
   , _downsampling( ( downsampling == 0 ) ? compute_automatic_downsampling() : downsampling )
   , rows( interpolation::n_volume( _dim, _lvl, _downsampling ) )
   , cols( polynomial::dimP( _dim, _q ) )
   , A( rows, cols )
   , Uh( cols, rows )
   , Si( cols )
   , V( cols, cols )
   , b( rows )
   , c( cols )
   {
      WALBERLA_ASSERT( rows >= cols, "Not enough sample points for given polynomial degree!" )
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
    * @brief Compute the maximum polynomial degree that should be used for approximation
    *    when using sample points determined by `level` and `downsampling` factor.
    *    Using a higher degree leads to an overdetermined system and should be avoided!
    */
   static constexpr uint8_t max_degree( uint_t lvl, uint_t downsampling = 1 )
   {
      downsampling = ( downsampling == 0 ) ? uint_t( 1 ) : downsampling;
      return interpolation::n_edge( lvl, downsampling ) - 1;
   }

   /**
    * @brief Constructs a LeastSquares object by setting up the Vandermonde matrix and computing an SVD.
    *
    * @param dim The spatial dimension of the domain T
    * @param degree The degree q of the polynomial space P_q(T)
    * @param lvl The level of refinement determining the mapping T→ℝ
    * @param downsampling The downsampling factor determining the number of sample points:
    *          0: choose automatically, 1: use all vertices, n: only use every n-th vertex in each direction
    *          ⇒ downsampling reduces the number of sample points by a factor of n^dim
    */
   LeastSquares( uint_t dim, uint_t degree, uint_t lvl, uint_t downsampling = 1 )
   : LeastSquares( dim, degree, lvl, downsampling, false, "" )
   {}

   /**
    * @brief Constructs a LeastSquares object using a precomputed SVD.
    *
    * @param path_to_svd Path to the directory where the files are stored.
    * @param dim The spatial dimension of the domain T
    * @param degree The degree q of the polynomial space P_q(T)
    * @param lvl The level of refinement determining the mapping T→ℝ
    * @param downsampling The downsampling factor determining the number of sample points:
    *          0: choose automatically, 1: use all vertices, n: only use every n-th vertex in each direction
    *          ⇒ downsampling reduces the number of sample points by a factor of n^dim
    */
   LeastSquares( const std::string& path_to_svd, uint_t dim, uint_t degree, uint_t lvl, uint_t downsampling = 1 )
   : LeastSquares( dim, degree, lvl, downsampling, true, path_to_svd )
   {}

   /**
    * @brief Store Vandermonde matrix and SVD, such that they can be used for `LeastSquares(path, lvl, downsampling)`.
    *
    * @param path Path to the directory where the files should be stored.
    */
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
   Iterator samplingIterator() const { return Iterator( _dim, _lvl, _downsampling ); }

   /**
    * @brief Set n-th coefficient of right hand side vector. Should be used together with `it = samplingIterator()`.
    *
    * @param n Index of the coefficient, ie, `n = it()`
    * @param f_xyz f( it.i(), it.j(), it.k() ), where f is the function to be approximated.
    */
   void setRHS( const idx_t n, const FLOAT f_xyz ) { b( n ) = f_xyz; }

   /**
    * @brief Set right hand side vector in one step.
    *    !Only use if you know what you are doing!
    * @param rhs Should have been set up using `samplingIterator()`.
    */
   void setRHS( const Vector& rhs ) { b = rhs; }

   /**
    * @brief Solve the least squares problem.
    *     !Use after setting all right-hand side coefficients!
    * @return The coefficients of the polynomial p minimizing ||p(S) - f(S)||.
    */
   const Vector& solve()
   {
      c = V * ( Si.cwiseProduct( Uh * b ) );
      return c;
   }

   /**
    * @brief Get the coefficients of the polynomial.
    *       ! call after `solve` !
    *
    * @return The coefficients of the polynomial.
    */
   const Vector& solution() const { return c; }

   /**
    * @brief Compute the residual of the least squares problem.
    *
    * @return discrete L^2 norm of the residual; ||r|| = sqrt(∑_{n=1}^N r_n^2 / N)
    */
   const FLOAT residual() const { return ( A * c - b ).norm() / std::sqrt( rows ); }

 private:
   // spacial dimension
   const uint_t _dim;
   // polynomial degree
   const uint_t _q;
   // level of refinement
   const uint_t _lvl;
   // downsampling factor
   const uint_t _downsampling;

 public:
   // number of interpolation points
   const int rows;
   // dimension of polynomial space
   const int cols;

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