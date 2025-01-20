/*
 * Copyright (c) 2025 Benjamin Mann.
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
#include <array>
#include <hyteg/Levelinfo.hpp>
#include <hyteg/eigen/EigenWrapper.hpp>
#include <hyteg/indexing/Common.hpp>
#include <hyteg/volumedofspace/CellDoFIndexing.hpp>
#include <hyteg/volumedofspace/FaceDoFIndexing.hpp>

#include "leastSquares.hpp"
#include "polynomial.hpp"

// various containers to store elementwise data such as surrogate polynomials

namespace hyteg {
namespace surrogate {

/**
 * @class LocalMatrixLike
 * @brief T^RxC matrix where R and C are the dimensions of the local element matrices for a given local polynomial degree
 * @tparam T The type of the elements stored in the matrix (e.g. polynomial::Polynomial).
 * @tparam DIM The spacial dimension of the PDE domain
 * @tparam SRC_DEGREE Polynomial degree of the local source space (domain of A_loc).
 * @tparam DST_DEGREE Polynomial degree of the local destination space (image of A_loc). (defaults to SRC_DEGREE).
 */
template < typename T, uint_t DIM, uint_t SRC_DEGREE, uint_t DST_DEGREE >
class LocalMatrixLike
{
 public:
   static constexpr idx_t R = polynomial::dimP( DIM, DST_DEGREE );
   static constexpr idx_t C = polynomial::dimP( DIM, SRC_DEGREE );

   LocalMatrixLike()
   : _data{}
   {}

   // get (i,j) entry of this
   inline T& operator()( idx_t i, idx_t j ) { return _data[static_cast< uint_t >( i )][static_cast< uint_t >( j )]; }
   // get (i,j) entry of this
   inline const T& operator()( idx_t i, idx_t j ) const { return _data[static_cast< uint_t >( i )][static_cast< uint_t >( j )]; }

   inline constexpr idx_t rows() const { return R; }
   inline constexpr idx_t cols() const { return C; }

 private:
   std::array< std::array< T, C >, R > _data;
};

template < typename T, uint_t DIM >
class ElementTypeWiseData
{
 public:
   using ElementType = typename std::conditional< ( DIM == 2 ), facedof::FaceType, celldof::CellType >::type;

   constexpr ElementTypeWiseData()
   : _data{}
   {}

   inline T&       operator[]( const ElementType& eType ) { return _data[static_cast< uint_t >( eType )]; }
   inline const T& operator[]( const ElementType& eType ) const { return _data[static_cast< uint_t >( eType )]; }

 private:
   std::array< T, ( DIM == 2 ) ? 2 : 6 > _data;
};

/**
 * @class LocalMatrixMap
 * @brief Map of local element matrices to store all local matrices of a given macro element.
 * @tparam FLOAT The type of the elements stored in the matrix (e.g. double, real_t, ...)
 * @tparam DIM The spacial dimension of the PDE domain
 * @tparam SRC_DEGREE Polynomial degree of the local source space (domain of A_loc).
 * @tparam DST_DEGREE Polynomial degree of the local destination space (image of A_loc). (defaults to SRC_DEGREE).
 */
template < typename FLOAT, uint_t DIM, uint_t SRC_DEGREE, uint_t DST_DEGREE >
class LocalMatrixMap
{
   static constexpr idx_t R_loc = polynomial::dimP( DIM, DST_DEGREE );
   static constexpr idx_t C_loc = polynomial::dimP( DIM, SRC_DEGREE );

 public:
   using local_matrix_t = Matrix< FLOAT, R_loc, C_loc >;
   using ElementType    = typename std::conditional< ( DIM == 2 ), facedof::FaceType, celldof::CellType >::type;

   LocalMatrixMap()
   : _data{}
   {}

   LocalMatrixMap( uint_t level )
   : LocalMatrixMap()
   {
      set_level( level );
   }

   // adjust size of the container such that all data of given `level` fits.
   void set_level( uint_t level )
   {
      _level = level;
      _data.resize( required_size( level ) );
   }

   // retrieve local element matrix for the element at `Index` `micro` of `FaceType`/`CellType` `eType`.
   inline local_matrix_t& operator()( const ElementType& eType, const indexing::Index& micro )
   {
      if constexpr ( DIM == 2 )
      {
         return _data[facedof::macroface::index( _level, micro.x(), micro.y(), eType )];
      }
      else
      {
         return _data[celldof::macrocell::index( _level, micro.x(), micro.y(), micro.z(), eType )];
      }
   }
   // retrieve local element matrix for the element at `Index` `micro` of `FaceType`/`CellType` `eType`.
   inline const local_matrix_t& operator()( const ElementType& eType, const indexing::Index& micro ) const
   {
      if constexpr ( DIM == 2 )
      {
         return _data[facedof::macroface::index( _level, micro.x(), micro.y(), eType )];
      }
      else
      {
         return _data[celldof::macrocell::index( _level, micro.x(), micro.y(), micro.z(), eType )];
      }
   }

 private:
   static constexpr uint_t required_size( uint_t level )
   {
      if constexpr ( DIM == 2 )
      {
         return levelinfo::num_microfaces_per_face( level );
      }
      else
      {
         return levelinfo::num_microcells_per_cell( level );
      }
   }

   uint_t                                                                    _level;
   std::vector< local_matrix_t, Eigen::aligned_allocator< local_matrix_t > > _data;
};

/**
 * @brief container for level wise data on each primitive
 *       usage: ElementWiseData[primitiveID][lvl] = T
 */
template < typename T >
struct ElementWiseData : public std::map< PrimitiveID, std::vector< T > >
{
   ElementWiseData( const std::shared_ptr< PrimitiveStorage >& storage, uint_t l_max )
   {
      uint_t dim = ( storage->hasGlobalCells() ) ? 3 : 2;
      if ( dim == 2 )
      {
         for ( auto& id : storage->getFaceIDs() )
         {
            ( *this )[id] = std::vector< T >( l_max + 1 );
         }
      }
      else
      {
         for ( auto& id : storage->getCellIDs() )
         {
            ( *this )[id] = std::vector< T >( l_max + 1 );
         }
      }
   }
};

// container for RHS vectors for lsq-fit
template < typename FLOAT, uint_t DIM, uint_t SRC_DEGREE, uint_t DST_DEGREE >
using RHS_matrix = LocalMatrixLike< typename LeastSquares< FLOAT >::Vector, DIM, SRC_DEGREE, DST_DEGREE >;
// container for precomputed element matrices
template < typename FLOAT, uint_t DIM, uint_t SRC_DEGREE, uint_t DST_DEGREE >
using PrecomputedData = ElementWiseData< LocalMatrixMap< FLOAT, DIM, SRC_DEGREE, DST_DEGREE > >;
// container for surrogates
template < typename FLOAT, uint_t DIM, uint_t SRC_DEGREE, uint_t DST_DEGREE >
using SurrogateData = ElementWiseData<
    ElementTypeWiseData< LocalMatrixLike< polynomial::Polynomial< FLOAT >, DIM, SRC_DEGREE, DST_DEGREE >, DIM > >;

} // namespace surrogate
} // namespace hyteg