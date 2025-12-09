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

#include "P1ElementwiseSurrogateOperator.hpp"

#include <hyteg/forms/P1RowSumForm.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp>
#include <hyteg/forms/form_hyteg_manual/SphericalElementFormMass.hpp>

#include "P1LocalOperations.hpp"

namespace hyteg {

template < class P1Form, uint8_t DEGREE, bool Symmetric >
P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::P1ElementwiseSurrogateOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel )
: P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >( storage, minLevel, maxLevel, P1Form() )
{}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::P1ElementwiseSurrogateOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel,
    const P1Form&                              form )
: Operator( storage, minLevel, maxLevel )
, form_( form )
, is_initialized_( false )
, lsq_( maxLevel + 1 )
, downsampling_( maxLevel + 1 )
, a_loc_2d_( storage, 0, std::min( maxLevel, min_lvl_for_surrogate - 1u ) )
, a_loc_3d_( storage, 0, std::min( maxLevel, min_lvl_for_surrogate - 1u ) )
, surrogate_2d_( storage, min_lvl_for_surrogate, maxLevel )
, surrogate_3d_( storage, min_lvl_for_surrogate, maxLevel )
, surrogate_cube_2d_( storage, min_lvl_for_surrogate, maxLevel )
, surrogate_cube_3d_( storage, min_lvl_for_surrogate, maxLevel )
{}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::init( size_t             downsampling,
                                                                        const std::string& path_to_svd,
                                                                        bool               needsInverseDiagEntries )
{
   uint_t dim = ( storage_->hasGlobalCells() ) ? 3 : 2;

   // precompute and store local stiffness matrices for level 1-3
   for ( uint_t level = 0; level < min_lvl_for_surrogate && level <= maxLevel_; ++level )
   {
      if ( dim == 2 )
      {
         precompute_local_stiffness_2d( level );
      }
      else
      {
         precompute_local_stiffness_3d( level );
      }
   }

   // approximate local stiffness matrices for level 4+ by polynomials
   for ( uint_t level = min_lvl_for_surrogate; level <= maxLevel_; ++level )
   {
      // initialize least squares approximation
      if ( lsq_[level] == nullptr || downsampling_ != downsampling )
      {
         /* In 2D, there is no 'blue' face for x=n.
            In 3D, no 'white-down' for x=n-1 and only 'white-up' for x=n.
            Therefore, we don't use these x-values as sample points.
         */
         auto offset = dim - 1;
         if ( path_to_svd == "" )
         {
            lsq_[level] = std::make_shared< LSQ >( dim, DEGREE, level, downsampling, offset );
         }
         else
         {
            lsq_[level] = std::make_shared< LSQ >( path_to_svd, dim, DEGREE, level, downsampling, offset );
         }
         downsampling_ = downsampling;
      }

      if ( dim == 2 )
      {
         compute_local_surrogates_2d( level );
      }
      else
      {
         compute_local_surrogates_3d( level );
      }
   }

   if ( needsInverseDiagEntries )
   {
      computeInverseDiagonalOperatorValues();
   }

   is_initialized_ = true;
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::store_svd( const std::string& path_to_svd )
{
   for ( uint_t level = min_lvl_for_surrogate; level <= maxLevel_; ++level )
   {
      if ( lsq_[level] != nullptr )
      {
         lsq_[level]->write_to_file( path_to_svd );
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::precompute_local_stiffness_2d( uint_t level )
{
   for ( const auto& [id, face] : storage_->getFaces() )
   {
      auto& a_loc = a_loc_2d_[id][level];

      a_loc.set_level( level );

      for ( const auto& fType : facedof::allFaceTypes )
      {
         for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
         {
            auto& elMat = a_loc( fType, micro );
            elMat.setZero();
            p1::assembleLocalElementMatrix2D( *face, level, micro, fType, form_, elMat );
         }
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::precompute_local_stiffness_3d( uint_t level )
{
   for ( const auto& [id, cell] : storage_->getCells() )
   {
      auto& a_loc = a_loc_3d_[id][level];
      a_loc.set_level( level );

      for ( const auto& cType : celldof::allCellTypes )
      {
         for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
         {
            auto& elMat = a_loc( cType, micro );
            elMat.setZero();
            p1::assembleLocalElementMatrix3D( *cell, level, micro, cType, form_, elMat );
         }
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::compute_local_surrogates_2d( uint_t level )
{
   auto& lsq = *lsq_[level];
   // initialize rhs vectors for lsq
   RHS_matrix< 2 > rhs;
   for ( idx_t i = 0; i < rhs.rows(); ++i )
   {
      for ( idx_t j = 0; j < rhs.cols(); ++j )
      {
         rhs( i, j ).resize( lsq.rows );
      }
   }

   for ( auto& [id, face] : storage_->getFaces() )
   {
      for ( const auto& fType : facedof::allFaceTypes )
      {
         // set up rhs vectors for each entry of the local stiffness matrix
         auto it = lsq.samplingIterator();
         while ( it != it.end() )
         {
            Matrix3r elMat( Matrix3r::Zero() );
            p1::assembleLocalElementMatrix2D( *face, level, it.ijk(), fType, form_, elMat );
            for ( idx_t i = 0; i < elMat.rows(); ++i )
            {
               for ( idx_t j = 0; j < elMat.cols(); ++j )
               {
                  rhs( i, j )[it()] = elMat( i, j );
               }
            }
            ++it;
         }
         // fit polynomials for each entry of the local stiffness matrix
         auto& surrogate = surrogate_2d_[id][level][fType];
         for ( idx_t i = 0; i < rhs.rows(); ++i )
         {
            for ( idx_t j = 0; j < rhs.cols(); ++j )
            {
               // apply least squares fit
               lsq.setRHS( rhs( i, j ) );
               auto& coeffs      = lsq.solve();
               surrogate( i, j ) = Poly< 2 >( coeffs );
            }
         }
      }

      // todo: directly approximate cube polys
      // construct cube polynomials
      auto& cube = surrogate_cube_2d_[id][level];
      auto& gray = surrogate_2d_[id][level][facedof::FaceType::GRAY];
      auto& blue = surrogate_2d_[id][level][facedof::FaceType::BLUE];
      for ( uint_t k = 0; k < cube( 0, 0 ).size(); ++k )
      {
         cube( 0, 0 )[k] = gray( 0, 0 )[k];
         cube( 0, 1 )[k] = gray( 0, 1 )[k];
         cube( 0, 2 )[k] = gray( 0, 2 )[k];
         // cube(0,3)[k] = 0;

         cube( 1, 0 )[k] = gray( 1, 0 )[k];
         cube( 1, 1 )[k] = gray( 1, 1 )[k] + blue( 0, 0 )[k];
         cube( 1, 2 )[k] = gray( 1, 2 )[k] + blue( 0, 2 )[k];
         cube( 1, 3 )[k] = /*             */ blue( 0, 1 )[k];

         cube( 2, 0 )[k] = gray( 2, 0 )[k];
         cube( 2, 1 )[k] = gray( 2, 1 )[k] + blue( 2, 0 )[k];
         cube( 2, 2 )[k] = gray( 2, 2 )[k] + blue( 2, 2 )[k];
         cube( 2, 3 )[k] = /*             */ blue( 2, 1 )[k];

         // cube( 3, 0 )[k] = 0;
         cube( 3, 1 )[k] = /*             */ blue( 1, 0 )[k];
         cube( 3, 2 )[k] = /*             */ blue( 1, 2 )[k];
         cube( 3, 3 )[k] = /*             */ blue( 1, 1 )[k];
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::compute_local_surrogates_3d( uint_t level )
{
   // todo: exploit symmetry
   // todo: remove computation of cType-wise polynomials

   auto& lsq = *lsq_[level];
   // initialize rhs vectors for lsq
   RHS_matrix< 3 > rhs;
   for ( idx_t i = 0; i < rhs.rows(); ++i )
   {
      for ( idx_t j = 0; j < rhs.cols(); ++j )
      {
         rhs( i, j ).resize( lsq.rows );
      }
   }
   RHS_cube< 3 > rhs_cube;
   for ( idx_t i = 0; i < rhs_cube.rows(); ++i )
   {
      for ( idx_t j = 0; j < rhs_cube.cols(); ++j )
      {
         rhs_cube( i, j ).resize( lsq.rows );
      }
   }

   for ( auto& [id, cell] : storage_->getCells() )
   {
      for ( const auto& cType : celldof::allCellTypes )
      {
         // set up rhs vectors for each entry of the local stiffness matrix
         auto it = lsq.samplingIterator();
         while ( it != it.end() )
         {
            Matrix4r elMat( Matrix4r::Zero() );
            p1::assembleLocalElementMatrix3D( *cell, level, it.ijk(), cType, form_, elMat );
            for ( idx_t i = 0; i < rhs.rows(); ++i )
            {
               for ( idx_t j = 0; j < rhs.cols(); ++j )
               {
                  rhs( i, j )[it()] = elMat( i, j );
               }
            }
            ++it;
         }
         // fit polynomials for each entry of the local stiffness matrix
         auto& surrogate = surrogate_3d_[id][level][cType];
         for ( idx_t i = 0; i < rhs.rows(); ++i )
         {
            for ( idx_t j = 0; j < rhs.cols(); ++j )
            {
               // apply least squares fit
               lsq.setRHS( rhs( i, j ) );
               auto& coeffs      = lsq.solve();
               surrogate( i, j ) = Poly< 3 >( coeffs );
            }
         }
      }

      // set up rhs vectors for each entry of the local cube matrix
      auto it = lsq.samplingIterator();
      while ( it != it.end() )
      {
         Matrix4r wu( Matrix4r::Zero() );
         Matrix4r bu( Matrix4r::Zero() );
         Matrix4r gu( Matrix4r::Zero() );
         Matrix4r wd( Matrix4r::Zero() );
         Matrix4r bd( Matrix4r::Zero() );
         Matrix4r gd( Matrix4r::Zero() );

         p1::assembleLocalElementMatrix3D( *cell, level, it.ijk(), celldof::CellType::WHITE_UP, form_, wu );
         p1::assembleLocalElementMatrix3D( *cell, level, it.ijk(), celldof::CellType::BLUE_UP, form_, bu );
         p1::assembleLocalElementMatrix3D( *cell, level, it.ijk(), celldof::CellType::GREEN_UP, form_, gu );
         p1::assembleLocalElementMatrix3D( *cell, level, it.ijk(), celldof::CellType::WHITE_DOWN, form_, wd );
         p1::assembleLocalElementMatrix3D( *cell, level, it.ijk(), celldof::CellType::BLUE_DOWN, form_, bd );
         p1::assembleLocalElementMatrix3D( *cell, level, it.ijk(), celldof::CellType::GREEN_DOWN, form_, gd );

         const auto k = it();

         rhs_cube( 0, 0 )[k] = wu( 0, 0 ) /*        */ /*        */ /*        */ /*        */ /*        */;
         rhs_cube( 0, 1 )[k] = wu( 0, 1 ) /*        */ /*        */ /*        */ /*        */ /*        */;
         rhs_cube( 0, 2 )[k] = wu( 0, 2 ) /*        */ /*        */ /*        */ /*        */ /*        */;
         // rhs_cube(0,3)[k] = 0;
         rhs_cube( 0, 4 )[k] = wu( 0, 3 ) /*        */ /*        */ /*        */ /*        */ /*        */;
         // rhs_cube(0,5)[k] = 0;
         // rhs_cube(0,6)[k] = 0;
         // rhs_cube(0,7)[k] = 0;

         rhs_cube( 1, 0 )[k] = wu( 1, 0 ) /*        */ /*        */ /*        */ /*        */ /*        */;
         rhs_cube( 1, 1 )[k] = wu( 1, 1 ) + bu( 0, 0 ) + gu( 0, 0 ) /*        */ /*        */ /*        */;
         rhs_cube( 1, 2 )[k] = wu( 1, 2 ) + bu( 0, 2 ) + gu( 0, 1 ) /*        */ /*        */ /*        */;
         rhs_cube( 1, 3 )[k] = /*        */ bu( 0, 1 ) /*        */ /*        */ /*        */ /*        */;
         rhs_cube( 1, 4 )[k] = wu( 1, 3 ) /*        */ + gu( 0, 3 ) /*        */ /*        */ /*        */;
         rhs_cube( 1, 5 )[k] = /*        */ bu( 0, 3 ) + gu( 0, 2 ) /*        */ /*        */ /*        */;
         // rhs_cube(1,6)[k] = 0;
         // rhs_cube(1,7)[k] = 0;

         rhs_cube( 2, 0 )[k] = wu( 2, 0 ) /*        */ /*        */ /*        */ /*        */ /*        */;
         rhs_cube( 2, 1 )[k] = wu( 2, 1 ) + bu( 2, 0 ) + gu( 1, 0 ) /*        */ /*        */ /*        */;
         rhs_cube( 2, 2 )[k] = wu( 2, 2 ) + bu( 2, 2 ) + gu( 1, 1 ) /*        */ + bd( 3, 3 ) + gd( 0, 0 );
         rhs_cube( 2, 3 )[k] = /*        */ bu( 2, 1 ) /*        */ /*        */ /*        */ + gd( 0, 1 );
         rhs_cube( 2, 4 )[k] = wu( 2, 3 ) /*        */ + gu( 1, 3 ) /*        */ + bd( 3, 2 ) /*        */;
         rhs_cube( 2, 5 )[k] = /*        */ bu( 2, 3 ) + gu( 1, 2 ) /*        */ + bd( 3, 0 ) + gd( 0, 2 );
         rhs_cube( 2, 6 )[k] = /*        */ /*        */ /*        */ /*        */ bd( 3, 1 ) + gd( 0, 3 );
         // rhs_cube(2,7)[k] = 0;

         // rhs_cube(3,0)[k] = 0;
         rhs_cube( 3, 1 )[k] = /*        */ bu( 1, 0 ) /*        */ /*        */ /*        */ /*        */;
         rhs_cube( 3, 2 )[k] = /*        */ bu( 1, 2 ) /*        */ /*        */ /*        */ + gd( 1, 0 );
         rhs_cube( 3, 3 )[k] = /*        */ bu( 1, 1 ) /*        */ + wd( 0, 0 ) /*        */ + gd( 1, 1 );
         // rhs_cube(3,4)[k] = 0;
         rhs_cube( 3, 5 )[k] = /*        */ bu( 1, 3 ) /*        */ + wd( 0, 3 ) /*        */ + gd( 1, 2 );
         rhs_cube( 3, 6 )[k] = /*        */ /*        */ /*        */ wd( 0, 2 ) /*        */ + gd( 1, 3 );
         rhs_cube( 3, 7 )[k] = /*        */ /*        */ /*        */ wd( 0, 1 ) /*        */ /*        */;

         rhs_cube( 4, 0 )[k] = wu( 3, 0 ) /*        */ /*        */ /*        */ /*        */ /*        */;
         rhs_cube( 4, 1 )[k] = wu( 3, 1 ) /*        */ + gu( 3, 0 ) /*        */ /*        */ /*        */;
         rhs_cube( 4, 2 )[k] = wu( 3, 2 ) /*        */ + gu( 3, 1 ) /*        */ + bd( 2, 3 ) /*        */;
         // rhs_cube(4,3)[k] = 0;
         rhs_cube( 4, 4 )[k] = wu( 3, 3 ) /*        */ + gu( 3, 3 ) /*        */ + bd( 2, 2 ) /*        */;
         rhs_cube( 4, 5 )[k] = /*        */ /*        */ gu( 3, 2 ) /*        */ + bd( 2, 0 ) /*        */;
         rhs_cube( 4, 6 )[k] = /*        */ /*        */ /*        */ /*        */ bd( 2, 1 ) /*        */;
         // rhs_cube(4,7)[k] = 0;

         // rhs_cube(5,0)[k] = 0;
         rhs_cube( 5, 1 )[k] = /*        */ bu( 3, 0 ) + gu( 2, 0 ) /*        */ /*        */ /*        */;
         rhs_cube( 5, 2 )[k] = /*        */ bu( 3, 2 ) + gu( 2, 1 ) /*        */ + bd( 0, 3 ) + gd( 2, 0 );
         rhs_cube( 5, 3 )[k] = /*        */ bu( 3, 1 ) /*        */ + wd( 3, 0 ) /*        */ + gd( 2, 1 );
         rhs_cube( 5, 4 )[k] = /*        */ /*        */ gu( 2, 3 ) /*        */ + bd( 0, 2 ) /*        */;
         rhs_cube( 5, 5 )[k] = /*        */ bu( 3, 3 ) + gu( 2, 2 ) + wd( 3, 3 ) + bd( 0, 0 ) + gd( 2, 2 );
         rhs_cube( 5, 6 )[k] = /*        */ /*        */ /*        */ wd( 3, 2 ) + bd( 0, 1 ) + gd( 2, 3 );
         rhs_cube( 5, 7 )[k] = /*        */ /*        */ /*        */ wd( 3, 1 ) /*        */ /*        */;

         // rhs_cube(6,0)[k] = 0;
         // rhs_cube(6,1)[k] = 0;
         rhs_cube( 6, 2 )[k] = /*        */ /*        */ /*        */ /*        */ bd( 1, 3 ) + gd( 3, 0 );
         rhs_cube( 6, 3 )[k] = /*        */ /*        */ /*        */ wd( 2, 0 ) /*        */ + gd( 3, 1 );
         rhs_cube( 6, 4 )[k] = /*        */ /*        */ /*        */ /*        */ bd( 1, 2 ) /*        */;
         rhs_cube( 6, 5 )[k] = /*        */ /*        */ /*        */ wd( 2, 3 ) + bd( 1, 0 ) + gd( 3, 2 );
         rhs_cube( 6, 6 )[k] = /*        */ /*        */ /*        */ wd( 2, 2 ) + bd( 1, 1 ) + gd( 3, 3 );
         rhs_cube( 6, 7 )[k] = /*        */ /*        */ /*        */ wd( 2, 1 ) /*        */ /*        */;

         // rhs_cube(7,0)[k] = 0;
         // rhs_cube(7,1)[k] = 0;
         // rhs_cube(7,2)[k] = 0;
         rhs_cube( 7, 3 )[k] = /*        */ /*        */ /*        */ wd( 1, 0 ) /*        */ /*        */;
         // rhs_cube(7,4)[k] = 0;
         rhs_cube( 7, 5 )[k] = /*        */ /*        */ /*        */ wd( 1, 3 ) /*        */ /*        */;
         rhs_cube( 7, 6 )[k] = /*        */ /*        */ /*        */ wd( 1, 2 ) /*        */ /*        */;
         rhs_cube( 7, 7 )[k] = /*        */ /*        */ /*        */ wd( 1, 1 ) /*        */ /*        */;

         ++it;
      }

      // these components always sum up to zero
      auto is_zero = []( idx_t i, idx_t j ) {
         if ( i == 0 )
         {
            return ( j == 3 || j == 5 || j == 6 || j == 7 );
         }
         if ( i == 1 )
         {
            return ( j == 6 || j == 7 );
         }
         if ( i == 2 )
         {
            return ( j == 7 );
         }
         if ( i == 3 )
         {
            return ( j == 0 || j == 4 );
         }
         if ( i == 4 )
         {
            return ( j == 3 || j == 7 );
         }
         if ( i == 5 )
         {
            return ( j == 0 );
         }
         if ( i == 6 )
         {
            return ( j == 0 || j == 1 );
         }
         if ( i == 7 )
         {
            return ( j == 0 || j == 1 || j == 2 || j == 4 );
         }
         return false;
      };

      // fit polynomials for each entry of the local cube matrix
      auto& cube = surrogate_cube_3d_[id][level];

      for ( idx_t i = 0; i < rhs_cube.rows(); ++i )
      {
         for ( idx_t j = 0; j < rhs_cube.cols(); ++j )
         {
            if ( !is_zero( i, j ) )
            {
               // apply least squares fit
               lsq.setRHS( rhs_cube( i, j ) );
               auto& coeffs = lsq.solve();
               cube( i, j ) = Poly< 3 >( coeffs );
            }
         }
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::apply( const P1Function< real_t >& src,
                                                                         const P1Function< real_t >& dst,
                                                                         size_t                      level,
                                                                         DoFType                     flag,
                                                                         UpdateType                  updateType ) const
{
   return gemv( real_c( 1 ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::gemv( const real_t&               alpha,
                                                                        const P1Function< real_t >& src,
                                                                        const real_t&               beta,
                                                                        const P1Function< real_t >& dst,
                                                                        size_t                      level,
                                                                        DoFType                     flag ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   this->startTiming( "gemv" );

   // Make sure that halos are up-to-date
   if ( this->storage_->hasGlobalCells() )
   {
      // Note that the order of communication is important, since the face -> cell communication may overwrite
      // parts of the halos that carry the macro-vertex and macro-edge unknowns.
      src.communicate< Face, Cell >( level );
      src.communicate< Edge, Cell >( level );
      src.communicate< Vertex, Cell >( level );
   }
   else
   {
      communication::syncFunctionBetweenPrimitives( src, level );
   }

   // Formerly updateType == Replace
   const bool betaIsZero = std::fpclassify( beta ) == FP_ZERO;
   // Formerly updateType == Add
   const bool betaIsOne = std::fpclassify( beta - real_c( 1.0 ) ) == FP_ZERO;

   if ( betaIsZero )
   {
      // We need to zero the destination array (including halos).
      // However, we must not zero out anything that is not flagged with the specified BCs.
      // Therefore we first zero out everything that flagged, and then, later,
      // the halos of the highest dim primitives.
      dst.interpolate( real_c( 0 ), level, flag );
   }
   else if ( !betaIsOne )
   {
      dst.assign( { beta }, { dst }, level, flag );
   }

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< real_t >, Cell > dstVertexDoFIdx = dst.getCellDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcVertexDoFIdx = src.getCellDataID();

         real_t* srcVertexData = cell.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = cell.getData( dstVertexDoFIdx )->getPointer( level );

         // Zero out dst halos only
         //
         // This is also necessary when using update type == Add.
         // During additive comm we then skip zeroing the data on the lower-dim primitives.
         for ( const auto& idx : vertexdof::macrocell::Iterator( level ) )
         {
            if ( !vertexdof::macrocell::isOnCellFace( idx, level ).empty() )
            {
               auto arrayIdx           = vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() );
               dstVertexData[arrayIdx] = real_c( 0 );
            }
         }

         // local mat-vec
         apply_3d( cell, level, srcVertexData, dstVertexData, alpha );
      }
      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      dst.communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      dst.communicateAdditively< Cell, Vertex >( level, DoFType::All ^ flag, *storage_, betaIsZero );
   }

   else
   {
      // we only perform computations on face primitives
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< real_t >, Face > dstVertexDoFIdx = dst.getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         // Zero out dst halos only
         //
         // This is also necessary when using update type == Add.
         // During additive comm we then skip zeroing the data on the lower-dim primitives.
         for ( const auto& idx : vertexdof::macroface::Iterator( level ) )
         {
            if ( vertexdof::macroface::isVertexOnBoundary( level, idx ) )
            {
               auto arrayIdx           = vertexdof::macroface::index( level, idx.x(), idx.y() );
               dstVertexData[arrayIdx] = real_c( 0 );
            }
         }

         // local mat-vec

         apply_2d( face, level, srcVertexData, dstVertexData, alpha );
      }
      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      dst.communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage_, betaIsZero );
   }

   this->stopTiming( "gemv" );
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::smooth_jac( const P1Function< real_t >& dst,
                                                                              const P1Function< real_t >& rhs,
                                                                              const P1Function< real_t >& src,
                                                                              real_t                      omega,
                                                                              size_t                      level,
                                                                              DoFType                     flag ) const
{
   this->startTiming( "smooth_jac" );

   // compute the current residual
   this->apply( src, dst, level, flag );
   dst.assign( { real_c( 1 ), real_c( -1 ) }, { rhs, dst }, level, flag );

   // perform Jacobi update step
   dst.multElementwise( { *getInverseDiagonalValues(), dst }, level, flag );
   dst.assign( { 1.0, omega }, { src, dst }, level, flag );

   this->stopTiming( "smooth_jac" );
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::apply_2d( const Face&         face,
                                                                            const uint_t        level,
                                                                            const real_t* const srcVertexData,
                                                                            real_t* const       dstVertexData,
                                                                            const real_t&       alpha ) const
{
   auto& id = face.getID();

   if ( level < min_lvl_for_surrogate )
   { // use precomputed local matrices
      auto& a_loc = a_loc_2d_.at( id )[level];
      for ( const auto& fType : facedof::allFaceTypes )
      {
         for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
         {
            p1::localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, dstVertexData, a_loc( fType, micro ), alpha );
         }
      }
   }
   else
   { // use surrogate polynomials to approximate local matrices
      // todo: use cube poly
      for ( const auto& fType : facedof::allFaceTypes )
      {
         auto& surrogate = surrogate_2d_.at( id )[level][fType];
         // domain of surrogates
         PolyDomain X( level );
         // local stiffness matrix
         Matrix3r elMat( Matrix3r::Zero() );
         // iterate over all micro-faces
         indexing::Index micro;
         const auto      row_end = idx_t( facedof::macroface::numFacesPerRowByType( level, fType ) );
         for ( micro.y() = 0; micro.y() < row_end; micro.y()++ )
         {
            // restrict to 1d polynomial
            auto y = X[micro.y()];
            for ( idx_t i = 0; i < surrogate.rows(); ++i )
            {
               for ( idx_t j = 0; j < surrogate.cols(); ++j )
               {
                  surrogate( i, j ).fix_y( y );
               }
            }

            for ( micro.x() = 0; micro.x() < row_end - micro.y(); micro.x()++ )
            {
               // evaluate p(x) using Horner's method
               auto x = X[micro.x()];
               for ( idx_t i = 0; i < surrogate.rows(); ++i )
               {
                  for ( idx_t j = 0; j < surrogate.cols(); ++j )
                  {
                     elMat( i, j ) = surrogate( i, j ).eval( x );
                  }
               }

               // local matvec
               p1::localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, dstVertexData, elMat, alpha );
            }
         }
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::apply_3d( const Cell&         cell,
                                                                            const uint_t        level,
                                                                            const real_t* const srcVertexData,
                                                                            real_t* const       dstVertexData,
                                                                            const real_t&       alpha ) const
{
   auto& id = cell.getID();

   if ( level < min_lvl_for_surrogate )
   { // use precomputed local matrices
      auto& a_loc = a_loc_3d_.at( id )[level];

      for ( const auto& cType : celldof::allCellTypes )
      {
         for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
         {
            p1::localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, dstVertexData, a_loc( cType, micro ), alpha );
         }
      }
   }
   else
   {
      if constexpr ( true ) //! cubes loop
      {
         auto& surrogate = surrogate_cube_3d_.at( id )[level];
         // monomial basis
         constexpr surrogate::polynomial::Basis< DEGREE > phi;
         // dimension of Pq
         constexpr auto dimP = surrogate::polynomial::dimP( 1, DEGREE );
         // domain
         const PolyDomain X( level );
         // restriction to 2d
         auto& pp00 = surrogate( 0, 0 ).get_2d_restriction();
         auto& pp01 = surrogate( 0, 1 ).get_2d_restriction();
         auto& pp02 = surrogate( 0, 2 ).get_2d_restriction();
         auto& pp03 = surrogate( 0, 3 ).get_2d_restriction();
         auto& pp04 = surrogate( 0, 4 ).get_2d_restriction();
         auto& pp05 = surrogate( 0, 5 ).get_2d_restriction();
         auto& pp06 = surrogate( 0, 6 ).get_2d_restriction();
         auto& pp07 = surrogate( 0, 7 ).get_2d_restriction();
         auto& pp10 = surrogate( 1, 0 ).get_2d_restriction();
         auto& pp11 = surrogate( 1, 1 ).get_2d_restriction();
         auto& pp12 = surrogate( 1, 2 ).get_2d_restriction();
         auto& pp13 = surrogate( 1, 3 ).get_2d_restriction();
         auto& pp14 = surrogate( 1, 4 ).get_2d_restriction();
         auto& pp15 = surrogate( 1, 5 ).get_2d_restriction();
         auto& pp16 = surrogate( 1, 6 ).get_2d_restriction();
         auto& pp17 = surrogate( 1, 7 ).get_2d_restriction();
         auto& pp20 = surrogate( 2, 0 ).get_2d_restriction();
         auto& pp21 = surrogate( 2, 1 ).get_2d_restriction();
         auto& pp22 = surrogate( 2, 2 ).get_2d_restriction();
         auto& pp23 = surrogate( 2, 3 ).get_2d_restriction();
         auto& pp24 = surrogate( 2, 4 ).get_2d_restriction();
         auto& pp25 = surrogate( 2, 5 ).get_2d_restriction();
         auto& pp26 = surrogate( 2, 6 ).get_2d_restriction();
         auto& pp27 = surrogate( 2, 7 ).get_2d_restriction();
         auto& pp30 = surrogate( 3, 0 ).get_2d_restriction();
         auto& pp31 = surrogate( 3, 1 ).get_2d_restriction();
         auto& pp32 = surrogate( 3, 2 ).get_2d_restriction();
         auto& pp33 = surrogate( 3, 3 ).get_2d_restriction();
         auto& pp34 = surrogate( 3, 4 ).get_2d_restriction();
         auto& pp35 = surrogate( 3, 5 ).get_2d_restriction();
         auto& pp36 = surrogate( 3, 6 ).get_2d_restriction();
         auto& pp37 = surrogate( 3, 7 ).get_2d_restriction();
         auto& pp40 = surrogate( 4, 0 ).get_2d_restriction();
         auto& pp41 = surrogate( 4, 1 ).get_2d_restriction();
         auto& pp42 = surrogate( 4, 2 ).get_2d_restriction();
         auto& pp43 = surrogate( 4, 3 ).get_2d_restriction();
         auto& pp44 = surrogate( 4, 4 ).get_2d_restriction();
         auto& pp45 = surrogate( 4, 5 ).get_2d_restriction();
         auto& pp46 = surrogate( 4, 6 ).get_2d_restriction();
         auto& pp47 = surrogate( 4, 7 ).get_2d_restriction();
         auto& pp50 = surrogate( 5, 0 ).get_2d_restriction();
         auto& pp51 = surrogate( 5, 1 ).get_2d_restriction();
         auto& pp52 = surrogate( 5, 2 ).get_2d_restriction();
         auto& pp53 = surrogate( 5, 3 ).get_2d_restriction();
         auto& pp54 = surrogate( 5, 4 ).get_2d_restriction();
         auto& pp55 = surrogate( 5, 5 ).get_2d_restriction();
         auto& pp56 = surrogate( 5, 6 ).get_2d_restriction();
         auto& pp57 = surrogate( 5, 7 ).get_2d_restriction();
         auto& pp60 = surrogate( 6, 0 ).get_2d_restriction();
         auto& pp61 = surrogate( 6, 1 ).get_2d_restriction();
         auto& pp62 = surrogate( 6, 2 ).get_2d_restriction();
         auto& pp63 = surrogate( 6, 3 ).get_2d_restriction();
         auto& pp64 = surrogate( 6, 4 ).get_2d_restriction();
         auto& pp65 = surrogate( 6, 5 ).get_2d_restriction();
         auto& pp66 = surrogate( 6, 6 ).get_2d_restriction();
         auto& pp67 = surrogate( 6, 7 ).get_2d_restriction();
         auto& pp70 = surrogate( 7, 0 ).get_2d_restriction();
         auto& pp71 = surrogate( 7, 1 ).get_2d_restriction();
         auto& pp72 = surrogate( 7, 2 ).get_2d_restriction();
         auto& pp73 = surrogate( 7, 3 ).get_2d_restriction();
         auto& pp74 = surrogate( 7, 4 ).get_2d_restriction();
         auto& pp75 = surrogate( 7, 5 ).get_2d_restriction();
         auto& pp76 = surrogate( 7, 6 ).get_2d_restriction();
         auto& pp77 = surrogate( 7, 7 ).get_2d_restriction();
         // restriction to 1d
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p00;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p01;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p02;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p03;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p04;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p05;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p06;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p07;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p10;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p11;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p12;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p13;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p14;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p15;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p16;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p17;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p20;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p21;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p22;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p23;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p24;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p25;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p26;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p27;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p30;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p31;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p32;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p33;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p34;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p35;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p36;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p37;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p40;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p41;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p42;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p43;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p44;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p45;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p46;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p47;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p50;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p51;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p52;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p53;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p54;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p55;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p56;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p57;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p60;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p61;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p62;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p63;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p64;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p65;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p66;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p67;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p70;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p71;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p72;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p73;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p74;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p75;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p76;
         surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p77;

         // iterate over all micro-cells
         indexing::Index micro;
         const auto      n = idx_t( celldof::macrocell::numCellsPerRowByType( level, celldof::CellType::WHITE_UP ) );
         for ( micro.z() = 0; micro.z() < n; micro.z()++ )
         {
            // fix z-coordinate
            const auto z = X[micro.z()];
            for ( auto& s_ij : surrogate )
            {
               s_ij.fix_z( z );
            }

            for ( micro.y() = 0; micro.y() < n - micro.z(); micro.y()++ )
            {
               // fix y-coordinate
               {
                  // code is copied from Polynomial::fix_coord to improve performance
                  const auto y = X[micro.y()];

                  // y^k for k=0,...,q
                  std::array< real_t, DEGREE + 1 > y_pow;
                  y_pow[0] = 1.0;
                  for ( uint_t k = 1; k <= DEGREE; ++k )
                  {
                     y_pow[k] = y_pow[k - 1] * y;
                  }

                  // first index where the 2d extension starts
                  auto jk = dimP;
                  // iterate over coefficients of 1d polynomial
                  for ( uint_t j = 0; j < dimP; ++j )
                  {
                     const auto max_k = uint_t( DEGREE - phi[j].degree() );
                     // k=0
                     p00[j] = pp00[j];
                     p01[j] = pp01[j];
                     p02[j] = pp02[j];
                     p03[j] = pp03[j];
                     p04[j] = pp04[j];
                     p05[j] = pp05[j];
                     p06[j] = pp06[j];
                     p07[j] = pp07[j];
                     p10[j] = pp10[j];
                     p11[j] = pp11[j];
                     p12[j] = pp12[j];
                     p13[j] = pp13[j];
                     p14[j] = pp14[j];
                     p15[j] = pp15[j];
                     p16[j] = pp16[j];
                     p17[j] = pp17[j];
                     p20[j] = pp20[j];
                     p21[j] = pp21[j];
                     p22[j] = pp22[j];
                     p23[j] = pp23[j];
                     p24[j] = pp24[j];
                     p25[j] = pp25[j];
                     p26[j] = pp26[j];
                     p27[j] = pp27[j];
                     p30[j] = pp30[j];
                     p31[j] = pp31[j];
                     p32[j] = pp32[j];
                     p33[j] = pp33[j];
                     p34[j] = pp34[j];
                     p35[j] = pp35[j];
                     p36[j] = pp36[j];
                     p37[j] = pp37[j];
                     p40[j] = pp40[j];
                     p41[j] = pp41[j];
                     p42[j] = pp42[j];
                     p43[j] = pp43[j];
                     p44[j] = pp44[j];
                     p45[j] = pp45[j];
                     p46[j] = pp46[j];
                     p47[j] = pp47[j];
                     p50[j] = pp50[j];
                     p51[j] = pp51[j];
                     p52[j] = pp52[j];
                     p53[j] = pp53[j];
                     p54[j] = pp54[j];
                     p55[j] = pp55[j];
                     p56[j] = pp56[j];
                     p57[j] = pp57[j];
                     p60[j] = pp60[j];
                     p61[j] = pp61[j];
                     p62[j] = pp62[j];
                     p63[j] = pp63[j];
                     p64[j] = pp64[j];
                     p65[j] = pp65[j];
                     p66[j] = pp66[j];
                     p67[j] = pp67[j];
                     p70[j] = pp70[j];
                     p71[j] = pp71[j];
                     p72[j] = pp72[j];
                     p73[j] = pp73[j];
                     p74[j] = pp74[j];
                     p75[j] = pp75[j];
                     p76[j] = pp76[j];
                     p77[j] = pp77[j];
                     for ( int k = 1; k <= max_k; ++k )
                     {
                        p00[j] += pp00[jk] * y_pow[k];
                        p01[j] += pp01[jk] * y_pow[k];
                        p02[j] += pp02[jk] * y_pow[k];
                        p03[j] += pp03[jk] * y_pow[k];
                        p04[j] += pp04[jk] * y_pow[k];
                        p05[j] += pp05[jk] * y_pow[k];
                        p06[j] += pp06[jk] * y_pow[k];
                        p07[j] += pp07[jk] * y_pow[k];
                        p10[j] += pp10[jk] * y_pow[k];
                        p11[j] += pp11[jk] * y_pow[k];
                        p12[j] += pp12[jk] * y_pow[k];
                        p13[j] += pp13[jk] * y_pow[k];
                        p14[j] += pp14[jk] * y_pow[k];
                        p15[j] += pp15[jk] * y_pow[k];
                        p16[j] += pp16[jk] * y_pow[k];
                        p17[j] += pp17[jk] * y_pow[k];
                        p20[j] += pp20[jk] * y_pow[k];
                        p21[j] += pp21[jk] * y_pow[k];
                        p22[j] += pp22[jk] * y_pow[k];
                        p23[j] += pp23[jk] * y_pow[k];
                        p24[j] += pp24[jk] * y_pow[k];
                        p25[j] += pp25[jk] * y_pow[k];
                        p26[j] += pp26[jk] * y_pow[k];
                        p27[j] += pp27[jk] * y_pow[k];
                        p30[j] += pp30[jk] * y_pow[k];
                        p31[j] += pp31[jk] * y_pow[k];
                        p32[j] += pp32[jk] * y_pow[k];
                        p33[j] += pp33[jk] * y_pow[k];
                        p34[j] += pp34[jk] * y_pow[k];
                        p35[j] += pp35[jk] * y_pow[k];
                        p36[j] += pp36[jk] * y_pow[k];
                        p37[j] += pp37[jk] * y_pow[k];
                        p40[j] += pp40[jk] * y_pow[k];
                        p41[j] += pp41[jk] * y_pow[k];
                        p42[j] += pp42[jk] * y_pow[k];
                        p43[j] += pp43[jk] * y_pow[k];
                        p44[j] += pp44[jk] * y_pow[k];
                        p45[j] += pp45[jk] * y_pow[k];
                        p46[j] += pp46[jk] * y_pow[k];
                        p47[j] += pp47[jk] * y_pow[k];
                        p50[j] += pp50[jk] * y_pow[k];
                        p51[j] += pp51[jk] * y_pow[k];
                        p52[j] += pp52[jk] * y_pow[k];
                        p53[j] += pp53[jk] * y_pow[k];
                        p54[j] += pp54[jk] * y_pow[k];
                        p55[j] += pp55[jk] * y_pow[k];
                        p56[j] += pp56[jk] * y_pow[k];
                        p57[j] += pp57[jk] * y_pow[k];
                        p60[j] += pp60[jk] * y_pow[k];
                        p61[j] += pp61[jk] * y_pow[k];
                        p62[j] += pp62[jk] * y_pow[k];
                        p63[j] += pp63[jk] * y_pow[k];
                        p64[j] += pp64[jk] * y_pow[k];
                        p65[j] += pp65[jk] * y_pow[k];
                        p66[j] += pp66[jk] * y_pow[k];
                        p67[j] += pp67[jk] * y_pow[k];
                        p70[j] += pp70[jk] * y_pow[k];
                        p71[j] += pp71[jk] * y_pow[k];
                        p72[j] += pp72[jk] * y_pow[k];
                        p73[j] += pp73[jk] * y_pow[k];
                        p74[j] += pp74[jk] * y_pow[k];
                        p75[j] += pp75[jk] * y_pow[k];
                        p76[j] += pp76[jk] * y_pow[k];
                        p77[j] += pp77[jk] * y_pow[k];
                        ++jk;
                     }
                  }
               }

               // global indices
               micro.x() = 0;
               std::array< uint_t, 8 > vertexDoFIndices{};
               p1::getGlobalCubeIndices3D( level, micro, vertexDoFIndices );

               /* constant part of local stiffness matrix
                  For some reason, accessing pij[0] within the x-loop costs performance
               */
               const real_t c00 = p00[0], c01 = p01[0], c02 = p02[0], c03 = p03[0], c04 = p04[0], c05 = p05[0], c06 = p06[0],
                            c07 = p07[0];
               const real_t c10 = p10[0], c11 = p11[0], c12 = p12[0], c13 = p13[0], c14 = p14[0], c15 = p15[0], c16 = p16[0],
                            c17 = p17[0];
               const real_t c20 = p20[0], c21 = p21[0], c22 = p22[0], c23 = p23[0], c24 = p24[0], c25 = p25[0], c26 = p26[0],
                            c27 = p27[0];
               const real_t c30 = p30[0], c31 = p31[0], c32 = p32[0], c33 = p33[0], c34 = p34[0], c35 = p35[0], c36 = p36[0],
                            c37 = p37[0];
               const real_t c40 = p40[0], c41 = p41[0], c42 = p42[0], c43 = p43[0], c44 = p44[0], c45 = p45[0], c46 = p46[0],
                            c47 = p47[0];
               const real_t c50 = p50[0], c51 = p51[0], c52 = p52[0], c53 = p53[0], c54 = p54[0], c55 = p55[0], c56 = p56[0],
                            c57 = p57[0];
               const real_t c60 = p60[0], c61 = p61[0], c62 = p62[0], c63 = p63[0], c64 = p64[0], c65 = p65[0], c66 = p66[0],
                            c67 = p67[0];
               const real_t c70 = p70[0], c71 = p71[0], c72 = p72[0], c73 = p73[0], c74 = p74[0], c75 = p75[0], c76 = p76[0],
                            c77 = p77[0];

               // loop over all full cubes in the row
               for ( micro.x() = 0; micro.x() < n - 2 - micro.z() - micro.y(); micro.x()++ )
               {
                  // global indices
                  const uint_t g0 = vertexDoFIndices[0] + micro.x();
                  const uint_t g1 = vertexDoFIndices[1] + micro.x();
                  const uint_t g2 = vertexDoFIndices[2] + micro.x();
                  const uint_t g3 = vertexDoFIndices[3] + micro.x();
                  const uint_t g4 = vertexDoFIndices[4] + micro.x();
                  const uint_t g5 = vertexDoFIndices[5] + micro.x();
                  const uint_t g6 = vertexDoFIndices[6] + micro.x();
                  const uint_t g7 = vertexDoFIndices[7] + micro.x();

                  // assemble local element vector v = alpha*src
                  const auto v0 = alpha * srcVertexData[g0];
                  const auto v1 = alpha * srcVertexData[g1];
                  const auto v2 = alpha * srcVertexData[g2];
                  const auto v3 = alpha * srcVertexData[g3];
                  const auto v4 = alpha * srcVertexData[g4];
                  const auto v5 = alpha * srcVertexData[g5];
                  const auto v6 = alpha * srcVertexData[g6];
                  const auto v7 = alpha * srcVertexData[g7];

                  // local stiffness matrix
                  real_t a00, a01, a02, a03, a04, a05, a06, a07;
                  real_t a10, a11, a12, a13, a14, a15, a16, a17;
                  real_t a20, a21, a22, a23, a24, a25, a26, a27;
                  real_t a30, a31, a32, a33, a34, a35, a36, a37;
                  real_t a40, a41, a42, a43, a44, a45, a46, a47;
                  real_t a50, a51, a52, a53, a54, a55, a56, a57;
                  real_t a60, a61, a62, a63, a64, a65, a66, a67;
                  real_t a70, a71, a72, a73, a74, a75, a76, a77;

                  // todo: exploit hard zeros

                  // constant part
                  if constexpr ( Symmetric )
                  {
                     a00 = c00;
                     a10 = c10, a11 = c11;
                     a20 = c20, a21 = c21, a22 = c22;
                     a30 = c30, a31 = c31, a32 = c32, a33 = c33;
                     a40 = c40, a41 = c41, a42 = c42, a43 = c43, a44 = c44;
                     a50 = c50, a51 = c51, a52 = c52, a53 = c53, a54 = c54, a55 = c55;
                     a60 = c60, a61 = c61, a62 = c62, a63 = c63, a64 = c64, a65 = c65, a66 = c66;
                     a70 = c70, a71 = c71, a72 = c72, a73 = c73, a74 = c74, a75 = c75, a76 = c76, a77 = c77;
                  }
                  else
                  {
                     a00 = c00, a01 = c01, a02 = c02, a03 = c03, a04 = c04, a05 = c05, a06 = c06, a07 = c07;
                     a10 = c10, a11 = c11, a12 = c12, a13 = c13, a14 = c14, a15 = c15, a16 = c16, a17 = c17;
                     a20 = c20, a21 = c21, a22 = c22, a23 = c23, a24 = c24, a25 = c25, a26 = c26, a27 = c27;
                     a30 = c30, a31 = c31, a32 = c32, a33 = c33, a34 = c34, a35 = c35, a36 = c36, a37 = c37;
                     a40 = c40, a41 = c41, a42 = c42, a43 = c43, a44 = c44, a45 = c45, a46 = c46, a47 = c47;
                     a50 = c50, a51 = c51, a52 = c52, a53 = c53, a54 = c54, a55 = c55, a56 = c56, a57 = c57;
                     a60 = c60, a61 = c61, a62 = c62, a63 = c63, a64 = c64, a65 = c65, a66 = c66, a67 = c67;
                     a70 = c70, a71 = c71, a72 = c72, a73 = c73, a74 = c74, a75 = c75, a76 = c76, a77 = c77;
                  }
                  // evaluate the 1d polynomials (variable part)
                  const auto x  = X[micro.x()];
                  auto       xk = real_t( 1.0 );
                  for ( uint_t k = 1; k < dimP; ++k )
                  {
                     xk *= x;
                     if constexpr ( Symmetric )
                     {
                        a00 += p00[k] * xk;
                        a10 += p10[k] * xk, a11 += p11[k] * xk;
                        a20 += p20[k] * xk, a21 += p21[k] * xk, a22 += p22[k] * xk;
                        a30 += p30[k] * xk, a31 += p31[k] * xk, a32 += p32[k] * xk, a33 += p33[k] * xk;
                        a40 += p40[k] * xk, a41 += p41[k] * xk, a42 += p42[k] * xk, a43 += p43[k] * xk, a44 += p44[k] * xk;
                        a50 += p50[k] * xk, a51 += p51[k] * xk, a52 += p52[k] * xk, a53 += p53[k] * xk, a54 += p54[k] * xk,
                            a55 += p55[k] * xk;
                        a60 += p60[k] * xk, a61 += p61[k] * xk, a62 += p62[k] * xk, a63 += p63[k] * xk, a64 += p64[k] * xk,
                            a65 += p65[k] * xk, a66 += p66[k] * xk;
                        a70 += p70[k] * xk, a71 += p71[k] * xk, a72 += p72[k] * xk, a73 += p73[k] * xk, a74 += p74[k] * xk,
                            a75 += p75[k] * xk, a76 += p76[k] * xk, a77 += p77[k] * xk;
                     }
                     else
                     {
                        a00 += p00[k] * xk, a01 += p01[k] * xk, a02 += p02[k] * xk, a03 += p03[k] * xk, a04 += p04[k] * xk,
                            a05 += p05[k] * xk, a06 += p06[k] * xk, a07 += p07[k] * xk;
                        a10 += p10[k] * xk, a11 += p11[k] * xk, a12 += p12[k] * xk, a13 += p13[k] * xk, a14 += p14[k] * xk,
                            a15 += p15[k] * xk, a16 += p16[k] * xk, a17 += p17[k] * xk;
                        a20 += p20[k] * xk, a21 += p21[k] * xk, a22 += p22[k] * xk, a23 += p23[k] * xk, a24 += p24[k] * xk,
                            a25 += p25[k] * xk, a26 += p26[k] * xk, a27 += p27[k] * xk;
                        a30 += p30[k] * xk, a31 += p31[k] * xk, a32 += p32[k] * xk, a33 += p33[k] * xk, a34 += p34[k] * xk,
                            a35 += p35[k] * xk, a36 += p36[k] * xk, a37 += p37[k] * xk;
                        a40 += p40[k] * xk, a41 += p41[k] * xk, a42 += p42[k] * xk, a43 += p43[k] * xk, a44 += p44[k] * xk,
                            a45 += p45[k] * xk, a46 += p46[k] * xk, a47 += p47[k] * xk;
                        a50 += p50[k] * xk, a51 += p51[k] * xk, a52 += p52[k] * xk, a53 += p53[k] * xk, a54 += p54[k] * xk,
                            a55 += p55[k] * xk, a56 += p56[k] * xk, a57 += p57[k] * xk;
                        a60 += p60[k] * xk, a61 += p61[k] * xk, a62 += p62[k] * xk, a63 += p63[k] * xk, a64 += p64[k] * xk,
                            a65 += p65[k] * xk, a66 += p66[k] * xk, a67 += p67[k] * xk;
                        a70 += p70[k] * xk, a71 += p71[k] * xk, a72 += p72[k] * xk, a73 += p73[k] * xk, a74 += p74[k] * xk,
                            a75 += p75[k] * xk, a76 += p76[k] * xk, a77 += p77[k] * xk;
                     }
                  }
                  if constexpr ( Symmetric )
                  {
                     a01 = a10, a02 = a20, a03 = a30, a04 = a40, a05 = a50, a06 = a60, a07 = a70;
                     /*      */ a12 = a21, a13 = a31, a14 = a41, a15 = a51, a16 = a61, a17 = a71;
                     /*                 */ a23 = a32, a24 = a42, a25 = a52, a26 = a62, a27 = a72;
                     /*                            */ a34 = a43, a35 = a53, a36 = a63, a37 = a73;
                     /*                                       */ a45 = a54, a46 = a64, a47 = a74;
                     /*                                                  */ a56 = a65, a57 = a75;
                     /*                                                             */ a67 = a76;
                  }

                  // local matvec w=Av
                  const auto w0 = a00 * v0 + a01 * v1 + a02 * v2 + a03 * v3 + a04 * v4 + a05 * v5 + a06 * v6 + a07 * v7;
                  const auto w1 = a10 * v0 + a11 * v1 + a12 * v2 + a13 * v3 + a14 * v4 + a15 * v5 + a16 * v6 + a17 * v7;
                  const auto w2 = a20 * v0 + a21 * v1 + a22 * v2 + a23 * v3 + a24 * v4 + a25 * v5 + a26 * v6 + a27 * v7;
                  const auto w3 = a30 * v0 + a31 * v1 + a32 * v2 + a33 * v3 + a34 * v4 + a35 * v5 + a36 * v6 + a37 * v7;
                  const auto w4 = a40 * v0 + a41 * v1 + a42 * v2 + a43 * v3 + a44 * v4 + a45 * v5 + a46 * v6 + a47 * v7;
                  const auto w5 = a50 * v0 + a51 * v1 + a52 * v2 + a53 * v3 + a54 * v4 + a55 * v5 + a56 * v6 + a57 * v7;
                  const auto w6 = a60 * v0 + a61 * v1 + a62 * v2 + a63 * v3 + a64 * v4 + a65 * v5 + a66 * v6 + a67 * v7;
                  const auto w7 = a70 * v0 + a71 * v1 + a72 * v2 + a73 * v3 + a74 * v4 + a75 * v5 + a76 * v6 + a77 * v7;

                  // write data back into global vector
                  dstVertexData[g0] += w0;
                  dstVertexData[g1] += w1;
                  dstVertexData[g2] += w2;
                  dstVertexData[g3] += w3;
                  dstVertexData[g4] += w4;
                  dstVertexData[g5] += w5;
                  dstVertexData[g6] += w6;
                  dstVertexData[g7] += w7;
               }

               // remainder: partial cube with missing white-down-element
               if ( micro.x() == n - 2 - micro.z() - micro.y() )
               {
                  // todo: use (scaled?) surrogates
                  for ( const auto& cType : celldof::allCellTypes )
                  {
                     if ( cType != celldof::CellType::WHITE_DOWN )
                     {
                        Matrix4r elMat( Matrix4r::Zero() );
                        p1::assembleLocalElementMatrix3D( cell, level, micro, cType, form_, elMat );
                        p1::localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, dstVertexData, elMat, alpha );
                     }
                  }
                  micro.x()++;
               }
               // remainder: single white-up-element at row end
               // if ( micro.x() == n - 1 - micro.z() - micro.y() ) // always true
               {
                  // todo: use (scaled?) surrogates
                  const auto cType = celldof::CellType::WHITE_UP;
                  Matrix4r   elMat( Matrix4r::Zero() );
                  p1::assembleLocalElementMatrix3D( cell, level, micro, cType, form_, elMat );
                  p1::localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, dstVertexData, elMat, alpha );
               }
            }
         }
      }
      else //! legacy: sawtooth. For now, we keep this for testing purposes. Should be removed when no longer required!
      {
         for ( const auto& cType : celldof::allCellTypes )
         {
            auto& surrogate = surrogate_3d_.at( id )[level][cType];
            // monomial basis
            constexpr surrogate::polynomial::Basis< DEGREE > phi;
            // dimension of Pq
            constexpr auto dimP = surrogate::polynomial::dimP( 1, DEGREE );
            // domain
            const PolyDomain X( level );
            // restriction to 2d
            auto& pp00 = surrogate( 0, 0 ).get_2d_restriction();
            auto& pp01 = surrogate( 0, 1 ).get_2d_restriction();
            auto& pp02 = surrogate( 0, 2 ).get_2d_restriction();
            auto& pp03 = surrogate( 0, 3 ).get_2d_restriction();
            auto& pp10 = surrogate( 1, 0 ).get_2d_restriction();
            auto& pp11 = surrogate( 1, 1 ).get_2d_restriction();
            auto& pp12 = surrogate( 1, 2 ).get_2d_restriction();
            auto& pp13 = surrogate( 1, 3 ).get_2d_restriction();
            auto& pp20 = surrogate( 2, 0 ).get_2d_restriction();
            auto& pp21 = surrogate( 2, 1 ).get_2d_restriction();
            auto& pp22 = surrogate( 2, 2 ).get_2d_restriction();
            auto& pp23 = surrogate( 2, 3 ).get_2d_restriction();
            auto& pp30 = surrogate( 3, 0 ).get_2d_restriction();
            auto& pp31 = surrogate( 3, 1 ).get_2d_restriction();
            auto& pp32 = surrogate( 3, 2 ).get_2d_restriction();
            auto& pp33 = surrogate( 3, 3 ).get_2d_restriction();
            // restriction to 1d
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p00;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p01;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p02;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p03;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p10;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p11;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p12;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p13;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p20;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p21;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p22;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p23;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p30;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p31;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p32;
            surrogate::polynomial::Polynomial< real_t, 1, DEGREE > p33;

            // iterate over all micro-cells
            indexing::Index micro;
            const auto      row_end = idx_t( celldof::macrocell::numCellsPerRowByType( level, cType ) );
            for ( micro.z() = 0; micro.z() < row_end; micro.z()++ )
            {
               // fix z-coordinate
               const auto z = X[micro.z()];
               for ( auto& s_ij : surrogate )
               {
                  s_ij.fix_z( z );
               }

               for ( micro.y() = 0; micro.y() < row_end - micro.z(); micro.y()++ )
               {
                  // fix y-coordinate
                  {
                     // code is copied from Polynomial::fix_coord to improve performance
                     const auto y = X[micro.y()];

                     // y^k for k=0,...,q
                     std::array< real_t, DEGREE + 1 > y_pow;
                     y_pow[0] = 1.0;
                     for ( uint_t k = 1; k <= DEGREE; ++k )
                     {
                        y_pow[k] = y_pow[k - 1] * y;
                     }

                     // first index where the 2d extension starts
                     auto jk = dimP;
                     // iterate over coefficients of 1d polynomial
                     for ( uint_t j = 0; j < dimP; ++j )
                     {
                        const auto max_k = uint_t( DEGREE - phi[j].degree() );
                        // k=0
                        p00[j] = pp00[j];
                        p01[j] = pp01[j];
                        p02[j] = pp02[j];
                        p03[j] = pp03[j];
                        p10[j] = pp10[j];
                        p11[j] = pp11[j];
                        p12[j] = pp12[j];
                        p13[j] = pp13[j];
                        p20[j] = pp20[j];
                        p21[j] = pp21[j];
                        p22[j] = pp22[j];
                        p23[j] = pp23[j];
                        p30[j] = pp30[j];
                        p31[j] = pp31[j];
                        p32[j] = pp32[j];
                        p33[j] = pp33[j];
                        for ( int k = 1; k <= max_k; ++k )
                        {
                           p00[j] += pp00[jk] * y_pow[k];
                           p01[j] += pp01[jk] * y_pow[k];
                           p02[j] += pp02[jk] * y_pow[k];
                           p03[j] += pp03[jk] * y_pow[k];
                           p10[j] += pp10[jk] * y_pow[k];
                           p11[j] += pp11[jk] * y_pow[k];
                           p12[j] += pp12[jk] * y_pow[k];
                           p13[j] += pp13[jk] * y_pow[k];
                           p20[j] += pp20[jk] * y_pow[k];
                           p21[j] += pp21[jk] * y_pow[k];
                           p22[j] += pp22[jk] * y_pow[k];
                           p23[j] += pp23[jk] * y_pow[k];
                           p30[j] += pp30[jk] * y_pow[k];
                           p31[j] += pp31[jk] * y_pow[k];
                           p32[j] += pp32[jk] * y_pow[k];
                           p33[j] += pp33[jk] * y_pow[k];
                           ++jk;
                        }
                     }
                  }

                  /* global dof indices
                        Computing global indices seems to be rather slow so we factor it
                        out of the x-loop.
                     */
                  micro.x() = 0;
                  std::array< uint_t, 4 > vertexDoFIndices{};
                  vertexdof::getVertexDoFDataIndicesFromMicroCell(micro, cType, level, vertexDoFIndices);

                  /* constant part of local stiffness matrix
                        For some reason, accessing pij[0] within the x-loop costs performance
                     */
                  const real_t c00 = p00[0], c01 = p01[0], c02 = p02[0], c03 = p03[0];
                  const real_t c10 = p10[0], c11 = p11[0], c12 = p12[0], c13 = p13[0];
                  const real_t c20 = p20[0], c21 = p21[0], c22 = p22[0], c23 = p23[0];
                  const real_t c30 = p30[0], c31 = p31[0], c32 = p32[0], c33 = p33[0];

                  for ( micro.x() = 0; micro.x() < row_end - micro.z() - micro.y(); micro.x()++ )
                  {
                     // global indices
                     const uint_t g0 = vertexDoFIndices[0] + micro.x();
                     const uint_t g1 = vertexDoFIndices[1] + micro.x();
                     const uint_t g2 = vertexDoFIndices[2] + micro.x();
                     const uint_t g3 = vertexDoFIndices[3] + micro.x();

                     // local stiffness matrix
                     real_t a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
                     // constant part
                     if constexpr ( Symmetric )
                     {
                        a00 = c00;
                        a10 = c10, a11 = c11;
                        a20 = c20, a21 = c21, a22 = c22;
                        a30 = c30, a31 = c31, a32 = c32, a33 = c33;
                     }
                     else
                     {
                        a00 = c00, a01 = c01, a02 = c02, a03 = c03;
                        a10 = c10, a11 = c11, a12 = c12, a13 = c13;
                        a20 = c20, a21 = c21, a22 = c22, a23 = c23;
                        a30 = c30, a31 = c31, a32 = c32, a33 = c33;
                     }
                     // evaluate the 1d polynomials (variable part)
                     const auto x  = X[micro.x()];
                     auto       xk = real_t( 1.0 );
                     for ( uint_t k = 1; k < dimP; ++k )
                     {
                        xk *= x;
                        if constexpr ( Symmetric )
                        {
                           a00 += p00[k] * xk;
                           a10 += p10[k] * xk, a11 += p11[k] * xk;
                           a20 += p20[k] * xk, a21 += p21[k] * xk, a22 += p22[k] * xk;
                           a30 += p30[k] * xk, a31 += p31[k] * xk, a32 += p32[k] * xk, a33 += p33[k] * xk;
                        }
                        else
                        {
                           a00 += p00[k] * xk, a01 += p01[k] * xk, a02 += p02[k] * xk, a03 += p03[k] * xk;
                           a10 += p10[k] * xk, a11 += p11[k] * xk, a12 += p12[k] * xk, a13 += p13[k] * xk;
                           a20 += p20[k] * xk, a21 += p21[k] * xk, a22 += p22[k] * xk, a23 += p23[k] * xk;
                           a30 += p30[k] * xk, a31 += p31[k] * xk, a32 += p32[k] * xk, a33 += p33[k] * xk;
                        }
                     }
                     if constexpr ( Symmetric )
                     {
                        a01 = a10, a02 = a20, a03 = a30;
                        /*      */ a12 = a21, a13 = a31;
                        /*                 */ a23 = a32;
                     }

                     // assemble local element vector
                     const auto src0 = alpha * srcVertexData[g0];
                     const auto src1 = alpha * srcVertexData[g1];
                     const auto src2 = alpha * srcVertexData[g2];
                     const auto src3 = alpha * srcVertexData[g3];

                     // local matvec
                     const auto dst0 = a00 * src0 + a01 * src1 + a02 * src2 + a03 * src3;
                     const auto dst1 = a10 * src0 + a11 * src1 + a12 * src2 + a13 * src3;
                     const auto dst2 = a20 * src0 + a21 * src1 + a22 * src2 + a23 * src3;
                     const auto dst3 = a30 * src0 + a31 * src1 + a32 * src2 + a33 * src3;

                     // write data back into global vector
                     dstVertexData[g0] += dst0;
                     dstVertexData[g1] += dst1;
                     dstVertexData[g2] += dst2;
                     dstVertexData[g3] += dst3;
                  }
               }
            }
         }
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::computeDiagonalOperatorValues( bool invert )
{
   std::shared_ptr< P1Function< real_t > > targetFunction;
   if ( invert )
   {
      if ( !inverseDiagonalValues_ )
      {
         inverseDiagonalValues_ =
             std::make_shared< P1Function< real_t > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = inverseDiagonalValues_;
   }
   else
   {
      if ( !diagonalValues_ )
      {
         diagonalValues_ = std::make_shared< P1Function< real_t > >( "diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = diagonalValues_;
   }

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      // Make sure that halos are up-to-date (can we improve communication here?)
      communication::syncFunctionBetweenPrimitives( *targetFunction, level );

      // Zero destination before performing additive computation
      targetFunction->setToZero( level );

      // For 3D we work on cells and for 2D on faces
      if ( storage_->hasGlobalCells() )
      {
         // we only perform computations on cell primitives
         for ( auto& macroIter : storage_->getCells() )
         {
            Cell& cell = *macroIter.second;

            // get hold of the actual numerical data
            auto    vertexDoFIdx = targetFunction->getCellDataID();
            real_t* vertexData   = cell.getData( vertexDoFIdx )->getPointer( level );

            // loop over micro-cells
            for ( const auto& cType : celldof::allCellTypes )
            {
               diagonal_contributions_3d( cell, level, cType, vertexData );
            }
         }

         // Push result to lower-dimensional primitives
         targetFunction->communicateAdditively< Cell, Face >( level );
         targetFunction->communicateAdditively< Cell, Edge >( level );
         targetFunction->communicateAdditively< Cell, Vertex >( level );
      }
      else
      {
         // we only perform computations on face primitives
         for ( auto& it : storage_->getFaces() )
         {
            Face& face = *it.second;

            // get hold of the actual numerical data in the two functions
            auto    vertexDoFIdx = targetFunction->getFaceDataID();
            real_t* vertexData   = face.getData( vertexDoFIdx )->getPointer( level );

            for ( const auto& fType : facedof::allFaceTypes )
            {
               diagonal_contributions_2d( face, level, fType, vertexData );
            }
         }

         // Push result to lower-dimensional primitives
         targetFunction->communicateAdditively< Face, Edge >( level );
         targetFunction->communicateAdditively< Face, Vertex >( level );
      }
      // Retrieve assembled data values
      targetFunction->communicate< Vertex, Edge >( level );
      targetFunction->communicate< Edge, Face >( level );
      targetFunction->communicate< Face, Cell >( level );

      // Invert values if desired (note: using false below means we only invert in the interior of the primitives,
      // the values in the halos are untouched; should be okay for using diagonalValue_ in smoothers)
      if ( invert )
      {
         targetFunction->invertElementwise( level, All, false );
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::diagonal_contributions_2d( const Face&             face,
                                                                                             const uint_t            level,
                                                                                             const facedof::FaceType fType,
                                                                                             real_t* const dstVertexData )
{
   // add local contributions of micro to global diagonal
   auto extract_diagonal_values = [&]( const indexing::Index& micro, const Matrix3r& elMat ) {
      // obtain data indices of dofs associated with micro-face
      std::array< uint_t, 3 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroFace( micro, fType, level, vertexDoFIndices );
      // extract matrix diagonal
      for ( idx_t k = 0; k < 3; ++k )
      {
         dstVertexData[vertexDoFIndices[uint_c( k )]] += elMat( k, k );
      }
   };

   auto& id = face.getID();

   if ( level < min_lvl_for_surrogate )
   { // use precomputed local matrices
      auto& a_loc = a_loc_2d_.at( id )[level];

      for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
      {
         extract_diagonal_values( micro, a_loc( fType, micro ) );
      }
   }
   else
   { // use surrogate polynomials to approximate local matrices
      auto& surrogate = surrogate_2d_.at( id )[level][fType];
      // domain of surrogates
      PolyDomain X( level );
      // local stiffness matrix
      Matrix3r elMat( Matrix3r::Zero() );
      // todo: use optimized polynomial evaluation
      for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
      {
         auto x = X( micro );

         for ( idx_t i = 0; i < surrogate.rows(); ++i )
         {
            elMat( i, i ) = surrogate( i, i ).eval_naive( x );
         }
         extract_diagonal_values( micro, elMat );
      }
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::diagonal_contributions_3d( const Cell&             cell,
                                                                                             const uint_t            level,
                                                                                             const celldof::CellType cType,
                                                                                             real_t* const dstVertexData )
{
   // add local contributions of micro to global diagonal
   auto extract_diagonal_values = [&]( const indexing::Index& micro, const Matrix4r& elMat ) {
      // obtain data indices of dofs associated with micro-cell
      std::array< uint_t, 4 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroCell( micro, cType, level, vertexDoFIndices );
      // extract matrix diagonal
      for ( idx_t k = 0; k < 4; ++k )
      {
         dstVertexData[vertexDoFIndices[uint_c( k )]] += elMat( k, k );
      }
   };

   auto& id = cell.getID();

   if ( level < min_lvl_for_surrogate )
   { // use precomputed local matrices
      auto& a_loc = a_loc_3d_.at( id )[level];

      for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
      {
         extract_diagonal_values( micro, a_loc( cType, micro ) );
      }
   }
   else
   { // use surrogate polynomials to approximate local matrices
      auto& surrogate = surrogate_3d_.at( id )[level][cType];
      // domain of surrogates
      PolyDomain X( level );
      // local stiffness matrix
      Matrix4r elMat( Matrix4r::Zero() );
      // todo: use optimized polynomial evaluation
      for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
      {
         auto x = X( micro );

         for ( idx_t i = 0; i < surrogate.rows(); ++i )
         {
            elMat( i, i ) = surrogate( i, i ).eval_naive( x );
         }
         extract_diagonal_values( micro, elMat );
      }
   }
}

// Assemble operator as sparse matrix
template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                            const P1Function< idx_t >&                  src,
                                                                            const P1Function< idx_t >&                  dst,
                                                                            uint_t                                      level,
                                                                            DoFType flag ) const
{
   // todo
   WALBERLA_UNUSED( mat );
   WALBERLA_UNUSED( src );
   WALBERLA_UNUSED( dst );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( flag );
   WALBERLA_ABORT( "Not implemented!" );
}
// // P1ElementwiseLaplaceOperator
// template class P1ElementwiseSurrogateOperator<
//     P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >;

// // P1ElementwisePolarLaplaceOperator
// template class P1ElementwiseSurrogateOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > >;

// // P1ElementwiseMassOperator
// template class P1ElementwiseSurrogateOperator<
//     P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >;

// // P1ElementwisePSPGOperator
// template class P1ElementwiseSurrogateOperator<
//     P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >;

// template class P1ElementwiseSurrogateOperator< P1LinearCombinationForm >;

// // P1ElementwiseBlendingMassOperator3D
// template class P1ElementwiseSurrogateOperator< forms::p1_mass_blending_q4 >;

// // P1ElementwiseBlendingLaplaceOperator
// template class P1ElementwiseSurrogateOperator< forms::p1_diffusion_blending_q3 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_diffusion_blending_q2 >;

// // Needed for P1Blending(Inverse)DiagonalOperator
// template class P1ElementwiseSurrogateOperator< P1RowSumForm >;

template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 0, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 1, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 2, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 3, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 4, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 5, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 6, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 7, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 8, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 9, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 10, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 11, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 12, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 13, true >;
template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3, 14, true >;

// template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_blending_q3 >;

// template class P1ElementwiseSurrogateOperator<
//     P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >;
// template class P1ElementwiseSurrogateOperator<
//     P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >;
// template class P1ElementwiseSurrogateOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > >;

// template class P1ElementwiseSurrogateOperator<
//     P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >;
// template class P1ElementwiseSurrogateOperator<
//     P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >;
// template class P1ElementwiseSurrogateOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > >;

// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_0_0_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_0_1_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_0_2_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_1_0_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_1_1_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_1_2_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_2_0_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_2_1_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsiloncc_2_2_affine_q2 >;

// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_0_0_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_0_1_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_0_2_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_1_0_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_1_1_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_1_2_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_2_0_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_2_1_affine_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_2_2_affine_q2 >;

// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_0_0_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_0_1_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_0_2_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_1_0_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_1_1_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_1_2_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_2_0_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_2_1_blending_q2 >;
// template class P1ElementwiseSurrogateOperator< forms::p1_epsilonvar_2_2_blending_q2 >;

// template class P1ElementwiseSurrogateOperator< forms::p1_k_mass_affine_q4 >;

// // This is a slight misuse of the P1ElementwiseOperator class, since the spherical
// // elements are not P1. However, the SphericalElementFunction, like the P1Function
// // is only an alias for the VertexDoFFunction, so we can re-use this operator.
// template class P1ElementwiseSurrogateOperator< SphericalElementFormMass >;

} // namespace hyteg
