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
, a_loc_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ) )
, a_loc_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ) )
, surrogate_2d_( storage, maxLevel )
, surrogate_3d_( storage, maxLevel )
, surrogate_cube_2d_( storage, maxLevel )
, surrogate_cube_3d_( storage, maxLevel )
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
   // todo correct interpolation at row end
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
   }
}

template < class P1Form, uint8_t DEGREE, bool Symmetric >
void P1ElementwiseSurrogateOperator< P1Form, DEGREE, Symmetric >::compute_local_surrogates_3d( uint_t level )
{
   // todo correct interpolation at row end
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
                  /* code is copied from Polynomial::fix_coord to improve performance
               */
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
               p1::getGlobalIndices3D( cType, level, micro, vertexDoFIndices );

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
                  const auto x    = X[micro.x()];
                  auto       xpow = real_t( 1.0 );
                  for ( uint_t k = 1; k < dimP; ++k )
                  {
                     xpow *= x;
                     if constexpr ( Symmetric )
                     {
                        a00 += p00[k] * xpow;
                        a10 += p10[k] * xpow, a11 += p11[k] * xpow;
                        a20 += p20[k] * xpow, a21 += p21[k] * xpow, a22 += p22[k] * xpow;
                        a30 += p30[k] * xpow, a31 += p31[k] * xpow, a32 += p32[k] * xpow, a33 += p33[k] * xpow;
                     }
                     else
                     {
                        a00 += p00[k] * xpow, a01 += p01[k] * xpow, a02 += p02[k] * xpow, a03 += p03[k] * xpow;
                        a10 += p10[k] * xpow, a11 += p11[k] * xpow, a12 += p12[k] * xpow, a13 += p13[k] * xpow;
                        a20 += p20[k] * xpow, a21 += p21[k] * xpow, a22 += p22[k] * xpow, a23 += p23[k] * xpow;
                        a30 += p30[k] * xpow, a31 += p31[k] * xpow, a32 += p32[k] * xpow, a33 += p33[k] * xpow;
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
