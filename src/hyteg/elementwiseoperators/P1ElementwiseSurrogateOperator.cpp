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

#include "P1ElementwiseOperator.hpp"

namespace hyteg {

template < class P1Form >
P1ElementwiseSurrogateOperator< P1Form >::P1ElementwiseSurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                          size_t                                     minLevel,
                                                                          size_t                                     maxLevel )
: P1ElementwiseSurrogateOperator< P1Form >( storage, minLevel, maxLevel, P1Form() )
{}

template < class P1Form >
P1ElementwiseSurrogateOperator< P1Form >::P1ElementwiseSurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                          size_t                                     minLevel,
                                                                          size_t                                     maxLevel,
                                                                          const P1Form&                              form )
: Operator( storage, minLevel, maxLevel )
, form_( form )
, is_initialized_( false )
, lsq_( maxLevel + 1 )
, downsampling_( maxLevel + 1 )
, poly_degree_( maxLevel + 1 )
, a_loc_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ) )
, a_loc_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ) )
, surrogate_2d_( storage, maxLevel )
, surrogate_3d_( storage, maxLevel )
{}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::init( uint8_t            poly_degree,
                                                     size_t             downsampling,
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
      // adjust downsampling for this level
      auto ds = downsampling;
      while ( LSQ::max_degree( level, ds ) < poly_degree && ds > 1 )
      {
         --ds;
      }
      // adjust polynomial degree for this level
      auto q = std::min( LSQ::max_degree( level, ds ), poly_degree );

      // initialize least squares approximation
      if ( lsq_[level] == nullptr || downsampling_[level] != ds || poly_degree_[level] != q )
      {
         if ( path_to_svd == "" )
         {
            lsq_[level] = std::make_shared< LSQ >( dim, q, level, ds );
         }
         else
         {
            lsq_[level] = std::make_shared< LSQ >( path_to_svd, dim, q, level, ds );
         }
         downsampling_[level] = ds;
         poly_degree_[level]  = q;
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

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::precompute_local_stiffness_2d( uint_t level )
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
            assembleLocalElementMatrix2D( *face, level, micro, fType, form_, elMat );
         }
      }
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::precompute_local_stiffness_3d( uint_t level )
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
            assembleLocalElementMatrix3D( *cell, level, micro, cType, form_, elMat );
         }
      }
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::compute_local_surrogates_2d( uint_t level )
{
   auto  q   = poly_degree_[level];
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
            assembleLocalElementMatrix2D( *face, level, it.ijk(), fType, form_, elMat );
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
               surrogate( i, j ) = Poly( 2, q, coeffs );
            }
         }
      }
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::compute_local_surrogates_3d( uint_t level )
{
   auto  q   = poly_degree_[level];
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
            assembleLocalElementMatrix3D( *cell, level, it.ijk(), cType, form_, elMat );
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
               surrogate( i, j ) = Poly( 3, q, coeffs );
            }
         }
      }
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::apply( const P1Function< real_t >& src,
                                                      const P1Function< real_t >& dst,
                                                      size_t                      level,
                                                      DoFType                     flag,
                                                      UpdateType                  updateType ) const
{
   return gemv( real_c( 1 ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::gemv( const real_t&               alpha,
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
         for ( const auto& cType : celldof::allCellTypes )
         {
            apply_3d( cell, level, cType, srcVertexData, dstVertexData, alpha );
         }
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
         for ( const auto& fType : facedof::allFaceTypes )
         {
            apply_2d( face, level, fType, srcVertexData, dstVertexData, alpha );
         }
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

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::smooth_jac( const P1Function< real_t >& dst,
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

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::apply_2d( const Face&             face,
                                                         const uint_t            level,
                                                         const facedof::FaceType fType,
                                                         const real_t* const     srcVertexData,
                                                         real_t* const           dstVertexData,
                                                         const real_t&           alpha ) const
{
   auto& id = face.getID();

   if ( level < min_lvl_for_surrogate )
   { // use precomputed local matrices
      auto& a_loc = a_loc_2d_.at( id )[level];

      for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
      {
         localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, dstVertexData, a_loc( fType, micro ), alpha );
      }
   }
   else
   { // use surrogate polynomials to approximate local matrices
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
            localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, dstVertexData, elMat, alpha );
         }
      }
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::apply_3d( const Cell&             cell,
                                                         const uint_t            level,
                                                         const celldof::CellType cType,
                                                         const real_t* const     srcVertexData,
                                                         real_t* const           dstVertexData,
                                                         const real_t&           alpha ) const
{
   auto& id = cell.getID();

   if ( level < min_lvl_for_surrogate )
   { // use precomputed local matrices
      auto& a_loc = a_loc_3d_.at( id )[level];

      for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
      {
         localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, dstVertexData, a_loc( cType, micro ), alpha );
      }
   }
   else
   {
      auto& surrogate = surrogate_3d_.at( id )[level][cType];
      // domain of surrogates
      PolyDomain X( level );
      // local stiffness matrix
      Matrix4r elMat( Matrix4r::Zero() );
      // iterate over all micro-cells
      indexing::Index micro;
      const auto      row_end = idx_t( celldof::macrocell::numCellsPerRowByType( level, cType ) );
      for ( micro.z() = 0; micro.z() < row_end; micro.z()++ )
      {
         // restrict to 2d polynomial
         auto z = X[micro.z()];
         for ( idx_t i = 0; i < surrogate.rows(); ++i )
         {
            for ( idx_t j = 0; j < surrogate.cols(); ++j )
            {
               surrogate( i, j ).fix_z( z );
            }
         }

         for ( micro.y() = 0; micro.y() < row_end - micro.z(); micro.y()++ )
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

            for ( micro.x() = 0; micro.x() < row_end - micro.z() - micro.y(); micro.x()++ )
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
               localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, dstVertexData, elMat, alpha );
            }
         }
      }
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::localMatrixVectorMultiply2D( const uint_t           level,
                                                                            const indexing::Index& microFace,
                                                                            facedof::FaceType      fType,
                                                                            const real_t* const    srcVertexData,
                                                                            real_t* const          dstVertexData,
                                                                            const Matrix3r&        elMat,
                                                                            const real_t&          alpha ) const
{
   WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

   // obtain data indices of dofs associated with micro-face
   std::array< uint_t, 3 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, fType, level, vertexDoFIndices );

   // assemble local element vector
   Point3D elVecOld, elVecNew;
   for ( idx_t k = 0; k < 3; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
   }

   // apply matrix (operator locally)
   elVecNew = alpha * ( elMat * elVecOld );

   // redistribute result from "local" to "global vector"
   for ( idx_t k = 0; k < 3; ++k )
   {
      dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::localMatrixVectorMultiply3D( const uint_t            level,
                                                                            const indexing::Index&  microCell,
                                                                            const celldof::CellType cType,
                                                                            const real_t* const     srcVertexData,
                                                                            real_t* const           dstVertexData,
                                                                            const Matrix4r&         elMat,
                                                                            const real_t&           alpha ) const
{
   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   // assemble local element vector
   Point4D elVecOld, elVecNew;
   for ( idx_t k = 0; k < 4; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
   }

   // apply matrix (operator locally)
   elVecNew = alpha * ( elMat * elVecOld );

   // redistribute result from "local" to "global vector"
   for ( idx_t k = 0; k < 4; ++k )
   {
      dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
   }
}

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::computeDiagonalOperatorValues( bool invert )
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

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::diagonal_contributions_2d( const Face&             face,
                                                                          const uint_t            level,
                                                                          const facedof::FaceType fType,
                                                                          real_t* const           dstVertexData )
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

template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::diagonal_contributions_3d( const Cell&             cell,
                                                                          const uint_t            level,
                                                                          const celldof::CellType cType,
                                                                          real_t* const           dstVertexData )
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
template < class P1Form >
void P1ElementwiseSurrogateOperator< P1Form >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                         const P1Function< idx_t >&                  src,
                                                         const P1Function< idx_t >&                  dst,
                                                         uint_t                                      level,
                                                         DoFType                                     flag ) const
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

template class P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3 >;
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
