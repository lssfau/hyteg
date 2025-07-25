/*
 * Copyright (c) 2017-2025 Marcus Mohr, Nils Kohl, Andreas Burkhart.
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

#include "P1ElementwiseOperator.hpp"

#include "hyteg/forms/P1RowSumForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_manual/SphericalElementFormMass.hpp"

#include "P1LocalOperations.hpp"

namespace hyteg {

template < class P1Form >
P1ElementwiseOperator< P1Form >::P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                        uint_t                                     minLevel,
                                                        uint_t                                     maxLevel )
: P1ElementwiseOperator< P1Form >( storage, minLevel, maxLevel, P1Form(), true )
{}

template < class P1Form >
P1ElementwiseOperator< P1Form >::P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                        uint_t                                     minLevel,
                                                        uint_t                                     maxLevel,
                                                        const P1Form&                              form )
: P1ElementwiseOperator< P1Form >( storage, minLevel, maxLevel, form, true )
{}

template < class P1Form >
P1ElementwiseOperator< P1Form >::P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                        uint_t                                     minLevel,
                                                        uint_t                                     maxLevel,
                                                        const P1Form&                              form,
                                                        bool                                       needsInverseDiagEntries )
: Operator( storage, minLevel, maxLevel )
, form_( form )
, localElementMatricesPrecomputed_( false )
{
   if ( needsInverseDiagEntries )
   {
      computeInverseDiagonalOperatorValues();
   }
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::apply( const P1Function< real_t >& src,
                                             const P1Function< real_t >& dst,
                                             uint_t                      level,
                                             DoFType                     flag,
                                             UpdateType                  updateType ) const
{
   return gemv( real_c( 1 ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::applyScaled( const real_t&               alpha,
                                                   const P1Function< real_t >& src,
                                                   const P1Function< real_t >& dst,
                                                   uint_t                      level,
                                                   DoFType                     flag,
                                                   UpdateType                  updateType ) const
{
   return gemv( alpha, src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::gemv( const real_t&               alpha,
                                            const P1Function< real_t >& src,
                                            const real_t&               beta,
                                            const P1Function< real_t >& dst,
                                            uint_t                      level,
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

         Matrix4r elMat( Matrix4r::Zero() );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               if ( localElementMatricesPrecomputed_ )
               {
                  elMat = localElementMatrix3D( cell, level, micro, cType );
               }
               else
               {
                  p1::assembleLocalElementMatrix3D( cell, level, micro, cType, form_, elMat );
               }

               p1::localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, dstVertexData, elMat, alpha );
            }
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

         Point3D         v0, v1, v2;
         indexing::Index nodeIdx;
         indexing::Index offset;

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

         Matrix3r elMat( Matrix3r::Zero() );

         // loop over micro-faces
         for ( const auto& fType : facedof::allFaceTypes )
         {
            for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
            {
               if ( localElementMatricesPrecomputed_ )
               {
                  elMat = localElementMatrix2D( face, level, micro, fType );
               }
               else
               {
                  p1::assembleLocalElementMatrix2D( face, level, micro, fType, form_, elMat );
               }

               p1::localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, dstVertexData, elMat, alpha );
            }
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
void P1ElementwiseOperator< P1Form >::smooth_jac_scaled( const real_t&               alpha,
                                                         const P1Function< real_t >& dst,
                                                         const P1Function< real_t >& rhs,
                                                         const P1Function< real_t >& src,
                                                         real_t                      omega,
                                                         size_t                      level,
                                                         DoFType                     flag ) const
{
   this->startTiming( "smooth_jac" );

   // compute the current residual
   this->applyScaled( alpha, src, dst, level, flag );
   dst.assign( { real_c( 1 ), real_c( -1 ) }, { rhs, dst }, level, flag );

   // perform Jacobi update step
   dst.multElementwise( { *getInverseDiagonalValues(), dst }, level, flag );
   dst.assign( { 1.0, omega }, { src, dst }, level, flag );

   this->stopTiming( "smooth_jac" );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::smooth_jac( const P1Function< real_t >& dst,
                                                  const P1Function< real_t >& rhs,
                                                  const P1Function< real_t >& src,
                                                  real_t                      omega,
                                                  size_t                      level,
                                                  DoFType                     flag ) const
{
   smooth_jac_scaled( real_c( 1 ), dst, rhs, src, omega, level, flag );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::computeDiagonalOperatorValues( bool invert, bool lumped, const real_t& alpha )
{
   std::shared_ptr< P1Function< real_t > > targetFunction;
   if ( invert )
   {
      if ( !lumped )
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
         if ( !lumpedInverseDiagonalValues_ )
         {
            lumpedInverseDiagonalValues_ =
                std::make_shared< P1Function< real_t > >( "lumped inverse diagonal entries", storage_, minLevel_, maxLevel_ );
         }
         targetFunction = lumpedInverseDiagonalValues_;
      }
   }
   else
   {
      if ( !lumped )
      {
         if ( !diagonalValues_ )
         {
            diagonalValues_ = std::make_shared< P1Function< real_t > >( "diagonal entries", storage_, minLevel_, maxLevel_ );
         }
         targetFunction = diagonalValues_;
      }
      else
      {
         if ( !lumpedDiagonalValues_ )
         {
            lumpedDiagonalValues_ =
                std::make_shared< P1Function< real_t > >( "lumped diagonal entries", storage_, minLevel_, maxLevel_ );
         }
         targetFunction = lumpedDiagonalValues_;
      }
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
            PrimitiveDataID< FunctionMemory< real_t >, Cell > diagVertexDoFIdx = targetFunction->getCellDataID();
            real_t* diagVertexData = cell.getData( diagVertexDoFIdx )->getPointer( level );

            // loop over micro-cells
            for ( const auto& cType : celldof::allCellTypes )
            {
               for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
               {
                  computeLocalDiagonalContributions3D( cell, level, micro, cType, diagVertexData, alpha, lumped );
               }
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

            uint_t          rowsize       = levelinfo::num_microvertices_per_edge( level );
            uint_t          inner_rowsize = rowsize;
            idx_t           xIdx, yIdx;
            Point3D         v0, v1, v2;
            indexing::Index nodeIdx;
            indexing::Index offset;

            // get hold of the actual numerical data in the two functions
            PrimitiveDataID< FunctionMemory< real_t >, Face > vertexDoFIdx = targetFunction->getFaceDataID();
            real_t*                                           vertexData   = face.getData( vertexDoFIdx )->getPointer( level );

            // now loop over micro-faces of macro-face
            for ( yIdx = 0; yIdx < idx_t( rowsize ) - 2; ++yIdx )
            {
               // loop over vertices in row with two associated triangles
               for ( xIdx = 1; xIdx < idx_t( inner_rowsize ) - 1; ++xIdx )
               {
                  // we associate two elements with current micro-vertex
                  computeLocalDiagonalContributions2D(
                      face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementN, vertexData, alpha, lumped );
                  computeLocalDiagonalContributions2D(
                      face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, vertexData, alpha, lumped );
               }
               --inner_rowsize;

               // final micro-vertex in row has only one associated micro-face
               computeLocalDiagonalContributions2D(
                   face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, vertexData, alpha, lumped );
            }

            // top north-west micro-element not treated, yet
            computeLocalDiagonalContributions2D(
                face, level, 1, yIdx, P1Elements::P1Elements2D::elementNW, vertexData, alpha, lumped );
         }

         // Push result to lower-dimensional primitives
         targetFunction->communicateAdditively< Face, Edge >( level );
         targetFunction->communicateAdditively< Face, Vertex >( level );

         // Retrieve assembled data values
         targetFunction->communicate< Vertex, Edge >( level );
         targetFunction->communicate< Edge, Face >( level );
      }

      // Invert values if desired (note: using false below means we only invert in the interior of the primitives,
      // the values in the halos are untouched; should be okay for using diagonalValue_ in smoothers)
      if ( invert )
      {
         targetFunction->invertElementwise( level, All, false );
      }
   }
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::computeAndStoreLocalElementMatrices()
{
   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      // For 3D we work on cells and for 2D on faces
      if ( storage_->hasGlobalCells() )
      {
         const uint_t numMicroCellsPerMacroCell = celldof::macrocell::numMicroCellsPerMacroCellTotal( level );

         for ( const auto& it : storage_->getCells() )
         {
            auto cellID = it.first;
            auto cell   = it.second;

            auto& elementMatrices = localElementMatrices3D_[cellID][level];

            if ( !localElementMatricesPrecomputed_ )
            {
               elementMatrices.resize( numMicroCellsPerMacroCell );
            }

            for ( const auto& cType : celldof::allCellTypes )
            {
               for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
               {
                  Matrix4r& elMat = localElementMatrix3D( *cell, level, micro, cType );
                  elMat.setZero();
                  p1::assembleLocalElementMatrix3D( *cell, level, micro, cType, form_, elMat );
               }
            }
         }
      }
      else
      {
         const uint_t numMicroFacesPerMacroFace = levelinfo::num_microfaces_per_face( level );

         for ( const auto& it : storage_->getFaces() )
         {
            auto faceID = it.first;
            auto face   = it.second;

            auto& elementMatrices = localElementMatrices2D_[faceID][level];

            if ( !localElementMatricesPrecomputed_ )
            {
               elementMatrices.resize( numMicroFacesPerMacroFace );
            }

            for ( const auto& fType : facedof::allFaceTypes )
            {
               for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
               {
                  Matrix3r& elMat = localElementMatrix2D( *face, level, micro, fType );
                  elMat.setZero();
                  p1::assembleLocalElementMatrix2D( *face, level, micro, fType, form_, elMat );
               }
            }
         }
      }
   }

   localElementMatricesPrecomputed_ = true;
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::computeLocalDiagonalContributions2D( const Face&                                face,
                                                                           const uint_t                               level,
                                                                           const idx_t                                xIdx,
                                                                           const idx_t                                yIdx,
                                                                           const P1Elements::P1Elements2D::P1Element& element,
                                                                           real_t* const dstVertexData,
                                                                           const real_t& alpha,
                                                                           bool          lumped )
{
   Matrix3r                elMat( Matrix3r::Zero() );
   indexing::Index         nodeIdx;
   indexing::Index         offset;
   Point3D                 v0, v1, v2;
   std::array< uint_t, 6 > dofDataIdx;
   P1Form                  form( form_ );

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( { v0, v1, v2 }, elMat );

   // get global indices for local dofs
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   // add local contributions to diagonal entries
   if ( !lumped )
   {
      dstVertexData[dofDataIdx[0]] += elMat( 0, 0 ) * alpha;
      dstVertexData[dofDataIdx[1]] += elMat( 1, 1 ) * alpha;
      dstVertexData[dofDataIdx[2]] += elMat( 2, 2 ) * alpha;
   }
   else
   {
      dstVertexData[dofDataIdx[0]] += elMat.colwise().sum()[0] * alpha;
      dstVertexData[dofDataIdx[1]] += elMat.colwise().sum()[1] * alpha;
      dstVertexData[dofDataIdx[2]] += elMat.colwise().sum()[2] * alpha;
   }
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::computeLocalDiagonalContributions3D( const Cell&             cell,
                                                                           const uint_t            level,
                                                                           const indexing::Index&  microCell,
                                                                           const celldof::CellType cType,
                                                                           real_t* const           vertexData,
                                                                           const real_t&           alpha,
                                                                           bool                    lumped )
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix4r elMat( Matrix4r::Zero() );
   P1Form   form( form_ );
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   // add contributions for central stencil weights
   for ( int k = 0; k < 4; ++k )
   {
      if ( !lumped )
      {
         vertexData[vertexDoFIndices[uint_c( k )]] += elMat( k, k ) * alpha;
      }
      else
      {
         vertexData[vertexDoFIndices[uint_c( k )]] += elMat.colwise().sum()[k] * alpha;
      }
   }
}

// Assemble operator as sparse matrix
template < class P1Form >
void P1ElementwiseOperator< P1Form >::toMatrixScaled( const real_t&                               alpha,
                                                      const std::shared_ptr< SparseMatrixProxy >& mat,
                                                      const P1Function< idx_t >&                  src,
                                                      const P1Function< idx_t >&                  dst,
                                                      uint_t                                      level,
                                                      DoFType                                     flag ) const
{
   // We currently ignore the flag provided!
   // WALBERLA_UNUSED( flag );
   if ( flag != All )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in P1ElementwiseOperator::assembleLocalMatrix(); using flag = All" );
   }

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data in the two indexing functions
         PrimitiveDataID< FunctionMemory< idx_t >, Cell > dstVertexDoFIdx = dst.getCellDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Cell > srcVertexDoFIdx = src.getCellDataID();

         idx_t* srcIdx = cell.getData( srcVertexDoFIdx )->getPointer( level );
         idx_t* dstIdx = cell.getData( dstVertexDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixAssembly3D( mat, cell, level, micro, cType, srcIdx, dstIdx, alpha );
            }
         }
      }
   }

   else
   {
      // we only perform computations on face primitives
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         uint_t          rowsize       = levelinfo::num_microvertices_per_edge( level );
         uint_t          inner_rowsize = rowsize;
         idx_t           xIdx, yIdx;
         Point3D         v0, v1, v2;
         indexing::Index nodeIdx;
         indexing::Index offset;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< idx_t >, Face > dstVertexDoFIdx = dst.getFaceDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcVertexDoFIdx = src.getFaceDataID();

         idx_t* srcIndices = face.getData( srcVertexDoFIdx )->getPointer( level );
         idx_t* dstIndices = face.getData( dstVertexDoFIdx )->getPointer( level );

         // the explicit uint_c cast prevents a segfault in intel compiler 2018.4
         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < idx_t( rowsize ) - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < idx_t( inner_rowsize ) - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               localMatrixAssembly2D(
                   mat, face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementN, srcIndices, dstIndices, alpha );
               localMatrixAssembly2D(
                   mat, face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, srcIndices, dstIndices, alpha );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixAssembly2D(
                mat, face, level, xIdx, yIdx, P1Elements::P1Elements2D::elementNW, srcIndices, dstIndices, alpha );
         }

         // top north-west micro-element not treated, yet
         localMatrixAssembly2D( mat, face, level, 1, yIdx, P1Elements::P1Elements2D::elementNW, srcIndices, dstIndices, alpha );
      }
   }
}

// Assemble operator as sparse matrix
template < class P1Form >
void P1ElementwiseOperator< P1Form >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                const P1Function< idx_t >&                  src,
                                                const P1Function< idx_t >&                  dst,
                                                uint_t                                      level,
                                                DoFType                                     flag ) const
{
   return toMatrixScaled( real_c( 1 ), mat, src, dst, level, flag );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                             const Face&                                 face,
                                                             const uint_t                                level,
                                                             const idx_t                                 xIdx,
                                                             const idx_t                                 yIdx,
                                                             const P1Elements::P1Elements2D::P1Element&  element,
                                                             const idx_t* const                          srcIdx,
                                                             const idx_t* const                          dstIdx,
                                                             const real_t&                               alpha ) const
{
   Matrix3r                elMat( Matrix3r::Zero() );
   indexing::Index         nodeIdx;
   indexing::Index         offset;
   Point3D                 v0, v1, v2;
   std::array< uint_t, 3 > dofDataIdx;
   P1Form                  form( form_ );

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( { v0, v1, v2 }, elMat );

   // determine global indices of our local DoFs
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   std::vector< uint_t > rowIdx( 3 );
   rowIdx[0] = uint_c( dstIdx[dofDataIdx[0]] );
   rowIdx[1] = uint_c( dstIdx[dofDataIdx[1]] );
   rowIdx[2] = uint_c( dstIdx[dofDataIdx[2]] );

   std::vector< uint_t > colIdx( 3 );
   colIdx[0] = uint_c( srcIdx[dofDataIdx[0]] );
   colIdx[1] = uint_c( srcIdx[dofDataIdx[1]] );
   colIdx[2] = uint_c( srcIdx[dofDataIdx[2]] );

   const uint_t          elMatSize = 9;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i] * alpha;
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

template < class P1Form >
void P1ElementwiseOperator< P1Form >::localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                             const Cell&                                 cell,
                                                             const uint_t                                level,
                                                             const indexing::Index&                      microCell,
                                                             const celldof::CellType                     cType,
                                                             const idx_t* const                          srcIdx,
                                                             const idx_t* const                          dstIdx,
                                                             const real_t&                               alpha ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix4r elMat( Matrix4r::Zero() );
   P1Form   form( form_ );
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFDataIdx;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFDataIdx );

   std::vector< uint_t > rowIdx( 4 );
   std::vector< uint_t > colIdx( 4 );
   for ( uint_t k = 0; k < 4; ++k )
   {
      rowIdx[k] = uint_c( dstIdx[vertexDoFDataIdx[k]] );
      colIdx[k] = uint_c( srcIdx[vertexDoFDataIdx[k]] );
   }

   const uint_t          elMatSize = 16;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i] * alpha;
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

// P1ElementwiseLaplaceOperator
template class P1ElementwiseOperator<
    P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >;

// P1ElementwisePolarLaplaceOperator
template class P1ElementwiseOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > >;

// P1ElementwiseMassOperator
template class P1ElementwiseOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >;

// P1ElementwisePSPGOperator
template class P1ElementwiseOperator<
    P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >;

// P1ElementwiseBlendingPSPGOperator
template class P1ElementwiseOperator< forms::p1_pspg_blending_q2 >;

template class P1ElementwiseOperator< P1LinearCombinationForm >;

// P1ElementwiseBlendingMassOperator3D
template class P1ElementwiseOperator< forms::p1_mass_blending_q4 >;

// P1ElementwiseBlendingLaplaceOperator
template class P1ElementwiseOperator< forms::p1_diffusion_blending_q3 >;
template class P1ElementwiseOperator< forms::p1_diffusion_blending_q2 >;

// Needed for P1Blending(Inverse)DiagonalOperator
template class P1ElementwiseOperator< P1RowSumForm >;

template class P1ElementwiseOperator< forms::p1_div_k_grad_affine_q3 >;
template class P1ElementwiseOperator< forms::p1_div_k_grad_blending_q3 >;

template class P1ElementwiseOperator<
    P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >;
template class P1ElementwiseOperator<
    P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >;
template class P1ElementwiseOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > >;

template class P1ElementwiseOperator<
    P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >;
template class P1ElementwiseOperator<
    P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >;
template class P1ElementwiseOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > >;

template class P1ElementwiseOperator< forms::p1_epsiloncc_0_0_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_0_1_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_0_2_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_1_0_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_1_1_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_1_2_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_2_0_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_2_1_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsiloncc_2_2_affine_q2 >;

template class P1ElementwiseOperator< forms::p1_epsilonvar_0_0_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_0_1_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_0_2_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_1_0_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_1_1_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_1_2_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_2_0_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_2_1_affine_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_2_2_affine_q2 >;

template class P1ElementwiseOperator< forms::p1_epsilonvar_0_0_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_0_1_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_0_2_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_1_0_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_1_1_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_1_2_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_2_0_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_2_1_blending_q2 >;
template class P1ElementwiseOperator< forms::p1_epsilonvar_2_2_blending_q2 >;

template class P1ElementwiseOperator< forms::p1_k_mass_affine_q4 >;

template class P1ElementwiseOperator< forms::p1_k_mass_centroid_blending_q4 >;

template class P1ElementwiseOperator< forms::p1_neighbour_form >;

// This is a slight misuse of the P1ElementwiseOperator class, since the spherical
// elements are not P1. However, the SphericalElementFunction, like the P1Function
// is only an alias for the VertexDoFFunction, so we can re-use this operator.
template class P1ElementwiseOperator< SphericalElementFormMass >;

} // namespace hyteg
