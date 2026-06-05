/*
 * Copyright (c) 2017-2026 Marcus Mohr, Andreas Burkhart.
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

#include "P3ElementwiseOperator.hpp"

#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

namespace hyteg {

template < class P3FormHyTeG >
P3ElementwiseOperator< P3FormHyTeG >::P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                             size_t                                     minLevel,
                                                             size_t                                     maxLevel )
: P3ElementwiseOperator< P3FormHyTeG >( storage, minLevel, maxLevel, P3FormHyTeG(), true )
{}

template < class P3FormHyTeG >
P3ElementwiseOperator< P3FormHyTeG >::P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                             size_t                                     minLevel,
                                                             size_t                                     maxLevel,
                                                             const P3FormHyTeG&                         form )
: P3ElementwiseOperator< P3FormHyTeG >( storage, minLevel, maxLevel, form, true )
{}

template < class P3FormHyTeG >
P3ElementwiseOperator< P3FormHyTeG >::P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                             size_t                                     minLevel,
                                                             size_t                                     maxLevel,
                                                             bool                                       needsInverseDiagEntries )
: P3ElementwiseOperator< P3FormHyTeG >( storage, minLevel, maxLevel, P3FormHyTeG(), needsInverseDiagEntries )
{}

template < class P3FormHyTeG >
P3ElementwiseOperator< P3FormHyTeG >::P3ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                             size_t                                     minLevel,
                                                             size_t                                     maxLevel,
                                                             const P3FormHyTeG&                         form,
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

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::gemv( const real_t&               alpha,
                                                 const P3Function< real_t >& src,
                                                 const real_t&               beta,
                                                 const P3Function< real_t >& dst,
                                                 uint_t                      level,
                                                 DoFType                     flag ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   this->startTiming( "apply" );

   // Make sure that halos are up-to-date
   this->storage_->getTimingTree()->start( "sync source communication" );

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
      // communication::syncFunctionBetweenPrimitives( src, level, communication::syncDirection_t::LOW2HIGH );
   }
   this->storage_->getTimingTree()->stop( "sync source communication" );

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
      WALBERLA_ABORT( "P3ElementwiseOperator::gemv() not implemented for 3D!" );
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
         PrimitiveDataID< FunctionMemory< real_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Face > dstEdgeDoFIdxBlue = dst.getEdgeDoFFunctionBlue().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcEdgeDoFIdxBlue = src.getEdgeDoFFunctionBlue().getFaceDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Face > dstEdgeDoFIdxRed = dst.getEdgeDoFFunctionRed().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcEdgeDoFIdxRed = src.getEdgeDoFFunctionRed().getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeDataBlue = face.getData( srcEdgeDoFIdxBlue )->getPointer( level );
         real_t* dstEdgeDataBlue = face.getData( dstEdgeDoFIdxBlue )->getPointer( level );

         real_t* srcEdgeDataRed = face.getData( srcEdgeDoFIdxRed )->getPointer( level );
         real_t* dstEdgeDataRed = face.getData( dstEdgeDoFIdxRed )->getPointer( level );

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

         for ( const auto& idx : edgedof::macroface::Iterator( level ) )
         {
            for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
            {
               if ( !edgedof::macroface::isInnerEdgeDoF( level, idx, orientation ) )
               {
                  auto arrayIdx             = edgedof::macroface::index( level, idx.x(), idx.y(), orientation );
                  dstEdgeDataBlue[arrayIdx] = real_c( 0 );
                  dstEdgeDataRed[arrayIdx]  = real_c( 0 );
               }
            }
         }

         Matrix10r elMat = Matrix10r::Zero();

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
                  assembleLocalElementMatrix2D( face, level, micro, fType, form_, elMat );
               }

               localMatrixVectorMultiply2D( level,
                                            micro,
                                            fType,
                                            srcVertexData,
                                            srcEdgeDataBlue,
                                            srcEdgeDataRed,
                                            src.getFaceDoFFunction().dofMemory( it.first, level ),
                                            dstVertexData,
                                            dstEdgeDataBlue,
                                            dstEdgeDataRed,
                                            dst.getFaceDoFFunction().dofMemory( it.first, level ),
                                            elMat,
                                            alpha );
            }
         }
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.getVertexDoFFunction().communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      dst.getVertexDoFFunction().communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      dst.getEdgeDoFFunctionBlue().communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      dst.getEdgeDoFFunctionRed().communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
   }

   this->stopTiming( "apply" );
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::apply( const P3Function< real_t >& src,
                                                  const P3Function< real_t >& dst,
                                                  uint_t                      level,
                                                  DoFType                     flag,
                                                  UpdateType                  updateType ) const
{
   return gemv( real_c( 1 ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::applyScaled( const real_t&               alpha,
                                                        const P3Function< real_t >& src,
                                                        const P3Function< real_t >& dst,
                                                        uint_t                      level,
                                                        DoFType                     flag,
                                                        UpdateType                  updateType ) const
{
   return gemv( alpha, src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::smooth_jac_scaled( const real_t&               alpha,
                                                              const P3Function< real_t >& dst,
                                                              const P3Function< real_t >& rhs,
                                                              const P3Function< real_t >& src,
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

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::smooth_jac( const P3Function< real_t >& dst,
                                                       const P3Function< real_t >& rhs,
                                                       const P3Function< real_t >& src,
                                                       real_t                      omega,
                                                       size_t                      level,
                                                       DoFType                     flag ) const
{
   smooth_jac_scaled( real_c( 1 ), dst, rhs, src, omega, level, flag );
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::localMatrixVectorMultiply2D( uint_t                 level,
                                                                        const indexing::Index& microFace,
                                                                        facedof::FaceType      fType,
                                                                        const real_t* const    srcVertexData,
                                                                        const real_t* const    srcEdgeDataBlue,
                                                                        const real_t* const    srcEdgeDataRed,
                                                                        const real_t* const    srcFaceData,
                                                                        real_t* const          dstVertexData,
                                                                        real_t* const          dstEdgeDataBlue,
                                                                        real_t* const          dstEdgeDataRed,
                                                                        real_t* const          dstFaceData,
                                                                        const Matrix10r&       elMat,
                                                                        const real_t&          alpha ) const
{
   // obtain data indices of dofs associated with micro-face
   std::array< uint_t, 3 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, fType, level, vertexDoFIndices );

   std::array< uint_t, 3 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering< true >( microFace, fType, level, edgeDoFIndices );

   uint_t faceDoFIdx = volumedofspace::indexing::index(
       microFace.x(), microFace.y(), fType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::AoS );

   // assemble local element vector
   Point10D elVecOld, elVecNew;
   for ( int k = 0; k < 3; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
   }
   for ( int k = 3; k < 6; ++k )
   {
      elVecOld[k] = srcEdgeDataBlue[edgeDoFIndices[uint_c( k - 3 )]];
   }
   for ( int k = 6; k < 9; ++k )
   {
      elVecOld[k] = srcEdgeDataRed[edgeDoFIndices[uint_c( k - 6 )]];
   }
   elVecOld[9] = srcFaceData[faceDoFIdx];

   // apply matrix (operator locally)
   elVecNew = alpha * ( elMat * elVecOld );

   // redistribute result from "local" to "global vector"
   for ( int k = 0; k < 3; ++k )
   {
      dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
   }
   for ( int k = 3; k < 6; ++k )
   {
      dstEdgeDataBlue[edgeDoFIndices[uint_c( k - 3 )]] += elVecNew[k];
   }
   for ( int k = 6; k < 9; ++k )
   {
      dstEdgeDataRed[edgeDoFIndices[uint_c( k - 6 )]] += elVecNew[k];
   }
   dstFaceData[faceDoFIdx] += elVecNew[9];
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::computeDiagonalOperatorValues( bool invert, bool lumped, const real_t& alpha )
{
   std::shared_ptr< P3Function< real_t > > targetFunction;

   if ( invert )
   {
      if ( !lumped )
      {
         if ( !inverseDiagonalValues_ )
         {
            inverseDiagonalValues_ =
                std::make_shared< P3Function< real_t > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
         }
         targetFunction = inverseDiagonalValues_;
      }
      else
      {
         if ( !lumpedInverseDiagonalValues_ )
         {
            lumpedInverseDiagonalValues_ =
                std::make_shared< P3Function< real_t > >( "lumped inverse diagonal entries", storage_, minLevel_, maxLevel_ );
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
            diagonalValues_ = std::make_shared< P3Function< real_t > >( "diagonal entries", storage_, minLevel_, maxLevel_ );
         }
         targetFunction = diagonalValues_;
      }
      else
      {
         if ( !lumpedDiagonalValues_ )
         {
            lumpedDiagonalValues_ =
                std::make_shared< P3Function< real_t > >( "lumped diagonal entries", storage_, minLevel_, maxLevel_ );
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
         WALBERLA_ABORT( "P3ElementwiseOperator::computeDiagonalOperatorValues() not implemented for 3D!" );
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
            PrimitiveDataID< FunctionMemory< real_t >, Face > vertexDoFIdx =
                targetFunction->getVertexDoFFunction().getFaceDataID();
            PrimitiveDataID< FunctionMemory< real_t >, Face > edgeDoFIdxBlue =
                targetFunction->getEdgeDoFFunctionBlue().getFaceDataID();
            PrimitiveDataID< FunctionMemory< real_t >, Face > edgeDoFIdxRed =
                targetFunction->getEdgeDoFFunctionRed().getFaceDataID();

            real_t* vertexData   = face.getData( vertexDoFIdx )->getPointer( level );
            real_t* edgeDataBlue = face.getData( edgeDoFIdxBlue )->getPointer( level );
            real_t* edgeDataRed  = face.getData( edgeDoFIdxRed )->getPointer( level );

            // now loop over micro-faces of macro-face
            for ( const auto& faceType : facedof::allFaceTypes )
            {
               for ( const auto& microFace : facedof::macroface::Iterator( level, faceType, 0 ) )
               {
                  computeLocalDiagonalContributions2D( face,
                                                       level,
                                                       microFace,
                                                       faceType,
                                                       vertexData,
                                                       edgeDataBlue,
                                                       edgeDataRed,
                                                       targetFunction->getFaceDoFFunction().dofMemory( it.first, level ),
                                                       alpha,
                                                       lumped );
               }
            }
         }

         // Push result to lower-dimensional primitives
         targetFunction->getVertexDoFFunction().communicateAdditively< Face, Edge >( level );
         targetFunction->getVertexDoFFunction().communicateAdditively< Face, Vertex >( level );
         targetFunction->getEdgeDoFFunctionBlue().communicateAdditively< Face, Edge >( level );
         targetFunction->getEdgeDoFFunctionRed().communicateAdditively< Face, Edge >( level );

         // Retrieve assembled data values
         targetFunction->getVertexDoFFunction().communicate< Vertex, Edge >( level );
         targetFunction->getVertexDoFFunction().communicate< Edge, Face >( level );
         targetFunction->getEdgeDoFFunctionBlue().communicate< Edge, Face >( level );
         targetFunction->getEdgeDoFFunctionRed().communicate< Edge, Face >( level );
      }

      // Invert values if desired (note: using false below means we only invert in the interior of the primitives,
      // the values in the halos are untouched; should be okay for using diagonalValue_ in smoothers)
      if ( invert )
      {
         targetFunction->invertElementwise( level, All, false );
      }
   }
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::computeAndStoreLocalElementMatrices()
{
   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      // For 3D we work on cells and for 2D on faces
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P3ElementwiseOperator::computeAndStoreLocalElementMatrices() not implemented for 3D!" );
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
                  Matrix10r& elMat = localElementMatrix2D( *face, level, micro, fType );
                  elMat.setZero();
                  assembleLocalElementMatrix2D( *face, level, micro, fType, form_, elMat );
               }
            }
         }
      }
   }

   localElementMatricesPrecomputed_ = true;
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::computeLocalDiagonalContributions2D( const Face&             face,
                                                                                const uint_t            level,
                                                                                const indexing::Index&  microFace,
                                                                                const facedof::FaceType faceType,
                                                                                real_t* const           dstVertexData,
                                                                                real_t* const           dstEdgeDataBlue,
                                                                                real_t* const           dstEdgeDataRed,
                                                                                real_t* const           dstFaceData,
                                                                                const real_t&           alpha,
                                                                                bool                    lumped )
{
   // obtain coordinates of vertices of given micro-face
   std::array< indexing::Index, 3 > microVertexIndex =
       facedof::macroface::getMicroVerticesFromMicroFace< true >( microFace, faceType );
   std::array< Point3D, 3 > coords;
   coords[0] = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndex[0] );
   coords[1] = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndex[1] );
   coords[2] = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndex[2] );

   // assemble local element matrix
   Matrix10r elMat = Matrix10r::Zero();
   form_.setGeometryMap( face.getGeometryMap() );
   form_.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-face
   std::array< uint_t, 3 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndices );

   std::array< uint_t, 3 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndices );

   uint_t faceDoFIdx = volumedofspace::indexing::index(
       microFace.x(), microFace.y(), faceType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::AoS );

   // add local contributions to diagonal entries
   if ( !lumped )
   {
      dstVertexData[vertexDoFIndices[0]] += elMat( 0, 0 ) * alpha;
      dstVertexData[vertexDoFIndices[1]] += elMat( 1, 1 ) * alpha;
      dstVertexData[vertexDoFIndices[2]] += elMat( 2, 2 ) * alpha;

      dstEdgeDataBlue[edgeDoFIndices[0]] += elMat( 3, 3 ) * alpha;
      dstEdgeDataBlue[edgeDoFIndices[1]] += elMat( 4, 4 ) * alpha;
      dstEdgeDataBlue[edgeDoFIndices[2]] += elMat( 5, 5 ) * alpha;

      dstEdgeDataRed[edgeDoFIndices[0]] += elMat( 6, 6 ) * alpha;
      dstEdgeDataRed[edgeDoFIndices[1]] += elMat( 7, 7 ) * alpha;
      dstEdgeDataRed[edgeDoFIndices[2]] += elMat( 8, 8 ) * alpha;

      dstFaceData[faceDoFIdx] += elMat( 9, 9 ) * alpha;
   }
   else
   {
      dstVertexData[vertexDoFIndices[0]] += elMat.colwise().sum()[0] * alpha;
      dstVertexData[vertexDoFIndices[1]] += elMat.colwise().sum()[1] * alpha;
      dstVertexData[vertexDoFIndices[2]] += elMat.colwise().sum()[2] * alpha;

      dstEdgeDataBlue[edgeDoFIndices[0]] += elMat.colwise().sum()[3] * alpha;
      dstEdgeDataBlue[edgeDoFIndices[1]] += elMat.colwise().sum()[4] * alpha;
      dstEdgeDataBlue[edgeDoFIndices[2]] += elMat.colwise().sum()[5] * alpha;

      dstEdgeDataRed[edgeDoFIndices[0]] += elMat.colwise().sum()[6] * alpha;
      dstEdgeDataRed[edgeDoFIndices[1]] += elMat.colwise().sum()[7] * alpha;
      dstEdgeDataRed[edgeDoFIndices[2]] += elMat.colwise().sum()[8] * alpha;

      dstFaceData[faceDoFIdx] += elMat.colwise().sum()[9] * alpha;
   }
}

template < class P3FormHyTeG >
P3FormHyTeG P3ElementwiseOperator< P3FormHyTeG >::getForm() const
{
   return form_;
}

// Assemble operator as scaled sparse matrix
template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::toMatrixScaled( const real_t&                               alpha,
                                                           const std::shared_ptr< SparseMatrixProxy >& mat,
                                                           const P3Function< idx_t >&                  src,
                                                           const P3Function< idx_t >&                  dst,
                                                           uint_t                                      level,
                                                           DoFType                                     flag ) const
{
   // We currently ignore the flag provided!
   // WALBERLA_UNUSED( flag );
   if ( flag != All )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in P3ElementwiseOperator::toMatrixScaled(); using flag = All" );
   }

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         WALBERLA_ABORT( "P3ElementwiseOperator::toMatrixScaled() not implemented for 3D!" );
      }
   }

   else
   {
      // we only perform computations on face primitives
      for ( auto& it : storage_->getFaces() )
      {
         const PrimitiveID& faceID = it.first;
         Face&              face   = *it.second;

         uint_t          rowsize       = levelinfo::num_microvertices_per_edge( level );
         uint_t          inner_rowsize = rowsize;
         idx_t           xIdx, yIdx;
         Point3D         v0, v1, v2;
         indexing::Index nodeIdx;
         indexing::Index offset;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< idx_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< idx_t >, Face > dstEdgeDoFIdxBlue = dst.getEdgeDoFFunctionBlue().getFaceDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcEdgeDoFIdxBlue = src.getEdgeDoFFunctionBlue().getFaceDataID();

         PrimitiveDataID< FunctionMemory< idx_t >, Face > dstEdgeDoFIdxRed = dst.getEdgeDoFFunctionRed().getFaceDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcEdgeDoFIdxRed = src.getEdgeDoFFunctionRed().getFaceDataID();

         // obtain pointers to linalg index data (ordering: src-vertex, src-edge-blue, src-edge-red, src-face,
         // dst-vertex, dst-edge-blue, dst-edge-red, dst-face)
         std::array< const idx_t*, 8 > indexDataPointer;

         indexDataPointer[0] = face.getData( srcVertexDoFIdx )->getPointer( level );
         indexDataPointer[1] = face.getData( srcEdgeDoFIdxBlue )->getPointer( level );
         indexDataPointer[2] = face.getData( srcEdgeDoFIdxRed )->getPointer( level );
         indexDataPointer[3] = src.getFaceDoFFunction().dofMemory( faceID, level ),

         indexDataPointer[4] = face.getData( dstVertexDoFIdx )->getPointer( level );
         indexDataPointer[5] = face.getData( dstEdgeDoFIdxBlue )->getPointer( level );
         indexDataPointer[6] = face.getData( dstEdgeDoFIdxRed )->getPointer( level );
         indexDataPointer[7] = dst.getFaceDoFFunction().dofMemory( faceID, level );

         // now loop over micro-faces of macro-face
         for ( const auto& faceType : facedof::allFaceTypes )
         {
            for ( const auto& microFace : facedof::macroface::Iterator( level, faceType, 0 ) )
            {
               localMatrixAssembly2D( mat, face, level, microFace, faceType, indexDataPointer, alpha );
            }
         }
      }
   }
}

// Assemble operator as sparse matrix
template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                     const P3Function< idx_t >&                  src,
                                                     const P3Function< idx_t >&                  dst,
                                                     uint_t                                      level,
                                                     DoFType                                     flag ) const
{
   return toMatrixScaled( real_c( 1 ), mat, src, dst, level, flag );
}

template < class P3FormHyTeG >
void P3ElementwiseOperator< P3FormHyTeG >::localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                  const Face&                                 face,
                                                                  const uint_t                                level,
                                                                  const indexing::Index&                      microFace,
                                                                  const facedof::FaceType                     faceType,
                                                                  std::array< const idx_t*, 8 >               indexDataPointer,
                                                                  const real_t&                               alpha ) const
{
   // indexDataPointer: pointers to linalg index data (ordering: src-vertex, src-edge-blue, src-edge-red, src-face,
   // dst-vertex, dst-edge-blue, dst-edge-red, dst-face)

   // obtain coordinates of vertices of given micro-face
   std::array< indexing::Index, 3 > microVertexIndex =
       facedof::macroface::getMicroVerticesFromMicroFace< true >( microFace, faceType );
   std::array< Point3D, 3 > coords;
   coords[0] = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndex[0] );
   coords[1] = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndex[1] );
   coords[2] = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndex[2] );

   // assemble local element matrix
   Matrix10r elMat = Matrix10r::Zero();
   form_.setGeometryMap( face.getGeometryMap() );
   form_.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-face
   std::array< uint_t, 3 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndices );

   std::array< uint_t, 3 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndices );

   uint_t faceDoFIdx = volumedofspace::indexing::index(
       microFace.x(), microFace.y(), faceType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::AoS );

   // determine row indices in global matrix
   std::vector< uint_t > rowIdx( 10 );

   rowIdx[0] = uint_c( indexDataPointer[0][vertexDoFIndices[0]] );
   rowIdx[1] = uint_c( indexDataPointer[0][vertexDoFIndices[1]] );
   rowIdx[2] = uint_c( indexDataPointer[0][vertexDoFIndices[2]] );

   rowIdx[3] = uint_c( indexDataPointer[1][edgeDoFIndices[0]] );
   rowIdx[4] = uint_c( indexDataPointer[1][edgeDoFIndices[1]] );
   rowIdx[5] = uint_c( indexDataPointer[1][edgeDoFIndices[2]] );

   rowIdx[6] = uint_c( indexDataPointer[2][edgeDoFIndices[0]] );
   rowIdx[7] = uint_c( indexDataPointer[2][edgeDoFIndices[1]] );
   rowIdx[8] = uint_c( indexDataPointer[2][edgeDoFIndices[2]] );

   rowIdx[9] = uint_c( indexDataPointer[3][faceDoFIdx] );

   // determine column indices in global matrix
   std::vector< uint_t > colIdx( 10 );

   colIdx[0] = uint_c( indexDataPointer[4][vertexDoFIndices[0]] );
   colIdx[1] = uint_c( indexDataPointer[4][vertexDoFIndices[1]] );
   colIdx[2] = uint_c( indexDataPointer[4][vertexDoFIndices[2]] );

   colIdx[3] = uint_c( indexDataPointer[5][edgeDoFIndices[0]] );
   colIdx[4] = uint_c( indexDataPointer[5][edgeDoFIndices[1]] );
   colIdx[5] = uint_c( indexDataPointer[5][edgeDoFIndices[2]] );

   colIdx[6] = uint_c( indexDataPointer[6][edgeDoFIndices[0]] );
   colIdx[7] = uint_c( indexDataPointer[6][edgeDoFIndices[1]] );
   colIdx[8] = uint_c( indexDataPointer[6][edgeDoFIndices[2]] );

   colIdx[9] = uint_c( indexDataPointer[7][faceDoFIdx] );

   // format conversion and (potential) scaling of element matrix
   const uint_t          elMatSize = 100;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i] * alpha;
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

// P3ElementwiseMassOperators
template class P3ElementwiseOperator< forms::p3_mass_affine_qe >;
template class P3ElementwiseOperator< forms::p3_mass_blending_q6 >;
template class P3ElementwiseOperator< forms::p3_diffusion_affine_q4 >;

} // namespace hyteg
