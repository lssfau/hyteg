/*
 * Copyright (c) 2017-2019 Marcus Mohr.
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

#include "P2ElementwiseOperator.hpp"

namespace hyteg {

template < class P2Form >
P2ElementwiseOperator< P2Form >::P2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                        size_t                                     minLevel,
                                                        size_t                                     maxLevel )
: P2ElementwiseOperator< P2Form >( storage, minLevel, maxLevel, P2Form(), true )
{}

template < class P2Form >
P2ElementwiseOperator< P2Form >::P2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                        size_t                                     minLevel,
                                                        size_t                                     maxLevel,
                                                        const P2Form&                              form )
: P2ElementwiseOperator< P2Form >( storage, minLevel, maxLevel, form, true )
{}

template < class P2Form >
P2ElementwiseOperator< P2Form >::P2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                        size_t                                     minLevel,
                                                        size_t                                     maxLevel,
                                                        const P2Form&                              form,
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

void localMatrixVectorMultiply3D( uint_t                 level,
                                  const indexing::Index& microCell,
                                  celldof::CellType      cType,
                                  const real_t* const    srcVertexData,
                                  const real_t* const    srcEdgeData,
                                  real_t* const          dstVertexData,
                                  real_t* const          dstEdgeData,
                                  const Matrix10r&       elMat )
{
   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // assemble local element vector
   Point10D elVecOld, elVecNew;
   for ( uint_t k = 0; k < 4; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[k]];
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      elVecOld[k] = srcEdgeData[edgeDoFIndices[k - 4]];
   }

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   for ( uint_t k = 0; k < 4; ++k )
   {
      dstVertexData[vertexDoFIndices[k]] += elVecNew[k];
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      dstEdgeData[edgeDoFIndices[k - 4]] += elVecNew[k];
   }
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::apply( const P2Function< real_t >& src,
                                             const P2Function< real_t >& dst,
                                             size_t                      level,
                                             DoFType                     flag,
                                             UpdateType                  updateType ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   this->startTiming( "apply" );

   this->storage_->getTimingTree()->start( "sync source communication" );
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
      communication::syncP2FunctionBetweenPrimitives( src, level );
   }
   this->storage_->getTimingTree()->stop( "sync source communication" );

   if ( updateType == Replace )
   {
      // We need to zero the destination array (including halos).
      // However, we must not zero out anything that is not flagged with the specified BCs.
      // Therefore we first zero out everything that flagged, and then, later,
      // the halos of the highest dim primitives.

      dst.interpolate( real_c( 0 ), level, flag );
   }

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< real_t >, Cell > dstVertexDoFIdx = dst.getVertexDoFFunction().getCellDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcVertexDoFIdx = src.getVertexDoFFunction().getCellDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Cell > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getCellDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcEdgeDoFIdx = src.getEdgeDoFFunction().getCellDataID();

         real_t* srcVertexData = cell.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = cell.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeData = cell.getData( srcEdgeDoFIdx )->getPointer( level );
         real_t* dstEdgeData = cell.getData( dstEdgeDoFIdx )->getPointer( level );

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

         for ( const auto& idx : edgedof::macrocell::Iterator( level ) )
         {
            for ( const auto& orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
            {
               if ( !edgedof::macrocell::isInnerEdgeDoF( level, idx, orientation ) )
               {
                  auto arrayIdx         = edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), orientation );
                  dstEdgeData[arrayIdx] = real_c( 0 );
               }
            }
         }

         Matrix10r elMat;

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
                  assembleLocalElementMatrix3D( cell, level, micro, cType, form_, elMat );
               }

               localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, srcEdgeData, dstVertexData, dstEdgeData, elMat );
            }
         }
      }

      this->storage_->getTimingTree()->start( "additive communication" );
      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.getVertexDoFFunction().communicateAdditively< Cell, Face >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getVertexDoFFunction().communicateAdditively< Cell, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getVertexDoFFunction().communicateAdditively< Cell, Vertex >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getEdgeDoFFunction().communicateAdditively< Cell, Face >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getEdgeDoFFunction().communicateAdditively< Cell, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      this->storage_->getTimingTree()->stop( "additive communication" );
   }

   else
   {
      // we only perform computations on face primitives
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         Point3D x0( face.coords[0] );
         Point3D x1( face.coords[1] );
         Point3D x2( face.coords[2] );

         Point3D                  v0, v1, v2;
         indexing::Index          nodeIdx;
         indexing::IndexIncrement offset;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< real_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Face > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcEdgeDoFIdx = src.getEdgeDoFFunction().getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeData = face.getData( srcEdgeDoFIdx )->getPointer( level );
         real_t* dstEdgeData = face.getData( dstEdgeDoFIdx )->getPointer( level );

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
                  auto arrayIdx         = edgedof::macroface::index( level, idx.x(), idx.y(), orientation );
                  dstEdgeData[arrayIdx] = real_c( 0 );
               }
            }
         }

         Matrix6r elMat;

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

               localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, srcEdgeData, dstVertexData, dstEdgeData, elMat );
            }
         }
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.getVertexDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getVertexDoFFunction().communicateAdditively< Face, Vertex >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.getEdgeDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
   }

   this->stopTiming( "apply" );
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::smooth_jac( const P2Function< real_t >& dst,
                                                  const P2Function< real_t >& rhs,
                                                  const P2Function< real_t >& src,
                                                  real_t                      omega,
                                                  size_t                      level,
                                                  DoFType                     flag ) const
{
   this->startTiming( "smooth_jac" );

   // compute the current residual
   this->apply( src, dst, level, flag );
   dst.assign( {real_c( 1 ), real_c( -1 )}, {rhs, dst}, level, flag );

   // perform Jacobi update step
   dst.multElementwise( {*getInverseDiagonalValues(), dst}, level, flag );
   dst.assign( {1.0, omega}, {src, dst}, level, flag );

   this->stopTiming( "smooth_jac" );
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::localMatrixVectorMultiply2D( uint_t                 level,
                                                                   const indexing::Index& microFace,
                                                                   facedof::FaceType      fType,
                                                                   const real_t* const    srcVertexData,
                                                                   const real_t* const    srcEdgeData,
                                                                   real_t* const          dstVertexData,
                                                                   real_t* const          dstEdgeData,
                                                                   const Matrix6r&        elMat ) const
{
   // obtain data indices of dofs associated with micro-face
   std::array< uint_t, 3 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, fType, level, vertexDoFIndices );

   std::array< uint_t, 3 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, fType, level, edgeDoFIndices );

   // assemble local element vector
   Point6D elVecOld, elVecNew;
   for ( uint_t k = 0; k < 3; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[k]];
   }
   for ( uint_t k = 3; k < 6; ++k )
   {
      elVecOld[k] = srcEdgeData[edgeDoFIndices[k - 3]];
   }

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   for ( uint_t k = 0; k < 3; ++k )
   {
      dstVertexData[vertexDoFIndices[k]] += elVecNew[k];
   }
   for ( uint_t k = 3; k < 6; ++k )
   {
      dstEdgeData[edgeDoFIndices[k - 3]] += elVecNew[k];
   }
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::computeDiagonalOperatorValues( bool invert )
{
   std::shared_ptr< P2Function< real_t > > targetFunction;
   if ( invert )
   {
      if ( !inverseDiagonalValues_ )
      {
         inverseDiagonalValues_ =
             std::make_shared< P2Function< real_t > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = inverseDiagonalValues_;
   }
   else
   {
      if ( !diagonalValues_ )
      {
         diagonalValues_ = std::make_shared< P2Function< real_t > >( "diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = diagonalValues_;
   }

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      // Make sure that halos are up-to-date (can we improve communication here?)
      communication::syncP2FunctionBetweenPrimitives( *targetFunction, level );

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
            PrimitiveDataID< FunctionMemory< real_t >, Cell > diagVertexDoFIdx =
                targetFunction->getVertexDoFFunction().getCellDataID();
            PrimitiveDataID< FunctionMemory< real_t >, Cell > diagEdgeDoFIdx =
                targetFunction->getEdgeDoFFunction().getCellDataID();

            real_t* diagVertexData = cell.getData( diagVertexDoFIdx )->getPointer( level );
            real_t* diagEdgeData   = cell.getData( diagEdgeDoFIdx )->getPointer( level );

            // loop over micro-cells
            for ( const auto& cType : celldof::allCellTypes )
            {
               for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
               {
                  computeLocalDiagonalContributions3D( cell, level, micro, cType, diagVertexData, diagEdgeData );
               }
            }
         }

         // Push result to lower-dimensional primitives
         targetFunction->getVertexDoFFunction().communicateAdditively< Cell, Face >( level );
         targetFunction->getVertexDoFFunction().communicateAdditively< Cell, Edge >( level );
         targetFunction->getVertexDoFFunction().communicateAdditively< Cell, Vertex >( level );
         targetFunction->getEdgeDoFFunction().communicateAdditively< Cell, Face >( level );
         targetFunction->getEdgeDoFFunction().communicateAdditively< Cell, Edge >( level );
      }

      else
      {
         // we only perform computations on face primitives
         for ( auto& it : storage_->getFaces() )
         {
            Face& face = *it.second;

            Point3D x0( face.coords[0] );
            Point3D x1( face.coords[1] );
            Point3D x2( face.coords[2] );

            uint_t                   rowsize       = levelinfo::num_microvertices_per_edge( level );
            uint_t                   inner_rowsize = rowsize;
            uint_t                   xIdx, yIdx;
            Point3D                  v0, v1, v2;
            indexing::Index          nodeIdx;
            indexing::IndexIncrement offset;

            // get hold of the actual numerical data in the two functions
            PrimitiveDataID< FunctionMemory< real_t >, Face > vertexDoFIdx =
                targetFunction->getVertexDoFFunction().getFaceDataID();
            PrimitiveDataID< FunctionMemory< real_t >, Face > edgeDoFIdx = targetFunction->getEdgeDoFFunction().getFaceDataID();

            real_t* vertexData = face.getData( vertexDoFIdx )->getPointer( level );
            real_t* edgeData   = face.getData( edgeDoFIdx )->getPointer( level );

            // now loop over micro-faces of macro-face
            for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx )
            {
               // loop over vertices in row with two associated triangles
               for ( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx )
               {
                  // we associate two elements with current micro-vertex
                  computeLocalDiagonalContributions2D(
                      face, level, xIdx, yIdx, P2Elements::P2Face::elementN, vertexData, edgeData );
                  computeLocalDiagonalContributions2D(
                      face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, vertexData, edgeData );
               }
               --inner_rowsize;

               // final micro-vertex in row has only one associated micro-face
               computeLocalDiagonalContributions2D(
                   face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, vertexData, edgeData );
            }

            // top north-west micro-element not treated, yet
            computeLocalDiagonalContributions2D( face, level, 1, yIdx, P2Elements::P2Face::elementNW, vertexData, edgeData );
         }

         // Push result to lower-dimensional primitives
         targetFunction->getVertexDoFFunction().communicateAdditively< Face, Edge >( level );
         targetFunction->getVertexDoFFunction().communicateAdditively< Face, Vertex >( level );
         targetFunction->getEdgeDoFFunction().communicateAdditively< Face, Edge >( level );

         // Retrieve assembled data values
         targetFunction->getVertexDoFFunction().communicate< Vertex, Edge >( level );
         targetFunction->getVertexDoFFunction().communicate< Edge, Face >( level );
         targetFunction->getEdgeDoFFunction().communicate< Edge, Face >( level );
      }

      // Invert values if desired (note: using false below means we only invert in the interior of the primitives,
      // the values in the halos are untouched; should be okay for using diagonalValue_ in smoothers)
      if ( invert )
      {
         targetFunction->invertElementwise( level, All, false );
      }
   }
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::computeAndStoreLocalElementMatrices()
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
                  Matrix10r& elMat = localElementMatrix3D( *cell, level, micro, cType );
                  elMat.setAll( 0 );
                  assembleLocalElementMatrix3D( *cell, level, micro, cType, form_, elMat );
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
                  Matrix6r& elMat = localElementMatrix2D( *face, level, micro, fType );
                  elMat.setAll( 0 );
                  assembleLocalElementMatrix2D( *face, level, micro, fType, form_, elMat );
               }
            }
         }
      }
   }

   localElementMatricesPrecomputed_ = true;
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::computeLocalDiagonalContributions2D( const Face&                  face,
                                                                           const uint_t                 level,
                                                                           const uint_t                 xIdx,
                                                                           const uint_t                 yIdx,
                                                                           const P2Elements::P2Element& element,
                                                                           real_t* const                dstVertexData,
                                                                           real_t* const                dstEdgeData )
{
   Matrix6r                 elMat;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;
   P2Form                   form( form_ );

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( {v0, v1, v2}, elMat );

   // get global indices for local dofs
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   dofDataIdx[3] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[4] );
   dofDataIdx[4] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[5] );
   dofDataIdx[5] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[3] );

   // add local contributions to diagonal entries
   dstVertexData[dofDataIdx[0]] += elMat( 0, 0 );
   dstVertexData[dofDataIdx[1]] += elMat( 1, 1 );
   dstVertexData[dofDataIdx[2]] += elMat( 2, 2 );

   dstEdgeData[dofDataIdx[3]] += elMat( 4, 4 );
   dstEdgeData[dofDataIdx[4]] += elMat( 5, 5 );
   dstEdgeData[dofDataIdx[5]] += elMat( 3, 3 );
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::computeLocalDiagonalContributions3D( const Cell&             cell,
                                                                           const uint_t            level,
                                                                           const indexing::Index&  microCell,
                                                                           const celldof::CellType cType,
                                                                           real_t* const           vertexData,
                                                                           real_t* const           edgeData )
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix10r elMat;
   P2Form    form( form_ );
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // add contributions for central stencil weights
   for ( uint_t k = 0; k < 4; ++k )
   {
      vertexData[vertexDoFIndices[k]] += elMat( k, k );
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      edgeData[edgeDoFIndices[k - 4]] += elMat( k, k );
   }
}
template < class P2Form >
P2Form P2ElementwiseOperator< P2Form >::getForm() const
{
   return form_;
}

// Assemble operator as sparse matrix
template < class P2Form >
void P2ElementwiseOperator< P2Form >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                const P2Function< idx_t >&                  src,
                                                const P2Function< idx_t >&                  dst,
                                                uint_t                                      level,
                                                DoFType                                     flag ) const
{
   // We currently ignore the flag provided!
   // WALBERLA_UNUSED( flag );
   if ( flag != All )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in P2ElementwiseOperator::assembleLocalMatrix(); using flag = All" );
   }

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data in the two indexing functions
         PrimitiveDataID< FunctionMemory< idx_t >, Cell > dstVertexDoFIdx = dst.getVertexDoFFunction().getCellDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Cell > srcVertexDoFIdx = src.getVertexDoFFunction().getCellDataID();

         PrimitiveDataID< FunctionMemory< idx_t >, Cell > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getCellDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Cell > srcEdgeDoFIdx = src.getEdgeDoFFunction().getCellDataID();

         idx_t* srcVertexIndices = cell.getData( srcVertexDoFIdx )->getPointer( level );
         idx_t* dstVertexIndices = cell.getData( dstVertexDoFIdx )->getPointer( level );

         idx_t* srcEdgeIndices = cell.getData( srcEdgeDoFIdx )->getPointer( level );
         idx_t* dstEdgeIndices = cell.getData( dstEdgeDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixAssembly3D(
                   mat, cell, level, micro, cType, srcVertexIndices, srcEdgeIndices, dstVertexIndices, dstEdgeIndices );
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

         Point3D x0( face.coords[0] );
         Point3D x1( face.coords[1] );
         Point3D x2( face.coords[2] );

         uint_t                   rowsize       = levelinfo::num_microvertices_per_edge( level );
         uint_t                   inner_rowsize = rowsize;
         uint_t                   xIdx, yIdx;
         Point3D                  v0, v1, v2;
         indexing::Index          nodeIdx;
         indexing::IndexIncrement offset;

         // get hold of the actual numerical data in the two functions
         PrimitiveDataID< FunctionMemory< idx_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< idx_t >, Face > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcEdgeDoFIdx = src.getEdgeDoFFunction().getFaceDataID();

         idx_t* srcVertexIndices = face.getData( srcVertexDoFIdx )->getPointer( level );
         idx_t* dstVertexIndices = face.getData( dstVertexDoFIdx )->getPointer( level );

         idx_t* srcEdgeIndices = face.getData( srcEdgeDoFIdx )->getPointer( level );
         idx_t* dstEdgeIndices = face.getData( dstEdgeDoFIdx )->getPointer( level );

         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               localMatrixAssembly2D( mat,
                                      face,
                                      level,
                                      xIdx,
                                      yIdx,
                                      P2Elements::P2Face::elementN,
                                      srcVertexIndices,
                                      srcEdgeIndices,
                                      dstVertexIndices,
                                      dstEdgeIndices );
               localMatrixAssembly2D( mat,
                                      face,
                                      level,
                                      xIdx,
                                      yIdx,
                                      P2Elements::P2Face::elementNW,
                                      srcVertexIndices,
                                      srcEdgeIndices,
                                      dstVertexIndices,
                                      dstEdgeIndices );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixAssembly2D( mat,
                                   face,
                                   level,
                                   xIdx,
                                   yIdx,
                                   P2Elements::P2Face::elementNW,
                                   srcVertexIndices,
                                   srcEdgeIndices,
                                   dstVertexIndices,
                                   dstEdgeIndices );
         }

         // top north-west micro-element not treated, yet
         localMatrixAssembly2D( mat,
                                face,
                                level,
                                1,
                                yIdx,
                                P2Elements::P2Face::elementNW,
                                srcVertexIndices,
                                srcEdgeIndices,
                                dstVertexIndices,
                                dstEdgeIndices );
      }
   }
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                             const Face&                                 face,
                                                             const uint_t                                level,
                                                             const uint_t                                xIdx,
                                                             const uint_t                                yIdx,
                                                             const P2Elements::P2Element&                element,
                                                             const idx_t* const                          srcVertexIdx,
                                                             const idx_t* const                          srcEdgeIdx,
                                                             const idx_t* const                          dstVertexIdx,
                                                             const idx_t* const                          dstEdgeIdx ) const
{
   Matrix6r                 elMat;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;
   P2Form                   form( form_ );

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( {v0, v1, v2}, elMat );

   // determine global indices of our local DoFs (note the tweaked ordering to go along with FEniCS indexing)
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   dofDataIdx[3] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[4] );
   dofDataIdx[4] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[5] );
   dofDataIdx[5] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[3] );

   std::vector< uint_t > rowIdx( 6 );
   rowIdx[0] = uint_c( dstVertexIdx[dofDataIdx[0]] );
   rowIdx[1] = uint_c( dstVertexIdx[dofDataIdx[1]] );
   rowIdx[2] = uint_c( dstVertexIdx[dofDataIdx[2]] );

   rowIdx[3] = uint_c( dstEdgeIdx[dofDataIdx[3]] );
   rowIdx[4] = uint_c( dstEdgeIdx[dofDataIdx[4]] );
   rowIdx[5] = uint_c( dstEdgeIdx[dofDataIdx[5]] );

   std::vector< uint_t > colIdx( 6 );
   colIdx[0] = uint_c( srcVertexIdx[dofDataIdx[0]] );
   colIdx[1] = uint_c( srcVertexIdx[dofDataIdx[1]] );
   colIdx[2] = uint_c( srcVertexIdx[dofDataIdx[2]] );

   colIdx[3] = uint_c( srcEdgeIdx[dofDataIdx[3]] );
   colIdx[4] = uint_c( srcEdgeIdx[dofDataIdx[4]] );
   colIdx[5] = uint_c( srcEdgeIdx[dofDataIdx[5]] );

   const uint_t          elMatSize = 36;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i];
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                             const Cell&                                 cell,
                                                             const uint_t                                level,
                                                             const indexing::Index&                      microCell,
                                                             const celldof::CellType                     cType,
                                                             const idx_t* const                          srcVertexIdx,
                                                             const idx_t* const                          srcEdgeIdx,
                                                             const idx_t* const                          dstVertexIdx,
                                                             const idx_t* const                          dstEdgeIdx ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   Matrix10r elMat;
   P2Form    form( form_ );
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   std::vector< uint_t > rowIdx( 10 );
   std::vector< uint_t > colIdx( 10 );

   for ( uint_t k = 0; k < 4; ++k )
   {
      rowIdx[k] = uint_c( dstVertexIdx[vertexDoFIndices[k]] );
      colIdx[k] = uint_c( srcVertexIdx[vertexDoFIndices[k]] );
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      rowIdx[k] = uint_c( dstEdgeIdx[edgeDoFIndices[k - 4]] );
      colIdx[k] = uint_c( srcEdgeIdx[edgeDoFIndices[k - 4]] );
   }

   const uint_t          elMatSize = 100;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i];
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

// P2ElementwiseLaplaceOperator
template class P2ElementwiseOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

// P2ElementwisePolarLaplaceOperator
template class P2ElementwiseOperator<
    P2FenicsForm< p2_polar_laplacian_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

// P2ElementwiseMassOperator
template class P2ElementwiseOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;

// P2ElementwiseDivKGradOperator
template class P2ElementwiseOperator< P2Form_divKgrad >;

// P2ElementwiseBlendingMassOperator
template class P2ElementwiseOperator< forms::p2_mass_blending_q5 >;

// P2ElementwiseBlendingLaplaceOperator
template class P2ElementwiseOperator< P2Form_laplace >;

// P2ElementwiseLinearCombinationOperator
template class P2ElementwiseOperator< P2LinearCombinationForm >;

// Needed for P2Blending(Inverse)DiagonalOperator
template class P2ElementwiseOperator< P2RowSumForm >;

template class P2ElementwiseOperator< forms::p2_div_k_grad_affine_q4 >;

template class P2ElementwiseOperator< forms::p2_epsiloncc_0_0_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsiloncc_0_1_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsiloncc_0_2_affine_q2 >;

template class P2ElementwiseOperator< forms::p2_epsiloncc_1_0_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsiloncc_1_1_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsiloncc_1_2_affine_q2 >;

template class P2ElementwiseOperator< forms::p2_epsiloncc_2_0_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsiloncc_2_1_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsiloncc_2_2_affine_q2 >;

template class P2ElementwiseOperator< forms::p2_epsilonvar_0_0_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_0_1_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_0_2_affine_q2 >;

template class P2ElementwiseOperator< forms::p2_epsilonvar_1_0_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_1_1_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_1_2_affine_q2 >;

template class P2ElementwiseOperator< forms::p2_epsilonvar_2_0_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_2_1_affine_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_2_2_affine_q2 >;

// Instantiations required for P2EpsilonOperator.hpp
template class P2ElementwiseOperator< forms::p2_epsilonvar_0_0_blending_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_0_1_blending_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_0_2_blending_q2 >;

template class P2ElementwiseOperator< forms::p2_epsilonvar_1_0_blending_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_1_1_blending_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_1_2_blending_q2 >;

template class P2ElementwiseOperator< forms::p2_epsilonvar_2_0_blending_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_2_1_blending_q2 >;
template class P2ElementwiseOperator< forms::p2_epsilonvar_2_2_blending_q2 >;

// Instantiations required for P2FullViscousOperator.hpp
template class P2ElementwiseOperator< forms::p2_full_stokesvar_0_0_blending_q3 >;
template class P2ElementwiseOperator< forms::p2_full_stokesvar_0_1_blending_q3 >;
template class P2ElementwiseOperator< forms::p2_full_stokesvar_0_2_blending_q3 >;

template class P2ElementwiseOperator< forms::p2_full_stokesvar_1_0_blending_q3 >;
template class P2ElementwiseOperator< forms::p2_full_stokesvar_1_1_blending_q3 >;
template class P2ElementwiseOperator< forms::p2_full_stokesvar_1_2_blending_q3 >;

template class P2ElementwiseOperator< forms::p2_full_stokesvar_2_0_blending_q3 >;
template class P2ElementwiseOperator< forms::p2_full_stokesvar_2_1_blending_q3 >;
template class P2ElementwiseOperator< forms::p2_full_stokesvar_2_2_blending_q3 >;

} // namespace hyteg
