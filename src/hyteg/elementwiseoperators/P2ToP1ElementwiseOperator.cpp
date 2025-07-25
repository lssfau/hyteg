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

#include "P2ToP1ElementwiseOperator.hpp"

namespace hyteg {

template < class P2toP1Form >
P2ToP1ElementwiseOperator< P2toP1Form >::P2ToP1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    size_t                                     minLevel,
                                                                    size_t                                     maxLevel )
: P2ToP1ElementwiseOperator( storage, minLevel, maxLevel, P2toP1Form() )
{}

template < class P2toP1Form >
P2ToP1ElementwiseOperator< P2toP1Form >::P2ToP1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    size_t                                     minLevel,
                                                                    size_t                                     maxLevel,
                                                                    const P2toP1Form&                          form )
: Operator( storage, minLevel, maxLevel )
, form_( form )
, localElementMatricesPrecomputed_( false )
{}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::computeAndStoreLocalElementMatrices()
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
                  Matrixr< 4, 10 >& elMat = localElementMatrix3D( *cell, level, micro, cType );
                  elMat.setZero();
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
                  Matrixr< 3, 6 >& elMat = localElementMatrix2D( *face, level, micro, fType );
                  elMat.setZero();
                  assembleLocalElementMatrix2D( *face, level, micro, fType, form_, elMat );
               }
            }
         }
      }
   }

   localElementMatricesPrecomputed_ = true;
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::gemv( const real_t&               alpha,
                                                    const P2Function< real_t >& src,
                                                    const real_t&               beta,
                                                    const P1Function< real_t >& dst,
                                                    uint_t                      level,
                                                    DoFType                     flag ) const
{
   this->startTiming( "apply" );

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
         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcVertexDoFIdx = src.getVertexDoFFunction().getCellDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Cell > srcEdgeDoFIdx = src.getEdgeDoFFunction().getCellDataID();

         real_t* srcVertexData = cell.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = cell.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeData = cell.getData( srcEdgeDoFIdx )->getPointer( level );

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

         Matrixr< 4, 10 > elMat = Matrixr< 4, 10 >::Zero();

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

               localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, srcEdgeData, dstVertexData, elMat, alpha );
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
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Face > srcEdgeDoFIdx = src.getEdgeDoFFunction().getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeData = face.getData( srcEdgeDoFIdx )->getPointer( level );

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

         Matrixr< 3, 6 > elMat = Matrixr< 3, 6 >::Zero();
         ;

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

               localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, srcEdgeData, dstVertexData, elMat, alpha );
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

   this->stopTiming( "apply" );
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::applyScaled( const real_t&               alpha,
                                                           const P2Function< real_t >& src,
                                                           const P1Function< real_t >& dst,
                                                           uint_t                      level,
                                                           DoFType                     flag,
                                                           UpdateType                  updateType ) const
{
   return gemv( real_c( alpha ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::apply( const P2Function< real_t >& src,
                                                     const P1Function< real_t >& dst,
                                                     size_t                      level,
                                                     DoFType                     flag,
                                                     UpdateType                  updateType ) const
{
   return gemv( real_c( 1 ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::localMatrixVectorMultiply2D( uint_t                 level,
                                                                           const indexing::Index& microFace,
                                                                           facedof::FaceType      fType,
                                                                           const real_t* const    srcVertexData,
                                                                           const real_t* const    srcEdgeData,
                                                                           real_t* const          dstVertexData,
                                                                           const Matrixr< 3, 6 >& elMat,
                                                                           const real_t&          alpha ) const
{
   // obtain data indices of dofs associated with micro-face
   std::array< uint_t, 3 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, fType, level, vertexDoFIndices );

   std::array< uint_t, 3 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, fType, level, edgeDoFIndices );

   // assemble local element vector
   Point6D elVecOld;
   Point3D elVecNew;
   for ( int k = 0; k < 3; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
   }
   for ( int k = 3; k < 6; ++k )
   {
      elVecOld[k] = srcEdgeData[edgeDoFIndices[uint_c( k - 3 )]];
   }

   // apply matrix (operator locally)
   elVecNew = alpha * ( elMat * elVecOld );

   // redistribute result from "local" to "global vector"
   for ( int k = 0; k < 3; ++k )
   {
      dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
   }
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::localMatrixVectorMultiply3D( const uint_t            level,
                                                                           const indexing::Index&  microCell,
                                                                           const celldof::CellType cType,
                                                                           const real_t* const     srcVertexData,
                                                                           const real_t* const     srcEdgeData,
                                                                           real_t* const           dstVertexData,
                                                                           const Matrixr< 4, 10 >& elMat,
                                                                           const real_t&           alpha ) const
{
   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // assemble local element vector
   Point10D elVecOld;
   Point4D  elVecNew;
   for ( int k = 0; k < 4; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
   }
   for ( int k = 4; k < 10; ++k )
   {
      elVecOld[k] = srcEdgeData[edgeDoFIndices[uint_c( k - 4 )]];
   }

   // apply matrix (operator locally)
   elVecNew = alpha * ( elMat * elVecOld );

   // redistribute result from "local" to "global vector"
   for ( int k = 0; k < 4; ++k )
   {
      dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
   }
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::assembleLocalElementMatrix2D( const Face&            face,
                                                                            uint_t                 level,
                                                                            const indexing::Index& microFace,
                                                                            facedof::FaceType      fType,
                                                                            P2toP1Form             form,
                                                                            Matrixr< 3, 6 >&       elMat ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 3 > verts = facedof::macroface::getMicroVerticesFromMicroFace( microFace, fType );
   std::array< Point3D, 3 >         coords;
   for ( uint_t k = 0; k < 3; ++k )
   {
      coords[k] = vertexdof::macroface::coordinateFromIndex( level, face, verts[k] );
   }

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( coords, elMat );
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::assembleLocalElementMatrix3D( const Cell&            cell,
                                                                            uint_t                 level,
                                                                            const indexing::Index& microCell,
                                                                            celldof::CellType      cType,
                                                                            P2toP1Form             form,
                                                                            Matrixr< 4, 10 >&      elMat ) const
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );
}

// Assemble operator as scaled sparse matrix
template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::toMatrixScaled( const real_t&                               alpha,
                                                              const std::shared_ptr< SparseMatrixProxy >& mat,
                                                              const P2Function< idx_t >&                  src,
                                                              const P1Function< idx_t >&                  dst,
                                                              uint_t                                      level,
                                                              DoFType                                     flag ) const
{
   // We currently ignore the flag provided!
   WALBERLA_UNUSED( flag );

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data in the two indexing functions
         PrimitiveDataID< FunctionMemory< idx_t >, Cell > dstVertexDoFIdx = dst.getCellDataID();
         PrimitiveDataID< FunctionMemory< idx_t >, Cell > srcVertexDoFIdx = src.getVertexDoFFunction().getCellDataID();

         PrimitiveDataID< FunctionMemory< idx_t >, Cell > srcEdgeDoFIdx = src.getEdgeDoFFunction().getCellDataID();

         idx_t* srcVertexIndices = cell.getData( srcVertexDoFIdx )->getPointer( level );
         idx_t* dstVertexIndices = cell.getData( dstVertexDoFIdx )->getPointer( level );

         idx_t* srcEdgeIndices = cell.getData( srcEdgeDoFIdx )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixAssembly3D( mat, cell, level, micro, cType, srcVertexIndices, srcEdgeIndices, dstVertexIndices, alpha );
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
         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< idx_t >, Face > srcEdgeDoFIdx = src.getEdgeDoFFunction().getFaceDataID();

         idx_t* srcVertexIndices = face.getData( srcVertexDoFIdx )->getPointer( level );
         idx_t* dstVertexIndices = face.getData( dstVertexDoFIdx )->getPointer( level );

         idx_t* srcEdgeIndices = face.getData( srcEdgeDoFIdx )->getPointer( level );

         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < idx_t( rowsize ) - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < idx_t( inner_rowsize ) - 1; ++xIdx )
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
                                      alpha );
               localMatrixAssembly2D( mat,
                                      face,
                                      level,
                                      xIdx,
                                      yIdx,
                                      P2Elements::P2Face::elementNW,
                                      srcVertexIndices,
                                      srcEdgeIndices,
                                      dstVertexIndices,
                                      alpha );
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
                                   alpha );
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
                                alpha );
      }
   }
}

// Assemble operator as sparse matrix
template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                        const P2Function< idx_t >&                  src,
                                                        const P1Function< idx_t >&                  dst,
                                                        uint_t                                      level,
                                                        DoFType                                     flag ) const
{
   return toMatrixScaled( real_c( 1 ), mat, src, dst, level, flag );
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                     const Face&                                 face,
                                                                     const uint_t                                level,
                                                                     const idx_t                                 xIdx,
                                                                     const idx_t                                 yIdx,
                                                                     const P2Elements::P2Element&                element,
                                                                     const idx_t* const                          srcVertexIdx,
                                                                     const idx_t* const                          srcEdgeIdx,
                                                                     const idx_t* const                          dstVertexIdx,
                                                                     const real_t&                               alpha ) const
{
   Matrixr< 3, 6 >         elMat = Matrixr< 3, 6 >::Zero();
   indexing::Index         nodeIdx;
   indexing::Index         offset;
   Point3D                 v0, v1, v2;
   std::array< uint_t, 6 > dofDataIdx;

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form_.setGeometryMap( face.getGeometryMap() );
   form_.integrateAll( { v0, v1, v2 }, elMat );

   // determine global indices of our local DoFs (note the tweaked ordering to go along with FEniCS indexing)
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   dofDataIdx[3] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[4] );
   dofDataIdx[4] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[5] );
   dofDataIdx[5] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[3] );

   std::vector< uint_t > rowIdx( 3 );
   rowIdx[0] = uint_c( dstVertexIdx[dofDataIdx[0]] );
   rowIdx[1] = uint_c( dstVertexIdx[dofDataIdx[1]] );
   rowIdx[2] = uint_c( dstVertexIdx[dofDataIdx[2]] );

   std::vector< uint_t > colIdx( 6 );
   colIdx[0] = uint_c( srcVertexIdx[dofDataIdx[0]] );
   colIdx[1] = uint_c( srcVertexIdx[dofDataIdx[1]] );
   colIdx[2] = uint_c( srcVertexIdx[dofDataIdx[2]] );

   colIdx[3] = uint_c( srcEdgeIdx[dofDataIdx[3]] );
   colIdx[4] = uint_c( srcEdgeIdx[dofDataIdx[4]] );
   colIdx[5] = uint_c( srcEdgeIdx[dofDataIdx[5]] );

   const uint_t          elMatSize = 3 * 6;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i] * alpha;
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

template < class P2toP1Form >
void P2ToP1ElementwiseOperator< P2toP1Form >::localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                     const Cell&                                 cell,
                                                                     const uint_t                                level,
                                                                     const indexing::Index&                      microCell,
                                                                     const celldof::CellType                     cType,
                                                                     const idx_t* const                          srcVertexIdx,
                                                                     const idx_t* const                          srcEdgeIdx,
                                                                     const idx_t* const                          dstVertexIdx,
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
   Matrixr< 4, 10 > elMat = Matrixr< 4, 10 >::Zero();

   form_.setGeometryMap( cell.getGeometryMap() );
   form_.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   std::vector< uint_t > rowIdx( 4 );
   std::vector< uint_t > colIdx( 10 );

   for ( uint_t k = 0; k < 4; ++k )
   {
      rowIdx[k] = uint_c( dstVertexIdx[vertexDoFIndices[k]] );
      colIdx[k] = uint_c( srcVertexIdx[vertexDoFIndices[k]] );
   }
   for ( uint_t k = 4; k < 10; ++k )
   {
      colIdx[k] = uint_c( srcEdgeIdx[edgeDoFIndices[k - 4]] );
   }

   const uint_t          elMatSize = 4 * 10;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i] * alpha;
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

template class P2ToP1ElementwiseOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;

template class P2ToP1ElementwiseOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;

template class P2ToP1ElementwiseOperator<
    P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

template class P2ToP1ElementwiseOperator< forms::p2_to_p1_div_0_blending_q2 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_div_1_blending_q2 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_div_2_blending_q2 >;

template class P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_0_affine_q2 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_1_affine_q2 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_2_affine_q2 >;

template class P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_0_blending_q3 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_1_blending_q3 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_2_blending_q3 >;

template class P2ToP1ElementwiseOperator< forms::p2_to_p1_k_mass_blending_q5 >;

template class P2ToP1ElementwiseOperator< forms::p2_to_p1_div_0_blending_q6 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_div_1_blending_q6 >;
template class P2ToP1ElementwiseOperator< forms::p2_to_p1_div_2_blending_q6 >;

} // namespace hyteg
