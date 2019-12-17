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
                                                        size_t                                     maxLevel,
                                                        bool                                       needsDiagEntries )
: Operator( storage, minLevel, maxLevel )
{
   if ( needsDiagEntries )
   {
      diagonalValues_ =
          std::unique_ptr< P2Function< real_t > >( new P2Function< real_t >( "diagonal entries", storage, minLevel, maxLevel ) );
      computeDiagonalOperatorValues( maxLevel, true );
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

   // Make sure that halos are up-to-date (can we improve communication here?)
   communication::syncP2FunctionBetweenPrimitives( src, level );

   if ( updateType == Add )
   {
      WALBERLA_ABORT( "P2ElementwiseOperator::apply does not support additive update!" );
      // Note: By zeroing only the dofs in the cells' or faces' halos of dst we could
      //       actually support an additive update type. That's incorrect, as the additive
      //       communication zeros before doing the accumulation.
      // 
      //       "Add" could be done, if we perform the computations on all primitives, which
      //       means that the additive communication is saved.
   }
   else
   {
      // We need to zero the destination array (including halos)
      dst.setToZero( level );
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

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               localMatrixVectorMultiply3D( cell, level, micro, cType, srcVertexData, srcEdgeData, dstVertexData, dstEdgeData );
            }
         }
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.getVertexDoFFunction().communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_ );
      dst.getVertexDoFFunction().communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_ );
      dst.getVertexDoFFunction().communicateAdditively< Cell, Vertex >( level, DoFType::All ^ flag, *storage_ );
      dst.getEdgeDoFFunction().communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_ );
      dst.getEdgeDoFFunction().communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_ );
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
         PrimitiveDataID< FunctionMemory< real_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getVertexDoFFunction().getFaceDataID();

         PrimitiveDataID< FunctionMemory< real_t >, Face > dstEdgeDoFIdx = dst.getEdgeDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > srcEdgeDoFIdx = src.getEdgeDoFFunction().getFaceDataID();

         real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );
         real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );

         real_t* srcEdgeData = face.getData( srcEdgeDoFIdx )->getPointer( level );
         real_t* dstEdgeData = face.getData( dstEdgeDoFIdx )->getPointer( level );

         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               localMatrixVectorMultiply2D( face,
                                            level,
                                            xIdx,
                                            yIdx,
                                            P2Elements::P2Face::elementN,
                                            srcVertexData,
                                            srcEdgeData,
                                            dstVertexData,
                                            dstEdgeData );
               localMatrixVectorMultiply2D( face,
                                            level,
                                            xIdx,
                                            yIdx,
                                            P2Elements::P2Face::elementNW,
                                            srcVertexData,
                                            srcEdgeData,
                                            dstVertexData,
                                            dstEdgeData );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            localMatrixVectorMultiply2D(
                face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, srcVertexData, srcEdgeData, dstVertexData, dstEdgeData );
         }

         // top north-west micro-element not treated, yet
         localMatrixVectorMultiply2D(
             face, level, 1, yIdx, P2Elements::P2Face::elementNW, srcVertexData, srcEdgeData, dstVertexData, dstEdgeData );
      }

      // Push result to lower-dimensional primitives
      //
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.getVertexDoFFunction().communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_ );
      dst.getVertexDoFFunction().communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage_ );
      dst.getEdgeDoFFunction().communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_ );
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
   dst.multElementwise( {( *diagonalValues_ ), dst}, level, flag );
   dst.assign( {1.0, omega}, {src, dst}, level, flag );

   this->stopTiming( "smooth_jac" );
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::localMatrixVectorMultiply2D( const Face&                  face,
                                                                   const uint_t                 level,
                                                                   const uint_t                 xIdx,
                                                                   const uint_t                 yIdx,
                                                                   const P2Elements::P2Element& element,
                                                                   const real_t* const          srcVertexData,
                                                                   const real_t* const          srcEdgeData,
                                                                   real_t* const                dstVertexData,
                                                                   real_t* const                dstEdgeData ) const
{
   WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );
   WALBERLA_ASSERT_UNEQUAL( srcEdgeData, dstEdgeData );

   Matrix6r                 elMat;
   Point6D                  elVecOld, elVecNew;
   indexing::Index          nodeIdx;
   indexing::IndexIncrement offset;
   Point3D                  v0, v1, v2;
   std::array< uint_t, 6 >  dofDataIdx;

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form_.integrateAll( {v0, v1, v2}, elMat );

   // assemble local element vector (note the tweaked ordering to go along with FEniCS indexing)
   dofDataIdx[0] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[0] );
   dofDataIdx[1] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[1] );
   dofDataIdx[2] = vertexdof::macroface::indexFromVertex( level, xIdx, yIdx, element[2] );

   dofDataIdx[3] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[4] );
   dofDataIdx[4] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[5] );
   dofDataIdx[5] = edgedof::macroface::indexFromVertex( level, xIdx, yIdx, element[3] );

   elVecOld[0] = srcVertexData[dofDataIdx[0]];
   elVecOld[1] = srcVertexData[dofDataIdx[1]];
   elVecOld[2] = srcVertexData[dofDataIdx[2]];

   elVecOld[3] = srcEdgeData[dofDataIdx[3]];
   elVecOld[4] = srcEdgeData[dofDataIdx[4]];
   elVecOld[5] = srcEdgeData[dofDataIdx[5]];

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   dstVertexData[dofDataIdx[0]] += elVecNew[0];
   dstVertexData[dofDataIdx[1]] += elVecNew[1];
   dstVertexData[dofDataIdx[2]] += elVecNew[2];

   dstEdgeData[dofDataIdx[3]] += elVecNew[3];
   dstEdgeData[dofDataIdx[4]] += elVecNew[4];
   dstEdgeData[dofDataIdx[5]] += elVecNew[5];
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::localMatrixVectorMultiply3D( const Cell&             cell,
                                                                   const uint_t            level,
                                                                   const indexing::Index&  microCell,
                                                                   const celldof::CellType cType,
                                                                   const real_t* const     srcVertexData,
                                                                   const real_t* const     srcEdgeData,
                                                                   real_t* const           dstVertexData,
                                                                   real_t* const           dstEdgeData ) const
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
   form_.integrateAll( coords, elMat );

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // assemble local element vector (we use FEniCS indexing? no we don't!)
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
void P2ElementwiseOperator< P2Form >::computeDiagonalOperatorValues( uint_t level, bool invert )
{
   WALBERLA_ASSERT_GREATER_EQUAL( level, minLevel_ );
   WALBERLA_ASSERT_LESS_EQUAL( level, maxLevel_ );
   WALBERLA_ASSERT_NOT_NULLPTR( diagonalValues_.get() );

   // Make sure that halos are up-to-date (can we improve communication here?)
   communication::syncP2FunctionBetweenPrimitives( *diagonalValues_, level );

   // Zero destination before performing additive computation
   diagonalValues_->setToZero( level );

   // For 3D we work on cells and for 2D on faces
   if ( storage_->hasGlobalCells() )
   {
      // WALBERLA_ABORT( "P2ElementwiseOperator::computeDiagonalOperatorValues() not implemented for cells, yet." );
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
         PrimitiveDataID< FunctionMemory< real_t >, Face > vertexDoFIdx = diagonalValues_->getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< real_t >, Face > edgeDoFIdx   = diagonalValues_->getEdgeDoFFunction().getFaceDataID();

         real_t* vertexData = face.getData( vertexDoFIdx )->getPointer( level );
         real_t* edgeData   = face.getData( edgeDoFIdx )->getPointer( level );

         // now loop over micro-faces of macro-face
         for ( yIdx = 0; yIdx < rowsize - 2; ++yIdx )
         {
            // loop over vertices in row with two associated triangles
            for ( xIdx = 1; xIdx < inner_rowsize - 1; ++xIdx )
            {
               // we associate two elements with current micro-vertex
               computeLocalDiagonalContributions2D( face, level, xIdx, yIdx, P2Elements::P2Face::elementN, vertexData, edgeData );
               computeLocalDiagonalContributions2D(
                   face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, vertexData, edgeData );
            }
            --inner_rowsize;

            // final micro-vertex in row has only one associated micro-face
            computeLocalDiagonalContributions2D( face, level, xIdx, yIdx, P2Elements::P2Face::elementNW, vertexData, edgeData );
         }

         // top north-west micro-element not treated, yet
         computeLocalDiagonalContributions2D( face, level, 1, yIdx, P2Elements::P2Face::elementNW, vertexData, edgeData );
      }

      // Push result to lower-dimensional primitives
      diagonalValues_->getVertexDoFFunction().communicateAdditively< Face, Edge >( level );
      diagonalValues_->getVertexDoFFunction().communicateAdditively< Face, Vertex >( level );
      diagonalValues_->getEdgeDoFFunction().communicateAdditively< Face, Edge >( level );

      // Retrieve assembled data values
      diagonalValues_->getVertexDoFFunction().communicate< Vertex, Edge >( level );
      diagonalValues_->getVertexDoFFunction().communicate< Edge, Face >( level );
      diagonalValues_->getEdgeDoFFunction().communicate< Edge, Face >( level );

      // Invert values if desired
      if ( invert )
      {
         diagonalValues_->invertElementwise( level );
      }
   }
}

template < class P2Form >
void P2ElementwiseOperator< P2Form >::computeLocalDiagonalContributions2D( Face&                        face,
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

   // determine vertices of micro-element
   nodeIdx = indexing::Index( xIdx, yIdx, 0 );
   v0      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[1] );
   v1      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );
   offset  = vertexdof::logicalIndexOffsetFromVertex( element[2] );
   v2      = vertexdof::macroface::coordinateFromIndex( level, face, nodeIdx + offset );

   // assemble local element matrix
   form_.integrateAll( {v0, v1, v2}, elMat );

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

// P2ElementwiseLaplaceOperator
template class P2ElementwiseOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

// P2ElementwiseMassOperator
template class P2ElementwiseOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;

} // namespace hyteg
