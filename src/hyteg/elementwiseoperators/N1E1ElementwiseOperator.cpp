/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "N1E1ElementwiseOperator.hpp"

#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"

namespace hyteg::n1e1 {

template < class N1E1FormType >
N1E1ElementwiseOperator< N1E1FormType >::N1E1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                  size_t                                     minLevel,
                                                                  size_t                                     maxLevel )
: N1E1ElementwiseOperator< N1E1FormType >( storage, minLevel, maxLevel, N1E1FormType() )
{}

template < class N1E1FormType >
N1E1ElementwiseOperator< N1E1FormType >::N1E1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                  size_t                                     minLevel,
                                                                  size_t                                     maxLevel,
                                                                  const N1E1FormType&                        form,
                                                                  const bool needsInverseDiagEntries )
: Operator( storage, minLevel, maxLevel )
, form_( form )
, edgeDirs_{ "edge directions", storage_, minLevel, maxLevel }
, localElementMatricesPrecomputed_( false )
{
   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      edgeDirs_.getDoFs()->interpolate( 1, level );
      edgeDirs_.communicate< Face, Cell >( level );
      edgeDirs_.communicate< Edge, Cell >( level );
   }

   if ( !storage_->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented for 2D." )
   }
   if ( needsInverseDiagEntries )
   {
      computeInverseDiagonalOperatorValues();
   }
}

void localMatrixVectorMultiply3D( uint_t                 level,
                                  const indexing::Index& microCell,
                                  celldof::CellType      cType,
                                  const real_t* const    srcEdgeData,
                                  real_t* const          dstEdgeData,
                                  const Matrix6r&        elMat )
{
   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 6 > edgeDoFIndices{};
   n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // assemble local element vector
   Point6D elVecOld, elVecNew;
   for ( int k = 0; k < 6; ++k )
   {
      elVecOld[k] = srcEdgeData[edgeDoFIndices[uint_c( k )]];
   }

   // apply matrix (operator locally)
   elVecNew = elMat * elVecOld;

   // redistribute result from "local" to "global vector"
   for ( int k = 0; k < 6; ++k )
   {
      dstEdgeData[edgeDoFIndices[uint_c( k )]] += elVecNew[k];
   }
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::apply( const N1E1VectorFunction< real_t >& src,
                                                     const N1E1VectorFunction< real_t >& dst,
                                                     const size_t                        level,
                                                     const DoFType                       flag,
                                                     const UpdateType                    updateType ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   this->startTiming( "apply" );

   this->storage_->getTimingTree()->start( "sync source communication" );
   // Make sure that halos are up-to-date
   //
   // Note: The order of communication is important, since the face -> cell communication may overwrite
   //       parts of the halos that carry the macro-edge unknowns.
   // Note: Edge directions can only be handled by the form.
   //       Therefore, the communication must be performed on the DoFs,
   //       otherwise signs would get flipped (again).
   src.getDoFs()->communicate< Face, Cell >( level );
   src.getDoFs()->communicate< Edge, Cell >( level );
   this->storage_->getTimingTree()->stop( "sync source communication" );

   if ( updateType == Replace )
   {
      // We need to zero the destination array (including halos).
      // However, we must not zero out anything that is not flagged with the specified BCs.
      // Therefore we first zero out everything that is flagged, and then, later,
      // the halos of the highest dim primitives.

      dst.interpolate( Point3D{ 0, 0, 0 }, level, flag );
   }

   // we only perform computations on cell primitives
   for ( auto& macroIter : storage_->getCells() )
   {
      Cell& cell = *macroIter.second;

      // get hold of the actual numerical data in the two functions
      PrimitiveDataID< FunctionMemory< real_t >, Cell > srcEdgeDoFIdx = src.getDoFs()->getCellDataID();
      PrimitiveDataID< FunctionMemory< real_t >, Cell > dstEdgeDoFIdx = dst.getDoFs()->getCellDataID();
      PrimitiveDataID< FunctionMemory< real_t >, Cell > edgeDirsId    = edgeDirs_.getDoFs()->getCellDataID();

      real_t* srcEdgeData  = cell.getData( srcEdgeDoFIdx )->getPointer( level );
      real_t* dstEdgeData  = cell.getData( dstEdgeDoFIdx )->getPointer( level );
      real_t* edgeDirsData = cell.getData( edgeDirsId )->getPointer( level );

      // Zero out dst halos only
      //
      // This is also necessary when using update type == Add.
      // During additive comm we then skip zeroing the data on the lower-dim primitives.
      edgedof::macrocell::setBoundaryToZero( level, cell, dstEdgeDoFIdx );

      Matrix6r elMat;

      // loop over micro-cells
      for ( const auto& cType : celldof::allCellTypes )
      {
         for ( const auto& micro : celldof::macrocell::Iterator( level, cType ) )
         {
            if ( localElementMatricesPrecomputed_ )
            {
               elMat = localElementMatrix3D( cell, level, micro, cType );
            }
            else
            {
               assembleLocalElementMatrix3D( cell, level, micro, cType, form_, edgeDirsData, elMat );
            }

            localMatrixVectorMultiply3D( level, micro, cType, srcEdgeData, dstEdgeData, elMat );
         }
      }
   }

   this->storage_->getTimingTree()->start( "additive communication" );
   // Push result to lower-dimensional primitives
   //
   // Note: Edge directions can only be handled by the form.
   //       Therefore, the communication must be performed on the DoFs,
   //       otherwise signs would get flipped (again).
   // Note: We could avoid communication here by implementing the apply() also for the respective
   //       lower dimensional primitives!
   dst.getDoFs()->communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
   dst.getDoFs()->communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
   this->storage_->getTimingTree()->stop( "additive communication" );

   this->stopTiming( "apply" );
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::computeAndStoreLocalElementMatrices()
{
   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      const uint_t numMicroCellsPerMacroCell = celldof::macrocell::numMicroCellsPerMacroCellTotal( level );

      for ( const auto& it : storage_->getCells() )
      {
         const PrimitiveID cellID = it.first;
         Cell&             cell   = *it.second;

         real_t* edgeDirsData = cell.getData( edgeDirs_.getDoFs()->getCellDataID() )->getPointer( level );

         auto& elementMatrices = localElementMatrices3D_[cellID][level];

         if ( !localElementMatricesPrecomputed_ )
         {
            elementMatrices.resize( numMicroCellsPerMacroCell );
         }

         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               Matrix6r& elMat = localElementMatrix3D( cell, level, micro, cType );
               elMat.setZero();
               assembleLocalElementMatrix3D( cell, level, micro, cType, form_, edgeDirsData, elMat );
            }
         }
      }
   }

   localElementMatricesPrecomputed_ = true;
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::computeDiagonalOperatorValues()
{
   computeDiagonalOperatorValues( false );
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::computeInverseDiagonalOperatorValues()
{
   computeDiagonalOperatorValues( true );
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::computeDiagonalOperatorValues( const bool invert )
{
   std::shared_ptr< N1E1VectorFunction< real_t > > targetFunction;
   if ( invert )
   {
      if ( inverseDiagonalValues_ == nullptr )
      {
         inverseDiagonalValues_ =
             std::make_shared< N1E1VectorFunction< real_t > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = inverseDiagonalValues_;
   }
   else
   {
      if ( diagonalValues_ == nullptr )
      {
         diagonalValues_ = std::make_shared< N1E1VectorFunction< real_t > >( "diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = diagonalValues_;
   }

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      targetFunction->setToZero( level );

      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data
         real_t* diagData     = cell.getData( targetFunction->getDoFs()->getCellDataID() )->getPointer( level );
         real_t* edgeDirsData = cell.getData( edgeDirs_.getDoFs()->getCellDataID() )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               computeLocalDiagonal( cell, level, micro, cType, edgeDirsData, diagData );
            }
         }
      }

      // Push result to lower-dimensional primitives.
      //
      // Note: Edge directions can only be handled by the form.
      //       Therefore, the communication must be performed on the DoFs,
      //       otherwise signs would get flipped (again).
      targetFunction->getDoFs()->communicateAdditively< Cell, Face >( level );
      targetFunction->getDoFs()->communicateAdditively< Cell, Edge >( level );

      if ( invert )
      {
         targetFunction->getDoFs()->invertElementwise( level );
      }
   }
}

template < class N1E1FormType >
std::shared_ptr< N1E1VectorFunction< real_t > > N1E1ElementwiseOperator< N1E1FormType >::getDiagonalValues() const
{
   WALBERLA_CHECK_NOT_NULLPTR(
       diagonalValues_, "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
   return diagonalValues_;
}

template < class N1E1FormType >
std::shared_ptr< N1E1VectorFunction< real_t > > N1E1ElementwiseOperator< N1E1FormType >::getInverseDiagonalValues() const
{
   WALBERLA_CHECK_NOT_NULLPTR(
       inverseDiagonalValues_,
       "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
   return inverseDiagonalValues_;
}

template < class N1E1FormType >
N1E1FormType N1E1ElementwiseOperator< N1E1FormType >::getForm() const
{
   return form_;
}

// Assemble operator as sparse matrix
template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                        const N1E1VectorFunction< idx_t >&          src,
                                                        const N1E1VectorFunction< idx_t >&          dst,
                                                        uint_t                                      level,
                                                        DoFType                                     flag ) const
{
   // We currently ignore the flag provided!
   // WALBERLA_UNUSED( flag );
   if ( flag != All )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in N1E1ElementwiseOperator::assembleLocalMatrix(); using flag = All" );
   }

   // we only perform computations on cell primitives
   for ( auto& macroIter : storage_->getCells() )
   {
      Cell& cell = *macroIter.second;

      // get hold of the actual numerical data in the two indexing functions
      PrimitiveDataID< FunctionMemory< idx_t >, Cell >  srcEdgeDoFId = src.getDoFs()->getCellDataID();
      PrimitiveDataID< FunctionMemory< idx_t >, Cell >  dstEdgeDoFId = dst.getDoFs()->getCellDataID();
      PrimitiveDataID< FunctionMemory< real_t >, Cell > edgeDirsId   = edgeDirs_.getDoFs()->getCellDataID();

      idx_t*  srcEdgeIndices = cell.getData( srcEdgeDoFId )->getPointer( level );
      idx_t*  dstEdgeIndices = cell.getData( dstEdgeDoFId )->getPointer( level );
      real_t* edgeDirsData   = cell.getData( edgeDirsId )->getPointer( level );

      // loop over micro-cells
      for ( const auto& cType : celldof::allCellTypes )
      {
         for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
         {
            localMatrixAssembly3D( mat, cell, level, micro, cType, srcEdgeIndices, dstEdgeIndices, edgeDirsData );
         }
      }
   }
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::computeLocalDiagonal( const Cell&             cell,
                                                                    const uint_t            level,
                                                                    const indexing::Index&  microCell,
                                                                    const celldof::CellType cType,
                                                                    const real_t* const     edgeDirsData,
                                                                    real_t* const           diagData )
{
   Matrix6r elMat;
   if ( localElementMatricesPrecomputed_ )
   {
      elMat = localElementMatrix3D( cell, level, microCell, cType );
   }
   else
   {
      assembleLocalElementMatrix3D( cell, level, microCell, cType, form_, edgeDirsData, elMat );
   }

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 6 > edgeDoFIndices{};
   n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // add contributions for central stencil weights
   for ( int k = 0; k < 6; ++k )
   {
      diagData[edgeDoFIndices[uint_c( k )]] += elMat( k, k );
   }
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                     const Cell&                                 cell,
                                                                     const uint_t                                level,
                                                                     const indexing::Index&                      microCell,
                                                                     const celldof::CellType                     cType,
                                                                     const idx_t* const                          srcEdgeIdx,
                                                                     const idx_t* const                          dstEdgeIdx,
                                                                     const real_t* const edgeDirsData ) const
{
   Matrix6r elMat;
   if ( localElementMatricesPrecomputed_ )
   {
      elMat = localElementMatrix3D( cell, level, microCell, cType );
   }
   else
   {
      assembleLocalElementMatrix3D( cell, level, microCell, cType, form_, edgeDirsData, elMat );
   }

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 6 > edgeDoFIndices{};
   n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   std::vector< uint_t > rowIdx( 6 );
   std::vector< uint_t > colIdx( 6 );

   for ( uint_t k = 0; k < 6; ++k )
   {
      rowIdx[k] = uint_c( dstEdgeIdx[edgeDoFIndices[k]] );
      colIdx[k] = uint_c( srcEdgeIdx[edgeDoFIndices[k]] );
   }

   const uint_t          elMatSize = 36;
   std::vector< real_t > blockMatData( elMatSize );
   for ( uint_t i = 0; i < elMatSize; i++ )
   {
      blockMatData[i] = elMat.data()[i];
   }

   // add local matrix into global matrix
   mat->addValues( rowIdx, colIdx, blockMatData );
}

// N1E1ElementwiseCurlCurlOperator
template class N1E1ElementwiseOperator< N1E1Form_curl_curl >;
// N1E1ElementwiseMassOperator
template class N1E1ElementwiseOperator< N1E1Form_mass >;
// N1E1ElementwiseLinearCombinationOperator
template class N1E1ElementwiseOperator< N1E1LinearCombinationForm >;
// N1E1ElementwiseLinearFormOperatorQ6
template class N1E1ElementwiseOperator< forms::n1e1_linear_form_affine_q6 >;

} // namespace hyteg::n1e1
