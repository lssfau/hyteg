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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/eigen/typeAliases.hpp"

namespace hyteg {
namespace n1e1 {

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
, localElementMatricesPrecomputed_( false )
{
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
   std::array< uint_t, 6 > edgeDoFIndices;
   n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // assemble local element vector
   Point6D elVecOld, elVecNew;
   for ( uint_t k = 0; k < 6; ++k )
   {
      elVecOld[k] = srcEdgeData[edgeDoFIndices[k]];
   }

   // apply matrix (operator locally)
   elVecNew = elMat.mul( elVecOld );

   // redistribute result from "local" to "global vector"
   for ( uint_t k = 0; k < 6; ++k )
   {
      dstEdgeData[edgeDoFIndices[k]] += elVecNew[k];
   }
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::apply( const N1E1VectorFunction< real_t >& src,
                                                     const N1E1VectorFunction< real_t >& dst,
                                                     size_t                              level,
                                                     DoFType                             flag,
                                                     UpdateType                          updateType ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   this->startTiming( "apply" );

   this->storage_->getTimingTree()->start( "sync source communication" );
   // Make sure that halos are up-to-date
   // Note that the order of communication is important, since the face -> cell communication may overwrite
   // parts of the halos that carry the macro-edge unknowns.

   src.communicate< Face, Cell >( level );
   src.communicate< Edge, Cell >( level );
   this->storage_->getTimingTree()->stop( "sync source communication" );

   if ( updateType == Replace )
   {
      // We need to zero the destination array (including halos).
      // However, we must not zero out anything that is not flagged with the specified BCs.
      // Therefore we first zero out everything that is flagged, and then, later,
      // the halos of the highest dim primitives.

      dst.interpolate( Eigen::Vector3r{ 0, 0, 0 }, level, flag );
   }

   // we only perform computations on cell primitives
   for ( auto& macroIter : storage_->getCells() )
   {
      Cell& cell = *macroIter.second;

      // get hold of the actual numerical data in the two functions
      PrimitiveDataID< FunctionMemory< real_t >, Cell > srcEdgeDoFIdx = src.getDoFs()->getCellDataID();
      PrimitiveDataID< FunctionMemory< real_t >, Cell > dstEdgeDoFIdx = dst.getDoFs()->getCellDataID();

      real_t* srcEdgeData = cell.getData( srcEdgeDoFIdx )->getPointer( level );
      real_t* dstEdgeData = cell.getData( dstEdgeDoFIdx )->getPointer( level );

      // Zero out dst halos only
      //
      // This is also necessary when using update type == Add.
      // During additive comm we then skip zeroing the data on the lower-dim primitives.

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

      Matrix6r elMat;

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

            localMatrixVectorMultiply3D( level, micro, cType, srcEdgeData, dstEdgeData, elMat );
         }
      }
   }

   this->storage_->getTimingTree()->start( "additive communication" );
   // Push result to lower-dimensional primitives
   //
   // Note: We could avoid communication here by implementing the apply() also for the respective
   //       lower dimensional primitives!
   dst.communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
   dst.communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
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
               Matrix6r& elMat = localElementMatrix3D( *cell, level, micro, cType );
               elMat.setAll( 0 );
               assembleLocalElementMatrix3D( *cell, level, micro, cType, form_, elMat );
            }
         }
      }
   }

   localElementMatricesPrecomputed_ = true;
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::computeInverseDiagonalOperatorValues()
{
   if ( inverseDiagonalValues_ == nullptr )
   {
      inverseDiagonalValues_ =
          std::make_shared< N1E1VectorFunction< real_t > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
   }
   else
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         inverseDiagonalValues_->setToZero( level );
      }
   }

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      // we only perform computations on cell primitives
      for ( auto& macroIter : storage_->getCells() )
      {
         Cell& cell = *macroIter.second;

         // get hold of the actual numerical data
         real_t* diagData = cell.getData( inverseDiagonalValues_->getDoFs()->getCellDataID() )->getPointer( level );

         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               computeLocalInverseDiagonal( cell, level, micro, cType, diagData );
            }
         }
      }

      // Push result to lower-dimensional primitives
      inverseDiagonalValues_->getDoFs()->communicateAdditively< Cell, Face >( level );
      inverseDiagonalValues_->getDoFs()->communicateAdditively< Cell, Edge >( level );

      inverseDiagonalValues_->getDoFs()->invertElementwise( level );
   }
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

   // TODO this is possible without allocating a new function
   N1E1VectorFunction< real_t > signs{ "signs", storage_, level, level };
   signs.getDoFs()->interpolate( 1, level );
   communication::syncFunctionBetweenPrimitives( signs, level );

   // we only perform computations on cell primitives
   for ( auto& macroIter : storage_->getCells() )
   {
      Cell& cell = *macroIter.second;

      // get hold of the actual numerical data in the two indexing functions
      PrimitiveDataID< FunctionMemory< idx_t >, Cell >  srcEdgeDoFId = src.getDoFs()->getCellDataID();
      PrimitiveDataID< FunctionMemory< idx_t >, Cell >  dstEdgeDoFId = dst.getDoFs()->getCellDataID();
      PrimitiveDataID< FunctionMemory< real_t >, Cell > signsId      = signs.getDoFs()->getCellDataID();

      idx_t*  srcEdgeIndices = cell.getData( srcEdgeDoFId )->getPointer( level );
      idx_t*  dstEdgeIndices = cell.getData( dstEdgeDoFId )->getPointer( level );
      real_t* signsData      = cell.getData( signsId )->getPointer( level );

      // loop over micro-cells
      for ( const auto& cType : celldof::allCellTypes )
      {
         for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
         {
            localMatrixAssembly3D( mat, cell, level, micro, cType, srcEdgeIndices, dstEdgeIndices, signsData );
         }
      }
   }
}

template < class N1E1FormType >
void N1E1ElementwiseOperator< N1E1FormType >::computeLocalInverseDiagonal( const Cell&             cell,
                                                                           const uint_t            level,
                                                                           const indexing::Index&  microCell,
                                                                           const celldof::CellType cType,
                                                                           real_t* const           diagData )
{
   Matrix6r elMat;
   if ( localElementMatricesPrecomputed_ )
   {
      elMat = localElementMatrix3D( cell, level, microCell, cType );
   }
   else
   {
      assembleLocalElementMatrix3D( cell, level, microCell, cType, form_, elMat );
   }

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 6 > edgeDoFIndices;
   n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   // add contributions for central stencil weights
   for ( uint_t k = 0; k < 6; ++k )
   {
      diagData[edgeDoFIndices[k]] += elMat( k, k );
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
                                                                     const real_t* const                         signs ) const
{
   Matrix6r elMat;
   if ( localElementMatricesPrecomputed_ )
   {
      elMat = localElementMatrix3D( cell, level, microCell, cType );
   }
   else
   {
      assembleLocalElementMatrix3D( cell, level, microCell, cType, form_, elMat );
   }

   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 6 > edgeDoFIndices;
   n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

   std::vector< uint_t > rowIdx( 6 );
   std::vector< uint_t > colIdx( 6 );

   for ( uint_t k = 0; k < 6; ++k )
   {
      rowIdx[k] = uint_c( dstEdgeIdx[edgeDoFIndices[k]] );
      colIdx[k] = uint_c( srcEdgeIdx[edgeDoFIndices[k]] );
   }

   for ( size_t i = 0; i < 6; ++i )
   {
      for ( size_t j = i + 1; j < 6; ++j )
      {
         const real_t sign = signs[edgeDoFIndices[i]] * signs[edgeDoFIndices[j]];
         elMat( i, j ) *= sign;
         elMat( j, i ) *= sign;
      }
   }

   std::vector< real_t > blockMatData( elMat.Size );
   for ( uint_t i = 0; i < elMat.Size; i++ )
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

} // namespace n1e1
} // namespace hyteg
