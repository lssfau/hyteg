/*
* Copyright (c) 2023-2025 Andreas Burkhart
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

// This file has been generated with the AHFC. If buggy try fixing the generator itself.

#pragma once

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/forms/form_hyteg_base/P1ToP2FormHyTeG.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

#include "ElementwiseDoFValueOperator.hpp"

namespace hyteg {

using walberla::real_t;

/// Elementwise DoF value operators are a variant of elementwise operators that allow you
/// to define a list of FEM functions via templating that the operator will communicate before
/// application and during usage the operator will fill an array called DoFValues_ with evaluations 
/// (in hyteg ordering) of the FEM function at the micro element nodes.
/// These blocks of evaluation values are put into the vector in the order of the given
/// FEM function list.
/// 
/// Technically an elementwise DoF value operator is both a form and an operator. During
/// form evaluation this allows the form in question access the DoFValues_ array.
/// 
/// Usually these operators are used by inheriting from the operator class and redefining the
/// "integrateAll" and "integrateRow0" functions of the form.
/// 
/// Elementwise DoF value operators are currently incompatible with an empty FEM function list.
/// In this case the templating breaks and you should use standard elementwise operators in that
/// case anyway.
template < typename... Functions >
class P1ToP2ElementwiseDoFValueOperator : public Operator< P1Function< real_t >, P2Function< real_t > >, public P1ToP2FormHyTeG
{
 public:
   P1ToP2ElementwiseDoFValueOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                      size_t                                     minLevel,
                                      size_t                                     maxLevel,
                                      Functions&... funcs )
   : Operator( storage, minLevel, maxLevel )
   , localElementMatricesPrecomputed_( false )
   , functuple_( funcs... )
   {
      static_assert( RestrictToFEMFunctions<>::value, "Functions... may contain only FEM Functions" );
      // resize the vector to the required size
      if ( storage_->hasGlobalCells() )
      {
         DoFValues_.resize( NumberOfLocalDoFsCell<>::value, 0 );
      }
      else
      {
         DoFValues_.resize( NumberOfLocalDoFsFace<>::value, 0 );
      }
   }

   void setDoFValues( const std::vector< real_t >& values ) { DoFValues_ = std::vector< real_t >( values ); }

   virtual void gemv( const real_t&               alpha,
                      const P1Function< real_t >& src,
                      const real_t&               beta,
                      const P2Function< real_t >& dst,
                      uint_t                      level,
                      DoFType                     flag ) const override
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

      communicateFunctions( functuple_, level );

      // For 3D we work on cells and for 2D on faces
      if ( storage_->hasGlobalCells() )
      {
         // we only perform computations on cell primitives
         for ( auto& macroIter : storage_->getCells() )
         {
            Cell& cell = *macroIter.second;

            // get hold of the actual data in the two functions
            PrimitiveDataID< FunctionMemory< real_t >, Cell > srcVertexDoFIdx = src.getCellDataID();
            real_t* srcVertexData = cell.getData( srcVertexDoFIdx )->getPointer( level );

            PrimitiveDataID< FunctionMemory< real_t >, Cell > dstVertexDoFIdx = dst.getVertexDoFFunction().getCellDataID();
            PrimitiveDataID< FunctionMemory< real_t >, Cell > dstEdgeDoFIdx   = dst.getEdgeDoFFunction().getCellDataID();
            real_t* dstVertexData = cell.getData( dstVertexDoFIdx )->getPointer( level );
            real_t* dstEdgeData   = cell.getData( dstEdgeDoFIdx )->getPointer( level );

            // Zero out dst halos only
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

            Matrixr< 10, 4 > elMat = Matrixr< 10, 4 >::Zero();

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
                     assembleLocalElementMatrix3D( cell, level, micro, cType, elMat );
                  }
                  localMatrixVectorMultiply3D( level, micro, cType, srcVertexData, dstVertexData, dstEdgeData, elMat, alpha );
               }
            }
         }

         // Push result to lower-dimensional primitives
         // Note: We could avoid communication here by implementing the apply() also for the respective
         //       lower dimensional primitives!
         dst.getVertexDoFFunction().communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_, betaIsZero );
         dst.getVertexDoFFunction().communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
         dst.getVertexDoFFunction().communicateAdditively< Cell, Vertex >( level, DoFType::All ^ flag, *storage_, betaIsZero );
         dst.getEdgeDoFFunction().communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_, betaIsZero );
         dst.getEdgeDoFFunction().communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      }
      else
      {
         // we only perform computations on face primitives
         for ( auto& it : storage_->getFaces() )
         {
            Face& face = *it.second;

            // get hold of the actual data in the two functions
            PrimitiveDataID< FunctionMemory< real_t >, Face > srcVertexDoFIdx = src.getFaceDataID();
            real_t* srcVertexData = face.getData( srcVertexDoFIdx )->getPointer( level );

            PrimitiveDataID< FunctionMemory< real_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
            PrimitiveDataID< FunctionMemory< real_t >, Face > dstEdgeDoFIdx   = dst.getEdgeDoFFunction().getFaceDataID();
            real_t* dstVertexData = face.getData( dstVertexDoFIdx )->getPointer( level );
            real_t* dstEdgeData   = face.getData( dstEdgeDoFIdx )->getPointer( level );

            // Zero out dst halos only
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

            Matrixr< 6, 3 > elMat = Matrixr< 6, 3 >::Zero();

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
                     assembleLocalElementMatrix2D( face, level, micro, fType, elMat );
                  }
                  localMatrixVectorMultiply2D( level, micro, fType, srcVertexData, dstVertexData, dstEdgeData, elMat, alpha );
               }
            }
         }

         // Push result to lower-dimensional primitives
         // Note: We could avoid communication here by implementing the apply() also for the respective
         //       lower dimensional primitives!
         dst.getVertexDoFFunction().communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
         dst.getVertexDoFFunction().communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage_, betaIsZero );
         dst.getEdgeDoFFunction().communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, betaIsZero );
      }

      this->stopTiming( "apply" );
   }

   virtual void applyScaled( const real_t&               alpha,
                             const P1Function< real_t >& src,
                             const P2Function< real_t >& dst,
                             uint_t                      level,
                             DoFType                     flag,
                             UpdateType                  updateType = Replace ) const override

   {
      return gemv( real_c( alpha ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
   }

   virtual void apply( const P1Function< real_t >& src,
                       const P2Function< real_t >& dst,
                       size_t                      level,
                       DoFType                     flag,
                       UpdateType                  updateType ) const override
   {
      return gemv( real_c( 1 ), src, ( updateType == Replace ? real_c( 0 ) : real_c( 1 ) ), dst, level, flag );
   }

   /// Assemble operator as sparse matrix with scaling
   ///
   /// \param alpha constant scaling of the matrix
   /// \param mat   a sparse matrix proxy
   /// \param src   P1Function for determining column indices
   /// \param dst   P2Function for determining row indices
   /// \param level le2el in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   virtual void toMatrixScaled( const real_t&                               alpha,
                                const std::shared_ptr< SparseMatrixProxy >& mat,
                                const P1Function< idx_t >&                  src,
                                const P2Function< idx_t >&                  dst,
                                uint_t                                      level,
                                DoFType                                     flag ) const override
   {
      // We currently ignore the flag provided!
      if ( flag != All )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in assembleLocalMatrix(); using flag = All" );
      }

      // For 3D we work on cells and for 2D on faces
      if ( storage_->hasGlobalCells() )
      {
         // we only perform computations on cell primitives
         for ( auto& macroIter : storage_->getCells() )
         {
            Cell& cell = *macroIter.second;

            // get hold of the actual indices in the two functions
            PrimitiveDataID< FunctionMemory< idx_t >, Cell > srcVertexDoFIdx = src.getCellDataID();
            idx_t* srcVertexIndices = cell.getData( srcVertexDoFIdx )->getPointer( level );

            PrimitiveDataID< FunctionMemory< idx_t >, Cell > dstVertexDoFIdx = dst.getVertexDoFFunction().getCellDataID();
            PrimitiveDataID< FunctionMemory< idx_t >, Cell > dstEdgeDoFIdx   = dst.getEdgeDoFFunction().getCellDataID();
            idx_t* dstVertexIndices = cell.getData( dstVertexDoFIdx )->getPointer( level );
            idx_t* dstEdgeIndices   = cell.getData( dstEdgeDoFIdx )->getPointer( level );

            // loop over micro-cells
            for ( const auto& cType : celldof::allCellTypes )
            {
               for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
               {
                  localMatrixAssembly3D(
                      mat, cell, level, micro, cType, srcVertexIndices, dstVertexIndices, dstEdgeIndices, alpha );
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

            // get hold of the actual indices in the two functions
            PrimitiveDataID< FunctionMemory< idx_t >, Face > srcVertexDoFIdx = src.getFaceDataID();
            idx_t* srcVertexIndices = face.getData( srcVertexDoFIdx )->getPointer( level );

            PrimitiveDataID< FunctionMemory< idx_t >, Face > dstVertexDoFIdx = dst.getVertexDoFFunction().getFaceDataID();
            PrimitiveDataID< FunctionMemory< idx_t >, Face > dstEdgeDoFIdx   = dst.getEdgeDoFFunction().getFaceDataID();
            idx_t* dstVertexIndices = face.getData( dstVertexDoFIdx )->getPointer( level );
            idx_t* dstEdgeIndices   = face.getData( dstEdgeDoFIdx )->getPointer( level );

            // loop over micro-faces
            for ( const auto& fType : facedof::allFaceTypes )
            {
               for ( const auto& micro : facedof::macroface::Iterator( level, fType, 0 ) )
               {
                  localMatrixAssembly2D(
                      mat, face, level, micro, fType, srcVertexIndices, dstVertexIndices, dstEdgeIndices, alpha );
               }
            }
         }
      }
   }

   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param src   Function for determining column indices
   /// \param dst   Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   virtual void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                          const P1Function< idx_t >&                  src,
                          const P2Function< idx_t >&                  dst,
                          uint_t                                      level,
                          DoFType                                     flag ) const override
   {
      return toMatrixScaled( real_c( 1 ), mat, src, dst, level, flag );
   }

   /// \brief Pre-computes the local stiffness matrices for each (micro-)element and stores them all in memory.
   ///
   /// If this method is called, all subsequent calls to apply() or smooth_*() use the stored element matrices.
   /// If the local element matrices need to be recomputed again, simply call this method again.
   virtual void computeAndStoreLocalElementMatrices()
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         // This is necessary, otherwise saving the local element matrices might change the result of applying the operator
         communicateFunctions( functuple_, level );

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
                     Matrixr< 10, 4 >& elMat = localElementMatrix3D( *cell, level, micro, cType );
                     elMat.setZero();
                     assembleLocalElementMatrix3D( *cell, level, micro, cType, elMat );
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
                     Matrixr< 6, 3 >& elMat = localElementMatrix2D( *face, level, micro, fType );
                     elMat.setZero();
                     assembleLocalElementMatrix2D( *face, level, micro, fType, elMat );
                  }
               }
            }
         }
      }

      localElementMatricesPrecomputed_ = true;
   }

 protected:
   // this allows us to determine the type of the I-th type in Functions...
   template < size_t I >
   using IRefType = typename std::tuple_element< I, std::tuple< Functions... > >::type;

   // this function iterates through all fem functions in functuple_ and updates their halos
   template < size_t I = 0, typename... Inp >
   void communicateFunctions( std::tuple< Inp... >& tup, uint_t level ) const
   {
      if constexpr ( sizeof...( Inp ) > 0 )
      {
         communicateDoFValue( std::get< I >( tup ).get(), level );
         if constexpr ( I + 1 != sizeof...( Inp ) )
         {
            communicateFunctions< I + 1 >( tup, level );
         }
      }
   }

   // this function iterates through all fem functions in functuple_ and retrieves the corresponding dofvalues on the microcell
   template < size_t I = 0, size_t offset = 0, typename... Inp >
   void getDoFsFromMicroCellFunctions( std::tuple< Inp... >&    tup,
                                       const Cell&              cell,
                                       const indexing::Index&   microCell,
                                       const celldof::CellType& cellType,
                                       const uint_t             level ) const
   {
      if constexpr ( sizeof...( Inp ) > 0 )
      {
         getDoFDataFromMicroCell( std::get< I >( tup ).get(), cell, microCell, cellType, level, DoFValues_, offset );
         if constexpr ( I + 1 != sizeof...( Inp ) )
         {
            getDoFsFromMicroCellFunctions< I + 1, offset + functionLocalDoFsCell< IRefType< I > >::value >(
                tup, cell, microCell, cellType, level );
         }
      }
   }

   // this function iterates through all fem functions in functuple_ and retrieves the corresponding dofvalues on the microface
   template < size_t I = 0, size_t offset = 0, typename... Inp >
   void getDoFsFromMicroFaceFunctions( std::tuple< Inp... >&    tup,
                                       const Face&              face,
                                       const indexing::Index&   microFace,
                                       const facedof::FaceType& fType,
                                       const uint_t             level ) const
   {
      if constexpr ( sizeof...( Inp ) > 0 )
      {
         getDoFDataFromMicroFace( std::get< I >( tup ).get(), face, microFace, fType, level, DoFValues_, offset );
         if constexpr ( I + 1 != sizeof...( Inp ) )
         {
            getDoFsFromMicroFaceFunctions< I + 1, offset + functionLocalDoFsFace< IRefType< I > >::value >(
                tup, face, microFace, fType, level );
         }
      }
   }

   // this determines the necessary length of the DoFValues_ vector, available at compile time
   template < int I = sizeof...( Functions ) - 1, bool Dummy = true >
   struct NumberOfLocalDoFsFace
   {
      static const size_t value = functionLocalDoFsFace< IRefType< I > >::value + NumberOfLocalDoFsFace< I - 1, Dummy >::value;
   };
   // default case
   template < bool Dummy >
   struct NumberOfLocalDoFsFace< -1, Dummy >
   {
      static const size_t value = 0;
   };

   template < int I = sizeof...( Functions ) - 1, bool Dummy = true >
   struct NumberOfLocalDoFsCell
   {
      static const size_t value = functionLocalDoFsCell< IRefType< I > >::value + NumberOfLocalDoFsCell< I - 1, Dummy >::value;
   };
   // default case
   template < bool Dummy >
   struct NumberOfLocalDoFsCell< -1, Dummy >
   {
      static const size_t value = 0;
   };

   // this determines if only allowed FEMFunctions types have been used as templates
   // is used for a static assert
   template < int I = sizeof...( Functions ) - 1, bool Dummy = true >
   struct RestrictToFEMFunctions
   {
      static const bool value = ( std::is_same< IRefType< I >, P1Function< real_t > >::value ||
                                  std::is_same< IRefType< I >, P2Function< real_t > >::value ) &&
                                RestrictToFEMFunctions< I - 1, Dummy >::value;
   };

   // default case
   template < bool Dummy >
   struct RestrictToFEMFunctions< -1, Dummy >
   {
      static const bool value = true;
   };

   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microFace      index associated with the current element = micro-face
   /// \param fType          type of micro-face (GRAY or BLUE)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param elMat          the element matrix to be multiplied
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply2D( const uint_t           level,
                                     const indexing::Index& microFace,
                                     facedof::FaceType      fType,
                                     const real_t* const    srcVertexData,
                                     real_t* const          dstVertexData,
                                     real_t* const          dstEdgeData,
                                     const Matrixr< 6, 3 >& elMat,
                                     const real_t&          alpha ) const
   {
      WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

      // obtain data indices of dofs associated with micro-face
      std::array< uint_t, 3 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, fType, level, vertexDoFIndices );

      std::array< uint_t, 3 > edgeDoFIndices;
      edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, fType, level, edgeDoFIndices );

      // assemble local element vector
      Point3D elVecOld;
      Point6D elVecNew;
      for ( int k = 0; k < 3; k++ )
      {
         elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
      }

      // apply matrix (operator locally)
      elVecNew = alpha * ( elMat * elVecOld );

      // redistribute result from "local" to "global vector"
      for ( int k = 0; k < 3; k++ )
      {
         dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
      }
      for ( int k = 3; k < 6; k++ )
      {
         dstEdgeData[edgeDoFIndices[uint_c( k - 3 )]] += elVecNew[k];
      }
   }

   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microCell      index associated with the current element = micro-cell
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param elMat          the element matrix to be multiplied
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply3D( const uint_t            level,
                                     const indexing::Index&  microCell,
                                     const celldof::CellType cType,
                                     const real_t* const     srcVertexData,
                                     real_t* const           dstVertexData,
                                     real_t* const           dstEdgeData,
                                     const Matrixr< 10, 4 >& elMat,
                                     const real_t&           alpha ) const
   {
      WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

      // obtain data indices of dofs associated with micro-cell
      std::array< uint_t, 4 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

      std::array< uint_t, 6 > edgeDoFIndices;
      edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

      // assemble local element vector
      Point4D  elVecOld;
      Point10D elVecNew;
      for ( int k = 0; k < 4; k++ )
      {
         elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
      }

      // apply matrix (operator locally)
      elVecNew = alpha * ( elMat * elVecOld );

      // redistribute result from "local" to "global vector"
      for ( int k = 0; k < 4; k++ )
      {
         dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
      }
      for ( int k = 4; k < 10; k++ )
      {
         dstEdgeData[edgeDoFIndices[uint_c( k - 4 )]] += elVecNew[k];
      }
   }

   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const indexing::Index&                      microFace,
                               const facedof::FaceType                     fType,
                               const idx_t* const                          srcVertexIdx,
                               const idx_t* const                          dstVertexIdx,
                               const idx_t* const                          dstEdgeIdx,
                               const real_t&                               alpha ) const
   {
      // determine coordinates of vertices of micro-element
      std::array< indexing::Index, 3 > verts = facedof::macroface::getMicroVerticesFromMicroFace( microFace, fType );
      std::array< Point3D, 3 >         coords;
      for ( uint_t k = 0; k < 3; k++ )
      {
         coords[k] = vertexdof::macroface::coordinateFromIndex( level, face, verts[k] );
      }

      // assemble local element matrix
      Matrixr< 6, 3 > elMat = Matrixr< 6, 3 >::Zero();

      getDoFsFromMicroFaceFunctions( functuple_, face, microFace, fType, level );

      setGeometryMap( face.getGeometryMap() );
      integrateAll( coords, elMat );

      // obtain data indices of dofs associated with micro-face
      std::array< uint_t, 3 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, fType, level, vertexDoFIndices );

      std::array< uint_t, 3 > edgeDoFIndices;
      edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, fType, level, edgeDoFIndices );

      std::vector< uint_t > rowIdx( 6 );
      std::vector< uint_t > colIdx( 3 );
      for ( uint_t k = 0; k < 3; k++ )
      {
         rowIdx[k] = uint_c( dstVertexIdx[vertexDoFIndices[k]] );
         colIdx[k] = uint_c( srcVertexIdx[vertexDoFIndices[k]] );
      }
      for ( uint_t k = 3; k < 6; k++ )
      {
         rowIdx[k] = uint_c( dstEdgeIdx[edgeDoFIndices[k - 3]] );
      }

      const uint_t          elMatSize = 18;
      std::vector< real_t > blockMatData( elMatSize );
      for ( uint_t i = 0; i < elMatSize; i++ )
      {
         blockMatData[i] = elMat.data()[i] * alpha;
      }

      // add local matrix into global matrix
      mat->addValues( rowIdx, colIdx, blockMatData );
   }

   void localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Cell&                                 cell,
                               const uint_t                                level,
                               const indexing::Index&                      microCell,
                               const celldof::CellType                     cType,
                               const idx_t* const                          srcVertexIdx,
                               const idx_t* const                          dstVertexIdx,
                               const idx_t* const                          dstEdgeIdx,
                               const real_t&                               alpha ) const
   {
      // determine coordinates of vertices of micro-element
      std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
      std::array< Point3D, 4 >         coords;
      for ( uint_t k = 0; k < 4; k++ )
      {
         coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
      }

      // assemble local element matrix
      Matrixr< 10, 4 > elMat = Matrixr< 10, 4 >::Zero();

      getDoFsFromMicroCellFunctions( functuple_, cell, microCell, cType, level );

      setGeometryMap( cell.getGeometryMap() );
      integrateAll( coords, elMat );

      // obtain data indices of dofs associated with micro-cell
      std::array< uint_t, 4 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cType, level, vertexDoFIndices );

      std::array< uint_t, 6 > edgeDoFIndices;
      edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cType, level, edgeDoFIndices );

      std::vector< uint_t > rowIdx( 10 );
      std::vector< uint_t > colIdx( 4 );
      for ( uint_t k = 0; k < 4; k++ )
      {
         rowIdx[k] = uint_c( dstVertexIdx[vertexDoFIndices[k]] );
         colIdx[k] = uint_c( srcVertexIdx[vertexDoFIndices[k]] );
      }
      for ( uint_t k = 4; k < 10; k++ )
      {
         rowIdx[k] = uint_c( dstEdgeIdx[edgeDoFIndices[k - 4]] );
      }

      const uint_t          elMatSize = 40;
      std::vector< real_t > blockMatData( elMatSize );
      for ( uint_t i = 0; i < elMatSize; i++ )
      {
         blockMatData[i] = elMat.data()[i] * alpha;
      }

      // add local matrix into global matrix
      mat->addValues( rowIdx, colIdx, blockMatData );
   }

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrixr< 6, 3 >&
       localElementMatrix2D( const Face& face, uint_t level, const indexing::Index& microFace, facedof::FaceType fType )
   {
      WALBERLA_ASSERT( !storage_->hasGlobalCells(), "Retriveing local element matrix for 2D in 3D run. Why?" )
      const auto idx = facedof::macroface::index( level, microFace.x(), microFace.y(), fType );
      WALBERLA_ASSERT( localElementMatrices2D_.count( face.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices2D_.at( face.getID() ).at( level ).empty() )
      return localElementMatrices2D_[face.getID()][level][idx];
   }

   /// \brief Returns a const reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   const Matrixr< 6, 3 >&
       localElementMatrix2D( const Face& face, uint_t level, const indexing::Index& microFace, facedof::FaceType fType ) const
   {
      WALBERLA_ASSERT( !storage_->hasGlobalCells(), "Retriveing local element matrix for 2D in 3D run. Why?" )
      const auto idx = facedof::macroface::index( level, microFace.x(), microFace.y(), fType );
      WALBERLA_ASSERT( localElementMatrices2D_.count( face.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices2D_.at( face.getID() ).at( level ).empty() )
      return localElementMatrices2D_.at( face.getID() ).at( level ).at( idx );
   }

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro cell.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrixr< 10, 4 >&
       localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType )
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells(), "Retriveing local element matrix for 3D in 2D run. Why?" )
      const auto idx = celldof::macrocell::index( level, microCell.x(), microCell.y(), microCell.z(), cType );
      WALBERLA_ASSERT( localElementMatrices3D_.count( cell.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices3D_.at( cell.getID() ).at( level ).empty() )
      return localElementMatrices3D_[cell.getID()][level][idx];
   }

   /// \brief Returns a const reference to the a precomputed element matrix of the specified micro cell.
   /// Probably crashes if local element matrices have not been precomputed.
   const Matrixr< 10, 4 >&
       localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells(), "Retriveing local element matrix for 3D in 2D run. Why?" )
      const auto idx = celldof::macrocell::index( level, microCell.x(), microCell.y(), microCell.z(), cType );
      WALBERLA_ASSERT( localElementMatrices3D_.count( cell.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices3D_.at( cell.getID() ).at( level ).empty() )
      return localElementMatrices3D_.at( cell.getID() ).at( level ).at( idx );
   }

   void assembleLocalElementMatrix2D( const Face&            face,
                                      uint_t                 level,
                                      const indexing::Index& microFace,
                                      facedof::FaceType      fType,
                                      Matrixr< 6, 3 >&       elMat ) const
   {
      // determine coordinates of vertices of micro-element
      std::array< indexing::Index, 3 > verts = facedof::macroface::getMicroVerticesFromMicroFace( microFace, fType );
      std::array< Point3D, 3 >         coords;
      for ( uint_t k = 0; k < 3; k++ )
      {
         coords[k] = vertexdof::macroface::coordinateFromIndex( level, face, verts[k] );
      }

      getDoFsFromMicroFaceFunctions( functuple_, face, microFace, fType, level );

      // assemble local element matrix
      setGeometryMap( face.getGeometryMap() );
      integrateAll( coords, elMat );
   }

   void assembleLocalElementMatrix3D( const Cell&            cell,
                                      uint_t                 level,
                                      const indexing::Index& microCell,
                                      celldof::CellType      cType,
                                      Matrixr< 10, 4 >&      elMat ) const
   {
      // determine coordinates of vertices of micro-element
      std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
      std::array< Point3D, 4 >         coords;
      for ( uint_t k = 0; k < 4; k++ )
      {
         coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
      }

      getDoFsFromMicroCellFunctions( functuple_, cell, microCell, cType, level );

      // assemble local element matrix
      setGeometryMap( cell.getGeometryMap() );
      integrateAll( coords, elMat );
   }

   mutable std::vector< real_t > DoFValues_;

   bool localElementMatricesPrecomputed_;

   mutable std::tuple< std::reference_wrapper< Functions >... > functuple_;

   /// Pre-computed local element matrices.
   /// localElementMatrices2D_[macroCellID][level][cellIdx] = Matrixr< 6, 3 >
   std::map< PrimitiveID, std::map< uint_t, std::vector< Matrixr< 6, 3 >, Eigen::aligned_allocator< Matrixr< 6, 3 > > > > >
       localElementMatrices2D_;

   /// Pre-computed local element matrices.
   /// localElementMatrices3D_[macroCellID][level][cellIdx] = Matrixr< 10, 4 >
   std::map< PrimitiveID, std::map< uint_t, std::vector< Matrixr< 10, 4 >, Eigen::aligned_allocator< Matrixr< 10, 4 > > > > >
       localElementMatrices3D_;
};
} // namespace hyteg