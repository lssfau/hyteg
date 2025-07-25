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
#pragma once

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p2_to_p1_div.h"
#include "hyteg/forms/form_fenics_generated/p2_to_p1_tet_div_tet.h"
#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_divT_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_divT_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_div_blending_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_div_blending_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_k_mass_blending_q5.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

using walberla::real_t;

template < class P2toP1Form >
class P2ToP1ElementwiseOperator : public Operator< P2Function< real_t >, P1Function< real_t > >
{
 public:
   P2ToP1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   P2ToP1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                              size_t                                     minLevel,
                              size_t                                     maxLevel,
                              const P2toP1Form&                          form );

   /// \brief Pre-computes the local stiffness matrices for each (micro-)element and stores them all in memory.
   ///
   /// If this method is called, all subsequent calls to apply() or smooth_*() use the stored element matrices.
   /// If the local element matrices need to be recomputed again, simply call this method again.
   void computeAndStoreLocalElementMatrices();

   void gemv( const real_t&               alpha,
              const P2Function< real_t >& src,
              const real_t&               beta,
              const P1Function< real_t >& dst,
              uint_t                      level,
              DoFType                     flag ) const override final;

   void applyScaled( const real_t&               alpha,
                     const P2Function< real_t >& src,
                     const P1Function< real_t >& dst,
                     uint_t                      level,
                     DoFType                     flag,
                     UpdateType                  updateType = Replace ) const override final;

   void apply( const P2Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final;

   /// Assemble operator as sparse matrix with scaling
   ///
   /// \param alpha constant scaling of the matrix
   /// \param mat   a sparse matrix proxy
   /// \param src   P2Function for determining column indices
   /// \param dst   P1Function for determining row indices
   /// \param level le2el in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void toMatrixScaled( const real_t&                               alpha,
                        const std::shared_ptr< SparseMatrixProxy >& mat,
                        const P2Function< idx_t >&                  src,
                        const P1Function< idx_t >&                  dst,
                        uint_t                                      level,
                        DoFType                                     flag ) const override;

   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param src   P2Function for determining column indices
   /// \param dst   P1Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  src,
                  const P1Function< idx_t >&                  dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override;

 private:
   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microFace      index associated with the current element = micro-face
   /// \param fType          type of micro-face (GRAY or BLUE)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param elMat          the 3x6 element matrix to be multiplied
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply2D( uint_t                 level,
                                     const indexing::Index& microFace,
                                     facedof::FaceType      fType,
                                     const real_t* const    srcVertexData,
                                     const real_t* const    srcEdgeData,
                                     real_t* const          dstVertexData,
                                     const Matrixr< 3, 6 >& elMat,
                                     const real_t&          alpha ) const;

   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microCell      index associated with the current element = micro-cell
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply3D( const uint_t            level,
                                     const indexing::Index&  microCell,
                                     const celldof::CellType cType,
                                     const real_t* const     srcVertexData,
                                     const real_t* const     srcEdgeData,
                                     real_t* const           dstVertexData,
                                     const Matrixr< 4, 10 >& elMat,
                                     const real_t&           alpha ) const;

   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const idx_t                                 xIdx,
                               const idx_t                                 yIdx,
                               const P2Elements::P2Element&                element,
                               const idx_t* const                          srcVertexIdx,
                               const idx_t* const                          srcEdgeIdx,
                               const idx_t* const                          dstVertexIdx,
                               const real_t&                               alpha ) const;

   void localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Cell&                                 cell,
                               const uint_t                                level,
                               const indexing::Index&                      microCell,
                               const celldof::CellType                     cType,
                               const idx_t* const                          srcVertexIdx,
                               const idx_t* const                          srcEdgeIdx,
                               const idx_t* const                          dstVertexIdx,
                               const real_t&                               alpha ) const;

   void assembleLocalElementMatrix2D( const Face&            face,
                                      uint_t                 level,
                                      const indexing::Index& microFace,
                                      facedof::FaceType      fType,
                                      P2toP1Form             form,
                                      Matrixr< 3, 6 >&       elMat ) const;

   void assembleLocalElementMatrix3D( const Cell&            cell,
                                      uint_t                 level,
                                      const indexing::Index& microCell,
                                      celldof::CellType      cType,
                                      P2toP1Form             form,
                                      Matrixr< 4, 10 >&      elMat ) const;

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrixr< 3, 6 >&
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
   const Matrixr< 3, 6 >&
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
   Matrixr< 4, 10 >&
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
   const Matrixr< 4, 10 >&
       localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells(), "Retriveing local element matrix for 3D in 2D run. Why?" )
      const auto idx = celldof::macrocell::index( level, microCell.x(), microCell.y(), microCell.z(), cType );
      WALBERLA_ASSERT( localElementMatrices3D_.count( cell.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices3D_.at( cell.getID() ).at( level ).empty() )
      return localElementMatrices3D_.at( cell.getID() ).at( level ).at( idx );
   }

   P2toP1Form form_;
   bool       localElementMatricesPrecomputed_;

   /// Pre-computed local element matrices.
   /// localElementMatrices_[macroCellID][level][faceIdx] = mat3x6
   std::map< PrimitiveID, std::map< uint_t, std::vector< Matrixr< 3, 6 >, Eigen::aligned_allocator< Matrixr< 3, 6 > > > > >
       localElementMatrices2D_;

   /// Pre-computed local element matrices.
   /// localElementMatrices_[macroCellID][level][cellIdx] = mat10x10
   std::map< PrimitiveID, std::map< uint_t, std::vector< Matrixr< 4, 10 >, Eigen::aligned_allocator< Matrixr< 4, 10 > > > > >
       localElementMatrices3D_;
};

typedef P2ToP1ElementwiseOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >
    P2ToP1ElementwiseDivxOperator;

typedef P2ToP1ElementwiseOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >
    P2ToP1ElementwiseDivyOperator;

typedef P2ToP1ElementwiseOperator< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >
    P2ToP1ElementwiseDivzOperator;

typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_div_0_blending_q2 > P2ToP1ElementwiseBlendingDivxOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_div_1_blending_q2 > P2ToP1ElementwiseBlendingDivyOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_div_2_blending_q2 > P2ToP1ElementwiseBlendingDivzOperator;

typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_0_affine_q2 > P2ToP1ElementwiseDivTxOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_1_affine_q2 > P2ToP1ElementwiseDivTyOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_2_affine_q2 > P2ToP1ElementwiseDivTzOperator;

typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_0_blending_q3 > P2ToP1ElementwiseBlendingDivTxOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_1_blending_q3 > P2ToP1ElementwiseBlendingDivTyOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_divT_2_blending_q3 > P2ToP1ElementwiseBlendingDivTzOperator;

typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_div_0_blending_q6 > P2ToP1ElementwiseBlendingDivxQFOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_div_1_blending_q6 > P2ToP1ElementwiseBlendingDivyQFOperator;
typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_div_2_blending_q6 > P2ToP1ElementwiseBlendingDivzQFOperator;

typedef P2ToP1ElementwiseOperator< forms::p2_to_p1_k_mass_blending_q5 > P2ToP1ElementwiseBlendingKMassOperator;

} // namespace hyteg
