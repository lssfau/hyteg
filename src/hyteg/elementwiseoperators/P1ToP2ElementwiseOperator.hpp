/*
 * Copyright (c) 2017-2020 Marcus Mohr, Nils Kohl.
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

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_fenics_base/P1ToP2FenicsForm.hpp"
#include "hyteg/forms/form_hyteg_manual/P1ToP2FormDivT.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormMass.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "hyteg/forms/form_fenics_generated/p1_to_p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_to_p2_tet_divt_tet.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {

using walberla::real_t;

template < class P1toP2Form >
class P1ToP2ElementwiseOperator : public Operator< P1Function< real_t >, P2Function< real_t > >
{
 public:
   P1ToP2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   /// \brief Pre-computes the local stiffness matrices for each (micro-)element and stores them all in memory.
   ///
   /// If this method is called, all subsequent calls to apply() or smooth_*() use the stored element matrices.
   /// If the local element matrices need to be recomputed again, simply call this method again.
   void computeAndStoreLocalElementMatrices();

   void apply( const P1Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const;

#ifdef HYTEG_BUILD_WITH_PETSC
   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param src   P1Function for determining column indices
   /// \param dst   P2Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void assembleLocalMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                             const P1Function< PetscInt >&               src,
                             const P2Function< PetscInt >&               dst,
                             uint_t                                      level,
                             DoFType                                     flag ) const;
#endif

 private:
   /// compute product of element local vector with element matrix
   ///
   /// \param face           face primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param xIdx           column index of vertex specifying micro-element
   /// \param yIdx           row index of vertex specifying micro-element
   /// \param element        element specification w.r.t. to micro-vertex
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply2D( const Face&                  face,
                                     const uint_t                 level,
                                     const uint_t                 xIdx,
                                     const uint_t                 yIdx,
                                     const P2Elements::P2Element& element,
                                     const real_t* const          srcVertexData,
                                     real_t* const                dstVertexData,
                                     real_t* const                dstEdgeData ) const;

   /// compute product of element local vector with element matrix
   ///
   /// \param cell           cell primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microCell      index associated with the current element = micro-cell
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply3D( const uint_t            level,
                                     const indexing::Index&  microCell,
                                     const celldof::CellType cType,
                                     const real_t* const     srcVertexData,
                                     real_t* const           dstVertexData,
                                     real_t* const           dstEdgeData,
                                     const Matrixr< 10, 4 >& elMat ) const;

#ifdef HYTEG_BUILD_WITH_PETSC
   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const uint_t                                xIdx,
                               const uint_t                                yIdx,
                               const P2Elements::P2Element&                element,
                               const PetscInt* const                       srcVertexIdx,
                               const PetscInt* const                       dstVertexIdx,
                               const PetscInt* const                       dstEdgeIdx ) const;

   void localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Cell&                                 cell,
                               const uint_t                                level,
                               const indexing::Index&                      microCell,
                               const celldof::CellType                     cType,
                               const PetscInt* const                       srcVertexIdx,
                               const PetscInt* const                       dstVertexIdx,
                               const PetscInt* const                       dstEdgeIdx ) const;
#endif

   void assembleLocalElementMatrix3D( const Cell&            cell,
                                      uint_t                 level,
                                      const indexing::Index& microCell,
                                      celldof::CellType      cType,
                                      P1toP2Form             form,
                                      Matrixr< 10, 4 >&      elMat ) const;

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro cell.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrixr< 10, 4 >&
       localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType )
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells(), "Retriveing local element matrix for 3D in 2D run. Why?" )
      const auto idx = celldof::macrocell::index( level, microCell.x(), microCell.y(), microCell.z(), cType );
      WALBERLA_ASSERT( localElementMatrices3D_.count( cell.getID().getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID().getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID().getID() ).at( level ).size() > 0 )
      return localElementMatrices3D_[cell.getID().getID()][level][idx];
   }

   /// \brief Returns a const reference to the a precomputed element matrix of the specified micro cell.
   /// Probably crashes if local element matrices have not been precomputed.
   const Matrixr< 10, 4 >&
       localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells(), "Retriveing local element matrix for 3D in 2D run. Why?" )
      const auto idx = celldof::macrocell::index( level, microCell.x(), microCell.y(), microCell.z(), cType );
      WALBERLA_ASSERT( localElementMatrices3D_.count( cell.getID().getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID().getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID().getID() ).at( level ).size() > 0 )
      return localElementMatrices3D_.at( cell.getID().getID() ).at( level ).at( idx );
   }

   bool localElementMatricesPrecomputed_;

   /// Pre-computed local element matrices.
   /// localElementMatrices_[macroCellID][level][cellIdx] = mat10x10
   std::map< PrimitiveID::IDType, std::map< uint_t, std::vector< Matrixr< 10, 4 > > > > localElementMatrices3D_;
};

typedef P1ToP2ElementwiseOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >
    P1ToP2ElementwiseDivTxOperator;

typedef P1ToP2ElementwiseOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >
    P1ToP2ElementwiseDivTyOperator;

typedef P1ToP2ElementwiseOperator< P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >
    P1ToP2ElementwiseDivTzOperator;

typedef P1ToP2ElementwiseOperator< P1ToP2Form_divt< 0 > > P1ToP2ElementwiseBlendingDivTxOperator;
typedef P1ToP2ElementwiseOperator< P1ToP2Form_divt< 1 > > P1ToP2ElementwiseBlendingDivTyOperator;
typedef P1ToP2ElementwiseOperator< P1ToP2Form_divt< 2 > > P1ToP2ElementwiseBlendingDivTzOperator;

} // namespace hyteg
