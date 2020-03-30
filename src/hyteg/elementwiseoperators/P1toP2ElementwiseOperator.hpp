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
#include "hyteg/forms/form_fenics_generated/p1_to_p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_to_p2_tet_divt_tet.h"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormMass.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

using walberla::real_t;

template < class P1toP2Form >
class P1toP2ElementwiseOperator : public Operator< P1Function< real_t >, P2Function< real_t > >
{
 public:
   P1toP2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   void apply( const P1Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const;

#ifdef HYTEG_BUILD_WITH_PETSC
   /// Assemble operator as sparse matrix for PETSc
   ///
   /// \param mat   PETSc's own matrix data structure
   /// \param src   P2Function for determining column indices
   /// \param dst   P2Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void assembleLocalMatrix( Mat&                          mat,
                             const P1Function< PetscInt >& src,
                             const P2Function< PetscInt >& dst,
                             uint_t                        level,
                             DoFType                       flag ) const;
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
   void localMatrixVectorMultiply3D( const Cell&             cell,
                                     const uint_t            level,
                                     const indexing::Index&  microCell,
                                     const celldof::CellType cType,
                                     const real_t* const     srcVertexData,
                                     real_t* const           dstVertexData,
                                     real_t* const           dstEdgeData ) const;

#ifdef HYTEG_BUILD_WITH_PETSC
   void localMatrixAssembly2D( Mat&                         mat,
                               const Face&                  face,
                               const uint_t                 level,
                               const uint_t                 xIdx,
                               const uint_t                 yIdx,
                               const P2Elements::P2Element& element,
                               const PetscInt* const        srcVertexIdx,
                               const PetscInt* const        dstVertexIdx,
                               const PetscInt* const        dstEdgeIdx ) const;

   void localMatrixAssembly3D( Mat&                    mat,
                               const Cell&             cell,
                               const uint_t            level,
                               const indexing::Index&  microCell,
                               const celldof::CellType cType,
                               const PetscInt* const   srcVertexIdx,
                               const PetscInt* const   dstVertexIdx,
                               const PetscInt* const   dstEdgeIdx ) const;
#endif
};

typedef P1toP2ElementwiseOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >
    P1toP2ElementwiseDivTxOperator;

typedef P1toP2ElementwiseOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >
    P1toP2ElementwiseDivTyOperator;

typedef P1toP2ElementwiseOperator< P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >
    P1toP2ElementwiseDivTzOperator;

} // namespace hyteg
