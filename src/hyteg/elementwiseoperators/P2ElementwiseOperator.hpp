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
#pragma once

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p2_polar_laplacian.h"
#include "hyteg/forms/form_hyteg_generated/P2FormDivKGradBlending.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormMass.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

using walberla::real_t;

template < class P2Form >
class P2ElementwiseOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   P2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          const P2Form&                              form );

   P2ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          const P2Form&                              form,
                          bool                                       needsInverseDiagEntries );

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false ); }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() { computeDiagonalOperatorValues( true ); }

   std::shared_ptr< P2Function< real_t > > getDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          diagonalValues_,
          "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return diagonalValues_;
   };

   std::shared_ptr< P2Function< real_t > > getInverseDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          inverseDiagonalValues_,
          "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return inverseDiagonalValues_;
   };

   void smooth_jac( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const P2Function< real_t >& src,
                    real_t                      omega,
                    size_t                      level,
                    DoFType                     flag ) const;

   void smooth_gs( const P2Function< real_t >&, const P2Function< real_t >&, size_t, DoFType ) const
   {
      WALBERLA_ABORT( "Gauss-Seidel not implemented for P2ElementwiseOperator." )
   }

   void smooth_sor( const P2Function< real_t >&, const P2Function< real_t >&, real_t, size_t, DoFType ) const
   {
      WALBERLA_ABORT( "SOR not implemented for P2ElementwiseOperator." )
   }

#ifdef HYTEG_BUILD_WITH_PETSC
   /// Assemble operator as sparse matrix.
   ///
   /// \param mat   a sparse matrix proxy
   /// \param src   P2Function for determining column indices
   /// \param dst   P2Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void assembleLocalMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                             const P2Function< PetscInt >&               src,
                             const P2Function< PetscInt >&               dst,
                             uint_t                                      level,
                             DoFType                                     flag ) const;
#endif

   P2Form getForm() const;

 private:
   /// compute product of element local vector with element matrix
   ///
   /// \param face           face primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param xIdx           column index of vertex specifying micro-element
   /// \param yIdx           row index of vertex specifying micro-element
   /// \param element        element specification w.r.t. to micro-vertex
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
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
                                     const real_t* const          srcEdgeData,
                                     real_t* const                dstVertexData,
                                     real_t* const                dstEdgeData ) const;

   /// Compute contributions to operator diagonal for given micro-face
   ///
   /// \param face           face primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param xIdx           column index of vertex specifying micro-element
   /// \param yIdx           row index of vertex specifying micro-element
   /// \param element        element specification w.r.t. to micro-vertex
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   void computeLocalDiagonalContributions2D( const Face&                  face,
                                             const uint_t                 level,
                                             const uint_t                 xIdx,
                                             const uint_t                 yIdx,
                                             const P2Elements::P2Element& element,
                                             real_t* const                dstVertexData,
                                             real_t* const                dstEdgeData );

   /// Compute contributions to operator diagonal for given micro-face
   ///
   /// \param cell           cell primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microCell      index associated with the current element = micro-cell
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param vertexData     pointer to DoF data for storing diagonal values for vertex DoF stencils
   /// \param edgeData       pointer to DoF data for storing diagonal values for edge DoF stencils
   void computeLocalDiagonalContributions3D( const Cell&             cell,
                                             const uint_t            level,
                                             const indexing::Index&  microCell,
                                             const celldof::CellType cType,
                                             real_t* const           vertexData,
                                             real_t* const           edgeData );

#ifdef HYTEG_BUILD_WITH_PETSC
   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const uint_t                                xIdx,
                               const uint_t                                yIdx,
                               const P2Elements::P2Element&                element,
                               const PetscInt* const                       srcVertexIdx,
                               const PetscInt* const                       srcEdgeIdx,
                               const PetscInt* const                       dstVertexIdx,
                               const PetscInt* const                       dstEdgeIdx ) const;

   void localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Cell&                                 cell,
                               const uint_t                                level,
                               const indexing::Index&                      microCell,
                               const celldof::CellType                     cType,
                               const PetscInt* const                       srcVertexIdx,
                               const PetscInt* const                       srcEdgeIdx,
                               const PetscInt* const                       dstVertexIdx,
                               const PetscInt* const                       dstEdgeIdx ) const;
#endif

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert );

   std::shared_ptr< P2Function< real_t > > diagonalValues_;
   std::shared_ptr< P2Function< real_t > > inverseDiagonalValues_;

   P2Form form_;
};

/// compute product of element local vector with element matrix
///
/// \param cell           cell primitive we operate on
/// \param level          level on which we operate in mesh hierarchy
/// \param microCell      index associated with the current element = micro-cell
/// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
/// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
/// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
/// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
/// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
///
/// \note The src and dst data arrays must not be identical.
template < class P2Form >
void localMatrixVectorMultiply3D( const Cell&            cell,
                                  uint_t                 level,
                                  const indexing::Index& microCell,
                                  celldof::CellType      cType,
                                  const real_t* const    srcVertexData,
                                  const real_t* const    srcEdgeData,
                                  real_t* const          dstVertexData,
                                  real_t* const          dstEdgeData,
                                  P2Form                 form );

typedef P2ElementwiseOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >
    P2ElementwiseLaplaceOperator;

typedef P2ElementwiseOperator<
    P2FenicsForm< p2_polar_laplacian_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >
    P2ElementwisePolarLaplaceOperator;

typedef P2ElementwiseOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
    P2ElementwiseMassOperator;

typedef P2ElementwiseOperator< P2Form_mass >    P2ElementwiseBlendingMassOperator;
typedef P2ElementwiseOperator< P2Form_laplace > P2ElementwiseBlendingLaplaceOperator;

typedef P2ElementwiseOperator< P2Form_divKgrad > P2ElementwiseDivKGradOperator;
typedef P2ElementwiseOperator< P2Form_divKgradBlending > P2ElementwiseDivKGradBlendingOperator;

typedef P2ElementwiseOperator< P2LinearCombinationForm > P2ElementwiseLinearCombinationOperator;

} // namespace hyteg
