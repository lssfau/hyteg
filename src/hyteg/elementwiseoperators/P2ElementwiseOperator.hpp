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
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4244) // should be fixed otherwise... 
#endif 
#include "hyteg/forms/form_fenics_generated/p2_polar_laplacian.h"
#ifdef _MSC_VER
#pragma warning( pop )
#endif
#include "hyteg/forms/form_hyteg_generated/p2/p2_div_k_grad_affine_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_full_stokes_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q5.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"



namespace hyteg {

using walberla::real_t;

template < class P2Form >
class P2ElementwiseOperator : public Operator< P2Function< real_t >, P2Function< real_t > >,
                              public WeightedJacobiSmoothable< P2Function< real_t > >,
                              public OperatorWithInverseDiagonal< P2Function< real_t > >
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
               UpdateType                  updateType = Replace ) const override final;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false ); }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() override final { computeDiagonalOperatorValues( true ); }

   std::shared_ptr< P2Function< real_t > > getDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          diagonalValues_,
          "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return diagonalValues_;
   };

   std::shared_ptr< P2Function< real_t > > getInverseDiagonalValues() const override
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          inverseDiagonalValues_,
          "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return inverseDiagonalValues_;
   };

   /// \brief Pre-computes the local stiffness matrices for each (micro-)element and stores them all in memory.
   ///
   /// If this method is called, all subsequent calls to apply() or smooth_*() use the stored element matrices.
   /// If the local element matrices need to be recomputed again, simply call this method again.
   void computeAndStoreLocalElementMatrices();

   void smooth_jac( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const P2Function< real_t >& src,
                    real_t                      omega,
                    size_t                      level,
                    DoFType                     flag ) const override;

   void smooth_gs( const P2Function< real_t >&, const P2Function< real_t >&, size_t, DoFType ) const
   {
      WALBERLA_ABORT( "Gauss-Seidel not implemented for P2ElementwiseOperator." )
   }

   void smooth_sor( const P2Function< real_t >&, const P2Function< real_t >&, real_t, size_t, DoFType ) const
   {
      WALBERLA_ABORT( "SOR not implemented for P2ElementwiseOperator." )
   }

   /// Assemble operator as sparse matrix.
   ///
   /// \param mat   a sparse matrix proxy
   /// \param src   P2Function for determining column indices
   /// \param dst   P2Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  src,
                  const P2Function< idx_t >&                  dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override;

   P2Form getForm() const;

 private:
   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microFace      index associated with the current element = micro-face
   /// \param fType          type of micro-face (GRAY or BLUE)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
   /// \param elMat          the 6x6 element matrix to be multiplied
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply2D( uint_t                 level,
                                     const indexing::Index& microFace,
                                     facedof::FaceType      fType,
                                     const real_t* const    srcVertexData,
                                     const real_t* const    srcEdgeData,
                                     real_t* const          dstVertexData,
                                     real_t* const          dstEdgeData,
                                     const Matrix6r&        elMat ) const;

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

   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const uint_t                                xIdx,
                               const uint_t                                yIdx,
                               const P2Elements::P2Element&                element,
                               const idx_t* const                          srcVertexIdx,
                               const idx_t* const                          srcEdgeIdx,
                               const idx_t* const                          dstVertexIdx,
                               const idx_t* const                          dstEdgeIdx ) const;

   void localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Cell&                                 cell,
                               const uint_t                                level,
                               const indexing::Index&                      microCell,
                               const celldof::CellType                     cType,
                               const idx_t* const                          srcVertexIdx,
                               const idx_t* const                          srcEdgeIdx,
                               const idx_t* const                          dstVertexIdx,
                               const idx_t* const                          dstEdgeIdx ) const;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert );

   std::shared_ptr< P2Function< real_t > > diagonalValues_;
   std::shared_ptr< P2Function< real_t > > inverseDiagonalValues_;

   P2Form form_;

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrix6r& localElementMatrix2D( const Face& face, uint_t level, const indexing::Index& microFace, facedof::FaceType fType )
   {
      WALBERLA_ASSERT( !storage_->hasGlobalCells(), "Retriveing local element matrix for 2D in 3D run. Why?" )
      const auto idx = facedof::macroface::index( level, microFace.x(), microFace.y(), fType );
      WALBERLA_ASSERT( localElementMatrices2D_.count( face.getID().getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID().getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID().getID() ).at( level ).size() > 0 )
      return localElementMatrices2D_[face.getID().getID()][level][idx];
   }

   /// \brief Returns a const reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   const Matrix6r&
       localElementMatrix2D( const Face& face, uint_t level, const indexing::Index& microFace, facedof::FaceType fType ) const
   {
      WALBERLA_ASSERT( !storage_->hasGlobalCells(), "Retriveing local element matrix for 2D in 3D run. Why?" )
      const auto idx = facedof::macroface::index( level, microFace.x(), microFace.y(), fType );
      WALBERLA_ASSERT( localElementMatrices2D_.count( face.getID().getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID().getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( localElementMatrices2D_.at( face.getID().getID() ).at( level ).size() > 0 )
      return localElementMatrices2D_.at( face.getID().getID() ).at( level ).at( idx );
   }

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro cell.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrix10r& localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType )
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
   const Matrix10r&
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
   /// localElementMatrices2D_[macroCellID][level][cellIdx] = mat6x6
   std::map< PrimitiveID::IDType, std::map< uint_t, std::vector< Matrix6r > > > localElementMatrices2D_;

   /// Pre-computed local element matrices.
   /// localElementMatrices3D_[macroCellID][level][cellIdx] = mat10x10
   std::map< PrimitiveID::IDType, std::map< uint_t, std::vector< Matrix10r > > > localElementMatrices3D_;
};

template < class P2Form >
void assembleLocalElementMatrix2D( const Face&            face,
                                   uint_t                 level,
                                   const indexing::Index& microFace,
                                   facedof::FaceType      fType,
                                   P2Form                 form,
                                   Matrix6r&              elMat )
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

template < class P2Form >
void assembleLocalElementMatrix3D( const Cell&            cell,
                                   uint_t                 level,
                                   const indexing::Index& microCell,
                                   celldof::CellType      cType,
                                   P2Form                 form,
                                   Matrix10r&             elMat )
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

/// compute product of element local vector with element matrix
///
/// \param level          level on which we operate in mesh hierarchy
/// \param microCell      index associated with the current element = micro-cell
/// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
/// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
/// \param srcEdgeData    pointer to DoF data on micro-edges (for reading data)
/// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
/// \param dstEdgeData    pointer to DoF data on micro-edges (for writing data)
/// \param elMat          the 10x10 element matrix to be multiplied
///
/// \note The src and dst data arrays must not be identical.
void localMatrixVectorMultiply3D( uint_t                 level,
                                  const indexing::Index& microCell,
                                  celldof::CellType      cType,
                                  const real_t* const    srcVertexData,
                                  const real_t* const    srcEdgeData,
                                  real_t* const          dstVertexData,
                                  real_t* const          dstEdgeData,
                                  const Matrix10r&       elMat );

typedef P2ElementwiseOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >
    P2ElementwiseLaplaceOperator;

typedef P2ElementwiseOperator<
    P2FenicsForm< p2_polar_laplacian_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >
    P2ElementwisePolarLaplaceOperator;

typedef P2ElementwiseOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
    P2ElementwiseMassOperator;

typedef P2ElementwiseOperator< forms::p2_mass_blending_q5 >    P2ElementwiseBlendingMassOperator;
typedef P2ElementwiseOperator< P2Form_laplace > P2ElementwiseBlendingLaplaceOperator;

typedef P2ElementwiseOperator< P2Form_divKgrad > P2ElementwiseDivKGradOperator;

typedef P2ElementwiseOperator< P2LinearCombinationForm > P2ElementwiseLinearCombinationOperator;

typedef P2ElementwiseOperator< forms::p2_div_k_grad_affine_q4 > P2ElementwiseAffineDivKGradOperator;

} // namespace hyteg
