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
#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_k_mass_affine_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_k_mass_centroid_blending_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_manifold_mass_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_pspg_blending_q2.hpp"
#include "hyteg/forms/form_hyteg_manual/p1_neighbour_form.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

using walberla::real_t;

template < class P1Form >
class P1ElementwiseOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                              public WeightedJacobiSmoothable< P1Function< real_t > >,
                              public OperatorWithInverseDiagonal< P1Function< real_t > >
{
 public:
   P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel );

   P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          uint_t                                     minLevel,
                          uint_t                                     maxLevel,
                          const P1Form&                              form );

   P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          uint_t                                     minLevel,
                          uint_t                                     maxLevel,
                          const P1Form&                              form,
                          bool                                       needsInverseDiagEntries );

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               uint_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final;

   void applyScaled( const real_t&               alpha,
                     const P1Function< real_t >& src,
                     const P1Function< real_t >& dst,
                     uint_t                      level,
                     DoFType                     flag,
                     UpdateType                  updateType = Replace ) const override final;

   void gemv( const real_t&               alpha,
              const P1Function< real_t >& src,
              const real_t&               beta,
              const P1Function< real_t >& dst,
              uint_t                      level,
              DoFType                     flag ) const override final;

   void smooth_jac_scaled( const real_t&               alpha,
                           const P1Function< real_t >& dst,
                           const P1Function< real_t >& rhs,
                           const P1Function< real_t >& src,
                           real_t                      omega,
                           uint_t                      level,
                           DoFType                     flag ) const override;

   void smooth_jac( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    const P1Function< real_t >& src,
                    real_t                      omega,
                    uint_t                      level,
                    DoFType                     flag ) const override;

   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param src   P1Function for determining column indices
   /// \param dst   P1Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< idx_t >&                  src,
                  const P1Function< idx_t >&                  dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override;

   /// Assemble operator as sparse matrix with scaling
   ///
   /// \param alpha constant scaling of the matrix
   /// \param mat   a sparse matrix proxy
   /// \param src   P1Function for determining column indices
   /// \param dst   P1Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  ignored
   ///
   /// \note src and dst are legal to and often will be the same function object
   void toMatrixScaled( const real_t&                               alpha,
                        const std::shared_ptr< SparseMatrixProxy >& mat,
                        const P1Function< idx_t >&                  src,
                        const P1Function< idx_t >&                  dst,
                        uint_t                                      level,
                        DoFType                                     flag ) const override;

   /// \brief Pre-computes the local stiffness matrices for each (micro-)element and stores them all in memory.
   ///
   /// If this method is called, all subsequent calls to apply() or smooth_*() use the stored element matrices.
   /// If the local element matrices need to be recomputed again, simply call this method again.
   void computeAndStoreLocalElementMatrices();

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false, false, real_c( 1 ) ); }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() override final { computeDiagonalOperatorValues( true, false, real_c( 1 ) ); }

   void computeLumpedDiagonalOperatorValues() { computeDiagonalOperatorValues( false, true, real_c( 1 ) ); }
   void computeLumpedInverseDiagonalOperatorValues() { computeDiagonalOperatorValues( true, true, real_c( 1 ) ); }

   /// Trigger (re)computation of scaled inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValuesScaled( const real_t& alpha ) override final
   {
      computeDiagonalOperatorValues( true, false, alpha );
   }

   void computeLumpedInverseDiagonalOperatorValuesScaled( const real_t& alpha )
   {
      computeDiagonalOperatorValues( true, true, alpha );
   }

   /// Trigger (re)computation of scaled diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValuesScaled( const real_t& alpha ) { computeDiagonalOperatorValues( false, false, alpha ); }

   void computeLumpedDiagonalOperatorValuesScaled( const real_t& alpha ) { computeDiagonalOperatorValues( false, true, alpha ); }

   std::shared_ptr< P1Function< real_t > > getDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          diagonalValues_,
          "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return diagonalValues_;
   };

   std::shared_ptr< P1Function< real_t > > getInverseDiagonalValues() const override
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          inverseDiagonalValues_,
          "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return inverseDiagonalValues_;
   };

   std::shared_ptr< P1Function< real_t > > getLumpedDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          lumpedDiagonalValues_,
          "Lumped diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return lumpedDiagonalValues_;
   };

   std::shared_ptr< P1Function< real_t > > getLumpedInverseDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          lumpedInverseDiagonalValues_,
          "Lumped inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return lumpedInverseDiagonalValues_;
   };

 private:
   /// Compute contributions to operator diagonal for given micro-face
   ///
   /// \param face           face primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param xIdx           column index of vertex specifying micro-element
   /// \param yIdx           row index of vertex specifying micro-element
   /// \param element        element specification w.r.t. to micro-vertex
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   void computeLocalDiagonalContributions2D( const Face&                                face,
                                             const uint_t                               level,
                                             const idx_t                                xIdx,
                                             const idx_t                                yIdx,
                                             const P1Elements::P1Elements2D::P1Element& element,
                                             real_t* const                              dstVertexData,
                                             const real_t&                              alpha,
                                             bool                                       lumped = false );

   /// Compute contributions to operator diagonal for given micro-face
   ///
   /// \param cell           cell primitive we operate on
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microCell      index associated with the current element = micro-cell
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param vertexData     pointer to DoF data for storing diagonal values for vertex DoF stencils
   void computeLocalDiagonalContributions3D( const Cell&             cell,
                                             const uint_t            level,
                                             const indexing::Index&  microCell,
                                             const celldof::CellType cType,
                                             real_t* const           vertexData,
                                             const real_t&           alpha,
                                             bool                    lumped = false );

   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const idx_t                                 xIdx,
                               const idx_t                                 yIdx,
                               const P1Elements::P1Elements2D::P1Element&  element,
                               const idx_t* const                          srcIdx,
                               const idx_t* const                          dstIdx,
                               const real_t&                               alpha ) const;

   void localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Cell&                                 cell,
                               const uint_t                                level,
                               const indexing::Index&                      microCell,
                               const celldof::CellType                     cType,
                               const idx_t* const                          srcIdx,
                               const idx_t* const                          dstIdx,
                               const real_t&                               alpha ) const;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert, bool lumped, const real_t& alpha );

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;
   std::shared_ptr< P1Function< real_t > > lumpedDiagonalValues_;
   std::shared_ptr< P1Function< real_t > > lumpedInverseDiagonalValues_;

   P1Form form_;

   /// \brief Returns a reference to the a precomputed element matrix of the specified micro face.
   /// Probably crashes if local element matrices have not been precomputed.
   Matrix3r& localElementMatrix2D( const Face& face, uint_t level, const indexing::Index& microFace, facedof::FaceType fType )
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
   const Matrix3r&
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
   Matrix4r& localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType )
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
   const Matrix4r&
       localElementMatrix3D( const Cell& cell, uint_t level, const indexing::Index& microCell, celldof::CellType cType ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells(), "Retriveing local element matrix for 3D in 2D run. Why?" )
      const auto idx = celldof::macrocell::index( level, microCell.x(), microCell.y(), microCell.z(), cType );
      WALBERLA_ASSERT( localElementMatrices3D_.count( cell.getID() ) > 0 )
      WALBERLA_ASSERT( localElementMatrices3D_.at( cell.getID() ).count( level ) > 0 )
      WALBERLA_ASSERT( !localElementMatrices3D_.at( cell.getID() ).at( level ).empty() )
      return localElementMatrices3D_.at( cell.getID() ).at( level ).at( idx );
   }

   bool localElementMatricesPrecomputed_;

   /// Pre-computed local element matrices.
   /// localElementMatrices2D_[macroCellID][level][cellIdx] = mat3x3
   std::map< PrimitiveID, std::map< uint_t, std::vector< Matrix3r, Eigen::aligned_allocator< Matrix3r > > > >
       localElementMatrices2D_;

   /// Pre-computed local element matrices.
   /// localElementMatrices3D_[macroCellID][level][cellIdx] = mat4x4
   std::map< PrimitiveID, std::map< uint_t, std::vector< Matrix4r, Eigen::aligned_allocator< Matrix4r > > > >
       localElementMatrices3D_;
};

typedef P1ElementwiseOperator<
    P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >
    P1ElementwiseLaplaceOperator;

typedef P1ElementwiseOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > > P1ElementwisePolarLaplaceOperator;

typedef P1ElementwiseOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >
    P1ElementwiseMassOperator;

typedef P1ElementwiseOperator< forms::p1_mass_blending_q4 > P1ElementwiseBlendingMassOperator;

typedef P1ElementwiseOperator< P1LinearCombinationForm > P1ElementwiseLinearCombinationOperator;

typedef P1ElementwiseOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >
    P1ElementwisePSPGOperator;

typedef P1ElementwiseOperator< forms::p1_diffusion_blending_q3 > P1ElementwiseBlendingLaplaceOperator;
typedef P1ElementwiseOperator< forms::p1_diffusion_blending_q2 > P1ElementwiseBlendingLaplaceOperatorQ2;

typedef P1ElementwiseOperator< P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >
    P1ElementwiseDivXOperator;
typedef P1ElementwiseOperator< P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >
    P1ElementwiseDivYOperator;
typedef P1ElementwiseOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > >
    P1ElementwiseDivZOperator;

typedef P1ElementwiseOperator< P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >
    P1ElementwiseDivTXOperator;
typedef P1ElementwiseOperator< P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >
    P1ElementwiseDivTYOperator;
typedef P1ElementwiseOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > >
    P1ElementwiseDivTZOperator;

typedef P1ElementwiseOperator< forms::p1_div_k_grad_affine_q3 >   P1ElementwiseAffineDivKGradOperator;
typedef P1ElementwiseOperator< forms::p1_div_k_grad_blending_q3 > P1ElementwiseBlendingDivKGradOperator;

typedef P1ElementwiseOperator< forms::p1_k_mass_affine_q4 > P1ElementwiseKMassOperator;

typedef P1ElementwiseOperator< forms::p1_k_mass_centroid_blending_q4 > P1ElementwiseBlendingKMassOperator_Centroid;

typedef P1ElementwiseOperator< forms::p1_pspg_blending_q2 > P1ElementwiseBlendingPSPGOperator;

typedef P1ElementwiseOperator< forms::p1_neighbour_form > P1ElementwiseNeighbourOperator;

typedef P1ElementwiseOperator< forms::p1_manifold_mass_blending_q3 > P1ElementwiseManifoldBlendingMassOperator;

} // namespace hyteg
