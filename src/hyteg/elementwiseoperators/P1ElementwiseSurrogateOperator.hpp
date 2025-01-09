/*
 * Copyright (c) 2017-2025 Marcus Mohr, Nils Kohl, Benjamin Mann.
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

#include <hyteg/polynomial/elementwise/leastSquares.hpp>
#include <hyteg/polynomial/elementwise/polynomial.hpp>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_k_mass_affine_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/polynomial/elementwise/leastSquares.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

using walberla::real_t;

template < class P1Form >
class P1ElementwiseSurrogateOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                                       public WeightedJacobiSmoothable< P1Function< real_t > >,
                                       public OperatorWithInverseDiagonal< P1Function< real_t > >
{
   static constexpr uint_t n_loc( uint_t dim ) { return ( D == 2 ) ? 3 : 4; }

   using Polynomial = surrogate::polynomial::Polynomial;
   template < size_t D >
   using P_matrix = surrogate::polynomial::Matrix< n_loc( D ) >;
   template < size_t D >
   using RHS_matrix = std::array< std::array< surrogate::LeastSquares::Vector, n_loc( D ) >, n_loc( D ) >;

 public:
   P1ElementwiseSurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   P1ElementwiseSurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                   size_t                                     minLevel,
                                   size_t                                     maxLevel,
                                   const P1Form&                              form );

   P1ElementwiseSurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                   size_t                                     minLevel,
                                   size_t                                     maxLevel,
                                   const P1Form&                              form,
                                   bool                                       needsInverseDiagEntries );

   void init( size_t             poly_degree,
              size_t             downsampling            = 1,
              const std::string& path_to_svd             = "",
              bool               needsInverseDiagEntries = true );

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final;

   void gemv( const real_t&               alpha,
              const P1Function< real_t >& src,
              const real_t&               beta,
              const P1Function< real_t >& dst,
              size_t                      level,
              DoFType                     flag ) const override final;

   void smooth_jac( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    const P1Function< real_t >& src,
                    real_t                      omega,
                    size_t                      level,
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

   /// \brief Pre-computes the local stiffness matrices for each (micro-)element and stores them all in memory.
   ///
   /// If this method is called, all subsequent calls to apply() or smooth_*() use the stored element matrices.
   /// If the local element matrices need to be recomputed again, simply call this method again.
   void computeAndStoreLocalElementMatrices();

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false ); }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() override final { computeDiagonalOperatorValues( true ); }

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

 private:
   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microFace      index associated with the current element = micro-face
   /// \param fType          type of micro-face (GRAY or BLUE)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param elMat          the 3x3 element matrix to be multiplied
   /// \param alpha          scaling factor that is applied to the local result vector
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply2D( const uint_t           level,
                                     const indexing::Index& microFace,
                                     facedof::FaceType      fType,
                                     const real_t* const    srcVertexData,
                                     real_t* const          dstVertexData,
                                     const Matrix3r&        elMat,
                                     const real_t&          alpha ) const;

   /// compute product of element local vector with element matrix
   ///
   /// \param level          level on which we operate in mesh hierarchy
   /// \param microCell      index associated with the current element = micro-cell
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param elMat          the 4x4 element matrix to be multiplied
   /// \param alpha          scaling factor that is applied to the local result vector
   ///
   /// \note The src and dst data arrays must not be identical.
   void localMatrixVectorMultiply3D( const uint_t            level,
                                     const indexing::Index&  microCell,
                                     const celldof::CellType cType,
                                     const real_t* const     srcVertexData,
                                     real_t* const           dstVertexData,
                                     const Matrix4r&         elMat,
                                     const real_t&           alpha ) const;

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
                                             real_t* const                              dstVertexData );

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
                                             real_t* const           vertexData );

   void precompute_local_stiffness_2d( uint_t lvl );
   void precompute_local_stiffness_3d( uint_t lvl );
   void compute_local_surrogates_2d( uint_t lvl );
   void compute_local_surrogates_3d( uint_t lvl );

   void localMatrixAssembly2D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Face&                                 face,
                               const uint_t                                level,
                               const idx_t                                 xIdx,
                               const idx_t                                 yIdx,
                               const P1Elements::P1Elements2D::P1Element&  element,
                               const idx_t* const                          srcIdx,
                               const idx_t* const                          dstIdx ) const;

   void localMatrixAssembly3D( const std::shared_ptr< SparseMatrixProxy >& mat,
                               const Cell&                                 cell,
                               const uint_t                                level,
                               const indexing::Index&                      microCell,
                               const celldof::CellType                     cType,
                               const idx_t* const                          srcIdx,
                               const idx_t* const                          dstIdx ) const;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert );

   bool is_initialized_;

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;

   P1Form form_;

   // least squares approximator for each level
   std::vector< std::shared_ptr< surrogate::LeastSquares > > lsq_;
   std::vector<uint_t> downsampling_;

   // polynomial degree for each level
   std::vector<uint_t> poly_degree_;


   // memory IDs of local stiffness matrices for level 1-3
   PrimitiveDataID< LevelWiseMemory< Matrix3r >, Face > stiffnessID_2d_;
   PrimitiveDataID< LevelWiseMemory< Matrix4r >, Cell > stiffnessID_2d_;
   // memory IDs of surrogates for level 4+ (one poly matrix for each element type)
   PrimitiveDataID< LevelWiseMemory< std::array< P_matrix< 2 >, 2 > >, Face > surrogateID_2d_;
   PrimitiveDataID< LevelWiseMemory< std::array< P_matrix< 2 >, 6 > >, Cell > surrogateID_3d_;
};

template < class P1Form >
void assembleLocalElementMatrix2D( const Face&            face,
                                   uint_t                 level,
                                   const indexing::Index& microFace,
                                   facedof::FaceType      fType,
                                   P1Form                 form,
                                   Matrix3r&              elMat )
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

template < class P1Form >
void assembleLocalElementMatrix3D( const Cell&            cell,
                                   uint_t                 level,
                                   const indexing::Index& microCell,
                                   celldof::CellType      cType,
                                   P1Form                 form,
                                   Matrix4r&              elMat )
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

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >
    P1ElementwiseLaplaceOperator;

typedef P1ElementwiseSurrogateOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > >
    P1ElementwisePolarLaplaceOperator;

typedef P1ElementwiseSurrogateOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >
    P1ElementwiseMassOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_mass_blending_q4 > P1ElementwiseBlendingMassOperator;

typedef P1ElementwiseSurrogateOperator< P1LinearCombinationForm > P1ElementwiseLinearCombinationOperator;

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >
    P1ElementwisePSPGOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_diffusion_blending_q3 > P1ElementwiseBlendingLaplaceOperator;
typedef P1ElementwiseSurrogateOperator< forms::p1_diffusion_blending_q2 > P1ElementwiseBlendingLaplaceOperatorQ2;

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >
    P1ElementwiseDivXOperator;
typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >
    P1ElementwiseDivYOperator;
typedef P1ElementwiseSurrogateOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > >
    P1ElementwiseDivZOperator;

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >
    P1ElementwiseDivTXOperator;
typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >
    P1ElementwiseDivTYOperator;
typedef P1ElementwiseSurrogateOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > >
    P1ElementwiseDivTZOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3 >   P1ElementwiseAffineDivKGradOperator;
typedef P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_blending_q3 > P1ElementwiseBlendingDivKGradOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_k_mass_affine_q4 > P1ElementwiseKMassOperator;

} // namespace hyteg
