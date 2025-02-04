/*
 * Copyright (c) 2025 Benjamin Mann.
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

#include <hyteg/communication/Syncing.hpp>
#include <hyteg/forms/P1LinearCombinationForm.hpp>
#include <hyteg/forms/form_fenics_base/P1FenicsForm.hpp>
#include <hyteg/forms/form_fenics_generated/p1_polar_laplacian.h>
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q2.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_k_mass_affine_q4.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp>
#include <hyteg/operators/Operator.hpp>
#include <hyteg/p1functionspace/P1Elements.hpp>
#include <hyteg/p1functionspace/P1Function.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>
#include <hyteg/polynomial/elementwise/data.hpp>
#include <hyteg/polynomial/elementwise/leastSquares.hpp>
#include <hyteg/polynomial/elementwise/polynomial.hpp>
#include <hyteg/solvers/Smoothables.hpp>
#include <hyteg/sparseassembly/SparseMatrixProxy.hpp>
#include <hyteg/volumedofspace/CellDoFIndexing.hpp>

namespace hyteg {

using walberla::real_t;

template < class P1Form >
class P1ElementwiseSurrogateOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                                       public WeightedJacobiSmoothable< P1Function< real_t > >,
                                       public OperatorWithInverseDiagonal< P1Function< real_t > >
{
   /* On lower levels, storing and evaluating polynomials is significantly less performant.
      Therefore, on levels 0-3 we precompute and store the system matrices, while we use
      surrogates for levels 4+
    */
   static constexpr uint_t min_lvl_for_surrogate = 4;

   /* Single precision LSQ leads to very poor accuracy of the resulting polynomials. Therefore,
      we use real_t for the polynomial evaluation only, while sticking to double precision LSQ.
    */
   using LSQ        = surrogate::LeastSquares< double >;
   using Poly       = surrogate::polynomial::Polynomial< real_t >;
   using PolyDomain = surrogate::polynomial::Domain< real_t >;

   template < uint_t DIM >
   using RHS_matrix = surrogate::RHS_matrix< real_t, DIM, 1, 1 >;
   template < uint_t DIM >
   using PrecomputedData = surrogate::PrecomputedData< real_t, DIM, 1, 1 >;
   template < uint_t DIM >
   using SurrogateData = surrogate::SurrogateData< real_t, DIM, 1, 1 >;

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

   /**
     * @brief Initializes the surrogate polynomials using an lsq-fit
     *
     * @param poly_degree The polynomial degree to be used.
     * @param downsampling The downsampling factor to be applied. Default is 0 (auto).
     * @param path_to_svd The file path to the SVD data. Default is an empty string (compute SVD on first call).
     * @param needsInverseDiagEntries Flag indicating whether inverse diagonal entries are needed. Default is true.
     */
   void init( uint8_t            poly_degree,
              size_t             downsampling            = 0,
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
   /// compute matrix vector product on all micro faces of given type
   ///
   /// \param face           macro face
   /// \param level          level on which we operate in mesh hierarchy
   /// \param fType          type of micro-face (GRAY or BLUE)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param alpha          scaling factor that is applied to the local result vector
   void apply_2d( const Face&             face,
                  const uint_t            level,
                  const facedof::FaceType ftype,
                  const real_t* const     srcVertexData,
                  real_t* const           dstVertexData,
                  const real_t&           alpha ) const;

   /// compute matrix vector product on all micro cells of given type
   ///
   /// \param cell           macro cell
   /// \param level          level on which we operate in mesh hierarchy
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   /// \param alpha          scaling factor that is applied to the local result vector
   void apply_3d( const Cell&             cell,
                  const uint_t            level,
                  const celldof::CellType cType,
                  const real_t* const     srcVertexData,
                  real_t* const           dstVertexData,
                  const real_t&           alpha ) const;

   /// compute contributions to diagonal of on all micro faces of given type
   ///
   /// \param face           macro face
   /// \param level          level on which we operate in mesh hierarchy
   /// \param fType          type of micro-face (GRAY or BLUE)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   void diagonal_contributions_2d( const Face&             face,
                                   const uint_t            level,
                                   const facedof::FaceType fType,
                                   real_t* const           dstVertexData );

   /// compute contributions to diagonal of on all micro cells of given type
   ///
   /// \param cell           macro cell
   /// \param level          level on which we operate in mesh hierarchy
   /// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
   /// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
   void diagonal_contributions_3d( const Cell&             cell,
                                   const uint_t            level,
                                   const celldof::CellType cType,
                                   real_t* const           dstVertexData );

   void precompute_local_stiffness_2d( uint_t lvl );
   void precompute_local_stiffness_3d( uint_t lvl );
   void compute_local_surrogates_2d( uint_t lvl );
   void compute_local_surrogates_3d( uint_t lvl );

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert );

   P1Form form_;

   bool is_initialized_;

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;

   // least squares approximator for each level
   std::vector< std::shared_ptr< LSQ > > lsq_;
   std::vector< uint_t >                 downsampling_;

   // polynomial degree for each level
   std::vector< uint8_t > poly_degree_;

   // precomputed local stiffness matrices for level 1-3
   PrecomputedData< 2 > a_loc_2d_;
   PrecomputedData< 3 > a_loc_3d_;
   // surrogates for level 4+ (one poly matrix for each element type)
   SurrogateData< 2 > surrogate_2d_;
   SurrogateData< 3 > surrogate_3d_;
};

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >
    P1ElementwiseSurrogateLaplaceOperator;

typedef P1ElementwiseSurrogateOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > >
    P1ElementwiseSurrogatePolarLaplaceOperator;

typedef P1ElementwiseSurrogateOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >
    P1ElementwiseSurrogateMassOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_mass_blending_q4 > P1ElementwiseSurrogateBlendingMassOperator;

typedef P1ElementwiseSurrogateOperator< P1LinearCombinationForm > P1ElementwiseSurrogateLinearCombinationOperator;

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >
    P1ElementwiseSurrogatePSPGOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_diffusion_blending_q3 > P1ElementwiseSurrogateBlendingLaplaceOperator;
typedef P1ElementwiseSurrogateOperator< forms::p1_diffusion_blending_q2 > P1ElementwiseSurrogateBlendingLaplaceOperatorQ2;

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >
    P1ElementwiseSurrogateDivXOperator;
typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >
    P1ElementwiseSurrogateDivYOperator;
typedef P1ElementwiseSurrogateOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > >
    P1ElementwiseSurrogateDivZOperator;

typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >
    P1ElementwiseSurrogateDivTXOperator;
typedef P1ElementwiseSurrogateOperator<
    P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >
    P1ElementwiseSurrogateDivTYOperator;
typedef P1ElementwiseSurrogateOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > >
    P1ElementwiseSurrogateDivTZOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_affine_q3 >   P1ElementwiseSurrogateAffineDivKGradOperator;
typedef P1ElementwiseSurrogateOperator< forms::p1_div_k_grad_blending_q3 > P1ElementwiseSurrogateBlendingDivKGradOperator;

typedef P1ElementwiseSurrogateOperator< forms::p1_k_mass_affine_q4 > P1ElementwiseSurrogateKMassOperator;

} // namespace hyteg
