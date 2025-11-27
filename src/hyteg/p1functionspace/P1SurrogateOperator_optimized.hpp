/*
 * Copyright (c) 2025 Benjamin Mann
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
#include <hyteg/Stencil.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp>
#include <hyteg/operators/Operator.hpp>
#include <hyteg/p1functionspace/P1Elements.hpp>
#include <hyteg/p1functionspace/P1Function.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>
#include <hyteg/p1functionspace/globalIndices.hpp>
#include <hyteg/polynomial/stencil/leastSquares.hpp>
#include <hyteg/polynomial/stencil/polynomial.hpp>
#include <hyteg/solvers/Smoothables.hpp>

#include "hyteg/primitives/PrimitiveID.hpp"

#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {

template < class P1Form, uint8_t DEGREE >
class P1SurrogateOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                            public GSSmoothable< P1Function< real_t > >,
                            public SORSmoothable< P1Function< real_t > >,
                            public WeightedJacobiSmoothable< P1Function< real_t > >,
                            public OperatorWithInverseDiagonal< P1Function< real_t > >
{
   /* On lower levels, storing and evaluating polynomials is significantly less performant.
      Therefore, on levels 0-3 we precompute and store the system matrices, while we use
      surrogates for levels 4+
    */
   static constexpr uint_t min_lvl_for_surrogate = 4;

   /* Single precision LSQ may lead to very poor accuracy of the resulting polynomials. Therefore,
      we use real_t for the polynomial evaluation only, while sticking to double precision LSQ.
    */
   using PolyDomain = surrogate::polynomial::Domain< real_t >;
   using LSQ        = p1::stencil::surrogate::LeastSquares< double >;

   // regular stencils, i.e., cell-dof (3d), face-dof (2d/3d), edge-dof (2d)
   template < uint8_t DIM >
   using Stencil = p1::stencil::StencilData< DIM >;
   // surrogate stencils
   template < uint8_t DIM_domain, uint8_t DIM_primitive >
   using PolyStencil = p1::stencil::surrogate::Polynomial< DIM_domain, DIM_primitive, DEGREE >;
   // irregular stencils, i.e., edge-dof (3d), vertex-dof (2d/3d)
   using VarStencil = std::vector< real_t >;
   // stencil data for LSQ fit
   template < uint8_t DIM >
   using LSQData = p1::stencil::StencilData< DIM, typename surrogate::LeastSquares< real_t >::Vector >;
   // containers for precomputed stencils. usage: map[id][lvl][k] -> stencil at k-th DoF of of el_id on lvl
   template < uint8_t DIM >
   using StencilMap    = surrogate::ElementWiseData< std::vector< Stencil< DIM > > >;
   using VarStencilMap = surrogate::ElementWiseData< std::vector< VarStencil > >;
   using VtxStencilMap = surrogate::ElementWiseData< VarStencil >;
   // container for surrogate stencils. usage: map[id][lvl] -> poly-stencil approximating the stencil of el_id on lvl
   template < uint8_t DIM_domain, uint8_t DIM_primitive >
   using PolyStencilMap = surrogate::ElementWiseData< PolyStencil< DIM_domain, DIM_primitive > >;
   // data structure for DoF indices
   using DofIdx = p1::stencil::StencilData< 3, walberla::uint_t >;

 public:
   // Ctor requiring form
   P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                        size_t                                     minLevel,
                        size_t                                     maxLevel,
                        const P1Form&                              form,
                        size_t                                     downsampling            = 1,
                        const std::string&                         path_to_svd             = "",
                        bool                                       needsInverseDiagEntries = false )
   : Operator< P1Function< real_t >, P1Function< real_t > >( storage, minLevel, maxLevel )
   , form_( form )
   , is_initialized_( false )
   , lsq_volume_( maxLevel + 1 )
   , lsq_interface_( maxLevel + 1 )
   , downsampling_( maxLevel + 1 )
   , stencil_vtx_( storage, maxLevel, 0 )
   , stencil_edge_3d_( storage, maxLevel, 1 )
   , stencil_edge_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 1 )
   , stencil_face_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 2 )
   , stencil_face_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 2 )
   , stencil_cell_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 3 )
   , surrogate_edge_2d_( storage, maxLevel, ( storage->hasGlobalCells() ) ? 99 : 1 ) // don't initialize if 3D
   , surrogate_face_2d_( storage, maxLevel, ( storage->hasGlobalCells() ) ? 99 : 2 ) // don't initialize if 3D
   , surrogate_face_3d_( storage, maxLevel, ( storage->hasGlobalCells() ) ? 2 : 99 ) // don't initialize if 2D
   , surrogate_cell_3d_( storage, maxLevel, 3 )
   {
      init( downsampling, path_to_svd, needsInverseDiagEntries );
   }

   // Ctor without form
   P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                        size_t                                     minLevel,
                        size_t                                     maxLevel,
                        size_t                                     downsampling            = 1,
                        const std::string&                         path_to_svd             = "",
                        bool                                       needsInverseDiagEntries = false )
   : P1SurrogateOperator( storage, minLevel, maxLevel, P1Form(), downsampling, path_to_svd, needsInverseDiagEntries )
   {}

   /* compute h^(d/2)*||A - Aq||_F with variable operator A and surrogate Aq
      @returns [h^(d/2)*||A - Aq||_F restricted to K] for all macro elements K
   */
   std::map< PrimitiveID, real_t > computeSurrogateError( uint_t level ) const
   {
      if ( storage_->hasGlobalCells() )
      {
         return computeSurrogateError3D( level );
      }
      else
      {
         return computeSurrogateError2D( level );
      }
   }

   void store_svd( const std::string& path_to_svd )
   {
      for ( uint_t level = min_lvl_for_surrogate; level <= maxLevel_; ++level )
      {
         if ( lsq_volume_[level] != nullptr )
         {
            lsq_volume_[level]->write_to_file( path_to_svd );
            lsq_interface_[level]->write_to_file( path_to_svd );
         }
      }
   }

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final;

   void smooth_gs( const P1Function< real_t >& dst,
                   const P1Function< real_t >& rhs,
                   size_t                      level,
                   DoFType                     flag ) const override final

   {
      smooth_sor( dst, rhs, 1.0, level, flag );
   }

   void smooth_sor( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override;

   void smooth_jac( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    const P1Function< real_t >& src,
                    const real_t                relax,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      this->startTiming( "smooth_jac" );

      // compute the current residual
      this->apply( src, dst, level, flag );
      dst.assign( { real_t( 1.0 ), real_t( -1.0 ) }, { rhs, dst }, level, flag );

      // perform Jacobi update step
      dst.multElementwise( { *getInverseDiagonalValues(), dst }, level, flag );
      dst.assign( { static_cast< real_t >( 1.0 ), relax }, { src, dst }, level, flag );

      this->stopTiming( "smooth_jac" );
   }

   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false ); }

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

   /**
    * @brief Returns string representation of the surrogate stencils for the given primitive and level.
    *
    * @param id The primitive ID for which to show polynomials.
    * @param lvl The level at which to retrieve the polynomial surrogates.
    * @return std::map<p1::stencil::Dir, std::string> Map from stencil direction to polynomial string.
    */
   std::map< p1::stencil::Dir, std::string > show_polynomials( const PrimitiveID& id, const uint_t lvl ) const;

 private:
   void init( size_t downsampling, const std::string& path_to_svd, bool needsInverseDiagEntries );

   std::map< PrimitiveID, real_t > computeSurrogateError2D( uint_t level ) const
   {
      std::map< PrimitiveID, real_t > errors;
      // todo
      return errors;
   }

   std::map< PrimitiveID, real_t > computeSurrogateError3D( uint_t level ) const
   {
      std::map< PrimitiveID, real_t > errors;
      // todo
      return errors;
   }

   void precompute_stencil_vtx_2d( uint_t lvl );
   void precompute_stencil_vtx_3d( uint_t lvl );
   void precompute_stencil_edge_2d( uint_t lvl );
   void precompute_stencil_edge_3d( uint_t lvl );
   void precompute_stencil_face_2d( uint_t lvl );
   void precompute_stencil_face_3d( uint_t lvl );
   void precompute_stencil_cell_3d( uint_t lvl );

   void compute_surrogates_edge_2d( uint_t lvl );
   void compute_surrogates_face_2d( uint_t lvl );
   void compute_surrogates_face_3d( uint_t lvl );
   void compute_surrogates_cell_3d( uint_t lvl );

   void compute_offsets_for_face_stencil_3d( uint_t lvl );

   void apply_edge_precomputed_2d( std::shared_ptr< hyteg::Edge > edge,
                                   uint_t                         lvl,
                                   const real_t* RESTRICT const   srcData,
                                   real_t* RESTRICT               dstData,
                                   UpdateType                     updateType ) const;

   void apply_edge_precomputed_3d( std::shared_ptr< hyteg::Edge > edge,
                                   uint_t                         lvl,
                                   const real_t* RESTRICT const   srcData,
                                   real_t* RESTRICT               dstData,
                                   UpdateType                     updateType ) const;

   void apply_face_precomputed_2d( std::shared_ptr< hyteg::Face > face,
                                   uint_t                         lvl,
                                   const real_t* RESTRICT const   srcData,
                                   real_t* RESTRICT               dstData,
                                   UpdateType                     updateType ) const;

   void apply_face_precomputed_3d( std::shared_ptr< hyteg::Face > face,
                                   uint_t                         lvl,
                                   const real_t* RESTRICT const   srcData,
                                   real_t* RESTRICT               dstData,
                                   UpdateType                     updateType ) const;

   void apply_cell_precomputed_3d( std::shared_ptr< hyteg::Cell > cell,
                                   uint_t                         lvl,
                                   const real_t* RESTRICT const   srcData,
                                   real_t* RESTRICT               dstData,
                                   UpdateType                     updateType ) const;

   void apply_edge_surrogate_2d( std::shared_ptr< hyteg::Edge > edge,
                                 uint_t                         lvl,
                                 const real_t* RESTRICT const   srcData,
                                 real_t* RESTRICT               dstData,
                                 UpdateType                     updateType ) const;

   void apply_face_surrogate_2d( std::shared_ptr< hyteg::Face > face,
                                 uint_t                         lvl,
                                 const real_t* RESTRICT const   srcData,
                                 real_t* RESTRICT               dstData,
                                 UpdateType                     updateType ) const;

   void apply_face_surrogate_3d( std::shared_ptr< hyteg::Face > face,
                                 uint_t                         lvl,
                                 const real_t* RESTRICT const   srcData,
                                 real_t* RESTRICT               dstData,
                                 UpdateType                     updateType ) const;

   void apply_cell_surrogate_3d( std::shared_ptr< hyteg::Cell > cell,
                                 uint_t                         lvl,
                                 const real_t* RESTRICT const   srcData,
                                 real_t* RESTRICT               dstData,
                                 UpdateType                     updateType ) const;

   void computeDiagonalOperatorValues( bool invert );

   void assemble_diagonalOperator_edge_precomputed_2d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_edge_precomputed_3d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_face_precomputed_2d( std::shared_ptr< hyteg::Face > face, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_face_precomputed_3d( std::shared_ptr< hyteg::Face > face, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_cell_precomputed_3d( std::shared_ptr< hyteg::Cell > cell, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_edge_surrogate_2d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_face_surrogate_2d( std::shared_ptr< hyteg::Face > face, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_face_surrogate_3d( std::shared_ptr< hyteg::Face > face, uint_t lvl, real_t* diagData );
   void assemble_diagonalOperator_cell_surrogate_3d( std::shared_ptr< hyteg::Cell > cell, uint_t lvl, real_t* diagData );

   P1Form form_;

   bool is_initialized_;

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;

   // least squares approximation for each level
   std::vector< std::shared_ptr< LSQ > > lsq_volume_;    // lsq for volume primitives (cells in 3d / faces in 2d)
   std::vector< std::shared_ptr< LSQ > > lsq_interface_; // lsq for interface primitives (faces in 3d / edges in 2d)
   uint_t                                downsampling_;

   // precomputed irregular stencils (vertices in 2D, vertices and edges in 3D
   VtxStencilMap stencil_vtx_;
   VarStencilMap stencil_edge_3d_;
   // precomputed regular stencils for level 1-3
   StencilMap< 2 > stencil_edge_2d_;
   StencilMap< 2 > stencil_face_2d_;
   StencilMap< 3 > stencil_face_3d_;
   StencilMap< 3 > stencil_cell_3d_;
   // surrogates for level 4+
   PolyStencilMap< 2, 1 > surrogate_edge_2d_;
   PolyStencilMap< 2, 2 > surrogate_face_2d_;
   PolyStencilMap< 3, 2 > surrogate_face_3d_;
   PolyStencilMap< 3, 3 > surrogate_cell_3d_;
   // logical offsets on 3d faces, i.e. stencil directions as {-1,0,1}^3 tuples
   std::map< PrimitiveID, p1::stencil::StencilData< 3, indexing::Index > > offsets_face_3d_;
};

template < uint8_t DEGREE >
using P1SurrogateBlendingMassOperator = P1SurrogateOperator< forms::p1_mass_blending_q4, DEGREE >;
template < uint8_t DEGREE >
using P1SurrogateBlendingLaplaceOperator = P1SurrogateOperator< forms::p1_diffusion_blending_q3, DEGREE >;
template < uint8_t DEGREE >
using P1SurrogateAffineDivKGradOperator = P1SurrogateOperator< forms::p1_div_k_grad_affine_q3, DEGREE >;
template < uint8_t DEGREE >
using P1SurrogateBlendingDivKGradOperator = P1SurrogateOperator< forms::p1_div_k_grad_blending_q3, DEGREE >;

} // namespace hyteg
