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
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q2.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_k_mass_affine_q4.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp>
#include <hyteg/operators/Operator.hpp>
#include <hyteg/p1functionspace/globalIndices.hpp>
#include <hyteg/polynomial/elementwise/data.hpp>
#include <hyteg/polynomial/elementwise/leastSquares.hpp>
#include <hyteg/polynomial/elementwise/polynomial.hpp>
#include <hyteg/solvers/Smoothables.hpp>

namespace hyteg {

template < class P1Form, uint8_t DEGREE, typename ValueType = real_t >
class P1SurrogateOperator : public Operator< P1Function< ValueType >, P1Function< ValueType > >,
                            public GSSmoothable< P1Function< ValueType > >,
                            public GSBackwardsSmoothable< P1Function< ValueType > >,
                            public SORSmoothable< P1Function< ValueType > >,
                            public SORBackwardsSmoothable< P1Function< ValueType > >,
                            public WeightedJacobiSmoothable< P1Function< ValueType > >,
                            public OperatorWithInverseDiagonal< P1Function< ValueType > >
{
   /* On lower levels, storing and evaluating polynomials is significantly less performant.
      Therefore, on levels 0-3 we precompute and store the system matrices, while we use
      surrogates for levels 4+
    */
   static constexpr uint_t min_lvl_for_surrogate = 4;

   /* Single precision LSQ may lead to very poor accuracy of the resulting polynomials. Therefore,
      we use real_t for the polynomial evaluation only, while sticking to double precision LSQ.
    */
   template < uint8_t DIM >
   using Poly       = surrogate::polynomial::Polynomial< real_t, DIM, DEGREE >;
   using PolyDomain = surrogate::polynomial::Domain< real_t >;
   using LSQ        = surrogate::LeastSquares< double >;

   // constant size stencils, i.e., cell-dof (3d), face-dof (2d/3d), edge-dof (2d)
   template < uint8_t DIM >
   using Stencil = p1::stencil::StencilData< DIM >;
   // surrogate stencils
   template < uint8_t DIM_domain, uint8_t DIM_primitive >
   using PolyStencil = p1::stencil::StencilData< DIM_domain, Poly< DIM_primitive > >;
   // variable size stencils, i.e., edge-dof (3d), vertex-dof (2d/3d)
   using VarStencil = std::vector< real_t >;
   // stencil data for LSQ fit
   template < uint8_t DIM >
   using LSQData = p1::stencil::StencilData< DIM, typename surrogate::LeastSquares< real_t >::Vector >;
   // containers for precomputed stencils. usage: map[id][lvl][k] -> stencil at k-th DoF of of el_id on lvl
   template < uint8_t DIM >
   using StencilMap = surrogate::ElementWiseData< std::vector< Stencil< DIM > > >;
   template < uint8_t DIM >
   using VarStencilMap = surrogate::ElementWiseData< std::vector< VarStencil > >;
   // container for surrogate stencils. usage: map[id][lvl] -> polystencil approximating the stencil of el_id on lvl
   template < uint8_t DIM_domain, uint8_t DIM_primitive >
   using PolyStencilMap = surrogate::ElementWiseData< PolyStencil< DIM_domain, DIM_primitive > >;

 public:
   P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : P1SurrogateOperator( storage, minLevel, maxLevel, P1Form() )

         P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                              size_t                                     minLevel,
                              size_t                                     maxLevel,
                              const P1Form&                              form )
   : Operator( storage, minLevel, maxLevel )
   , form_( form )
   , is_initialized_( false )
   , lsq_( maxLevel + 1 )
   , downsampling_( maxLevel + 1 )
   , stencil_vtx_( storage, maxLevel, 0 )
   , stencil_edge_3d_( storage, maxLevel, 1 )
   , stencil_edge_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 1 )
   , stencil_face_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 2 )
   , stencil_face_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 2 )
   , stencil_cell_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 3 )
   , surrogate_edge_2d_( storage, maxLevel, 1 )
   , surrogate_face_2d_( storage, maxLevel, 2 )
   , surrogate_face_3d_( storage, maxLevel, 2 )
   , surrogate_cell_3d_( storage, maxLevel, 3 )
   , varStencil_2d_( storage, maxLevel )
   , varStencil_3d_( storage, maxLevel )
   , surrogate_2d_( storage, maxLevel )
   , surrogate_3d_( storage, maxLevel )
   {}

   // todo: continue copying stuff from elementwiseSurrogate

   void interpolateStencils( uint_t polyDegree, uint_t interpolationLevel )
   {
      // compute polynomial coefficients
      // todo perform QR only once
      if ( storage_->hasGlobalCells() )
      {
         interpolate3D( polyDegree, interpolationLevel );
      }
      else
      {
         interpolate2D( polyDegree, interpolationLevel );
      }
   }

   /* compute h^(d/2)*||A - Aq||_F with variable operator A and surrogate Aq
      @returns [h^(d/2)*||A - Aq||_F restricted to K] for all macro elements K
   */
   std::vector< real_t > computeSurrogateError( uint_t level ) const
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

   const PrimitiveDataID< LevelWiseMemory< StencilPoly_face >, Face > getFacePolyID()
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_LOG_WARNING( "P1SurrogateOperator::getFacePolyID() called for 3D mesh!" );
      }
      return facePolyID_;
   }

   const PrimitiveDataID< LevelWiseMemory< StencilPoly_cell >, Cell > getCellPolyID()
   {
      if ( !storage_->hasGlobalCells() )
      {
         WALBERLA_LOG_WARNING( "P1SurrogateOperator::getCellPolyID() called for 2D mesh!" );
      }
      return cellPolyID_;
   }

 protected:
   static const uint_t faceStencilSize2D = 9;

   std::vector< real_t > computeSurrogateError2D( uint_t level ) const
   {
      uint_t stencilSize = faceStencilSize2D;
      uint_t rowsizeY    = levelinfo::num_microvertices_per_edge( level );

      std::vector< real_t > err;

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         uint_t rowsize       = rowsizeY;
         uint_t inner_rowsize = rowsize;

         real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );

         assemble_variableStencil_face_init( face, level );
         assemble_stencil_face_init( face, level );

         std::vector< real_t > variableStencil( stencilSize_ );

         real_t normF2 = real_c( 0 );

         for ( uint_t j = 1; j < rowsize - 2; ++j )
         {
            assemble_stencil_face_init_y( j );

            for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
            {
               assemble_variableStencil_face( variableStencil.data(), i, j );
               assemble_stencil_face( opr_data, i, j );

               for ( uint_t s = 0; s < stencilSize; ++s )
               {
                  real_t err_ij = opr_data[s] - variableStencil[s];
                  normF2 += err_ij * err_ij;
               }
            }
            --inner_rowsize;
         }
         err.push_back( std::sqrt( h_ * h_ * normF2 ) );
      }

      return err;
   }

   std::vector< real_t > computeSurrogateError3D( uint_t level ) const
   {
      typedef stencilDirection sd;

      const uint_t rowsizeZ = levelinfo::num_microvertices_per_edge( level );

      std::vector< real_t >              err;
      vertexdof::macrocell::StencilMap_T variableStencil;

      for ( auto& it : storage_->getCells() )
      {
         Cell& cell = *it.second;

         auto& operatorData = cell.getData( cellStencilID_ )->getData( level );

         assemble_variableStencil_cell_init( cell, level );
         assemble_stencil_cell_init( cell, level );

         real_t normF2 = 0;

         uint_t rowsizeY, rowsizeX;

         // skip level 0 (no interior points)
         if ( level > 0 )
         {
            for ( uint_t k = 1; k < rowsizeZ - 3; ++k )
            {
               assemble_stencil_cell_init_z( k );

               rowsizeY = rowsizeZ - k;

               for ( uint_t j = 1; j < rowsizeY - 2; ++j )
               {
                  assemble_stencil_cell_init_y( j );

                  rowsizeX = rowsizeY - j;

                  for ( uint_t i = 1; i < rowsizeX - 1; ++i )
                  {
                     assemble_variableStencil_cell( variableStencil, i, j, k );
                     assemble_stencil_cell( operatorData, i, j, k );

                     for ( const auto& neighbor : vertexdof::macrocell::neighborsWithCenter )
                     {
                        real_t sur     = operatorData[vertexdof::logicalIndexOffsetFromVertex( neighbor )];
                        real_t var     = variableStencil[vertexdof::logicalIndexOffsetFromVertex( neighbor )];
                        real_t err_ijk = sur - var;

                        normF2 += err_ijk * err_ijk;
                     }
                  }
               }
            }
         }
         err.push_back( std::sqrt( h_ * h_ * h_ * normF2 ) );
      }

      return err;
   }

   /* interpolate polynomials
   */
   void interpolate2D( uint_t polyDegree, uint_t maxInterpolationLevel )
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
         const uint_t interpolationLevel = std::min( level, maxInterpolationLevel );
         const uint_t lvlDiff            = level - interpolationLevel;
         const uint_t rowsizeY           = levelinfo::num_microvertices_per_edge( interpolationLevel );
         const real_t h_il               = 1.0 / ( walberla::real_c( rowsizeY - 1 ) );

         for ( auto& it : storage_->getFaces() )
         {
            Face& face          = *it.second;
            auto  stencilMemory = face.getData( faceStencilID_ )->getPointer( level );
            auto& stencilPoly   = face.getData( facePolyID_ )->getData( level );

            Interpolator_face interpolator;

            assemble_variableStencil_face_init( face, level );

            // initialize polynomials
            // auto stencilSize   = face.getData(faceStencilID_)->getSize(level); // always returns 27!
            auto stencilSize = faceStencilSize2D;

            for ( uint_t c = 0; c < stencilSize; ++c )
            {
               stencilPoly.push_back( Poly2D( polyDegree ) );
               interpolator.push_back( Interpolator2D( polyDegree, interpolationLevel, level < maxLevel_ ) );
            }

            // add sample points
            uint_t  rowsizeX;
            Point2D x;

            for ( uint_t j = 1; j < rowsizeY - 2; ++j )
            {
               x[1]     = j * h_il;
               rowsizeX = rowsizeY - j;

               for ( uint_t i = 1; i < rowsizeX - 1; ++i )
               {
                  x[0] = i * h_il;

                  assemble_variableStencil_face( stencilMemory, i << lvlDiff, j << lvlDiff );

                  for ( uint_t c = 0; c < stencilSize; ++c )
                  {
                     interpolator[c].addInterpolationPoint( x, stencilMemory[c] );
                  }
               }
            }

            // find polynomials by L2 fit
            for ( uint_t c = 0; c < stencilSize; ++c )
            {
               interpolator[c].interpolate( stencilPoly[c] );
            }
         }
      }

      // initialize polynomial evaluator
      // for ( auto& it : storage_->getFaces() )
      // {
      // auto stencilSize   = it.second->getData(faceStencilID_)->getSize(maxLevel_); // always returns 27!
      auto stencilSize = faceStencilSize2D;

      for ( uint_t c = 0; c < stencilSize; ++c )
      {
         facePolyEvaluator_.push_back( Polynomial2DEvaluator( polyDegree ) );
      }

      //    break; // we use the same evaluator for all faces
      // }
   }

   /* interpolate polynomials
   */
   void interpolate3D( uint_t polyDegree, uint_t maxInterpolationLevel )
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
         // skip level 0 (no interior points)
         if ( level == 0 )
            continue;

         const uint_t interpolationLevel = std::min( level, maxInterpolationLevel );
         const uint_t lvlDiff            = level - interpolationLevel;
         const uint_t rowsizeZ           = levelinfo::num_microvertices_per_edge( interpolationLevel );
         const real_t h_il               = 1.0 / ( walberla::real_c( rowsizeZ - 1 ) );

         for ( const auto& it : storage_->getCells() )
         {
            Cell& cell          = *it.second;
            auto& stencilMemory = cell.getData( cellStencilID_ )->getData( level );
            auto& stencilPoly   = cell.getData( cellPolyID_ )->getData( level );

            Interpolator_cell interpolator;

            assemble_variableStencil_cell_init( cell, level );

            // initialize polynomials
            assemble_variableStencil_cell( stencilMemory, 1, 1, 1 );

            for ( auto& [idx, val] : stencilMemory )
            {
               stencilPoly.insert_or_assign( idx, Poly3D( polyDegree ) );
               interpolator.insert_or_assign( idx, Interpolator3D( polyDegree, interpolationLevel, level < maxLevel_ ) );
            }

            // add sample points
            uint_t  rowsizeY, rowsizeX;
            Point3D x;

            for ( uint_t k = 1; k < rowsizeZ - 3; ++k )
            {
               x[2]     = k * h_il;
               rowsizeY = rowsizeZ - k;

               for ( uint_t j = 1; j < rowsizeY - 2; ++j )
               {
                  x[1]     = j * h_il;
                  rowsizeX = rowsizeY - j;

                  for ( uint_t i = 1; i < rowsizeX - 1; ++i )
                  {
                     x[0] = i * h_il;

                     assemble_variableStencil_cell( stencilMemory, i << lvlDiff, j << lvlDiff, k << lvlDiff );

                     for ( auto& [idx, val] : stencilMemory )
                     {
                        interpolator[idx].addInterpolationPoint( x, val );
                     }
                  }
               }
            }

            // find polynomials by L2 fit
            for ( auto& [idx, interp] : interpolator )
            {
               interp.interpolate( stencilPoly[idx] );
            }
         }
      }

      // initialize polynomial evaluator
      for ( auto& it : storage_->getCells() )
      {
         auto& stencilMemory = it.second->getData( cellStencilID_ )->getData( maxLevel_ );

         for ( auto& [idx, val] : stencilMemory )
         {
            cellPolyEvaluator_.insert_or_assign( idx, Polynomial3DEvaluator( polyDegree ) );
         }

         break; // we use the same evaluator for all cells
      }
   }

   /// stencil assembly ///////////

   /* Initialize assembly of variable edge stencil.
      Will be called before iterating over edge whenever the stencil is applied.
   */
   inline void assemble_stencil_edge_init( Edge& edge, const uint_t level ) const
   {
      assemble_variableStencil_edge_init( edge, level );
   }

   /* Assembly of edge stencil.
      Will be called before stencil is applied to a particuar edge-DoF.
   */
   inline void assemble_stencil_edge( real_t* edge_stencil, const uint_t i ) const
   {
      assemble_variableStencil_edge( edge_stencil, i );
   }

   /* Initialize assembly of face stencil.
      Will be called before iterating over face whenever the stencil is applied.
   */
   inline void assemble_stencil_face_init( Face& face, const uint_t level ) const
   {
      if ( storage_->hasGlobalCells() )
      {
         assemble_variableStencil_face_init( face, level );
      }
      else
      {
         h_ = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );

         auto& stencilPoly = face.getData( facePolyID_ )->getData( level );

         for ( uint_t c = 0; c < facePolyEvaluator_.size(); ++c )
         {
            facePolyEvaluator_[c].setPolynomial( stencilPoly[c] );
         }
      }
   }

   inline void assemble_stencil_face_init_y( const uint_t j ) const
   {
      if ( !( storage_->hasGlobalCells() ) )
      {
         real_t y = h_ * j;

         for ( auto& evaluator : facePolyEvaluator_ )
         {
            evaluator.setY( y );
            if constexpr ( USE_INCREMENTAL_EVAL )
            {
               evaluator.setStartX( 0, h_ );
            }
         }
      }
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 2d domain.
   */
   inline void assemble_stencil_face( real_t* face_stencil, const uint_t i, const uint_t j ) const
   {
      real_t x = h_ * i;

      for ( uint_t c = 0; c < facePolyEvaluator_.size(); ++c )
      {
         if constexpr ( USE_INCREMENTAL_EVAL )
         {
            face_stencil[c] = facePolyEvaluator_[c].incrementEval();
         }
         else
         {
            face_stencil[c] = facePolyEvaluator_[c].evalX( x );
         }
      }
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 3D domain.
   */
   inline void assemble_stencil_face3D( vertexdof::macroface::StencilMap_T& face_stencil, const uint_t i, const uint_t j ) const
   {
      assemble_variableStencil_face3D( face_stencil, i, j );
   }

   /* Initialize assembly of cell stencil.
      Will be called before iterating over cell whenever the stencil is applied.
   */
   inline void assemble_stencil_cell_init( Cell& cell, const uint_t level ) const
   {
      h_ = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );

      auto& stencilPoly = cell.getData( cellPolyID_ )->getData( level );

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         evaluator.setPolynomial( stencilPoly[idx] );
      }
   }

   inline void assemble_stencil_cell_init_z( const uint_t k ) const
   {
      real_t z = h_ * k;

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         evaluator.setZ( z );
      }
   }

   inline void assemble_stencil_cell_init_y( const uint_t j ) const
   {
      real_t y = h_ * j;

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         evaluator.setY( y );
         if constexpr ( USE_INCREMENTAL_EVAL )
         {
            evaluator.setStartX( 0, h_ );
         }
      }
   }

   /* Assembly of cell stencil.
      Will be called before stencil is applied to a particuar cell-DoF.
   */
   inline void assemble_stencil_cell( vertexdof::macrocell::StencilMap_T& cell_stencil,
                                      const uint_t                        i,
                                      const uint_t                        j,
                                      const uint_t                        k ) const
   {
      real_t x = h_ * i;

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         if constexpr ( USE_INCREMENTAL_EVAL )
         {
            cell_stencil[idx] = evaluator.incrementEval();
         }
         else
         {
            cell_stencil[idx] = evaluator.evalX( x );
         }
      }
   }

   inline bool backwards_sor_available() const { return false; }
   inline bool variableStencil() const { return true; }

   P1Form form_;

   bool is_initialized_;

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;

   // least squares approximator for each level
   std::vector< std::shared_ptr< LSQ > > lsq_;
   uint_t                                downsampling_;

   // precomputed stencils for vertices and edges (3d)
   VarStencilMap stencil_vtx_;
   VarStencilMap stencil_edge_3d_;
   // precomputed stencils for level 1-3
   StencilMap< 2 > stencil_edge_2d_;
   StencilMap< 2 > stencil_face_2d_;
   StencilMap< 3 > stencil_face_3d_;
   StencilMap< 3 > stencil_cell_3d_;
   // surrogates for level 4+
   PolyStencilMap< 2, 1 > surrogate_edge_2d_;
   PolyStencilMap< 2, 2 > surrogate_face_2d_;
   PolyStencilMap< 3, 2 > surrogate_face_3d_;
   PolyStencilMap< 3, 3 > surrogate_cell_3d_;
};

// todo test other forms

// typedef P1SurrogateOperator< forms::p1_diffusion_blending_q1 >  P1SurrogateLaplaceOperator;
// typedef P1SurrogateOperator< forms::p1_mass_blending_q4 >       P1SurrogateMassOperator;
// typedef P1SurrogateOperator< forms::p1_div_k_grad_blending_q3 > P1SurrogateDivkGradOperator;
// typedef P1SurrogateOperator< forms::p1_div_k_grad_affine_q3 >   P1SurrogateAffineDivkGradOperator;

} // namespace hyteg
