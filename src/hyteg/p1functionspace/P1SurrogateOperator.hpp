/*
 * Copyright (c) 2021 Benjamin Mann
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

#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/p1functionspace/P1Operator.hpp"
#include "hyteg/polynomial/LSQPInterpolator.hpp"
#include "hyteg/polynomial/PolynomialEvaluator.hpp"

namespace hyteg {

template < class P1Form, bool USE_INCREMENTAL_EVAL = true >
class P1SurrogateOperator : public P1Operator< P1Form >
{
   using Poly2D = Polynomial2D< MonomialBasis2D >;
   using Poly3D = Polynomial3D< MonomialBasis3D >;

   // todo add polynomials for macro-faces in 3D

   using StencilPoly_cell = std::map< indexing::IndexIncrement, Poly3D >;
   // using StencilPoly_face3D = std::map< uint_t, StencilPoly_cell > ;
   using StencilPoly_face = std::vector< Poly2D >;

   using Interpolator3D    = LSQPInterpolator3D< MonomialBasis3D, LSQPType::VERTEX >;
   using Interpolator2D    = LSQPInterpolator< MonomialBasis2D, LSQPType::VERTEX >;
   using Interpolator_cell = std::map< indexing::IndexIncrement, Interpolator3D >;
   // using Interpolator_face3D = std::map< uint_t, Interpolator_cell > ;
   using Interpolator_face = std::vector< Interpolator2D >;

   using Evaluator_cell = std::map< indexing::IndexIncrement, Polynomial3DEvaluator >;
   // using Evaluator_face3D = std::map< uint_t, Evaluator_cell > ;
   using Evaluator_face = std::vector< Polynomial2DEvaluator >;

   using P1Operator< P1Form >::P1Operator;
   using P1Operator< P1Form >::storage_;
   using P1Operator< P1Form >::h_;
   using P1Operator< P1Form >::stencilSize_;
   using P1Operator< P1Form >::minLevel_;
   using P1Operator< P1Form >::maxLevel_;
   using P1Operator< P1Form >::vertexStencilID_;
   using P1Operator< P1Form >::edgeStencilID_;
   using P1Operator< P1Form >::faceStencilID_;
   using P1Operator< P1Form >::edgeStencil3DID_;
   using P1Operator< P1Form >::faceStencil3DID_;
   using P1Operator< P1Form >::cellStencilID_;
   using P1Operator< P1Form >::assemble_variableStencil_edge_init;
   using P1Operator< P1Form >::assemble_variableStencil_face_init;
   using P1Operator< P1Form >::assemble_variableStencil_cell_init;
   using P1Operator< P1Form >::assemble_variableStencil_edge;
   using P1Operator< P1Form >::assemble_variableStencil_edge3D;
   using P1Operator< P1Form >::assemble_variableStencil_face;
   using P1Operator< P1Form >::assemble_variableStencil_face3D;
   using P1Operator< P1Form >::assemble_variableStencil_cell;

 public:
   P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : P1SurrogateOperator( storage, minLevel, maxLevel, P1Form() )
   {}

   P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form )
   : P1Operator< P1Form >( storage, minLevel, maxLevel, form )
   {
      // todo add polynomials for macro-bounaries
      auto cellDataHandling =
          std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< StencilPoly_cell >, Cell > >( minLevel_, maxLevel_ );

      auto faceDataHandling =
          std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< StencilPoly_face >, Face > >( minLevel_, maxLevel_ );

      storage->addCellData( cellPolyID_, cellDataHandling, "P1OperatorCellPolynomial" );
      storage->addFaceData( facePolyID_, faceDataHandling, "P1OperatorFacePolynomial" );
   }

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
   std::vector< double > computeSurrogateError( uint_t level ) const
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

   std::vector< double > computeSurrogateError2D( uint_t level ) const
   {
      uint_t stencilSize = faceStencilSize2D;
      uint_t rowsizeY    = levelinfo::num_microvertices_per_edge( level );

      std::vector< double > err;

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         uint_t rowsize       = rowsizeY;
         uint_t inner_rowsize = rowsize;

         real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );

         assemble_variableStencil_face_init( face, level );
         assemble_stencil_face_init( face, level );

         std::vector< double > variableStencil( stencilSize_ );

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
                  double err_ij = opr_data[s] - variableStencil[s];
                  normF2 += err_ij * err_ij;
               }
            }
            --inner_rowsize;
         }
         err.push_back( std::sqrt( h_ * h_ * normF2 ) );
      }

      return err;
   }

   std::vector< double > computeSurrogateError3D( uint_t level ) const
   {
      typedef stencilDirection sd;

      const uint_t rowsizeZ = levelinfo::num_microvertices_per_edge( level );

      std::vector< double >              err;
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
                        double sur     = operatorData[vertexdof::logicalIndexOffsetFromVertex( neighbor )];
                        double var     = variableStencil[vertexdof::logicalIndexOffsetFromVertex( neighbor )];
                        double err_ijk = sur - var;

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

   // uint_t degree_;
   // uint_t interpolationLevel_;
   PrimitiveDataID< LevelWiseMemory< StencilPoly_face >, Face > facePolyID_;
   PrimitiveDataID< LevelWiseMemory< StencilPoly_cell >, Cell > cellPolyID_;

   mutable Evaluator_cell cellPolyEvaluator_;
   mutable Evaluator_face facePolyEvaluator_;
};

// todo test other forms

typedef P1SurrogateOperator< forms::p1_diffusion_blending_q1 >  P1SurrogateLaplaceOperator;
typedef P1SurrogateOperator< forms::p1_mass_blending_q4 >       P1SurrogateMassOperator;
typedef P1SurrogateOperator< forms::p1_div_k_grad_blending_q3 > P1SurrogateDivkGradOperator;
typedef P1SurrogateOperator< forms::p1_div_k_grad_affine_q3 >   P1SurrogateAffineDivkGradOperator;

} // namespace hyteg
