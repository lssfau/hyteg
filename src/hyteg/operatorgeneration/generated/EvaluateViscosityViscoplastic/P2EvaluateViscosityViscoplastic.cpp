/*
* Copyright (c) 2017-2024 Nils Kohl, Daniel Bauer, Fabian Böhm.
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

/*
* The entire file was generated with the HyTeG Operator Generator.
*
* Avoid modifying this file. If buggy, consider fixing the generator itself.
*/

// Unfortunately, the inverse diagonal kernel wrapper triggers a GCC bug (maybe
// (related to) https://gcc.gnu.org/bugzilla/show_bug.cgi?id=107087) causing a
// warning in an internal standard library header (bits/stl_algobase.h). As a
// workaround, we disable the warning and include this header indirectly through
// a public header.
#include <waLBerlaDefinitions.h>
#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnonnull"
#endif
#include <cmath>
#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#include "P2EvaluateViscosityViscoplastic.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P2EvaluateViscosityViscoplastic::P2EvaluateViscosityViscoplastic( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                  size_t                                     minLevel,
                                                                  size_t                                     maxLevel,
                                                                  const P1Function< walberla::float64 >&     _mu_lin,
                                                                  const P1Function< walberla::float64 >&     _ux,
                                                                  const P1Function< walberla::float64 >&     _uy,
                                                                  walberla::float64 mu_star_P2EvaluateViscosityViscoplastic,
                                                                  walberla::float64 sigma_y_P2EvaluateViscosityViscoplastic )
: Operator( storage, minLevel, maxLevel )
, mu_lin( _mu_lin )
, ux( _ux )
, uy( _uy )
, mu_star_P2EvaluateViscosityViscoplastic_( mu_star_P2EvaluateViscosityViscoplastic )
, sigma_y_P2EvaluateViscosityViscoplastic_( sigma_y_P2EvaluateViscosityViscoplastic )
{}

void P2EvaluateViscosityViscoplastic::apply( const P1Function< walberla::float64 >& src,
                                             const P1Function< walberla::float64 >& dst,
                                             uint_t                                 level,
                                             DoFType                                flag,
                                             UpdateType                             updateType ) const
{
   this->startTiming( "apply" );

   // Make sure that halos are up-to-date
   this->timingTree_->start( "pre-communication" );
   if ( this->storage_->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented." );
   }
   else
   {
      communication::syncFunctionBetweenPrimitives( src, level, communication::syncDirection_t::LOW2HIGH );
      communication::syncFunctionBetweenPrimitives( mu_lin, level, communication::syncDirection_t::LOW2HIGH );
      communication::syncFunctionBetweenPrimitives( ux, level, communication::syncDirection_t::LOW2HIGH );
      communication::syncFunctionBetweenPrimitives( uy, level, communication::syncDirection_t::LOW2HIGH );
   }
   this->timingTree_->stop( "pre-communication" );

   if ( updateType == Replace )
   {
      // We need to zero the destination array (including halos).
      // However, we must not zero out anything that is not flagged with the specified BCs.
      // Therefore, we first zero out everything that flagged, and then, later,
      // the halos of the highest dim primitives.
      dst.interpolate( walberla::numeric_cast< walberla::float64 >( 0 ), level, flag );
   }

   if ( storage_->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented." );
   }
   else
   {
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         // get hold of the actual numerical data in the functions
         walberla::float64* _data_src    = face.getData( src.getFaceDataID() )->getPointer( level );
         walberla::float64* _data_dst    = face.getData( dst.getFaceDataID() )->getPointer( level );
         walberla::float64* _data_mu_lin = face.getData( mu_lin.getFaceDataID() )->getPointer( level );
         walberla::float64* _data_ux     = face.getData( ux.getFaceDataID() )->getPointer( level );
         walberla::float64* _data_uy     = face.getData( uy.getFaceDataID() )->getPointer( level );

         // Zero out dst halos only
         //
         // This is also necessary when using update type == Add.
         // During additive comm we then skip zeroing the data on the lower-dim primitives.
         for ( const auto& idx : vertexdof::macroface::Iterator( level ) )
         {
            if ( vertexdof::macroface::isVertexOnBoundary( level, idx ) )
            {
               auto arrayIdx       = vertexdof::macroface::index( level, idx.x(), idx.y() );
               _data_dst[arrayIdx] = walberla::float64( 0 );
            }
         }

         const auto micro_edges_per_macro_edge                = (int64_t) levelinfo::num_microedges_per_edge( level );
         const auto micro_edges_per_macro_edge_float          = (walberla::float64) levelinfo::num_microedges_per_edge( level );
         const walberla::float64 macro_vertex_coord_id_0comp0 = (walberla::float64) face.getCoordinates()[0][0];
         const walberla::float64 macro_vertex_coord_id_0comp1 = (walberla::float64) face.getCoordinates()[0][1];
         const walberla::float64 macro_vertex_coord_id_1comp0 = (walberla::float64) face.getCoordinates()[1][0];
         const walberla::float64 macro_vertex_coord_id_1comp1 = (walberla::float64) face.getCoordinates()[1][1];
         const walberla::float64 macro_vertex_coord_id_2comp0 = (walberla::float64) face.getCoordinates()[2][0];
         const walberla::float64 macro_vertex_coord_id_2comp1 = (walberla::float64) face.getCoordinates()[2][1];

         this->timingTree_->start( "kernel" );

         apply_P2EvaluateViscosityViscoplastic_macro_2D(

             _data_dst,
             _data_mu_lin,
             _data_src,
             _data_ux,
             _data_uy,
             macro_vertex_coord_id_0comp0,
             macro_vertex_coord_id_0comp1,
             macro_vertex_coord_id_1comp0,
             macro_vertex_coord_id_1comp1,
             macro_vertex_coord_id_2comp0,
             macro_vertex_coord_id_2comp1,
             micro_edges_per_macro_edge,
             micro_edges_per_macro_edge_float,
             mu_star_P2EvaluateViscosityViscoplastic_,
             sigma_y_P2EvaluateViscosityViscoplastic_ );

         this->timingTree_->stop( "kernel" );
      }

      // Push result to lower-dimensional primitives
      //
      this->timingTree_->start( "post-communication" );
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      this->timingTree_->stop( "post-communication" );
   }

   this->stopTiming( "apply" );
}
void P2EvaluateViscosityViscoplastic::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                const P1Function< idx_t >&                  src,
                                                const P1Function< idx_t >&                  dst,
                                                uint_t                                      level,
                                                DoFType                                     flag ) const
{
   this->startTiming( "toMatrix" );

   // We currently ignore the flag provided!
   if ( flag != All )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in toMatrix; using flag = All" );
   }

   if ( storage_->hasGlobalCells() )
   {
      this->timingTree_->start( "pre-communication" );
      mu_lin.communicate< Face, Cell >( level );
      mu_lin.communicate< Edge, Cell >( level );
      mu_lin.communicate< Vertex, Cell >( level );
      ux.communicate< Face, Cell >( level );
      ux.communicate< Edge, Cell >( level );
      ux.communicate< Vertex, Cell >( level );
      uy.communicate< Face, Cell >( level );
      uy.communicate< Edge, Cell >( level );
      uy.communicate< Vertex, Cell >( level );
      this->timingTree_->stop( "pre-communication" );

      WALBERLA_ABORT( "Not implemented." );
   }
   else
   {
      this->timingTree_->start( "pre-communication" );
      communication::syncFunctionBetweenPrimitives( mu_lin, level, communication::syncDirection_t::LOW2HIGH );
      communication::syncFunctionBetweenPrimitives( ux, level, communication::syncDirection_t::LOW2HIGH );
      communication::syncFunctionBetweenPrimitives( uy, level, communication::syncDirection_t::LOW2HIGH );
      this->timingTree_->stop( "pre-communication" );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         // get hold of the actual numerical data
         idx_t*             _data_src    = face.getData( src.getFaceDataID() )->getPointer( level );
         idx_t*             _data_dst    = face.getData( dst.getFaceDataID() )->getPointer( level );
         walberla::float64* _data_mu_lin = face.getData( mu_lin.getFaceDataID() )->getPointer( level );
         walberla::float64* _data_ux     = face.getData( ux.getFaceDataID() )->getPointer( level );
         walberla::float64* _data_uy     = face.getData( uy.getFaceDataID() )->getPointer( level );

         const auto micro_edges_per_macro_edge                = (int64_t) levelinfo::num_microedges_per_edge( level );
         const auto micro_edges_per_macro_edge_float          = (walberla::float64) levelinfo::num_microedges_per_edge( level );
         const walberla::float64 macro_vertex_coord_id_0comp0 = (walberla::float64) face.getCoordinates()[0][0];
         const walberla::float64 macro_vertex_coord_id_0comp1 = (walberla::float64) face.getCoordinates()[0][1];
         const walberla::float64 macro_vertex_coord_id_1comp0 = (walberla::float64) face.getCoordinates()[1][0];
         const walberla::float64 macro_vertex_coord_id_1comp1 = (walberla::float64) face.getCoordinates()[1][1];
         const walberla::float64 macro_vertex_coord_id_2comp0 = (walberla::float64) face.getCoordinates()[2][0];
         const walberla::float64 macro_vertex_coord_id_2comp1 = (walberla::float64) face.getCoordinates()[2][1];

         this->timingTree_->start( "kernel" );

         toMatrix_P2EvaluateViscosityViscoplastic_macro_2D(

             _data_dst,
             _data_mu_lin,
             _data_src,
             _data_ux,
             _data_uy,
             macro_vertex_coord_id_0comp0,
             macro_vertex_coord_id_0comp1,
             macro_vertex_coord_id_1comp0,
             macro_vertex_coord_id_1comp1,
             macro_vertex_coord_id_2comp0,
             macro_vertex_coord_id_2comp1,
             mat,
             micro_edges_per_macro_edge,
             micro_edges_per_macro_edge_float,
             mu_star_P2EvaluateViscosityViscoplastic_,
             sigma_y_P2EvaluateViscosityViscoplastic_ );

         this->timingTree_->stop( "kernel" );
      }
   }
   this->stopTiming( "toMatrix" );
}
void P2EvaluateViscosityViscoplastic::computeInverseDiagonalOperatorValues()
{
   this->startTiming( "computeInverseDiagonalOperatorValues" );

   if ( invDiag_ == nullptr )
   {
      invDiag_ =
          std::make_shared< P1Function< walberla::float64 > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
   }

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      invDiag_->setToZero( level );

      if ( storage_->hasGlobalCells() )
      {
         this->timingTree_->start( "pre-communication" );
         mu_lin.communicate< Face, Cell >( level );
         mu_lin.communicate< Edge, Cell >( level );
         mu_lin.communicate< Vertex, Cell >( level );
         ux.communicate< Face, Cell >( level );
         ux.communicate< Edge, Cell >( level );
         ux.communicate< Vertex, Cell >( level );
         uy.communicate< Face, Cell >( level );
         uy.communicate< Edge, Cell >( level );
         uy.communicate< Vertex, Cell >( level );
         this->timingTree_->stop( "pre-communication" );

         WALBERLA_ABORT( "Not implemented." );
         ( *invDiag_ ).invertElementwise( level );
      }
      else
      {
         this->timingTree_->start( "pre-communication" );
         communication::syncFunctionBetweenPrimitives( mu_lin, level, communication::syncDirection_t::LOW2HIGH );
         communication::syncFunctionBetweenPrimitives( ux, level, communication::syncDirection_t::LOW2HIGH );
         communication::syncFunctionBetweenPrimitives( uy, level, communication::syncDirection_t::LOW2HIGH );
         this->timingTree_->stop( "pre-communication" );

         for ( auto& it : storage_->getFaces() )
         {
            Face& face = *it.second;

            // get hold of the actual numerical data
            walberla::float64* _data_invDiag_ = face.getData( ( *invDiag_ ).getFaceDataID() )->getPointer( level );
            walberla::float64* _data_mu_lin   = face.getData( mu_lin.getFaceDataID() )->getPointer( level );
            walberla::float64* _data_ux       = face.getData( ux.getFaceDataID() )->getPointer( level );
            walberla::float64* _data_uy       = face.getData( uy.getFaceDataID() )->getPointer( level );

            const auto micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
            const auto micro_edges_per_macro_edge_float = (walberla::float64) levelinfo::num_microedges_per_edge( level );
            const walberla::float64 macro_vertex_coord_id_0comp0 = (walberla::float64) face.getCoordinates()[0][0];
            const walberla::float64 macro_vertex_coord_id_0comp1 = (walberla::float64) face.getCoordinates()[0][1];
            const walberla::float64 macro_vertex_coord_id_1comp0 = (walberla::float64) face.getCoordinates()[1][0];
            const walberla::float64 macro_vertex_coord_id_1comp1 = (walberla::float64) face.getCoordinates()[1][1];
            const walberla::float64 macro_vertex_coord_id_2comp0 = (walberla::float64) face.getCoordinates()[2][0];
            const walberla::float64 macro_vertex_coord_id_2comp1 = (walberla::float64) face.getCoordinates()[2][1];

            this->timingTree_->start( "kernel" );

            computeInverseDiagonalOperatorValues_P2EvaluateViscosityViscoplastic_macro_2D(

                _data_invDiag_,
                _data_mu_lin,
                _data_ux,
                _data_uy,
                macro_vertex_coord_id_0comp0,
                macro_vertex_coord_id_0comp1,
                macro_vertex_coord_id_1comp0,
                macro_vertex_coord_id_1comp1,
                macro_vertex_coord_id_2comp0,
                macro_vertex_coord_id_2comp1,
                micro_edges_per_macro_edge,
                micro_edges_per_macro_edge_float,
                mu_star_P2EvaluateViscosityViscoplastic_,
                sigma_y_P2EvaluateViscosityViscoplastic_ );

            this->timingTree_->stop( "kernel" );
         }

         // Push result to lower-dimensional primitives
         //
         this->timingTree_->start( "post-communication" );
         // Note: We could avoid communication here by implementing the apply() also for the respective
         //       lower dimensional primitives!
         ( *invDiag_ ).communicateAdditively< Face, Edge >( level );
         ( *invDiag_ ).communicateAdditively< Face, Vertex >( level );
         this->timingTree_->stop( "post-communication" );
         ( *invDiag_ ).invertElementwise( level );
      }
   }

   this->stopTiming( "computeInverseDiagonalOperatorValues" );
}
std::shared_ptr< P1Function< walberla::float64 > > P2EvaluateViscosityViscoplastic::getInverseDiagonalValues() const
{
   return invDiag_;
}
void P2EvaluateViscosityViscoplastic::apply_P2EvaluateViscosityViscoplastic_macro_2D(
    walberla::float64* RESTRICT _data_dst,
    walberla::float64* RESTRICT _data_mu_lin,
    walberla::float64* RESTRICT _data_src,
    walberla::float64* RESTRICT _data_ux,
    walberla::float64* RESTRICT _data_uy,
    walberla::float64           macro_vertex_coord_id_0comp0,
    walberla::float64           macro_vertex_coord_id_0comp1,
    walberla::float64           macro_vertex_coord_id_1comp0,
    walberla::float64           macro_vertex_coord_id_1comp1,
    walberla::float64           macro_vertex_coord_id_2comp0,
    walberla::float64           macro_vertex_coord_id_2comp1,
    int64_t                     micro_edges_per_macro_edge,
    walberla::float64           micro_edges_per_macro_edge_float,
    walberla::float64           mu_star,
    walberla::float64           sigma_y ) const
{
   {
      {
         /* FaceType.GRAY */
         const walberla::float64 tmp_coords_jac_0_GRAY   = 1.0 / ( micro_edges_per_macro_edge_float ) * 1.0;
         const walberla::float64 p_affine_const_0_0_GRAY = macro_vertex_coord_id_0comp0;
         const walberla::float64 p_affine_const_0_1_GRAY = macro_vertex_coord_id_0comp1;
         const walberla::float64 p_affine_const_1_0_GRAY =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 );
         const walberla::float64 p_affine_const_1_1_GRAY =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 );
         const walberla::float64 p_affine_const_2_0_GRAY =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 );
         const walberla::float64 p_affine_const_2_1_GRAY =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 );
         const walberla::float64 jac_affine_0_0_GRAY = -p_affine_const_0_0_GRAY + p_affine_const_1_0_GRAY;
         const walberla::float64 jac_affine_0_1_GRAY = -p_affine_const_0_0_GRAY + p_affine_const_2_0_GRAY;
         const walberla::float64 jac_affine_1_0_GRAY = -p_affine_const_0_1_GRAY + p_affine_const_1_1_GRAY;
         const walberla::float64 jac_affine_1_1_GRAY = -p_affine_const_0_1_GRAY + p_affine_const_2_1_GRAY;
         const walberla::float64 tmp_coords_jac_1_GRAY =
             1.0 / ( jac_affine_0_0_GRAY * jac_affine_1_1_GRAY - jac_affine_0_1_GRAY * jac_affine_1_0_GRAY );
         const walberla::float64 jac_affine_inv_0_0_GRAY = jac_affine_1_1_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_0_1_GRAY = -jac_affine_0_1_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_1_0_GRAY = -jac_affine_1_0_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_1_1_GRAY = jac_affine_0_0_GRAY * tmp_coords_jac_1_GRAY;
         for ( int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1 )
            for ( int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge; ctr_0 += 1 )
            {
               const walberla::float64 p_affine_0_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_0_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 src_dof_0 =
                   _data_src[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_1 =
                   _data_src[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_2 = _data_src[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                             ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 mu_lin_dof_0 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 mu_lin_dof_1 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 mu_lin_dof_2 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_0 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_1 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_2 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_0 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_1 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_2        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 tmp_kernel_op_0 = -ux_dof_0 + ux_dof_1;
               const walberla::float64 tmp_kernel_op_1 = -ux_dof_0 + ux_dof_2;
               const walberla::float64 tmp_kernel_op_2 = -uy_dof_0 + uy_dof_1;
               const walberla::float64 tmp_kernel_op_3 = -uy_dof_0 + uy_dof_2;
               const walberla::float64 tmp_kernel_op_4 = jac_affine_inv_0_0_GRAY * tmp_kernel_op_2;
               const walberla::float64 tmp_kernel_op_5 = jac_affine_inv_0_1_GRAY * tmp_kernel_op_0;
               const walberla::float64 tmp_kernel_op_6 = jac_affine_inv_1_0_GRAY * tmp_kernel_op_3;
               const walberla::float64 tmp_kernel_op_7 = jac_affine_inv_1_1_GRAY * tmp_kernel_op_1;
               const walberla::float64 tmp_kernel_op_8 =
                   1.0 /
                   ( mu_star +
                     sigma_y *
                         pow( ( ( jac_affine_inv_0_0_GRAY * tmp_kernel_op_0 + jac_affine_inv_1_0_GRAY * tmp_kernel_op_1 ) *
                                ( jac_affine_inv_0_0_GRAY * tmp_kernel_op_0 + jac_affine_inv_1_0_GRAY * tmp_kernel_op_1 ) ) +
                                  ( ( jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 + jac_affine_inv_1_1_GRAY * tmp_kernel_op_3 ) *
                                    ( jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 + jac_affine_inv_1_1_GRAY * tmp_kernel_op_3 ) ) +
                                  ( tmp_kernel_op_4 + tmp_kernel_op_5 + tmp_kernel_op_6 + tmp_kernel_op_7 ) *
                                      ( tmp_kernel_op_4 * 0.5 + tmp_kernel_op_5 * 0.5 + tmp_kernel_op_6 * 0.5 +
                                        tmp_kernel_op_7 * 0.5 ),
                              -0.50000000000000000 ) ) *
                   1.0;
               const walberla::float64 elMatVec_0 = src_dof_0 * 1.0 / ( 1.0 / ( mu_lin_dof_0 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatVec_1 = src_dof_1 * 1.0 / ( 1.0 / ( mu_lin_dof_1 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatVec_2 = src_dof_2 * 1.0 / ( 1.0 / ( mu_lin_dof_2 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatVec_0 +
                   _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatVec_1 +
                   _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_2 + _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
            }
      }
      {
         /* FaceType.BLUE */
         const walberla::float64 tmp_coords_jac_0_BLUE = 1.0 / ( micro_edges_per_macro_edge_float ) * 1.0;
         const walberla::float64 tmp_coords_jac_1_BLUE =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 );
         const walberla::float64 tmp_coords_jac_2_BLUE =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 );
         const walberla::float64 tmp_coords_jac_3_BLUE =
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 );
         const walberla::float64 tmp_coords_jac_4_BLUE =
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 );
         const walberla::float64 p_affine_const_0_0_BLUE = tmp_coords_jac_1_BLUE;
         const walberla::float64 p_affine_const_0_1_BLUE = tmp_coords_jac_2_BLUE;
         const walberla::float64 p_affine_const_1_0_BLUE = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3_BLUE;
         const walberla::float64 p_affine_const_1_1_BLUE = macro_vertex_coord_id_0comp1 + tmp_coords_jac_4_BLUE;
         const walberla::float64 p_affine_const_2_0_BLUE = tmp_coords_jac_1_BLUE + tmp_coords_jac_3_BLUE;
         const walberla::float64 p_affine_const_2_1_BLUE = tmp_coords_jac_2_BLUE + tmp_coords_jac_4_BLUE;
         const walberla::float64 jac_affine_0_0_BLUE     = -p_affine_const_0_0_BLUE + p_affine_const_1_0_BLUE;
         const walberla::float64 jac_affine_0_1_BLUE     = -p_affine_const_0_0_BLUE + p_affine_const_2_0_BLUE;
         const walberla::float64 jac_affine_1_0_BLUE     = -p_affine_const_0_1_BLUE + p_affine_const_1_1_BLUE;
         const walberla::float64 jac_affine_1_1_BLUE     = -p_affine_const_0_1_BLUE + p_affine_const_2_1_BLUE;
         const walberla::float64 tmp_coords_jac_5_BLUE =
             1.0 / ( jac_affine_0_0_BLUE * jac_affine_1_1_BLUE - jac_affine_0_1_BLUE * jac_affine_1_0_BLUE );
         const walberla::float64 jac_affine_inv_0_0_BLUE = jac_affine_1_1_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_0_1_BLUE = -jac_affine_0_1_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_1_0_BLUE = -jac_affine_1_0_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_1_1_BLUE = jac_affine_0_0_BLUE * tmp_coords_jac_5_BLUE;
         for ( int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1 )
            for ( int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1; ctr_0 += 1 )
            {
               const walberla::float64 p_affine_0_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_0_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_1_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 src_dof_0 =
                   _data_src[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_1 = _data_src[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                             ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_2 = _data_src[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                             ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 mu_lin_dof_0 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 mu_lin_dof_1 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 mu_lin_dof_2 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_0 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_1 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_2 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_0 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_1        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_2        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 tmp_kernel_op_0 = -ux_dof_0 + ux_dof_1;
               const walberla::float64 tmp_kernel_op_1 = -ux_dof_0 + ux_dof_2;
               const walberla::float64 tmp_kernel_op_2 = -uy_dof_0 + uy_dof_1;
               const walberla::float64 tmp_kernel_op_3 = -uy_dof_0 + uy_dof_2;
               const walberla::float64 tmp_kernel_op_4 = jac_affine_inv_0_0_BLUE * tmp_kernel_op_2;
               const walberla::float64 tmp_kernel_op_5 = jac_affine_inv_0_1_BLUE * tmp_kernel_op_0;
               const walberla::float64 tmp_kernel_op_6 = jac_affine_inv_1_0_BLUE * tmp_kernel_op_3;
               const walberla::float64 tmp_kernel_op_7 = jac_affine_inv_1_1_BLUE * tmp_kernel_op_1;
               const walberla::float64 tmp_kernel_op_8 =
                   1.0 /
                   ( mu_star +
                     sigma_y *
                         pow( ( ( jac_affine_inv_0_0_BLUE * tmp_kernel_op_0 + jac_affine_inv_1_0_BLUE * tmp_kernel_op_1 ) *
                                ( jac_affine_inv_0_0_BLUE * tmp_kernel_op_0 + jac_affine_inv_1_0_BLUE * tmp_kernel_op_1 ) ) +
                                  ( ( jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 + jac_affine_inv_1_1_BLUE * tmp_kernel_op_3 ) *
                                    ( jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 + jac_affine_inv_1_1_BLUE * tmp_kernel_op_3 ) ) +
                                  ( tmp_kernel_op_4 + tmp_kernel_op_5 + tmp_kernel_op_6 + tmp_kernel_op_7 ) *
                                      ( tmp_kernel_op_4 * 0.5 + tmp_kernel_op_5 * 0.5 + tmp_kernel_op_6 * 0.5 +
                                        tmp_kernel_op_7 * 0.5 ),
                              -0.50000000000000000 ) ) *
                   1.0;
               const walberla::float64 elMatVec_0 = src_dof_0 * 1.0 / ( 1.0 / ( mu_lin_dof_0 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatVec_1 = src_dof_1 * 1.0 / ( 1.0 / ( mu_lin_dof_1 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatVec_2 = src_dof_2 * 1.0 / ( 1.0 / ( mu_lin_dof_2 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatVec_0 +
                   _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_1 + _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] =
                   elMatVec_2 + _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
            }
      }
   }
}
void P2EvaluateViscosityViscoplastic::toMatrix_P2EvaluateViscosityViscoplastic_macro_2D(
    idx_t* RESTRICT                      _data_dst,
    walberla::float64* RESTRICT          _data_mu_lin,
    idx_t* RESTRICT                      _data_src,
    walberla::float64* RESTRICT          _data_ux,
    walberla::float64* RESTRICT          _data_uy,
    walberla::float64                    macro_vertex_coord_id_0comp0,
    walberla::float64                    macro_vertex_coord_id_0comp1,
    walberla::float64                    macro_vertex_coord_id_1comp0,
    walberla::float64                    macro_vertex_coord_id_1comp1,
    walberla::float64                    macro_vertex_coord_id_2comp0,
    walberla::float64                    macro_vertex_coord_id_2comp1,
    std::shared_ptr< SparseMatrixProxy > mat,
    int64_t                              micro_edges_per_macro_edge,
    walberla::float64                    micro_edges_per_macro_edge_float,
    walberla::float64                    mu_star,
    walberla::float64                    sigma_y ) const
{
   {
      {
         /* FaceType.GRAY */
         const walberla::float64 tmp_coords_jac_0_GRAY   = 1.0 / ( micro_edges_per_macro_edge_float ) * 1.0;
         const walberla::float64 p_affine_const_0_0_GRAY = macro_vertex_coord_id_0comp0;
         const walberla::float64 p_affine_const_0_1_GRAY = macro_vertex_coord_id_0comp1;
         const walberla::float64 p_affine_const_1_0_GRAY =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 );
         const walberla::float64 p_affine_const_1_1_GRAY =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 );
         const walberla::float64 p_affine_const_2_0_GRAY =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 );
         const walberla::float64 p_affine_const_2_1_GRAY =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 );
         const walberla::float64 jac_affine_0_0_GRAY = -p_affine_const_0_0_GRAY + p_affine_const_1_0_GRAY;
         const walberla::float64 jac_affine_0_1_GRAY = -p_affine_const_0_0_GRAY + p_affine_const_2_0_GRAY;
         const walberla::float64 jac_affine_1_0_GRAY = -p_affine_const_0_1_GRAY + p_affine_const_1_1_GRAY;
         const walberla::float64 jac_affine_1_1_GRAY = -p_affine_const_0_1_GRAY + p_affine_const_2_1_GRAY;
         const walberla::float64 tmp_coords_jac_1_GRAY =
             1.0 / ( jac_affine_0_0_GRAY * jac_affine_1_1_GRAY - jac_affine_0_1_GRAY * jac_affine_1_0_GRAY );
         const walberla::float64 jac_affine_inv_0_0_GRAY = jac_affine_1_1_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_0_1_GRAY = -jac_affine_0_1_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_1_0_GRAY = -jac_affine_1_0_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_1_1_GRAY = jac_affine_0_0_GRAY * tmp_coords_jac_1_GRAY;
         for ( int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1 )
            for ( int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge; ctr_0 += 1 )
            {
               const walberla::float64 p_affine_0_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_0_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 mu_lin_dof_0 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 mu_lin_dof_1 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 mu_lin_dof_2 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_0 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_1 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_2 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_0 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_1 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_2        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 tmp_kernel_op_0 = -ux_dof_0 + ux_dof_1;
               const walberla::float64 tmp_kernel_op_1 = -ux_dof_0 + ux_dof_2;
               const walberla::float64 tmp_kernel_op_2 = -uy_dof_0 + uy_dof_1;
               const walberla::float64 tmp_kernel_op_3 = -uy_dof_0 + uy_dof_2;
               const walberla::float64 tmp_kernel_op_4 = jac_affine_inv_0_0_GRAY * tmp_kernel_op_2;
               const walberla::float64 tmp_kernel_op_5 = jac_affine_inv_0_1_GRAY * tmp_kernel_op_0;
               const walberla::float64 tmp_kernel_op_6 = jac_affine_inv_1_0_GRAY * tmp_kernel_op_3;
               const walberla::float64 tmp_kernel_op_7 = jac_affine_inv_1_1_GRAY * tmp_kernel_op_1;
               const walberla::float64 tmp_kernel_op_8 =
                   1.0 /
                   ( mu_star +
                     sigma_y *
                         pow( ( ( jac_affine_inv_0_0_GRAY * tmp_kernel_op_0 + jac_affine_inv_1_0_GRAY * tmp_kernel_op_1 ) *
                                ( jac_affine_inv_0_0_GRAY * tmp_kernel_op_0 + jac_affine_inv_1_0_GRAY * tmp_kernel_op_1 ) ) +
                                  ( ( jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 + jac_affine_inv_1_1_GRAY * tmp_kernel_op_3 ) *
                                    ( jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 + jac_affine_inv_1_1_GRAY * tmp_kernel_op_3 ) ) +
                                  ( tmp_kernel_op_4 + tmp_kernel_op_5 + tmp_kernel_op_6 + tmp_kernel_op_7 ) *
                                      ( tmp_kernel_op_4 * 0.5 + tmp_kernel_op_5 * 0.5 + tmp_kernel_op_6 * 0.5 +
                                        tmp_kernel_op_7 * 0.5 ),
                              -0.50000000000000000 ) ) *
                   1.0;
               const walberla::float64 elMat_0_0 = 1.0 / ( 1.0 / ( mu_lin_dof_0 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const int64_t           elMat_0_1 = 0;
               const int64_t           elMat_0_2 = 0;
               const int64_t           elMat_1_0 = 0;
               const walberla::float64 elMat_1_1 = 1.0 / ( 1.0 / ( mu_lin_dof_1 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const int64_t           elMat_1_2 = 0;
               const int64_t           elMat_2_0 = 0;
               const int64_t           elMat_2_1 = 0;
               const walberla::float64 elMat_2_2 = 1.0 / ( 1.0 / ( mu_lin_dof_2 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;

               std::vector< uint_t > _data_rowIdx( 3 );
               std::vector< uint_t > _data_colIdx( 3 );
               std::vector< real_t > _data_mat( 9 );

               _data_rowIdx[0] = ( (uint64_t) ( _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_rowIdx[1] = ( (uint64_t) ( _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[2] = ( (uint64_t) ( _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[0] = ( (uint64_t) ( _data_src[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_colIdx[1] = ( (uint64_t) ( _data_src[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[2] = ( (uint64_t) ( _data_src[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );

               /* Apply basis transformation */

               _data_mat[0] = ( (real_t) ( elMat_0_0 ) );
               _data_mat[1] = ( (real_t) ( elMat_0_1 ) );
               _data_mat[2] = ( (real_t) ( elMat_0_2 ) );
               _data_mat[3] = ( (real_t) ( elMat_1_0 ) );
               _data_mat[4] = ( (real_t) ( elMat_1_1 ) );
               _data_mat[5] = ( (real_t) ( elMat_1_2 ) );
               _data_mat[6] = ( (real_t) ( elMat_2_0 ) );
               _data_mat[7] = ( (real_t) ( elMat_2_1 ) );
               _data_mat[8] = ( (real_t) ( elMat_2_2 ) );

               mat->addValues( _data_rowIdx, _data_colIdx, _data_mat );
            }
      }
      {
         /* FaceType.BLUE */
         const walberla::float64 tmp_coords_jac_0_BLUE = 1.0 / ( micro_edges_per_macro_edge_float ) * 1.0;
         const walberla::float64 tmp_coords_jac_1_BLUE =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 );
         const walberla::float64 tmp_coords_jac_2_BLUE =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 );
         const walberla::float64 tmp_coords_jac_3_BLUE =
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 );
         const walberla::float64 tmp_coords_jac_4_BLUE =
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 );
         const walberla::float64 p_affine_const_0_0_BLUE = tmp_coords_jac_1_BLUE;
         const walberla::float64 p_affine_const_0_1_BLUE = tmp_coords_jac_2_BLUE;
         const walberla::float64 p_affine_const_1_0_BLUE = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3_BLUE;
         const walberla::float64 p_affine_const_1_1_BLUE = macro_vertex_coord_id_0comp1 + tmp_coords_jac_4_BLUE;
         const walberla::float64 p_affine_const_2_0_BLUE = tmp_coords_jac_1_BLUE + tmp_coords_jac_3_BLUE;
         const walberla::float64 p_affine_const_2_1_BLUE = tmp_coords_jac_2_BLUE + tmp_coords_jac_4_BLUE;
         const walberla::float64 jac_affine_0_0_BLUE     = -p_affine_const_0_0_BLUE + p_affine_const_1_0_BLUE;
         const walberla::float64 jac_affine_0_1_BLUE     = -p_affine_const_0_0_BLUE + p_affine_const_2_0_BLUE;
         const walberla::float64 jac_affine_1_0_BLUE     = -p_affine_const_0_1_BLUE + p_affine_const_1_1_BLUE;
         const walberla::float64 jac_affine_1_1_BLUE     = -p_affine_const_0_1_BLUE + p_affine_const_2_1_BLUE;
         const walberla::float64 tmp_coords_jac_5_BLUE =
             1.0 / ( jac_affine_0_0_BLUE * jac_affine_1_1_BLUE - jac_affine_0_1_BLUE * jac_affine_1_0_BLUE );
         const walberla::float64 jac_affine_inv_0_0_BLUE = jac_affine_1_1_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_0_1_BLUE = -jac_affine_0_1_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_1_0_BLUE = -jac_affine_1_0_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_1_1_BLUE = jac_affine_0_0_BLUE * tmp_coords_jac_5_BLUE;
         for ( int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1 )
            for ( int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1; ctr_0 += 1 )
            {
               const walberla::float64 p_affine_0_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_0_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_1_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 mu_lin_dof_0 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 mu_lin_dof_1 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 mu_lin_dof_2 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_0 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_1 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_2 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_0 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_1        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_2        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 tmp_kernel_op_0 = -ux_dof_0 + ux_dof_1;
               const walberla::float64 tmp_kernel_op_1 = -ux_dof_0 + ux_dof_2;
               const walberla::float64 tmp_kernel_op_2 = -uy_dof_0 + uy_dof_1;
               const walberla::float64 tmp_kernel_op_3 = -uy_dof_0 + uy_dof_2;
               const walberla::float64 tmp_kernel_op_4 = jac_affine_inv_0_0_BLUE * tmp_kernel_op_2;
               const walberla::float64 tmp_kernel_op_5 = jac_affine_inv_0_1_BLUE * tmp_kernel_op_0;
               const walberla::float64 tmp_kernel_op_6 = jac_affine_inv_1_0_BLUE * tmp_kernel_op_3;
               const walberla::float64 tmp_kernel_op_7 = jac_affine_inv_1_1_BLUE * tmp_kernel_op_1;
               const walberla::float64 tmp_kernel_op_8 =
                   1.0 /
                   ( mu_star +
                     sigma_y *
                         pow( ( ( jac_affine_inv_0_0_BLUE * tmp_kernel_op_0 + jac_affine_inv_1_0_BLUE * tmp_kernel_op_1 ) *
                                ( jac_affine_inv_0_0_BLUE * tmp_kernel_op_0 + jac_affine_inv_1_0_BLUE * tmp_kernel_op_1 ) ) +
                                  ( ( jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 + jac_affine_inv_1_1_BLUE * tmp_kernel_op_3 ) *
                                    ( jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 + jac_affine_inv_1_1_BLUE * tmp_kernel_op_3 ) ) +
                                  ( tmp_kernel_op_4 + tmp_kernel_op_5 + tmp_kernel_op_6 + tmp_kernel_op_7 ) *
                                      ( tmp_kernel_op_4 * 0.5 + tmp_kernel_op_5 * 0.5 + tmp_kernel_op_6 * 0.5 +
                                        tmp_kernel_op_7 * 0.5 ),
                              -0.50000000000000000 ) ) *
                   1.0;
               const walberla::float64 elMat_0_0 = 1.0 / ( 1.0 / ( mu_lin_dof_0 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const int64_t           elMat_0_1 = 0;
               const int64_t           elMat_0_2 = 0;
               const int64_t           elMat_1_0 = 0;
               const walberla::float64 elMat_1_1 = 1.0 / ( 1.0 / ( mu_lin_dof_1 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const int64_t           elMat_1_2 = 0;
               const int64_t           elMat_2_0 = 0;
               const int64_t           elMat_2_1 = 0;
               const walberla::float64 elMat_2_2 = 1.0 / ( 1.0 / ( mu_lin_dof_2 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;

               std::vector< uint_t > _data_rowIdx( 3 );
               std::vector< uint_t > _data_colIdx( 3 );
               std::vector< real_t > _data_mat( 9 );

               _data_rowIdx[0] = ( (uint64_t) ( _data_dst[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[1] = ( (uint64_t) ( _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_rowIdx[2] = ( (uint64_t) ( _data_dst[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[0] = ( (uint64_t) ( _data_src[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[1] = ( (uint64_t) ( _data_src[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[2] = ( (uint64_t) ( _data_src[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                          ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] ) );

               /* Apply basis transformation */

               _data_mat[0] = ( (real_t) ( elMat_0_0 ) );
               _data_mat[1] = ( (real_t) ( elMat_0_1 ) );
               _data_mat[2] = ( (real_t) ( elMat_0_2 ) );
               _data_mat[3] = ( (real_t) ( elMat_1_0 ) );
               _data_mat[4] = ( (real_t) ( elMat_1_1 ) );
               _data_mat[5] = ( (real_t) ( elMat_1_2 ) );
               _data_mat[6] = ( (real_t) ( elMat_2_0 ) );
               _data_mat[7] = ( (real_t) ( elMat_2_1 ) );
               _data_mat[8] = ( (real_t) ( elMat_2_2 ) );

               mat->addValues( _data_rowIdx, _data_colIdx, _data_mat );
            }
      }
   }
}
void P2EvaluateViscosityViscoplastic::computeInverseDiagonalOperatorValues_P2EvaluateViscosityViscoplastic_macro_2D(
    walberla::float64* RESTRICT _data_invDiag_,
    walberla::float64* RESTRICT _data_mu_lin,
    walberla::float64* RESTRICT _data_ux,
    walberla::float64* RESTRICT _data_uy,
    walberla::float64           macro_vertex_coord_id_0comp0,
    walberla::float64           macro_vertex_coord_id_0comp1,
    walberla::float64           macro_vertex_coord_id_1comp0,
    walberla::float64           macro_vertex_coord_id_1comp1,
    walberla::float64           macro_vertex_coord_id_2comp0,
    walberla::float64           macro_vertex_coord_id_2comp1,
    int64_t                     micro_edges_per_macro_edge,
    walberla::float64           micro_edges_per_macro_edge_float,
    walberla::float64           mu_star,
    walberla::float64           sigma_y ) const
{
   {
      {
         /* FaceType.GRAY */
         const walberla::float64 tmp_coords_jac_0_GRAY   = 1.0 / ( micro_edges_per_macro_edge_float ) * 1.0;
         const walberla::float64 p_affine_const_0_0_GRAY = macro_vertex_coord_id_0comp0;
         const walberla::float64 p_affine_const_0_1_GRAY = macro_vertex_coord_id_0comp1;
         const walberla::float64 p_affine_const_1_0_GRAY =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 );
         const walberla::float64 p_affine_const_1_1_GRAY =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 );
         const walberla::float64 p_affine_const_2_0_GRAY =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 );
         const walberla::float64 p_affine_const_2_1_GRAY =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_GRAY * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 );
         const walberla::float64 jac_affine_0_0_GRAY = -p_affine_const_0_0_GRAY + p_affine_const_1_0_GRAY;
         const walberla::float64 jac_affine_0_1_GRAY = -p_affine_const_0_0_GRAY + p_affine_const_2_0_GRAY;
         const walberla::float64 jac_affine_1_0_GRAY = -p_affine_const_0_1_GRAY + p_affine_const_1_1_GRAY;
         const walberla::float64 jac_affine_1_1_GRAY = -p_affine_const_0_1_GRAY + p_affine_const_2_1_GRAY;
         const walberla::float64 tmp_coords_jac_1_GRAY =
             1.0 / ( jac_affine_0_0_GRAY * jac_affine_1_1_GRAY - jac_affine_0_1_GRAY * jac_affine_1_0_GRAY );
         const walberla::float64 jac_affine_inv_0_0_GRAY = jac_affine_1_1_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_0_1_GRAY = -jac_affine_0_1_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_1_0_GRAY = -jac_affine_1_0_GRAY * tmp_coords_jac_1_GRAY;
         const walberla::float64 jac_affine_inv_1_1_GRAY = jac_affine_0_0_GRAY * tmp_coords_jac_1_GRAY;
         for ( int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1 )
            for ( int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge; ctr_0 += 1 )
            {
               const walberla::float64 p_affine_0_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_0_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 mu_lin_dof_0 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 mu_lin_dof_1 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 mu_lin_dof_2 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_0 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_1 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_2 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_0 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_1 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_2        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 tmp_kernel_op_0 = -ux_dof_0 + ux_dof_1;
               const walberla::float64 tmp_kernel_op_1 = -ux_dof_0 + ux_dof_2;
               const walberla::float64 tmp_kernel_op_2 = -uy_dof_0 + uy_dof_1;
               const walberla::float64 tmp_kernel_op_3 = -uy_dof_0 + uy_dof_2;
               const walberla::float64 tmp_kernel_op_4 = jac_affine_inv_0_0_GRAY * tmp_kernel_op_2;
               const walberla::float64 tmp_kernel_op_5 = jac_affine_inv_0_1_GRAY * tmp_kernel_op_0;
               const walberla::float64 tmp_kernel_op_6 = jac_affine_inv_1_0_GRAY * tmp_kernel_op_3;
               const walberla::float64 tmp_kernel_op_7 = jac_affine_inv_1_1_GRAY * tmp_kernel_op_1;
               const walberla::float64 tmp_kernel_op_8 =
                   1.0 /
                   ( mu_star +
                     sigma_y *
                         pow( ( ( jac_affine_inv_0_0_GRAY * tmp_kernel_op_0 + jac_affine_inv_1_0_GRAY * tmp_kernel_op_1 ) *
                                ( jac_affine_inv_0_0_GRAY * tmp_kernel_op_0 + jac_affine_inv_1_0_GRAY * tmp_kernel_op_1 ) ) +
                                  ( ( jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 + jac_affine_inv_1_1_GRAY * tmp_kernel_op_3 ) *
                                    ( jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 + jac_affine_inv_1_1_GRAY * tmp_kernel_op_3 ) ) +
                                  ( tmp_kernel_op_4 + tmp_kernel_op_5 + tmp_kernel_op_6 + tmp_kernel_op_7 ) *
                                      ( tmp_kernel_op_4 * 0.5 + tmp_kernel_op_5 * 0.5 + tmp_kernel_op_6 * 0.5 +
                                        tmp_kernel_op_7 * 0.5 ),
                              -0.50000000000000000 ) ) *
                   1.0;
               const walberla::float64 elMatDiag_0 = 1.0 / ( 1.0 / ( mu_lin_dof_0 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatDiag_1 = 1.0 / ( 1.0 / ( mu_lin_dof_1 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatDiag_2 = 1.0 / ( 1.0 / ( mu_lin_dof_2 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               _data_invDiag_[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatDiag_0 +
                   _data_invDiag_[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_invDiag_[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_1 +
                   _data_invDiag_[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag_[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                              ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_2 + _data_invDiag_[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
            }
      }
      {
         /* FaceType.BLUE */
         const walberla::float64 tmp_coords_jac_0_BLUE = 1.0 / ( micro_edges_per_macro_edge_float ) * 1.0;
         const walberla::float64 tmp_coords_jac_1_BLUE =
             macro_vertex_coord_id_0comp0 +
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 );
         const walberla::float64 tmp_coords_jac_2_BLUE =
             macro_vertex_coord_id_0comp1 +
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 );
         const walberla::float64 tmp_coords_jac_3_BLUE =
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 );
         const walberla::float64 tmp_coords_jac_4_BLUE =
             tmp_coords_jac_0_BLUE * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 );
         const walberla::float64 p_affine_const_0_0_BLUE = tmp_coords_jac_1_BLUE;
         const walberla::float64 p_affine_const_0_1_BLUE = tmp_coords_jac_2_BLUE;
         const walberla::float64 p_affine_const_1_0_BLUE = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3_BLUE;
         const walberla::float64 p_affine_const_1_1_BLUE = macro_vertex_coord_id_0comp1 + tmp_coords_jac_4_BLUE;
         const walberla::float64 p_affine_const_2_0_BLUE = tmp_coords_jac_1_BLUE + tmp_coords_jac_3_BLUE;
         const walberla::float64 p_affine_const_2_1_BLUE = tmp_coords_jac_2_BLUE + tmp_coords_jac_4_BLUE;
         const walberla::float64 jac_affine_0_0_BLUE     = -p_affine_const_0_0_BLUE + p_affine_const_1_0_BLUE;
         const walberla::float64 jac_affine_0_1_BLUE     = -p_affine_const_0_0_BLUE + p_affine_const_2_0_BLUE;
         const walberla::float64 jac_affine_1_0_BLUE     = -p_affine_const_0_1_BLUE + p_affine_const_1_1_BLUE;
         const walberla::float64 jac_affine_1_1_BLUE     = -p_affine_const_0_1_BLUE + p_affine_const_2_1_BLUE;
         const walberla::float64 tmp_coords_jac_5_BLUE =
             1.0 / ( jac_affine_0_0_BLUE * jac_affine_1_1_BLUE - jac_affine_0_1_BLUE * jac_affine_1_0_BLUE );
         const walberla::float64 jac_affine_inv_0_0_BLUE = jac_affine_1_1_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_0_1_BLUE = -jac_affine_0_1_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_1_0_BLUE = -jac_affine_1_0_BLUE * tmp_coords_jac_5_BLUE;
         const walberla::float64 jac_affine_inv_1_1_BLUE = jac_affine_0_0_BLUE * tmp_coords_jac_5_BLUE;
         for ( int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1 )
            for ( int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1; ctr_0 += 1 )
            {
               const walberla::float64 p_affine_0_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_0_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 ) );
               const walberla::float64 p_affine_1_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_1_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_0 + 1 ) ) +
                   1.0 / ( micro_edges_per_macro_edge_float ) * ( -macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1 ) *
                       1.0 * ( ( walberla::float64 )( ctr_1 + 1 ) );
               const walberla::float64 mu_lin_dof_0 =
                   _data_mu_lin[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 mu_lin_dof_1 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 mu_lin_dof_2 = _data_mu_lin[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_0 =
                   _data_ux[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 ux_dof_1 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 ux_dof_2 = _data_ux[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_0 =
                   _data_uy[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 uy_dof_1        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_2        = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 tmp_kernel_op_0 = -ux_dof_0 + ux_dof_1;
               const walberla::float64 tmp_kernel_op_1 = -ux_dof_0 + ux_dof_2;
               const walberla::float64 tmp_kernel_op_2 = -uy_dof_0 + uy_dof_1;
               const walberla::float64 tmp_kernel_op_3 = -uy_dof_0 + uy_dof_2;
               const walberla::float64 tmp_kernel_op_4 = jac_affine_inv_0_0_BLUE * tmp_kernel_op_2;
               const walberla::float64 tmp_kernel_op_5 = jac_affine_inv_0_1_BLUE * tmp_kernel_op_0;
               const walberla::float64 tmp_kernel_op_6 = jac_affine_inv_1_0_BLUE * tmp_kernel_op_3;
               const walberla::float64 tmp_kernel_op_7 = jac_affine_inv_1_1_BLUE * tmp_kernel_op_1;
               const walberla::float64 tmp_kernel_op_8 =
                   1.0 /
                   ( mu_star +
                     sigma_y *
                         pow( ( ( jac_affine_inv_0_0_BLUE * tmp_kernel_op_0 + jac_affine_inv_1_0_BLUE * tmp_kernel_op_1 ) *
                                ( jac_affine_inv_0_0_BLUE * tmp_kernel_op_0 + jac_affine_inv_1_0_BLUE * tmp_kernel_op_1 ) ) +
                                  ( ( jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 + jac_affine_inv_1_1_BLUE * tmp_kernel_op_3 ) *
                                    ( jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 + jac_affine_inv_1_1_BLUE * tmp_kernel_op_3 ) ) +
                                  ( tmp_kernel_op_4 + tmp_kernel_op_5 + tmp_kernel_op_6 + tmp_kernel_op_7 ) *
                                      ( tmp_kernel_op_4 * 0.5 + tmp_kernel_op_5 * 0.5 + tmp_kernel_op_6 * 0.5 +
                                        tmp_kernel_op_7 * 0.5 ),
                              -0.50000000000000000 ) ) *
                   1.0;
               const walberla::float64 elMatDiag_0 = 1.0 / ( 1.0 / ( mu_lin_dof_0 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatDiag_1 = 1.0 / ( 1.0 / ( mu_lin_dof_1 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               const walberla::float64 elMatDiag_2 = 1.0 / ( 1.0 / ( mu_lin_dof_2 ) * 1.0 + tmp_kernel_op_8 ) * 2.0;
               _data_invDiag_[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_0 +
                   _data_invDiag_[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag_[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                              ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_1 + _data_invDiag_[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_invDiag_[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                              ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_2 + _data_invDiag_[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
            }
      }
   }
}

} // namespace operatorgeneration

} // namespace hyteg
