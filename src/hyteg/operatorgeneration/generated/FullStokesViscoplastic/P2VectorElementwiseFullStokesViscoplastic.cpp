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

#include "P2VectorElementwiseFullStokesViscoplastic.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P2VectorElementwiseFullStokesViscoplastic::P2VectorElementwiseFullStokesViscoplastic(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel,
    const P1Function< walberla::float64 >&     _mu_lin,
    const P1Function< walberla::float64 >&     _ux,
    const P1Function< walberla::float64 >&     _uy,
    walberla::float64                          mu_star_P2VectorElementwiseFullStokesViscoplastic,
    walberla::float64                          sigma_y_P2VectorElementwiseFullStokesViscoplastic )
: Operator( storage, minLevel, maxLevel )
, mu_lin( _mu_lin )
, ux( _ux )
, uy( _uy )
, mu_star_P2VectorElementwiseFullStokesViscoplastic_( mu_star_P2VectorElementwiseFullStokesViscoplastic )
, sigma_y_P2VectorElementwiseFullStokesViscoplastic_( sigma_y_P2VectorElementwiseFullStokesViscoplastic )
{}

void P2VectorElementwiseFullStokesViscoplastic::apply( const P2VectorFunction< walberla::float64 >& src,
                                                       const P2VectorFunction< walberla::float64 >& dst,
                                                       uint_t                                       level,
                                                       DoFType                                      flag,
                                                       UpdateType                                   updateType ) const
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
      communication::syncVectorFunctionBetweenPrimitives( src, level, communication::syncDirection_t::LOW2HIGH );
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
         walberla::float64* _data_src_vertex_0 =
             face.getData( src[0].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         walberla::float64* _data_src_edge_0 = face.getData( src[0].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         walberla::float64* _data_src_vertex_1 =
             face.getData( src[1].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         walberla::float64* _data_src_edge_1 = face.getData( src[1].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

         walberla::float64* _data_dst_vertex_0 =
             face.getData( dst[0].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         walberla::float64* _data_dst_edge_0 = face.getData( dst[0].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         walberla::float64* _data_dst_vertex_1 =
             face.getData( dst[1].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         walberla::float64* _data_dst_edge_1 = face.getData( dst[1].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

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
               auto arrayIdx                = vertexdof::macroface::index( level, idx.x(), idx.y() );
               _data_dst_vertex_0[arrayIdx] = walberla::numeric_cast< walberla::float64 >( 0 );
               _data_dst_vertex_1[arrayIdx] = walberla::numeric_cast< walberla::float64 >( 0 );
            }
         }
         for ( const auto& idx : edgedof::macroface::Iterator( level ) )
         {
            for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
            {
               if ( !edgedof::macroface::isInnerEdgeDoF( level, idx, orientation ) )
               {
                  auto arrayIdx              = edgedof::macroface::index( level, idx.x(), idx.y(), orientation );
                  _data_dst_edge_0[arrayIdx] = walberla::numeric_cast< walberla::float64 >( 0 );
                  _data_dst_edge_1[arrayIdx] = walberla::numeric_cast< walberla::float64 >( 0 );
               }
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

         apply_P2VectorElementwiseFullStokesViscoplastic_macro_2D(

             _data_dst_edge_0,
             _data_dst_edge_1,
             _data_dst_vertex_0,
             _data_dst_vertex_1,
             _data_mu_lin,
             _data_src_edge_0,
             _data_src_edge_1,
             _data_src_vertex_0,
             _data_src_vertex_1,
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
             mu_star_P2VectorElementwiseFullStokesViscoplastic_,
             sigma_y_P2VectorElementwiseFullStokesViscoplastic_ );

         this->timingTree_->stop( "kernel" );
      }

      // Push result to lower-dimensional primitives
      //
      this->timingTree_->start( "post-communication" );
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst[0].getVertexDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst[0].getVertexDoFFunction().communicateAdditively< Face, Vertex >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst[0].getEdgeDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst[1].getVertexDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst[1].getVertexDoFFunction().communicateAdditively< Face, Vertex >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst[1].getEdgeDoFFunction().communicateAdditively< Face, Edge >(
          level, DoFType::All ^ flag, *storage_, updateType == Replace );
      this->timingTree_->stop( "post-communication" );
   }

   this->stopTiming( "apply" );
}
void P2VectorElementwiseFullStokesViscoplastic::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                          const P2VectorFunction< idx_t >&            src,
                                                          const P2VectorFunction< idx_t >&            dst,
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
         idx_t* _data_src_vertex_0 = face.getData( src[0].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         idx_t* _data_src_edge_0   = face.getData( src[0].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         idx_t* _data_src_vertex_1 = face.getData( src[1].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         idx_t* _data_src_edge_1   = face.getData( src[1].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

         idx_t* _data_dst_vertex_0 = face.getData( dst[0].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         idx_t* _data_dst_edge_0   = face.getData( dst[0].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         idx_t* _data_dst_vertex_1 = face.getData( dst[1].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         idx_t* _data_dst_edge_1   = face.getData( dst[1].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

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

         toMatrix_P2VectorElementwiseFullStokesViscoplastic_macro_2D(

             _data_dst_edge_0,
             _data_dst_edge_1,
             _data_dst_vertex_0,
             _data_dst_vertex_1,
             _data_mu_lin,
             _data_src_edge_0,
             _data_src_edge_1,
             _data_src_vertex_0,
             _data_src_vertex_1,
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
             mu_star_P2VectorElementwiseFullStokesViscoplastic_,
             sigma_y_P2VectorElementwiseFullStokesViscoplastic_ );

         this->timingTree_->stop( "kernel" );
      }
   }
   this->stopTiming( "toMatrix" );
}
void P2VectorElementwiseFullStokesViscoplastic::computeInverseDiagonalOperatorValues()
{
   this->startTiming( "computeInverseDiagonalOperatorValues" );

   if ( invDiag_ == nullptr )
   {
      invDiag_ =
          std::make_shared< P2VectorFunction< walberla::float64 > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
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
         ( *invDiag_ )[0].invertElementwise( level );
         ( *invDiag_ )[1].invertElementwise( level );
         ( *invDiag_ )[2].invertElementwise( level );
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
            walberla::float64* _data_invDiag__vertex_0 =
                face.getData( ( *invDiag_ )[0].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
            walberla::float64* _data_invDiag__edge_0 =
                face.getData( ( *invDiag_ )[0].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
            walberla::float64* _data_invDiag__vertex_1 =
                face.getData( ( *invDiag_ )[1].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
            walberla::float64* _data_invDiag__edge_1 =
                face.getData( ( *invDiag_ )[1].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

            walberla::float64* _data_mu_lin = face.getData( mu_lin.getFaceDataID() )->getPointer( level );
            walberla::float64* _data_ux     = face.getData( ux.getFaceDataID() )->getPointer( level );
            walberla::float64* _data_uy     = face.getData( uy.getFaceDataID() )->getPointer( level );

            const auto micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
            const auto micro_edges_per_macro_edge_float = (walberla::float64) levelinfo::num_microedges_per_edge( level );
            const walberla::float64 macro_vertex_coord_id_0comp0 = (walberla::float64) face.getCoordinates()[0][0];
            const walberla::float64 macro_vertex_coord_id_0comp1 = (walberla::float64) face.getCoordinates()[0][1];
            const walberla::float64 macro_vertex_coord_id_1comp0 = (walberla::float64) face.getCoordinates()[1][0];
            const walberla::float64 macro_vertex_coord_id_1comp1 = (walberla::float64) face.getCoordinates()[1][1];
            const walberla::float64 macro_vertex_coord_id_2comp0 = (walberla::float64) face.getCoordinates()[2][0];
            const walberla::float64 macro_vertex_coord_id_2comp1 = (walberla::float64) face.getCoordinates()[2][1];

            this->timingTree_->start( "kernel" );

            computeInverseDiagonalOperatorValues_P2VectorElementwiseFullStokesViscoplastic_macro_2D(

                _data_invDiag__edge_0,
                _data_invDiag__edge_1,
                _data_invDiag__vertex_0,
                _data_invDiag__vertex_1,
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
                mu_star_P2VectorElementwiseFullStokesViscoplastic_,
                sigma_y_P2VectorElementwiseFullStokesViscoplastic_ );

            this->timingTree_->stop( "kernel" );
         }

         // Push result to lower-dimensional primitives
         //
         this->timingTree_->start( "post-communication" );
         // Note: We could avoid communication here by implementing the apply() also for the respective
         //       lower dimensional primitives!
         ( *invDiag_ )[0].getVertexDoFFunction().communicateAdditively< Face, Edge >( level );
         ( *invDiag_ )[0].getVertexDoFFunction().communicateAdditively< Face, Vertex >( level );
         ( *invDiag_ )[0].getEdgeDoFFunction().communicateAdditively< Face, Edge >( level );
         ( *invDiag_ )[1].getVertexDoFFunction().communicateAdditively< Face, Edge >( level );
         ( *invDiag_ )[1].getVertexDoFFunction().communicateAdditively< Face, Vertex >( level );
         ( *invDiag_ )[1].getEdgeDoFFunction().communicateAdditively< Face, Edge >( level );
         this->timingTree_->stop( "post-communication" );
         ( *invDiag_ )[0].invertElementwise( level );
         ( *invDiag_ )[1].invertElementwise( level );
      }
   }

   this->stopTiming( "computeInverseDiagonalOperatorValues" );
}
std::shared_ptr< P2VectorFunction< walberla::float64 > >
    P2VectorElementwiseFullStokesViscoplastic::getInverseDiagonalValues() const
{
   return invDiag_;
}
void P2VectorElementwiseFullStokesViscoplastic::apply_P2VectorElementwiseFullStokesViscoplastic_macro_2D(
    walberla::float64* RESTRICT _data_dst_edge_0,
    walberla::float64* RESTRICT _data_dst_edge_1,
    walberla::float64* RESTRICT _data_dst_vertex_0,
    walberla::float64* RESTRICT _data_dst_vertex_1,
    walberla::float64* RESTRICT _data_mu_lin,
    walberla::float64* RESTRICT _data_src_edge_0,
    walberla::float64* RESTRICT _data_src_edge_1,
    walberla::float64* RESTRICT _data_src_vertex_0,
    walberla::float64* RESTRICT _data_src_vertex_1,
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
      const walberla::float64 _data_q_w[] = { 0.16666666666666666, 0.16666666666666666, 0.16666666666666666 };

      const walberla::float64 _data_q_p_0[] = { 0.16666666666666666, 0.66666666666666663, 0.16666666666666666 };

      const walberla::float64 _data_q_p_1[] = { 0.66666666666666663, 0.16666666666666666, 0.16666666666666666 };

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
             jac_affine_0_0_GRAY * jac_affine_1_1_GRAY - jac_affine_0_1_GRAY * jac_affine_1_0_GRAY;
         const walberla::float64 tmp_coords_jac_2_GRAY   = 1.0 / ( tmp_coords_jac_1_GRAY );
         const walberla::float64 jac_affine_inv_0_0_GRAY = jac_affine_1_1_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_0_1_GRAY = -jac_affine_0_1_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_1_0_GRAY = -jac_affine_1_0_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_1_1_GRAY = jac_affine_0_0_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 abs_det_jac_affine_GRAY = abs( tmp_coords_jac_1_GRAY );
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
                   _data_src_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_1 = _data_src_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_2 = _data_src_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_3 =
                   _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_4 =
                   _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_5 =
                   _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_6 =
                   _data_src_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_7 = _data_src_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_8 = _data_src_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_9 =
                   _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_10 =
                   _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_11 =
                   _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
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
               const walberla::float64 uy_dof_2    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               walberla::float64       q_acc_0_0   = 0.0;
               walberla::float64       q_acc_0_1   = 0.0;
               walberla::float64       q_acc_0_2   = 0.0;
               walberla::float64       q_acc_0_3   = 0.0;
               walberla::float64       q_acc_0_4   = 0.0;
               walberla::float64       q_acc_0_5   = 0.0;
               walberla::float64       q_acc_0_6   = 0.0;
               walberla::float64       q_acc_0_7   = 0.0;
               walberla::float64       q_acc_0_8   = 0.0;
               walberla::float64       q_acc_0_9   = 0.0;
               walberla::float64       q_acc_0_10  = 0.0;
               walberla::float64       q_acc_0_11  = 0.0;
               walberla::float64       q_acc_1_1   = 0.0;
               walberla::float64       q_acc_1_2   = 0.0;
               walberla::float64       q_acc_1_3   = 0.0;
               walberla::float64       q_acc_1_4   = 0.0;
               walberla::float64       q_acc_1_5   = 0.0;
               walberla::float64       q_acc_1_6   = 0.0;
               walberla::float64       q_acc_1_7   = 0.0;
               walberla::float64       q_acc_1_8   = 0.0;
               walberla::float64       q_acc_1_9   = 0.0;
               walberla::float64       q_acc_1_10  = 0.0;
               walberla::float64       q_acc_1_11  = 0.0;
               walberla::float64       q_acc_2_2   = 0.0;
               walberla::float64       q_acc_2_3   = 0.0;
               walberla::float64       q_acc_2_4   = 0.0;
               walberla::float64       q_acc_2_5   = 0.0;
               walberla::float64       q_acc_2_6   = 0.0;
               walberla::float64       q_acc_2_7   = 0.0;
               walberla::float64       q_acc_2_8   = 0.0;
               walberla::float64       q_acc_2_9   = 0.0;
               walberla::float64       q_acc_2_10  = 0.0;
               walberla::float64       q_acc_2_11  = 0.0;
               walberla::float64       q_acc_3_3   = 0.0;
               walberla::float64       q_acc_3_4   = 0.0;
               walberla::float64       q_acc_3_5   = 0.0;
               walberla::float64       q_acc_3_6   = 0.0;
               walberla::float64       q_acc_3_7   = 0.0;
               walberla::float64       q_acc_3_8   = 0.0;
               walberla::float64       q_acc_3_9   = 0.0;
               walberla::float64       q_acc_3_10  = 0.0;
               walberla::float64       q_acc_3_11  = 0.0;
               walberla::float64       q_acc_4_4   = 0.0;
               walberla::float64       q_acc_4_5   = 0.0;
               walberla::float64       q_acc_4_6   = 0.0;
               walberla::float64       q_acc_4_7   = 0.0;
               walberla::float64       q_acc_4_8   = 0.0;
               walberla::float64       q_acc_4_9   = 0.0;
               walberla::float64       q_acc_4_10  = 0.0;
               walberla::float64       q_acc_4_11  = 0.0;
               walberla::float64       q_acc_5_5   = 0.0;
               walberla::float64       q_acc_5_6   = 0.0;
               walberla::float64       q_acc_5_7   = 0.0;
               walberla::float64       q_acc_5_8   = 0.0;
               walberla::float64       q_acc_5_9   = 0.0;
               walberla::float64       q_acc_5_10  = 0.0;
               walberla::float64       q_acc_5_11  = 0.0;
               walberla::float64       q_acc_6_6   = 0.0;
               walberla::float64       q_acc_6_7   = 0.0;
               walberla::float64       q_acc_6_8   = 0.0;
               walberla::float64       q_acc_6_9   = 0.0;
               walberla::float64       q_acc_6_10  = 0.0;
               walberla::float64       q_acc_6_11  = 0.0;
               walberla::float64       q_acc_7_7   = 0.0;
               walberla::float64       q_acc_7_8   = 0.0;
               walberla::float64       q_acc_7_9   = 0.0;
               walberla::float64       q_acc_7_10  = 0.0;
               walberla::float64       q_acc_7_11  = 0.0;
               walberla::float64       q_acc_8_8   = 0.0;
               walberla::float64       q_acc_8_9   = 0.0;
               walberla::float64       q_acc_8_10  = 0.0;
               walberla::float64       q_acc_8_11  = 0.0;
               walberla::float64       q_acc_9_9   = 0.0;
               walberla::float64       q_acc_9_10  = 0.0;
               walberla::float64       q_acc_9_11  = 0.0;
               walberla::float64       q_acc_10_10 = 0.0;
               walberla::float64       q_acc_10_11 = 0.0;
               walberla::float64       q_acc_11_11 = 0.0;
               for ( int64_t q = 0; q < 3; q += 1 )
               {
                  const walberla::float64 tmp_qloop_0  = 4.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_1  = 4.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_2  = tmp_qloop_0 + tmp_qloop_1 - 3.0;
                  const walberla::float64 tmp_qloop_3  = jac_affine_inv_0_0_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_4  = jac_affine_inv_1_0_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_5  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_6  = abs_det_jac_affine_GRAY * tmp_qloop_5;
                  const walberla::float64 tmp_qloop_7  = jac_affine_inv_0_1_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_8  = jac_affine_inv_1_1_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_9  = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_10 = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_11 = -ux_dof_0 + ux_dof_1;
                  const walberla::float64 tmp_qloop_12 = -ux_dof_0 + ux_dof_2;
                  const walberla::float64 tmp_qloop_13 = -uy_dof_0 + uy_dof_1;
                  const walberla::float64 tmp_qloop_14 = -uy_dof_0 + uy_dof_2;
                  const walberla::float64 tmp_qloop_15 = jac_affine_inv_0_0_GRAY * tmp_qloop_13;
                  const walberla::float64 tmp_qloop_16 = jac_affine_inv_0_1_GRAY * tmp_qloop_11;
                  const walberla::float64 tmp_qloop_17 = jac_affine_inv_1_0_GRAY * tmp_qloop_14;
                  const walberla::float64 tmp_qloop_18 = jac_affine_inv_1_1_GRAY * tmp_qloop_12;
                  const walberla::float64 tmp_qloop_19 =
                      1.0 /
                      ( 1.0 /
                            ( mu_star +
                              sigma_y * 1.0 /
                                  ( pow( ( ( jac_affine_inv_0_0_GRAY * tmp_qloop_11 + jac_affine_inv_1_0_GRAY * tmp_qloop_12 ) *
                                           ( jac_affine_inv_0_0_GRAY * tmp_qloop_11 + jac_affine_inv_1_0_GRAY * tmp_qloop_12 ) ) +
                                             ( ( jac_affine_inv_0_1_GRAY * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_GRAY * tmp_qloop_14 ) *
                                               ( jac_affine_inv_0_1_GRAY * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_GRAY * tmp_qloop_14 ) ) +
                                             ( tmp_qloop_15 + tmp_qloop_16 + tmp_qloop_17 + tmp_qloop_18 ) *
                                                 ( tmp_qloop_15 * 0.5 + tmp_qloop_16 * 0.5 + tmp_qloop_17 * 0.5 +
                                                   tmp_qloop_18 * 0.5 ),
                                         0.50000000000000000 ) +
                                    9.9999999999999995e-21 ) ) *
                            1.0 +
                        1.0 /
                            ( mu_lin_dof_0 * ( 1.0 - _data_q_p_0[q] - _data_q_p_1[q] ) + mu_lin_dof_1 * _data_q_p_0[q] +
                              mu_lin_dof_2 * _data_q_p_1[q] ) *
                            1.0 ) *
                      _data_q_w[q];
                  const walberla::float64 tmp_qloop_20 = tmp_qloop_19 * 2.0;
                  const walberla::float64 tmp_qloop_21 = tmp_qloop_0 - 1.0;
                  const walberla::float64 tmp_qloop_22 = jac_affine_inv_0_0_GRAY * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_23 = tmp_qloop_6 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_24 = tmp_qloop_5 * 2.0;
                  const walberla::float64 tmp_qloop_25 = jac_affine_inv_0_1_GRAY * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_26 = tmp_qloop_1 - 1.0;
                  const walberla::float64 tmp_qloop_27 = jac_affine_inv_1_0_GRAY * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_28 = jac_affine_inv_1_1_GRAY * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_29 = 2.6666666666666665 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_30 = jac_affine_inv_1_0_GRAY * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_31 = 2.6666666666666665 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_32 = jac_affine_inv_0_0_GRAY * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_33 = tmp_qloop_30 + tmp_qloop_32;
                  const walberla::float64 tmp_qloop_34 = jac_affine_inv_1_0_GRAY * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_35 = jac_affine_inv_0_0_GRAY * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_36 = tmp_qloop_34 + tmp_qloop_35;
                  const walberla::float64 tmp_qloop_37 = jac_affine_inv_1_1_GRAY * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_38 = jac_affine_inv_0_1_GRAY * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_39 = tmp_qloop_37 + tmp_qloop_38;
                  const walberla::float64 tmp_qloop_40 = -tmp_qloop_0 - 8.0 * _data_q_p_1[q] + 4.0;
                  const walberla::float64 tmp_qloop_41 =
                      jac_affine_inv_1_0_GRAY * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_32;
                  const walberla::float64 tmp_qloop_42 = jac_affine_inv_1_0_GRAY * tmp_qloop_40 - tmp_qloop_35;
                  const walberla::float64 tmp_qloop_43 = jac_affine_inv_1_1_GRAY * tmp_qloop_40 - tmp_qloop_38;
                  const walberla::float64 tmp_qloop_44 = -tmp_qloop_1 - 8.0 * _data_q_p_0[q] + 4.0;
                  const walberla::float64 tmp_qloop_45 =
                      jac_affine_inv_0_0_GRAY * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_30;
                  const walberla::float64 tmp_qloop_46 = jac_affine_inv_0_0_GRAY * tmp_qloop_44 - tmp_qloop_34;
                  const walberla::float64 tmp_qloop_47 = jac_affine_inv_0_1_GRAY * tmp_qloop_44 - tmp_qloop_37;
                  const walberla::float64 tmp_qloop_48 = tmp_qloop_7 * 0.66666666666666663 + tmp_qloop_8 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_49 = abs_det_jac_affine_GRAY * tmp_qloop_22;
                  const walberla::float64 tmp_qloop_50 = abs_det_jac_affine_GRAY * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_51 = jac_affine_inv_1_1_GRAY * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_52 = jac_affine_inv_0_1_GRAY * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_53 = tmp_qloop_51 + tmp_qloop_52;
                  const walberla::float64 tmp_qloop_54 = abs_det_jac_affine_GRAY * tmp_qloop_36;
                  const walberla::float64 tmp_qloop_55 =
                      jac_affine_inv_1_1_GRAY * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_52;
                  const walberla::float64 tmp_qloop_56 = abs_det_jac_affine_GRAY * tmp_qloop_42;
                  const walberla::float64 tmp_qloop_57 =
                      jac_affine_inv_0_1_GRAY * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_51;
                  const walberla::float64 tmp_qloop_58  = abs_det_jac_affine_GRAY * tmp_qloop_46;
                  const walberla::float64 tmp_qloop_59  = ( jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY );
                  const walberla::float64 tmp_qloop_60  = ( tmp_qloop_21 * tmp_qloop_21 );
                  const walberla::float64 tmp_qloop_61  = abs_det_jac_affine_GRAY * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_62  = tmp_qloop_60 * tmp_qloop_61;
                  const walberla::float64 tmp_qloop_63  = tmp_qloop_60 * 2.0;
                  const walberla::float64 tmp_qloop_64  = ( jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY );
                  const walberla::float64 tmp_qloop_65  = tmp_qloop_60 * 1.0;
                  const walberla::float64 tmp_qloop_66  = tmp_qloop_22 * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_67  = tmp_qloop_25 * 1.0;
                  const walberla::float64 tmp_qloop_68  = tmp_qloop_36 * 2.0;
                  const walberla::float64 tmp_qloop_69  = tmp_qloop_42 * 2.0;
                  const walberla::float64 tmp_qloop_70  = tmp_qloop_46 * 2.0;
                  const walberla::float64 tmp_qloop_71  = abs_det_jac_affine_GRAY * tmp_qloop_19 * 0.66666666666666674;
                  const walberla::float64 tmp_qloop_72  = abs_det_jac_affine_GRAY * tmp_qloop_57;
                  const walberla::float64 tmp_qloop_73  = ( jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY );
                  const walberla::float64 tmp_qloop_74  = ( tmp_qloop_26 * tmp_qloop_26 );
                  const walberla::float64 tmp_qloop_75  = tmp_qloop_61 * tmp_qloop_74;
                  const walberla::float64 tmp_qloop_76  = tmp_qloop_74 * 2.0;
                  const walberla::float64 tmp_qloop_77  = ( jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY );
                  const walberla::float64 tmp_qloop_78  = tmp_qloop_74 * 1.0;
                  const walberla::float64 tmp_qloop_79  = tmp_qloop_28 * 1.0;
                  const walberla::float64 tmp_qloop_80  = 2.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_81  = jac_affine_inv_1_1_GRAY * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_82  = 2.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_83  = jac_affine_inv_0_1_GRAY * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_84  = tmp_qloop_81 + tmp_qloop_83;
                  const walberla::float64 tmp_qloop_85  = tmp_qloop_84 * 2.0;
                  const walberla::float64 tmp_qloop_86  = tmp_qloop_54 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_87  = jac_affine_inv_1_1_GRAY * tmp_qloop_40 * 0.5 - tmp_qloop_83;
                  const walberla::float64 tmp_qloop_88  = tmp_qloop_87 * 2.0;
                  const walberla::float64 tmp_qloop_89  = tmp_qloop_56 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_90  = jac_affine_inv_0_1_GRAY * tmp_qloop_44 * 0.5 - tmp_qloop_81;
                  const walberla::float64 tmp_qloop_91  = tmp_qloop_90 * 2.0;
                  const walberla::float64 tmp_qloop_92  = tmp_qloop_58 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_93  = abs_det_jac_affine_GRAY * tmp_qloop_9;
                  const walberla::float64 tmp_qloop_94  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_95  = tmp_qloop_93 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_96  = tmp_qloop_9 * 2.0;
                  const walberla::float64 tmp_qloop_97  = tmp_qloop_25 * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_98  = tmp_qloop_22 * 1.0;
                  const walberla::float64 tmp_qloop_99  = abs_det_jac_affine_GRAY * tmp_qloop_25;
                  const walberla::float64 tmp_qloop_100 = tmp_qloop_39 * 2.0;
                  const walberla::float64 tmp_qloop_101 = tmp_qloop_43 * 2.0;
                  const walberla::float64 tmp_qloop_102 = tmp_qloop_47 * 2.0;
                  const walberla::float64 tmp_qloop_103 = abs_det_jac_affine_GRAY * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_104 = tmp_qloop_27 * 1.0;
                  const walberla::float64 tmp_qloop_105 = abs_det_jac_affine_GRAY * tmp_qloop_39;
                  const walberla::float64 tmp_qloop_106 = jac_affine_inv_1_0_GRAY * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_107 = jac_affine_inv_0_0_GRAY * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_108 = tmp_qloop_106 * 2.0 + tmp_qloop_107 * 2.0;
                  const walberla::float64 tmp_qloop_109 = abs_det_jac_affine_GRAY * tmp_qloop_43;
                  const walberla::float64 tmp_qloop_110 = jac_affine_inv_1_0_GRAY * tmp_qloop_40 + tmp_qloop_107 * -2.0;
                  const walberla::float64 q_tmp_0_0 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_9 + ( tmp_qloop_5 * tmp_qloop_5 ) * 2.0 ) -
                        tmp_qloop_6 * ( tmp_qloop_3 * 0.66666666666666663 + tmp_qloop_4 * 0.66666666666666663 ) );
                  const walberla::float64 q_tmp_0_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_25 + tmp_qloop_22 * tmp_qloop_24 ) -
                                       tmp_qloop_22 * tmp_qloop_23 );
                  const walberla::float64 q_tmp_0_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_28 + tmp_qloop_24 * tmp_qloop_27 ) -
                                       tmp_qloop_23 * tmp_qloop_27 );
                  const walberla::float64 q_tmp_0_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_39 + tmp_qloop_24 * tmp_qloop_36 ) -
                                       tmp_qloop_33 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_43 + tmp_qloop_24 * tmp_qloop_42 ) -
                                       tmp_qloop_41 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_47 + tmp_qloop_24 * tmp_qloop_46 ) -
                                       tmp_qloop_45 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_6 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_6 - tmp_qloop_48 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_7 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_49 - tmp_qloop_23 * tmp_qloop_25 );
                  const walberla::float64 q_tmp_0_8 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_50 - tmp_qloop_23 * tmp_qloop_28 );
                  const walberla::float64 q_tmp_0_9 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_54 - tmp_qloop_53 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_10 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_56 - tmp_qloop_55 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_11 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_58 - tmp_qloop_57 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_1_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_59 * tmp_qloop_63 + tmp_qloop_64 * tmp_qloop_65 ) -
                                       tmp_qloop_59 * tmp_qloop_62 );
                  const walberla::float64 q_tmp_1_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_28 * tmp_qloop_67 + tmp_qloop_66 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_66 );
                  const walberla::float64 q_tmp_1_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_67 ) -
                                       tmp_qloop_33 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_67 ) -
                                       tmp_qloop_41 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_67 ) -
                                       tmp_qloop_45 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_7 =
                      jac_affine_inv_0_0_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_60 * tmp_qloop_71;
                  const walberla::float64 q_tmp_1_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * jac_affine_inv_1_0_GRAY *
                                           tmp_qloop_21 * tmp_qloop_26 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_28 * tmp_qloop_61 );
                  const walberla::float64 q_tmp_1_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_1_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_1_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_2_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_73 * tmp_qloop_76 + tmp_qloop_77 * tmp_qloop_78 ) -
                                       tmp_qloop_73 * tmp_qloop_75 );
                  const walberla::float64 q_tmp_2_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_79 ) -
                                       tmp_qloop_33 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_79 ) -
                                       tmp_qloop_41 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_79 ) -
                                       tmp_qloop_45 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_27 * tmp_qloop_61 + tmp_qloop_49 * tmp_qloop_79 );
                  const walberla::float64 q_tmp_2_8 =
                      jac_affine_inv_1_0_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_71 * tmp_qloop_74;
                  const walberla::float64 q_tmp_2_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_2_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_2_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_27 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_3_3 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( ( tmp_qloop_36 * tmp_qloop_36 ) * 2.0 + tmp_qloop_39 * tmp_qloop_85 ) -
                        tmp_qloop_33 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_42 * tmp_qloop_68 + tmp_qloop_43 * tmp_qloop_85 ) -
                                       tmp_qloop_41 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_46 * tmp_qloop_68 + tmp_qloop_47 * tmp_qloop_85 ) -
                                       tmp_qloop_45 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * tmp_qloop_5 * tmp_qloop_84 * 2.0 - tmp_qloop_48 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_86 + tmp_qloop_49 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_86 + tmp_qloop_50 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_54 + tmp_qloop_54 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_42 * tmp_qloop_84 * 2.0 - tmp_qloop_54 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_3_11 =
                      tmp_qloop_20 * ( -tmp_qloop_54 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_4_4 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( ( tmp_qloop_42 * tmp_qloop_42 ) * 2.0 + tmp_qloop_43 * tmp_qloop_88 ) -
                        tmp_qloop_41 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_46 * tmp_qloop_69 + tmp_qloop_47 * tmp_qloop_88 ) -
                                       tmp_qloop_45 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * tmp_qloop_5 * tmp_qloop_87 * 2.0 - tmp_qloop_48 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_89 + tmp_qloop_49 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_89 + tmp_qloop_50 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_56 + tmp_qloop_54 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_42 * tmp_qloop_87 * 2.0 - tmp_qloop_55 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_11 =
                      tmp_qloop_20 * ( -tmp_qloop_56 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_5_5 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( ( tmp_qloop_46 * tmp_qloop_46 ) * 2.0 + tmp_qloop_47 * tmp_qloop_91 ) -
                        tmp_qloop_45 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * tmp_qloop_5 * tmp_qloop_90 * 2.0 - tmp_qloop_48 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_92 + tmp_qloop_49 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_92 + tmp_qloop_50 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_36 * tmp_qloop_90 * 2.0 - tmp_qloop_53 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_42 * tmp_qloop_90 * 2.0 - tmp_qloop_55 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_46 * tmp_qloop_90 * 2.0 - tmp_qloop_57 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_6_6 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_5 * tmp_qloop_94 + ( tmp_qloop_9 * tmp_qloop_9 ) * 2.0 ) -
                        tmp_qloop_48 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_94 + tmp_qloop_25 * tmp_qloop_96 ) -
                                       tmp_qloop_25 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_94 + tmp_qloop_28 * tmp_qloop_96 ) -
                                       tmp_qloop_28 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_36 * tmp_qloop_94 + tmp_qloop_39 * tmp_qloop_96 ) -
                                       tmp_qloop_53 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_42 * tmp_qloop_94 + tmp_qloop_43 * tmp_qloop_96 ) -
                                       tmp_qloop_55 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_46 * tmp_qloop_94 + tmp_qloop_47 * tmp_qloop_96 ) -
                                       tmp_qloop_57 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_7_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_59 * tmp_qloop_65 + tmp_qloop_63 * tmp_qloop_64 ) -
                                       tmp_qloop_62 * tmp_qloop_64 );
                  const walberla::float64 q_tmp_7_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_98 + tmp_qloop_97 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_97 );
                  const walberla::float64 q_tmp_7_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_25 + tmp_qloop_36 * tmp_qloop_98 ) -
                                       tmp_qloop_53 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_101 * tmp_qloop_25 + tmp_qloop_42 * tmp_qloop_98 ) -
                                       tmp_qloop_55 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_102 * tmp_qloop_25 + tmp_qloop_46 * tmp_qloop_98 ) -
                                       tmp_qloop_25 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_8_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_73 * tmp_qloop_78 + tmp_qloop_76 * tmp_qloop_77 ) -
                                       tmp_qloop_75 * tmp_qloop_77 );
                  const walberla::float64 q_tmp_8_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_36 ) -
                                       tmp_qloop_103 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_8_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_101 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_42 ) -
                                       tmp_qloop_103 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_8_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_102 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_46 ) -
                                       tmp_qloop_28 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_9_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_108 * tmp_qloop_36 + ( tmp_qloop_39 * tmp_qloop_39 ) * 2.0 ) -
                        tmp_qloop_105 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_9_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_43 + tmp_qloop_108 * tmp_qloop_42 ) -
                                       tmp_qloop_105 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_9_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_47 + tmp_qloop_108 * tmp_qloop_46 ) -
                                       tmp_qloop_105 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_10_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_110 * tmp_qloop_42 + ( tmp_qloop_43 * tmp_qloop_43 ) * 2.0 ) -
                        tmp_qloop_109 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_10_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_101 * tmp_qloop_47 + tmp_qloop_110 * tmp_qloop_46 ) -
                                       tmp_qloop_109 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_11_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY *
                            ( tmp_qloop_46 * ( jac_affine_inv_0_0_GRAY * tmp_qloop_44 * 0.5 - tmp_qloop_106 ) * 2.0 +
                              ( tmp_qloop_47 * tmp_qloop_47 ) * 2.0 ) -
                        tmp_qloop_47 * tmp_qloop_72 );
                  q_acc_0_0   = q_acc_0_0 + q_tmp_0_0;
                  q_acc_0_1   = q_acc_0_1 + q_tmp_0_1;
                  q_acc_0_2   = q_acc_0_2 + q_tmp_0_2;
                  q_acc_0_3   = q_acc_0_3 + q_tmp_0_3;
                  q_acc_0_4   = q_acc_0_4 + q_tmp_0_4;
                  q_acc_0_5   = q_acc_0_5 + q_tmp_0_5;
                  q_acc_0_6   = q_acc_0_6 + q_tmp_0_6;
                  q_acc_0_7   = q_acc_0_7 + q_tmp_0_7;
                  q_acc_0_8   = q_acc_0_8 + q_tmp_0_8;
                  q_acc_0_9   = q_acc_0_9 + q_tmp_0_9;
                  q_acc_0_10  = q_acc_0_10 + q_tmp_0_10;
                  q_acc_0_11  = q_acc_0_11 + q_tmp_0_11;
                  q_acc_1_1   = q_acc_1_1 + q_tmp_1_1;
                  q_acc_1_2   = q_acc_1_2 + q_tmp_1_2;
                  q_acc_1_3   = q_acc_1_3 + q_tmp_1_3;
                  q_acc_1_4   = q_acc_1_4 + q_tmp_1_4;
                  q_acc_1_5   = q_acc_1_5 + q_tmp_1_5;
                  q_acc_1_6   = q_acc_1_6 + q_tmp_1_6;
                  q_acc_1_7   = q_acc_1_7 + q_tmp_1_7;
                  q_acc_1_8   = q_acc_1_8 + q_tmp_1_8;
                  q_acc_1_9   = q_acc_1_9 + q_tmp_1_9;
                  q_acc_1_10  = q_acc_1_10 + q_tmp_1_10;
                  q_acc_1_11  = q_acc_1_11 + q_tmp_1_11;
                  q_acc_2_2   = q_acc_2_2 + q_tmp_2_2;
                  q_acc_2_3   = q_acc_2_3 + q_tmp_2_3;
                  q_acc_2_4   = q_acc_2_4 + q_tmp_2_4;
                  q_acc_2_5   = q_acc_2_5 + q_tmp_2_5;
                  q_acc_2_6   = q_acc_2_6 + q_tmp_2_6;
                  q_acc_2_7   = q_acc_2_7 + q_tmp_2_7;
                  q_acc_2_8   = q_acc_2_8 + q_tmp_2_8;
                  q_acc_2_9   = q_acc_2_9 + q_tmp_2_9;
                  q_acc_2_10  = q_acc_2_10 + q_tmp_2_10;
                  q_acc_2_11  = q_acc_2_11 + q_tmp_2_11;
                  q_acc_3_3   = q_acc_3_3 + q_tmp_3_3;
                  q_acc_3_4   = q_acc_3_4 + q_tmp_3_4;
                  q_acc_3_5   = q_acc_3_5 + q_tmp_3_5;
                  q_acc_3_6   = q_acc_3_6 + q_tmp_3_6;
                  q_acc_3_7   = q_acc_3_7 + q_tmp_3_7;
                  q_acc_3_8   = q_acc_3_8 + q_tmp_3_8;
                  q_acc_3_9   = q_acc_3_9 + q_tmp_3_9;
                  q_acc_3_10  = q_acc_3_10 + q_tmp_3_10;
                  q_acc_3_11  = q_acc_3_11 + q_tmp_3_11;
                  q_acc_4_4   = q_acc_4_4 + q_tmp_4_4;
                  q_acc_4_5   = q_acc_4_5 + q_tmp_4_5;
                  q_acc_4_6   = q_acc_4_6 + q_tmp_4_6;
                  q_acc_4_7   = q_acc_4_7 + q_tmp_4_7;
                  q_acc_4_8   = q_acc_4_8 + q_tmp_4_8;
                  q_acc_4_9   = q_acc_4_9 + q_tmp_4_9;
                  q_acc_4_10  = q_acc_4_10 + q_tmp_4_10;
                  q_acc_4_11  = q_acc_4_11 + q_tmp_4_11;
                  q_acc_5_5   = q_acc_5_5 + q_tmp_5_5;
                  q_acc_5_6   = q_acc_5_6 + q_tmp_5_6;
                  q_acc_5_7   = q_acc_5_7 + q_tmp_5_7;
                  q_acc_5_8   = q_acc_5_8 + q_tmp_5_8;
                  q_acc_5_9   = q_acc_5_9 + q_tmp_5_9;
                  q_acc_5_10  = q_acc_5_10 + q_tmp_5_10;
                  q_acc_5_11  = q_acc_5_11 + q_tmp_5_11;
                  q_acc_6_6   = q_acc_6_6 + q_tmp_6_6;
                  q_acc_6_7   = q_acc_6_7 + q_tmp_6_7;
                  q_acc_6_8   = q_acc_6_8 + q_tmp_6_8;
                  q_acc_6_9   = q_acc_6_9 + q_tmp_6_9;
                  q_acc_6_10  = q_acc_6_10 + q_tmp_6_10;
                  q_acc_6_11  = q_acc_6_11 + q_tmp_6_11;
                  q_acc_7_7   = q_acc_7_7 + q_tmp_7_7;
                  q_acc_7_8   = q_acc_7_8 + q_tmp_7_8;
                  q_acc_7_9   = q_acc_7_9 + q_tmp_7_9;
                  q_acc_7_10  = q_acc_7_10 + q_tmp_7_10;
                  q_acc_7_11  = q_acc_7_11 + q_tmp_7_11;
                  q_acc_8_8   = q_acc_8_8 + q_tmp_8_8;
                  q_acc_8_9   = q_acc_8_9 + q_tmp_8_9;
                  q_acc_8_10  = q_acc_8_10 + q_tmp_8_10;
                  q_acc_8_11  = q_acc_8_11 + q_tmp_8_11;
                  q_acc_9_9   = q_acc_9_9 + q_tmp_9_9;
                  q_acc_9_10  = q_acc_9_10 + q_tmp_9_10;
                  q_acc_9_11  = q_acc_9_11 + q_tmp_9_11;
                  q_acc_10_10 = q_acc_10_10 + q_tmp_10_10;
                  q_acc_10_11 = q_acc_10_11 + q_tmp_10_11;
                  q_acc_11_11 = q_acc_11_11 + q_tmp_11_11;
               }
               const walberla::float64 elMatVec_0 = q_acc_0_0 * src_dof_0 + q_acc_0_1 * src_dof_1 + q_acc_0_10 * src_dof_10 +
                                                    q_acc_0_11 * src_dof_11 + q_acc_0_2 * src_dof_2 + q_acc_0_3 * src_dof_3 +
                                                    q_acc_0_4 * src_dof_4 + q_acc_0_5 * src_dof_5 + q_acc_0_6 * src_dof_6 +
                                                    q_acc_0_7 * src_dof_7 + q_acc_0_8 * src_dof_8 + q_acc_0_9 * src_dof_9;
               const walberla::float64 elMatVec_1 = q_acc_0_1 * src_dof_0 + q_acc_1_1 * src_dof_1 + q_acc_1_10 * src_dof_10 +
                                                    q_acc_1_11 * src_dof_11 + q_acc_1_2 * src_dof_2 + q_acc_1_3 * src_dof_3 +
                                                    q_acc_1_4 * src_dof_4 + q_acc_1_5 * src_dof_5 + q_acc_1_6 * src_dof_6 +
                                                    q_acc_1_7 * src_dof_7 + q_acc_1_8 * src_dof_8 + q_acc_1_9 * src_dof_9;
               const walberla::float64 elMatVec_2 = q_acc_0_2 * src_dof_0 + q_acc_1_2 * src_dof_1 + q_acc_2_10 * src_dof_10 +
                                                    q_acc_2_11 * src_dof_11 + q_acc_2_2 * src_dof_2 + q_acc_2_3 * src_dof_3 +
                                                    q_acc_2_4 * src_dof_4 + q_acc_2_5 * src_dof_5 + q_acc_2_6 * src_dof_6 +
                                                    q_acc_2_7 * src_dof_7 + q_acc_2_8 * src_dof_8 + q_acc_2_9 * src_dof_9;
               const walberla::float64 elMatVec_3 = q_acc_0_3 * src_dof_0 + q_acc_1_3 * src_dof_1 + q_acc_2_3 * src_dof_2 +
                                                    q_acc_3_10 * src_dof_10 + q_acc_3_11 * src_dof_11 + q_acc_3_3 * src_dof_3 +
                                                    q_acc_3_4 * src_dof_4 + q_acc_3_5 * src_dof_5 + q_acc_3_6 * src_dof_6 +
                                                    q_acc_3_7 * src_dof_7 + q_acc_3_8 * src_dof_8 + q_acc_3_9 * src_dof_9;
               const walberla::float64 elMatVec_4 = q_acc_0_4 * src_dof_0 + q_acc_1_4 * src_dof_1 + q_acc_2_4 * src_dof_2 +
                                                    q_acc_3_4 * src_dof_3 + q_acc_4_10 * src_dof_10 + q_acc_4_11 * src_dof_11 +
                                                    q_acc_4_4 * src_dof_4 + q_acc_4_5 * src_dof_5 + q_acc_4_6 * src_dof_6 +
                                                    q_acc_4_7 * src_dof_7 + q_acc_4_8 * src_dof_8 + q_acc_4_9 * src_dof_9;
               const walberla::float64 elMatVec_5 = q_acc_0_5 * src_dof_0 + q_acc_1_5 * src_dof_1 + q_acc_2_5 * src_dof_2 +
                                                    q_acc_3_5 * src_dof_3 + q_acc_4_5 * src_dof_4 + q_acc_5_10 * src_dof_10 +
                                                    q_acc_5_11 * src_dof_11 + q_acc_5_5 * src_dof_5 + q_acc_5_6 * src_dof_6 +
                                                    q_acc_5_7 * src_dof_7 + q_acc_5_8 * src_dof_8 + q_acc_5_9 * src_dof_9;
               const walberla::float64 elMatVec_6 = q_acc_0_6 * src_dof_0 + q_acc_1_6 * src_dof_1 + q_acc_2_6 * src_dof_2 +
                                                    q_acc_3_6 * src_dof_3 + q_acc_4_6 * src_dof_4 + q_acc_5_6 * src_dof_5 +
                                                    q_acc_6_10 * src_dof_10 + q_acc_6_11 * src_dof_11 + q_acc_6_6 * src_dof_6 +
                                                    q_acc_6_7 * src_dof_7 + q_acc_6_8 * src_dof_8 + q_acc_6_9 * src_dof_9;
               const walberla::float64 elMatVec_7 = q_acc_0_7 * src_dof_0 + q_acc_1_7 * src_dof_1 + q_acc_2_7 * src_dof_2 +
                                                    q_acc_3_7 * src_dof_3 + q_acc_4_7 * src_dof_4 + q_acc_5_7 * src_dof_5 +
                                                    q_acc_6_7 * src_dof_6 + q_acc_7_10 * src_dof_10 + q_acc_7_11 * src_dof_11 +
                                                    q_acc_7_7 * src_dof_7 + q_acc_7_8 * src_dof_8 + q_acc_7_9 * src_dof_9;
               const walberla::float64 elMatVec_8 = q_acc_0_8 * src_dof_0 + q_acc_1_8 * src_dof_1 + q_acc_2_8 * src_dof_2 +
                                                    q_acc_3_8 * src_dof_3 + q_acc_4_8 * src_dof_4 + q_acc_5_8 * src_dof_5 +
                                                    q_acc_6_8 * src_dof_6 + q_acc_7_8 * src_dof_7 + q_acc_8_10 * src_dof_10 +
                                                    q_acc_8_11 * src_dof_11 + q_acc_8_8 * src_dof_8 + q_acc_8_9 * src_dof_9;
               const walberla::float64 elMatVec_9 = q_acc_0_9 * src_dof_0 + q_acc_1_9 * src_dof_1 + q_acc_2_9 * src_dof_2 +
                                                    q_acc_3_9 * src_dof_3 + q_acc_4_9 * src_dof_4 + q_acc_5_9 * src_dof_5 +
                                                    q_acc_6_9 * src_dof_6 + q_acc_7_9 * src_dof_7 + q_acc_8_9 * src_dof_8 +
                                                    q_acc_9_10 * src_dof_10 + q_acc_9_11 * src_dof_11 + q_acc_9_9 * src_dof_9;
               const walberla::float64 elMatVec_10 =
                   q_acc_0_10 * src_dof_0 + q_acc_10_10 * src_dof_10 + q_acc_10_11 * src_dof_11 + q_acc_1_10 * src_dof_1 +
                   q_acc_2_10 * src_dof_2 + q_acc_3_10 * src_dof_3 + q_acc_4_10 * src_dof_4 + q_acc_5_10 * src_dof_5 +
                   q_acc_6_10 * src_dof_6 + q_acc_7_10 * src_dof_7 + q_acc_8_10 * src_dof_8 + q_acc_9_10 * src_dof_9;
               const walberla::float64 elMatVec_11 =
                   q_acc_0_11 * src_dof_0 + q_acc_10_11 * src_dof_10 + q_acc_11_11 * src_dof_11 + q_acc_1_11 * src_dof_1 +
                   q_acc_2_11 * src_dof_2 + q_acc_3_11 * src_dof_3 + q_acc_4_11 * src_dof_4 + q_acc_5_11 * src_dof_5 +
                   q_acc_6_11 * src_dof_6 + q_acc_7_11 * src_dof_7 + q_acc_8_11 * src_dof_8 + q_acc_9_11 * src_dof_9;
               _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatVec_0 +
                   _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                  1] = elMatVec_1 + _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                  ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_2 + _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatVec_3 +
                   _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatVec_4 +
                   _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatVec_5 +
                   _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatVec_6 +
                   _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                  1] = elMatVec_7 + _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                  ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_8 + _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatVec_9 +
                   _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatVec_10 +
                   _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatVec_11 +
                   _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
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
             jac_affine_0_0_BLUE * jac_affine_1_1_BLUE - jac_affine_0_1_BLUE * jac_affine_1_0_BLUE;
         const walberla::float64 tmp_coords_jac_6_BLUE   = 1.0 / ( tmp_coords_jac_5_BLUE );
         const walberla::float64 jac_affine_inv_0_0_BLUE = jac_affine_1_1_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_0_1_BLUE = -jac_affine_0_1_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_1_0_BLUE = -jac_affine_1_0_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_1_1_BLUE = jac_affine_0_0_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 abs_det_jac_affine_BLUE = abs( tmp_coords_jac_5_BLUE );
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
               const walberla::float64 src_dof_0 = _data_src_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_1 = _data_src_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_2 = _data_src_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_3 = _data_src_edge_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                                    ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_4 =
                   _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_5 =
                   _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_6 = _data_src_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_7 = _data_src_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_8 = _data_src_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                      ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_9 = _data_src_edge_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                                    ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 src_dof_10 =
                   _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1];
               const walberla::float64 src_dof_11 =
                   _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
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
               const walberla::float64 uy_dof_1    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_2    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               walberla::float64       q_acc_0_0   = 0.0;
               walberla::float64       q_acc_0_1   = 0.0;
               walberla::float64       q_acc_0_2   = 0.0;
               walberla::float64       q_acc_0_3   = 0.0;
               walberla::float64       q_acc_0_4   = 0.0;
               walberla::float64       q_acc_0_5   = 0.0;
               walberla::float64       q_acc_0_6   = 0.0;
               walberla::float64       q_acc_0_7   = 0.0;
               walberla::float64       q_acc_0_8   = 0.0;
               walberla::float64       q_acc_0_9   = 0.0;
               walberla::float64       q_acc_0_10  = 0.0;
               walberla::float64       q_acc_0_11  = 0.0;
               walberla::float64       q_acc_1_1   = 0.0;
               walberla::float64       q_acc_1_2   = 0.0;
               walberla::float64       q_acc_1_3   = 0.0;
               walberla::float64       q_acc_1_4   = 0.0;
               walberla::float64       q_acc_1_5   = 0.0;
               walberla::float64       q_acc_1_6   = 0.0;
               walberla::float64       q_acc_1_7   = 0.0;
               walberla::float64       q_acc_1_8   = 0.0;
               walberla::float64       q_acc_1_9   = 0.0;
               walberla::float64       q_acc_1_10  = 0.0;
               walberla::float64       q_acc_1_11  = 0.0;
               walberla::float64       q_acc_2_2   = 0.0;
               walberla::float64       q_acc_2_3   = 0.0;
               walberla::float64       q_acc_2_4   = 0.0;
               walberla::float64       q_acc_2_5   = 0.0;
               walberla::float64       q_acc_2_6   = 0.0;
               walberla::float64       q_acc_2_7   = 0.0;
               walberla::float64       q_acc_2_8   = 0.0;
               walberla::float64       q_acc_2_9   = 0.0;
               walberla::float64       q_acc_2_10  = 0.0;
               walberla::float64       q_acc_2_11  = 0.0;
               walberla::float64       q_acc_3_3   = 0.0;
               walberla::float64       q_acc_3_4   = 0.0;
               walberla::float64       q_acc_3_5   = 0.0;
               walberla::float64       q_acc_3_6   = 0.0;
               walberla::float64       q_acc_3_7   = 0.0;
               walberla::float64       q_acc_3_8   = 0.0;
               walberla::float64       q_acc_3_9   = 0.0;
               walberla::float64       q_acc_3_10  = 0.0;
               walberla::float64       q_acc_3_11  = 0.0;
               walberla::float64       q_acc_4_4   = 0.0;
               walberla::float64       q_acc_4_5   = 0.0;
               walberla::float64       q_acc_4_6   = 0.0;
               walberla::float64       q_acc_4_7   = 0.0;
               walberla::float64       q_acc_4_8   = 0.0;
               walberla::float64       q_acc_4_9   = 0.0;
               walberla::float64       q_acc_4_10  = 0.0;
               walberla::float64       q_acc_4_11  = 0.0;
               walberla::float64       q_acc_5_5   = 0.0;
               walberla::float64       q_acc_5_6   = 0.0;
               walberla::float64       q_acc_5_7   = 0.0;
               walberla::float64       q_acc_5_8   = 0.0;
               walberla::float64       q_acc_5_9   = 0.0;
               walberla::float64       q_acc_5_10  = 0.0;
               walberla::float64       q_acc_5_11  = 0.0;
               walberla::float64       q_acc_6_6   = 0.0;
               walberla::float64       q_acc_6_7   = 0.0;
               walberla::float64       q_acc_6_8   = 0.0;
               walberla::float64       q_acc_6_9   = 0.0;
               walberla::float64       q_acc_6_10  = 0.0;
               walberla::float64       q_acc_6_11  = 0.0;
               walberla::float64       q_acc_7_7   = 0.0;
               walberla::float64       q_acc_7_8   = 0.0;
               walberla::float64       q_acc_7_9   = 0.0;
               walberla::float64       q_acc_7_10  = 0.0;
               walberla::float64       q_acc_7_11  = 0.0;
               walberla::float64       q_acc_8_8   = 0.0;
               walberla::float64       q_acc_8_9   = 0.0;
               walberla::float64       q_acc_8_10  = 0.0;
               walberla::float64       q_acc_8_11  = 0.0;
               walberla::float64       q_acc_9_9   = 0.0;
               walberla::float64       q_acc_9_10  = 0.0;
               walberla::float64       q_acc_9_11  = 0.0;
               walberla::float64       q_acc_10_10 = 0.0;
               walberla::float64       q_acc_10_11 = 0.0;
               walberla::float64       q_acc_11_11 = 0.0;
               for ( int64_t q = 0; q < 3; q += 1 )
               {
                  const walberla::float64 tmp_qloop_0  = 4.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_1  = 4.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_2  = tmp_qloop_0 + tmp_qloop_1 - 3.0;
                  const walberla::float64 tmp_qloop_3  = jac_affine_inv_0_0_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_4  = jac_affine_inv_1_0_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_5  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_6  = abs_det_jac_affine_BLUE * tmp_qloop_5;
                  const walberla::float64 tmp_qloop_7  = jac_affine_inv_0_1_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_8  = jac_affine_inv_1_1_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_9  = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_10 = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_11 = -ux_dof_0 + ux_dof_1;
                  const walberla::float64 tmp_qloop_12 = -ux_dof_0 + ux_dof_2;
                  const walberla::float64 tmp_qloop_13 = -uy_dof_0 + uy_dof_1;
                  const walberla::float64 tmp_qloop_14 = -uy_dof_0 + uy_dof_2;
                  const walberla::float64 tmp_qloop_15 = jac_affine_inv_0_0_BLUE * tmp_qloop_13;
                  const walberla::float64 tmp_qloop_16 = jac_affine_inv_0_1_BLUE * tmp_qloop_11;
                  const walberla::float64 tmp_qloop_17 = jac_affine_inv_1_0_BLUE * tmp_qloop_14;
                  const walberla::float64 tmp_qloop_18 = jac_affine_inv_1_1_BLUE * tmp_qloop_12;
                  const walberla::float64 tmp_qloop_19 =
                      1.0 /
                      ( 1.0 /
                            ( mu_star +
                              sigma_y * 1.0 /
                                  ( pow( ( ( jac_affine_inv_0_0_BLUE * tmp_qloop_11 + jac_affine_inv_1_0_BLUE * tmp_qloop_12 ) *
                                           ( jac_affine_inv_0_0_BLUE * tmp_qloop_11 + jac_affine_inv_1_0_BLUE * tmp_qloop_12 ) ) +
                                             ( ( jac_affine_inv_0_1_BLUE * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_BLUE * tmp_qloop_14 ) *
                                               ( jac_affine_inv_0_1_BLUE * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_BLUE * tmp_qloop_14 ) ) +
                                             ( tmp_qloop_15 + tmp_qloop_16 + tmp_qloop_17 + tmp_qloop_18 ) *
                                                 ( tmp_qloop_15 * 0.5 + tmp_qloop_16 * 0.5 + tmp_qloop_17 * 0.5 +
                                                   tmp_qloop_18 * 0.5 ),
                                         0.50000000000000000 ) +
                                    9.9999999999999995e-21 ) ) *
                            1.0 +
                        1.0 /
                            ( mu_lin_dof_0 * ( 1.0 - _data_q_p_0[q] - _data_q_p_1[q] ) + mu_lin_dof_1 * _data_q_p_0[q] +
                              mu_lin_dof_2 * _data_q_p_1[q] ) *
                            1.0 ) *
                      _data_q_w[q];
                  const walberla::float64 tmp_qloop_20 = tmp_qloop_19 * 2.0;
                  const walberla::float64 tmp_qloop_21 = tmp_qloop_0 - 1.0;
                  const walberla::float64 tmp_qloop_22 = jac_affine_inv_0_0_BLUE * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_23 = tmp_qloop_6 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_24 = tmp_qloop_5 * 2.0;
                  const walberla::float64 tmp_qloop_25 = jac_affine_inv_0_1_BLUE * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_26 = tmp_qloop_1 - 1.0;
                  const walberla::float64 tmp_qloop_27 = jac_affine_inv_1_0_BLUE * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_28 = jac_affine_inv_1_1_BLUE * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_29 = 2.6666666666666665 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_30 = jac_affine_inv_1_0_BLUE * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_31 = 2.6666666666666665 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_32 = jac_affine_inv_0_0_BLUE * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_33 = tmp_qloop_30 + tmp_qloop_32;
                  const walberla::float64 tmp_qloop_34 = jac_affine_inv_1_0_BLUE * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_35 = jac_affine_inv_0_0_BLUE * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_36 = tmp_qloop_34 + tmp_qloop_35;
                  const walberla::float64 tmp_qloop_37 = jac_affine_inv_1_1_BLUE * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_38 = jac_affine_inv_0_1_BLUE * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_39 = tmp_qloop_37 + tmp_qloop_38;
                  const walberla::float64 tmp_qloop_40 = -tmp_qloop_0 - 8.0 * _data_q_p_1[q] + 4.0;
                  const walberla::float64 tmp_qloop_41 =
                      jac_affine_inv_1_0_BLUE * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_32;
                  const walberla::float64 tmp_qloop_42 = jac_affine_inv_1_0_BLUE * tmp_qloop_40 - tmp_qloop_35;
                  const walberla::float64 tmp_qloop_43 = jac_affine_inv_1_1_BLUE * tmp_qloop_40 - tmp_qloop_38;
                  const walberla::float64 tmp_qloop_44 = -tmp_qloop_1 - 8.0 * _data_q_p_0[q] + 4.0;
                  const walberla::float64 tmp_qloop_45 =
                      jac_affine_inv_0_0_BLUE * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_30;
                  const walberla::float64 tmp_qloop_46 = jac_affine_inv_0_0_BLUE * tmp_qloop_44 - tmp_qloop_34;
                  const walberla::float64 tmp_qloop_47 = jac_affine_inv_0_1_BLUE * tmp_qloop_44 - tmp_qloop_37;
                  const walberla::float64 tmp_qloop_48 = tmp_qloop_7 * 0.66666666666666663 + tmp_qloop_8 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_49 = abs_det_jac_affine_BLUE * tmp_qloop_22;
                  const walberla::float64 tmp_qloop_50 = abs_det_jac_affine_BLUE * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_51 = jac_affine_inv_1_1_BLUE * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_52 = jac_affine_inv_0_1_BLUE * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_53 = tmp_qloop_51 + tmp_qloop_52;
                  const walberla::float64 tmp_qloop_54 = abs_det_jac_affine_BLUE * tmp_qloop_36;
                  const walberla::float64 tmp_qloop_55 =
                      jac_affine_inv_1_1_BLUE * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_52;
                  const walberla::float64 tmp_qloop_56 = abs_det_jac_affine_BLUE * tmp_qloop_42;
                  const walberla::float64 tmp_qloop_57 =
                      jac_affine_inv_0_1_BLUE * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_51;
                  const walberla::float64 tmp_qloop_58  = abs_det_jac_affine_BLUE * tmp_qloop_46;
                  const walberla::float64 tmp_qloop_59  = ( jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE );
                  const walberla::float64 tmp_qloop_60  = ( tmp_qloop_21 * tmp_qloop_21 );
                  const walberla::float64 tmp_qloop_61  = abs_det_jac_affine_BLUE * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_62  = tmp_qloop_60 * tmp_qloop_61;
                  const walberla::float64 tmp_qloop_63  = tmp_qloop_60 * 2.0;
                  const walberla::float64 tmp_qloop_64  = ( jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE );
                  const walberla::float64 tmp_qloop_65  = tmp_qloop_60 * 1.0;
                  const walberla::float64 tmp_qloop_66  = tmp_qloop_22 * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_67  = tmp_qloop_25 * 1.0;
                  const walberla::float64 tmp_qloop_68  = tmp_qloop_36 * 2.0;
                  const walberla::float64 tmp_qloop_69  = tmp_qloop_42 * 2.0;
                  const walberla::float64 tmp_qloop_70  = tmp_qloop_46 * 2.0;
                  const walberla::float64 tmp_qloop_71  = abs_det_jac_affine_BLUE * tmp_qloop_19 * 0.66666666666666674;
                  const walberla::float64 tmp_qloop_72  = abs_det_jac_affine_BLUE * tmp_qloop_57;
                  const walberla::float64 tmp_qloop_73  = ( jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE );
                  const walberla::float64 tmp_qloop_74  = ( tmp_qloop_26 * tmp_qloop_26 );
                  const walberla::float64 tmp_qloop_75  = tmp_qloop_61 * tmp_qloop_74;
                  const walberla::float64 tmp_qloop_76  = tmp_qloop_74 * 2.0;
                  const walberla::float64 tmp_qloop_77  = ( jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE );
                  const walberla::float64 tmp_qloop_78  = tmp_qloop_74 * 1.0;
                  const walberla::float64 tmp_qloop_79  = tmp_qloop_28 * 1.0;
                  const walberla::float64 tmp_qloop_80  = 2.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_81  = jac_affine_inv_1_1_BLUE * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_82  = 2.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_83  = jac_affine_inv_0_1_BLUE * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_84  = tmp_qloop_81 + tmp_qloop_83;
                  const walberla::float64 tmp_qloop_85  = tmp_qloop_84 * 2.0;
                  const walberla::float64 tmp_qloop_86  = tmp_qloop_54 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_87  = jac_affine_inv_1_1_BLUE * tmp_qloop_40 * 0.5 - tmp_qloop_83;
                  const walberla::float64 tmp_qloop_88  = tmp_qloop_87 * 2.0;
                  const walberla::float64 tmp_qloop_89  = tmp_qloop_56 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_90  = jac_affine_inv_0_1_BLUE * tmp_qloop_44 * 0.5 - tmp_qloop_81;
                  const walberla::float64 tmp_qloop_91  = tmp_qloop_90 * 2.0;
                  const walberla::float64 tmp_qloop_92  = tmp_qloop_58 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_93  = abs_det_jac_affine_BLUE * tmp_qloop_9;
                  const walberla::float64 tmp_qloop_94  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_95  = tmp_qloop_93 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_96  = tmp_qloop_9 * 2.0;
                  const walberla::float64 tmp_qloop_97  = tmp_qloop_25 * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_98  = tmp_qloop_22 * 1.0;
                  const walberla::float64 tmp_qloop_99  = abs_det_jac_affine_BLUE * tmp_qloop_25;
                  const walberla::float64 tmp_qloop_100 = tmp_qloop_39 * 2.0;
                  const walberla::float64 tmp_qloop_101 = tmp_qloop_43 * 2.0;
                  const walberla::float64 tmp_qloop_102 = tmp_qloop_47 * 2.0;
                  const walberla::float64 tmp_qloop_103 = abs_det_jac_affine_BLUE * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_104 = tmp_qloop_27 * 1.0;
                  const walberla::float64 tmp_qloop_105 = abs_det_jac_affine_BLUE * tmp_qloop_39;
                  const walberla::float64 tmp_qloop_106 = jac_affine_inv_1_0_BLUE * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_107 = jac_affine_inv_0_0_BLUE * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_108 = tmp_qloop_106 * 2.0 + tmp_qloop_107 * 2.0;
                  const walberla::float64 tmp_qloop_109 = abs_det_jac_affine_BLUE * tmp_qloop_43;
                  const walberla::float64 tmp_qloop_110 = jac_affine_inv_1_0_BLUE * tmp_qloop_40 + tmp_qloop_107 * -2.0;
                  const walberla::float64 q_tmp_0_0 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_9 + ( tmp_qloop_5 * tmp_qloop_5 ) * 2.0 ) -
                        tmp_qloop_6 * ( tmp_qloop_3 * 0.66666666666666663 + tmp_qloop_4 * 0.66666666666666663 ) );
                  const walberla::float64 q_tmp_0_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_25 + tmp_qloop_22 * tmp_qloop_24 ) -
                                       tmp_qloop_22 * tmp_qloop_23 );
                  const walberla::float64 q_tmp_0_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_28 + tmp_qloop_24 * tmp_qloop_27 ) -
                                       tmp_qloop_23 * tmp_qloop_27 );
                  const walberla::float64 q_tmp_0_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_39 + tmp_qloop_24 * tmp_qloop_36 ) -
                                       tmp_qloop_33 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_43 + tmp_qloop_24 * tmp_qloop_42 ) -
                                       tmp_qloop_41 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_47 + tmp_qloop_24 * tmp_qloop_46 ) -
                                       tmp_qloop_45 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_6 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_6 - tmp_qloop_48 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_7 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_49 - tmp_qloop_23 * tmp_qloop_25 );
                  const walberla::float64 q_tmp_0_8 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_50 - tmp_qloop_23 * tmp_qloop_28 );
                  const walberla::float64 q_tmp_0_9 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_54 - tmp_qloop_53 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_10 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_56 - tmp_qloop_55 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_11 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_58 - tmp_qloop_57 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_1_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_59 * tmp_qloop_63 + tmp_qloop_64 * tmp_qloop_65 ) -
                                       tmp_qloop_59 * tmp_qloop_62 );
                  const walberla::float64 q_tmp_1_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_28 * tmp_qloop_67 + tmp_qloop_66 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_66 );
                  const walberla::float64 q_tmp_1_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_67 ) -
                                       tmp_qloop_33 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_67 ) -
                                       tmp_qloop_41 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_67 ) -
                                       tmp_qloop_45 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_7 =
                      jac_affine_inv_0_0_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_60 * tmp_qloop_71;
                  const walberla::float64 q_tmp_1_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * jac_affine_inv_1_0_BLUE *
                                           tmp_qloop_21 * tmp_qloop_26 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_28 * tmp_qloop_61 );
                  const walberla::float64 q_tmp_1_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_1_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_1_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_2_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_73 * tmp_qloop_76 + tmp_qloop_77 * tmp_qloop_78 ) -
                                       tmp_qloop_73 * tmp_qloop_75 );
                  const walberla::float64 q_tmp_2_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_79 ) -
                                       tmp_qloop_33 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_79 ) -
                                       tmp_qloop_41 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_79 ) -
                                       tmp_qloop_45 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_27 * tmp_qloop_61 + tmp_qloop_49 * tmp_qloop_79 );
                  const walberla::float64 q_tmp_2_8 =
                      jac_affine_inv_1_0_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_71 * tmp_qloop_74;
                  const walberla::float64 q_tmp_2_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_2_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_2_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_27 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_3_3 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( ( tmp_qloop_36 * tmp_qloop_36 ) * 2.0 + tmp_qloop_39 * tmp_qloop_85 ) -
                        tmp_qloop_33 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_42 * tmp_qloop_68 + tmp_qloop_43 * tmp_qloop_85 ) -
                                       tmp_qloop_41 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_46 * tmp_qloop_68 + tmp_qloop_47 * tmp_qloop_85 ) -
                                       tmp_qloop_45 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * tmp_qloop_5 * tmp_qloop_84 * 2.0 - tmp_qloop_48 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_86 + tmp_qloop_49 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_86 + tmp_qloop_50 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_54 + tmp_qloop_54 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_42 * tmp_qloop_84 * 2.0 - tmp_qloop_54 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_3_11 =
                      tmp_qloop_20 * ( -tmp_qloop_54 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_4_4 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( ( tmp_qloop_42 * tmp_qloop_42 ) * 2.0 + tmp_qloop_43 * tmp_qloop_88 ) -
                        tmp_qloop_41 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_46 * tmp_qloop_69 + tmp_qloop_47 * tmp_qloop_88 ) -
                                       tmp_qloop_45 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * tmp_qloop_5 * tmp_qloop_87 * 2.0 - tmp_qloop_48 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_89 + tmp_qloop_49 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_89 + tmp_qloop_50 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_56 + tmp_qloop_54 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_42 * tmp_qloop_87 * 2.0 - tmp_qloop_55 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_11 =
                      tmp_qloop_20 * ( -tmp_qloop_56 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_5_5 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( ( tmp_qloop_46 * tmp_qloop_46 ) * 2.0 + tmp_qloop_47 * tmp_qloop_91 ) -
                        tmp_qloop_45 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * tmp_qloop_5 * tmp_qloop_90 * 2.0 - tmp_qloop_48 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_92 + tmp_qloop_49 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_92 + tmp_qloop_50 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_36 * tmp_qloop_90 * 2.0 - tmp_qloop_53 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_42 * tmp_qloop_90 * 2.0 - tmp_qloop_55 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_46 * tmp_qloop_90 * 2.0 - tmp_qloop_57 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_6_6 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_5 * tmp_qloop_94 + ( tmp_qloop_9 * tmp_qloop_9 ) * 2.0 ) -
                        tmp_qloop_48 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_94 + tmp_qloop_25 * tmp_qloop_96 ) -
                                       tmp_qloop_25 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_94 + tmp_qloop_28 * tmp_qloop_96 ) -
                                       tmp_qloop_28 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_36 * tmp_qloop_94 + tmp_qloop_39 * tmp_qloop_96 ) -
                                       tmp_qloop_53 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_42 * tmp_qloop_94 + tmp_qloop_43 * tmp_qloop_96 ) -
                                       tmp_qloop_55 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_46 * tmp_qloop_94 + tmp_qloop_47 * tmp_qloop_96 ) -
                                       tmp_qloop_57 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_7_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_59 * tmp_qloop_65 + tmp_qloop_63 * tmp_qloop_64 ) -
                                       tmp_qloop_62 * tmp_qloop_64 );
                  const walberla::float64 q_tmp_7_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_98 + tmp_qloop_97 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_97 );
                  const walberla::float64 q_tmp_7_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_25 + tmp_qloop_36 * tmp_qloop_98 ) -
                                       tmp_qloop_53 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_101 * tmp_qloop_25 + tmp_qloop_42 * tmp_qloop_98 ) -
                                       tmp_qloop_55 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_102 * tmp_qloop_25 + tmp_qloop_46 * tmp_qloop_98 ) -
                                       tmp_qloop_25 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_8_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_73 * tmp_qloop_78 + tmp_qloop_76 * tmp_qloop_77 ) -
                                       tmp_qloop_75 * tmp_qloop_77 );
                  const walberla::float64 q_tmp_8_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_36 ) -
                                       tmp_qloop_103 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_8_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_101 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_42 ) -
                                       tmp_qloop_103 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_8_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_102 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_46 ) -
                                       tmp_qloop_28 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_9_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_108 * tmp_qloop_36 + ( tmp_qloop_39 * tmp_qloop_39 ) * 2.0 ) -
                        tmp_qloop_105 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_9_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_43 + tmp_qloop_108 * tmp_qloop_42 ) -
                                       tmp_qloop_105 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_9_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_47 + tmp_qloop_108 * tmp_qloop_46 ) -
                                       tmp_qloop_105 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_10_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_110 * tmp_qloop_42 + ( tmp_qloop_43 * tmp_qloop_43 ) * 2.0 ) -
                        tmp_qloop_109 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_10_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_101 * tmp_qloop_47 + tmp_qloop_110 * tmp_qloop_46 ) -
                                       tmp_qloop_109 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_11_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE *
                            ( tmp_qloop_46 * ( jac_affine_inv_0_0_BLUE * tmp_qloop_44 * 0.5 - tmp_qloop_106 ) * 2.0 +
                              ( tmp_qloop_47 * tmp_qloop_47 ) * 2.0 ) -
                        tmp_qloop_47 * tmp_qloop_72 );
                  q_acc_0_0   = q_acc_0_0 + q_tmp_0_0;
                  q_acc_0_1   = q_acc_0_1 + q_tmp_0_1;
                  q_acc_0_2   = q_acc_0_2 + q_tmp_0_2;
                  q_acc_0_3   = q_acc_0_3 + q_tmp_0_3;
                  q_acc_0_4   = q_acc_0_4 + q_tmp_0_4;
                  q_acc_0_5   = q_acc_0_5 + q_tmp_0_5;
                  q_acc_0_6   = q_acc_0_6 + q_tmp_0_6;
                  q_acc_0_7   = q_acc_0_7 + q_tmp_0_7;
                  q_acc_0_8   = q_acc_0_8 + q_tmp_0_8;
                  q_acc_0_9   = q_acc_0_9 + q_tmp_0_9;
                  q_acc_0_10  = q_acc_0_10 + q_tmp_0_10;
                  q_acc_0_11  = q_acc_0_11 + q_tmp_0_11;
                  q_acc_1_1   = q_acc_1_1 + q_tmp_1_1;
                  q_acc_1_2   = q_acc_1_2 + q_tmp_1_2;
                  q_acc_1_3   = q_acc_1_3 + q_tmp_1_3;
                  q_acc_1_4   = q_acc_1_4 + q_tmp_1_4;
                  q_acc_1_5   = q_acc_1_5 + q_tmp_1_5;
                  q_acc_1_6   = q_acc_1_6 + q_tmp_1_6;
                  q_acc_1_7   = q_acc_1_7 + q_tmp_1_7;
                  q_acc_1_8   = q_acc_1_8 + q_tmp_1_8;
                  q_acc_1_9   = q_acc_1_9 + q_tmp_1_9;
                  q_acc_1_10  = q_acc_1_10 + q_tmp_1_10;
                  q_acc_1_11  = q_acc_1_11 + q_tmp_1_11;
                  q_acc_2_2   = q_acc_2_2 + q_tmp_2_2;
                  q_acc_2_3   = q_acc_2_3 + q_tmp_2_3;
                  q_acc_2_4   = q_acc_2_4 + q_tmp_2_4;
                  q_acc_2_5   = q_acc_2_5 + q_tmp_2_5;
                  q_acc_2_6   = q_acc_2_6 + q_tmp_2_6;
                  q_acc_2_7   = q_acc_2_7 + q_tmp_2_7;
                  q_acc_2_8   = q_acc_2_8 + q_tmp_2_8;
                  q_acc_2_9   = q_acc_2_9 + q_tmp_2_9;
                  q_acc_2_10  = q_acc_2_10 + q_tmp_2_10;
                  q_acc_2_11  = q_acc_2_11 + q_tmp_2_11;
                  q_acc_3_3   = q_acc_3_3 + q_tmp_3_3;
                  q_acc_3_4   = q_acc_3_4 + q_tmp_3_4;
                  q_acc_3_5   = q_acc_3_5 + q_tmp_3_5;
                  q_acc_3_6   = q_acc_3_6 + q_tmp_3_6;
                  q_acc_3_7   = q_acc_3_7 + q_tmp_3_7;
                  q_acc_3_8   = q_acc_3_8 + q_tmp_3_8;
                  q_acc_3_9   = q_acc_3_9 + q_tmp_3_9;
                  q_acc_3_10  = q_acc_3_10 + q_tmp_3_10;
                  q_acc_3_11  = q_acc_3_11 + q_tmp_3_11;
                  q_acc_4_4   = q_acc_4_4 + q_tmp_4_4;
                  q_acc_4_5   = q_acc_4_5 + q_tmp_4_5;
                  q_acc_4_6   = q_acc_4_6 + q_tmp_4_6;
                  q_acc_4_7   = q_acc_4_7 + q_tmp_4_7;
                  q_acc_4_8   = q_acc_4_8 + q_tmp_4_8;
                  q_acc_4_9   = q_acc_4_9 + q_tmp_4_9;
                  q_acc_4_10  = q_acc_4_10 + q_tmp_4_10;
                  q_acc_4_11  = q_acc_4_11 + q_tmp_4_11;
                  q_acc_5_5   = q_acc_5_5 + q_tmp_5_5;
                  q_acc_5_6   = q_acc_5_6 + q_tmp_5_6;
                  q_acc_5_7   = q_acc_5_7 + q_tmp_5_7;
                  q_acc_5_8   = q_acc_5_8 + q_tmp_5_8;
                  q_acc_5_9   = q_acc_5_9 + q_tmp_5_9;
                  q_acc_5_10  = q_acc_5_10 + q_tmp_5_10;
                  q_acc_5_11  = q_acc_5_11 + q_tmp_5_11;
                  q_acc_6_6   = q_acc_6_6 + q_tmp_6_6;
                  q_acc_6_7   = q_acc_6_7 + q_tmp_6_7;
                  q_acc_6_8   = q_acc_6_8 + q_tmp_6_8;
                  q_acc_6_9   = q_acc_6_9 + q_tmp_6_9;
                  q_acc_6_10  = q_acc_6_10 + q_tmp_6_10;
                  q_acc_6_11  = q_acc_6_11 + q_tmp_6_11;
                  q_acc_7_7   = q_acc_7_7 + q_tmp_7_7;
                  q_acc_7_8   = q_acc_7_8 + q_tmp_7_8;
                  q_acc_7_9   = q_acc_7_9 + q_tmp_7_9;
                  q_acc_7_10  = q_acc_7_10 + q_tmp_7_10;
                  q_acc_7_11  = q_acc_7_11 + q_tmp_7_11;
                  q_acc_8_8   = q_acc_8_8 + q_tmp_8_8;
                  q_acc_8_9   = q_acc_8_9 + q_tmp_8_9;
                  q_acc_8_10  = q_acc_8_10 + q_tmp_8_10;
                  q_acc_8_11  = q_acc_8_11 + q_tmp_8_11;
                  q_acc_9_9   = q_acc_9_9 + q_tmp_9_9;
                  q_acc_9_10  = q_acc_9_10 + q_tmp_9_10;
                  q_acc_9_11  = q_acc_9_11 + q_tmp_9_11;
                  q_acc_10_10 = q_acc_10_10 + q_tmp_10_10;
                  q_acc_10_11 = q_acc_10_11 + q_tmp_10_11;
                  q_acc_11_11 = q_acc_11_11 + q_tmp_11_11;
               }
               const walberla::float64 elMatVec_0 = q_acc_0_0 * src_dof_0 + q_acc_0_1 * src_dof_1 + q_acc_0_10 * src_dof_10 +
                                                    q_acc_0_11 * src_dof_11 + q_acc_0_2 * src_dof_2 + q_acc_0_3 * src_dof_3 +
                                                    q_acc_0_4 * src_dof_4 + q_acc_0_5 * src_dof_5 + q_acc_0_6 * src_dof_6 +
                                                    q_acc_0_7 * src_dof_7 + q_acc_0_8 * src_dof_8 + q_acc_0_9 * src_dof_9;
               const walberla::float64 elMatVec_1 = q_acc_0_1 * src_dof_0 + q_acc_1_1 * src_dof_1 + q_acc_1_10 * src_dof_10 +
                                                    q_acc_1_11 * src_dof_11 + q_acc_1_2 * src_dof_2 + q_acc_1_3 * src_dof_3 +
                                                    q_acc_1_4 * src_dof_4 + q_acc_1_5 * src_dof_5 + q_acc_1_6 * src_dof_6 +
                                                    q_acc_1_7 * src_dof_7 + q_acc_1_8 * src_dof_8 + q_acc_1_9 * src_dof_9;
               const walberla::float64 elMatVec_2 = q_acc_0_2 * src_dof_0 + q_acc_1_2 * src_dof_1 + q_acc_2_10 * src_dof_10 +
                                                    q_acc_2_11 * src_dof_11 + q_acc_2_2 * src_dof_2 + q_acc_2_3 * src_dof_3 +
                                                    q_acc_2_4 * src_dof_4 + q_acc_2_5 * src_dof_5 + q_acc_2_6 * src_dof_6 +
                                                    q_acc_2_7 * src_dof_7 + q_acc_2_8 * src_dof_8 + q_acc_2_9 * src_dof_9;
               const walberla::float64 elMatVec_3 = q_acc_0_3 * src_dof_0 + q_acc_1_3 * src_dof_1 + q_acc_2_3 * src_dof_2 +
                                                    q_acc_3_10 * src_dof_10 + q_acc_3_11 * src_dof_11 + q_acc_3_3 * src_dof_3 +
                                                    q_acc_3_4 * src_dof_4 + q_acc_3_5 * src_dof_5 + q_acc_3_6 * src_dof_6 +
                                                    q_acc_3_7 * src_dof_7 + q_acc_3_8 * src_dof_8 + q_acc_3_9 * src_dof_9;
               const walberla::float64 elMatVec_4 = q_acc_0_4 * src_dof_0 + q_acc_1_4 * src_dof_1 + q_acc_2_4 * src_dof_2 +
                                                    q_acc_3_4 * src_dof_3 + q_acc_4_10 * src_dof_10 + q_acc_4_11 * src_dof_11 +
                                                    q_acc_4_4 * src_dof_4 + q_acc_4_5 * src_dof_5 + q_acc_4_6 * src_dof_6 +
                                                    q_acc_4_7 * src_dof_7 + q_acc_4_8 * src_dof_8 + q_acc_4_9 * src_dof_9;
               const walberla::float64 elMatVec_5 = q_acc_0_5 * src_dof_0 + q_acc_1_5 * src_dof_1 + q_acc_2_5 * src_dof_2 +
                                                    q_acc_3_5 * src_dof_3 + q_acc_4_5 * src_dof_4 + q_acc_5_10 * src_dof_10 +
                                                    q_acc_5_11 * src_dof_11 + q_acc_5_5 * src_dof_5 + q_acc_5_6 * src_dof_6 +
                                                    q_acc_5_7 * src_dof_7 + q_acc_5_8 * src_dof_8 + q_acc_5_9 * src_dof_9;
               const walberla::float64 elMatVec_6 = q_acc_0_6 * src_dof_0 + q_acc_1_6 * src_dof_1 + q_acc_2_6 * src_dof_2 +
                                                    q_acc_3_6 * src_dof_3 + q_acc_4_6 * src_dof_4 + q_acc_5_6 * src_dof_5 +
                                                    q_acc_6_10 * src_dof_10 + q_acc_6_11 * src_dof_11 + q_acc_6_6 * src_dof_6 +
                                                    q_acc_6_7 * src_dof_7 + q_acc_6_8 * src_dof_8 + q_acc_6_9 * src_dof_9;
               const walberla::float64 elMatVec_7 = q_acc_0_7 * src_dof_0 + q_acc_1_7 * src_dof_1 + q_acc_2_7 * src_dof_2 +
                                                    q_acc_3_7 * src_dof_3 + q_acc_4_7 * src_dof_4 + q_acc_5_7 * src_dof_5 +
                                                    q_acc_6_7 * src_dof_6 + q_acc_7_10 * src_dof_10 + q_acc_7_11 * src_dof_11 +
                                                    q_acc_7_7 * src_dof_7 + q_acc_7_8 * src_dof_8 + q_acc_7_9 * src_dof_9;
               const walberla::float64 elMatVec_8 = q_acc_0_8 * src_dof_0 + q_acc_1_8 * src_dof_1 + q_acc_2_8 * src_dof_2 +
                                                    q_acc_3_8 * src_dof_3 + q_acc_4_8 * src_dof_4 + q_acc_5_8 * src_dof_5 +
                                                    q_acc_6_8 * src_dof_6 + q_acc_7_8 * src_dof_7 + q_acc_8_10 * src_dof_10 +
                                                    q_acc_8_11 * src_dof_11 + q_acc_8_8 * src_dof_8 + q_acc_8_9 * src_dof_9;
               const walberla::float64 elMatVec_9 = q_acc_0_9 * src_dof_0 + q_acc_1_9 * src_dof_1 + q_acc_2_9 * src_dof_2 +
                                                    q_acc_3_9 * src_dof_3 + q_acc_4_9 * src_dof_4 + q_acc_5_9 * src_dof_5 +
                                                    q_acc_6_9 * src_dof_6 + q_acc_7_9 * src_dof_7 + q_acc_8_9 * src_dof_8 +
                                                    q_acc_9_10 * src_dof_10 + q_acc_9_11 * src_dof_11 + q_acc_9_9 * src_dof_9;
               const walberla::float64 elMatVec_10 =
                   q_acc_0_10 * src_dof_0 + q_acc_10_10 * src_dof_10 + q_acc_10_11 * src_dof_11 + q_acc_1_10 * src_dof_1 +
                   q_acc_2_10 * src_dof_2 + q_acc_3_10 * src_dof_3 + q_acc_4_10 * src_dof_4 + q_acc_5_10 * src_dof_5 +
                   q_acc_6_10 * src_dof_6 + q_acc_7_10 * src_dof_7 + q_acc_8_10 * src_dof_8 + q_acc_9_10 * src_dof_9;
               const walberla::float64 elMatVec_11 =
                   q_acc_0_11 * src_dof_0 + q_acc_10_11 * src_dof_10 + q_acc_11_11 * src_dof_11 + q_acc_1_11 * src_dof_1 +
                   q_acc_2_11 * src_dof_2 + q_acc_3_11 * src_dof_3 + q_acc_4_11 * src_dof_4 + q_acc_5_11 * src_dof_5 +
                   q_acc_6_11 * src_dof_6 + q_acc_7_11 * src_dof_7 + q_acc_8_11 * src_dof_8 + q_acc_9_11 * src_dof_9;
               _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                  1] = elMatVec_0 + _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                  ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_1 + _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                  ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] =
                   elMatVec_2 + _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               _data_dst_edge_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_3 + _data_dst_edge_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                 ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] =
                   elMatVec_4 +
                   _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1];
               _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatVec_5 +
                   _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                  1] = elMatVec_6 + _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                  ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_7 + _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                  ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] =
                   elMatVec_8 + _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               _data_dst_edge_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatVec_9 + _data_dst_edge_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                 ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] =
                   elMatVec_10 +
                   _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1];
               _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatVec_11 +
                   _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
            }
      }
   }
}
void P2VectorElementwiseFullStokesViscoplastic::toMatrix_P2VectorElementwiseFullStokesViscoplastic_macro_2D(
    idx_t* RESTRICT                      _data_dst_edge_0,
    idx_t* RESTRICT                      _data_dst_edge_1,
    idx_t* RESTRICT                      _data_dst_vertex_0,
    idx_t* RESTRICT                      _data_dst_vertex_1,
    walberla::float64* RESTRICT          _data_mu_lin,
    idx_t* RESTRICT                      _data_src_edge_0,
    idx_t* RESTRICT                      _data_src_edge_1,
    idx_t* RESTRICT                      _data_src_vertex_0,
    idx_t* RESTRICT                      _data_src_vertex_1,
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
      const walberla::float64 _data_q_w[] = { 0.16666666666666666, 0.16666666666666666, 0.16666666666666666 };

      const walberla::float64 _data_q_p_0[] = { 0.16666666666666666, 0.66666666666666663, 0.16666666666666666 };

      const walberla::float64 _data_q_p_1[] = { 0.66666666666666663, 0.16666666666666666, 0.16666666666666666 };

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
             jac_affine_0_0_GRAY * jac_affine_1_1_GRAY - jac_affine_0_1_GRAY * jac_affine_1_0_GRAY;
         const walberla::float64 tmp_coords_jac_2_GRAY   = 1.0 / ( tmp_coords_jac_1_GRAY );
         const walberla::float64 jac_affine_inv_0_0_GRAY = jac_affine_1_1_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_0_1_GRAY = -jac_affine_0_1_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_1_0_GRAY = -jac_affine_1_0_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_1_1_GRAY = jac_affine_0_0_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 abs_det_jac_affine_GRAY = abs( tmp_coords_jac_1_GRAY );
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
               const walberla::float64 uy_dof_2    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               walberla::float64       q_acc_0_0   = 0.0;
               walberla::float64       q_acc_0_1   = 0.0;
               walberla::float64       q_acc_0_2   = 0.0;
               walberla::float64       q_acc_0_3   = 0.0;
               walberla::float64       q_acc_0_4   = 0.0;
               walberla::float64       q_acc_0_5   = 0.0;
               walberla::float64       q_acc_0_6   = 0.0;
               walberla::float64       q_acc_0_7   = 0.0;
               walberla::float64       q_acc_0_8   = 0.0;
               walberla::float64       q_acc_0_9   = 0.0;
               walberla::float64       q_acc_0_10  = 0.0;
               walberla::float64       q_acc_0_11  = 0.0;
               walberla::float64       q_acc_1_1   = 0.0;
               walberla::float64       q_acc_1_2   = 0.0;
               walberla::float64       q_acc_1_3   = 0.0;
               walberla::float64       q_acc_1_4   = 0.0;
               walberla::float64       q_acc_1_5   = 0.0;
               walberla::float64       q_acc_1_6   = 0.0;
               walberla::float64       q_acc_1_7   = 0.0;
               walberla::float64       q_acc_1_8   = 0.0;
               walberla::float64       q_acc_1_9   = 0.0;
               walberla::float64       q_acc_1_10  = 0.0;
               walberla::float64       q_acc_1_11  = 0.0;
               walberla::float64       q_acc_2_2   = 0.0;
               walberla::float64       q_acc_2_3   = 0.0;
               walberla::float64       q_acc_2_4   = 0.0;
               walberla::float64       q_acc_2_5   = 0.0;
               walberla::float64       q_acc_2_6   = 0.0;
               walberla::float64       q_acc_2_7   = 0.0;
               walberla::float64       q_acc_2_8   = 0.0;
               walberla::float64       q_acc_2_9   = 0.0;
               walberla::float64       q_acc_2_10  = 0.0;
               walberla::float64       q_acc_2_11  = 0.0;
               walberla::float64       q_acc_3_3   = 0.0;
               walberla::float64       q_acc_3_4   = 0.0;
               walberla::float64       q_acc_3_5   = 0.0;
               walberla::float64       q_acc_3_6   = 0.0;
               walberla::float64       q_acc_3_7   = 0.0;
               walberla::float64       q_acc_3_8   = 0.0;
               walberla::float64       q_acc_3_9   = 0.0;
               walberla::float64       q_acc_3_10  = 0.0;
               walberla::float64       q_acc_3_11  = 0.0;
               walberla::float64       q_acc_4_4   = 0.0;
               walberla::float64       q_acc_4_5   = 0.0;
               walberla::float64       q_acc_4_6   = 0.0;
               walberla::float64       q_acc_4_7   = 0.0;
               walberla::float64       q_acc_4_8   = 0.0;
               walberla::float64       q_acc_4_9   = 0.0;
               walberla::float64       q_acc_4_10  = 0.0;
               walberla::float64       q_acc_4_11  = 0.0;
               walberla::float64       q_acc_5_5   = 0.0;
               walberla::float64       q_acc_5_6   = 0.0;
               walberla::float64       q_acc_5_7   = 0.0;
               walberla::float64       q_acc_5_8   = 0.0;
               walberla::float64       q_acc_5_9   = 0.0;
               walberla::float64       q_acc_5_10  = 0.0;
               walberla::float64       q_acc_5_11  = 0.0;
               walberla::float64       q_acc_6_6   = 0.0;
               walberla::float64       q_acc_6_7   = 0.0;
               walberla::float64       q_acc_6_8   = 0.0;
               walberla::float64       q_acc_6_9   = 0.0;
               walberla::float64       q_acc_6_10  = 0.0;
               walberla::float64       q_acc_6_11  = 0.0;
               walberla::float64       q_acc_7_7   = 0.0;
               walberla::float64       q_acc_7_8   = 0.0;
               walberla::float64       q_acc_7_9   = 0.0;
               walberla::float64       q_acc_7_10  = 0.0;
               walberla::float64       q_acc_7_11  = 0.0;
               walberla::float64       q_acc_8_8   = 0.0;
               walberla::float64       q_acc_8_9   = 0.0;
               walberla::float64       q_acc_8_10  = 0.0;
               walberla::float64       q_acc_8_11  = 0.0;
               walberla::float64       q_acc_9_9   = 0.0;
               walberla::float64       q_acc_9_10  = 0.0;
               walberla::float64       q_acc_9_11  = 0.0;
               walberla::float64       q_acc_10_10 = 0.0;
               walberla::float64       q_acc_10_11 = 0.0;
               walberla::float64       q_acc_11_11 = 0.0;
               for ( int64_t q = 0; q < 3; q += 1 )
               {
                  const walberla::float64 tmp_qloop_0  = 4.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_1  = 4.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_2  = tmp_qloop_0 + tmp_qloop_1 - 3.0;
                  const walberla::float64 tmp_qloop_3  = jac_affine_inv_0_0_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_4  = jac_affine_inv_1_0_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_5  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_6  = abs_det_jac_affine_GRAY * tmp_qloop_5;
                  const walberla::float64 tmp_qloop_7  = jac_affine_inv_0_1_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_8  = jac_affine_inv_1_1_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_9  = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_10 = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_11 = -ux_dof_0 + ux_dof_1;
                  const walberla::float64 tmp_qloop_12 = -ux_dof_0 + ux_dof_2;
                  const walberla::float64 tmp_qloop_13 = -uy_dof_0 + uy_dof_1;
                  const walberla::float64 tmp_qloop_14 = -uy_dof_0 + uy_dof_2;
                  const walberla::float64 tmp_qloop_15 = jac_affine_inv_0_0_GRAY * tmp_qloop_13;
                  const walberla::float64 tmp_qloop_16 = jac_affine_inv_0_1_GRAY * tmp_qloop_11;
                  const walberla::float64 tmp_qloop_17 = jac_affine_inv_1_0_GRAY * tmp_qloop_14;
                  const walberla::float64 tmp_qloop_18 = jac_affine_inv_1_1_GRAY * tmp_qloop_12;
                  const walberla::float64 tmp_qloop_19 =
                      1.0 /
                      ( 1.0 /
                            ( mu_star +
                              sigma_y * 1.0 /
                                  ( pow( ( ( jac_affine_inv_0_0_GRAY * tmp_qloop_11 + jac_affine_inv_1_0_GRAY * tmp_qloop_12 ) *
                                           ( jac_affine_inv_0_0_GRAY * tmp_qloop_11 + jac_affine_inv_1_0_GRAY * tmp_qloop_12 ) ) +
                                             ( ( jac_affine_inv_0_1_GRAY * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_GRAY * tmp_qloop_14 ) *
                                               ( jac_affine_inv_0_1_GRAY * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_GRAY * tmp_qloop_14 ) ) +
                                             ( tmp_qloop_15 + tmp_qloop_16 + tmp_qloop_17 + tmp_qloop_18 ) *
                                                 ( tmp_qloop_15 * 0.5 + tmp_qloop_16 * 0.5 + tmp_qloop_17 * 0.5 +
                                                   tmp_qloop_18 * 0.5 ),
                                         0.50000000000000000 ) +
                                    9.9999999999999995e-21 ) ) *
                            1.0 +
                        1.0 /
                            ( mu_lin_dof_0 * ( 1.0 - _data_q_p_0[q] - _data_q_p_1[q] ) + mu_lin_dof_1 * _data_q_p_0[q] +
                              mu_lin_dof_2 * _data_q_p_1[q] ) *
                            1.0 ) *
                      _data_q_w[q];
                  const walberla::float64 tmp_qloop_20 = tmp_qloop_19 * 2.0;
                  const walberla::float64 tmp_qloop_21 = tmp_qloop_0 - 1.0;
                  const walberla::float64 tmp_qloop_22 = jac_affine_inv_0_0_GRAY * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_23 = tmp_qloop_6 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_24 = tmp_qloop_5 * 2.0;
                  const walberla::float64 tmp_qloop_25 = jac_affine_inv_0_1_GRAY * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_26 = tmp_qloop_1 - 1.0;
                  const walberla::float64 tmp_qloop_27 = jac_affine_inv_1_0_GRAY * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_28 = jac_affine_inv_1_1_GRAY * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_29 = 2.6666666666666665 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_30 = jac_affine_inv_1_0_GRAY * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_31 = 2.6666666666666665 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_32 = jac_affine_inv_0_0_GRAY * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_33 = tmp_qloop_30 + tmp_qloop_32;
                  const walberla::float64 tmp_qloop_34 = jac_affine_inv_1_0_GRAY * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_35 = jac_affine_inv_0_0_GRAY * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_36 = tmp_qloop_34 + tmp_qloop_35;
                  const walberla::float64 tmp_qloop_37 = jac_affine_inv_1_1_GRAY * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_38 = jac_affine_inv_0_1_GRAY * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_39 = tmp_qloop_37 + tmp_qloop_38;
                  const walberla::float64 tmp_qloop_40 = -tmp_qloop_0 - 8.0 * _data_q_p_1[q] + 4.0;
                  const walberla::float64 tmp_qloop_41 =
                      jac_affine_inv_1_0_GRAY * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_32;
                  const walberla::float64 tmp_qloop_42 = jac_affine_inv_1_0_GRAY * tmp_qloop_40 - tmp_qloop_35;
                  const walberla::float64 tmp_qloop_43 = jac_affine_inv_1_1_GRAY * tmp_qloop_40 - tmp_qloop_38;
                  const walberla::float64 tmp_qloop_44 = -tmp_qloop_1 - 8.0 * _data_q_p_0[q] + 4.0;
                  const walberla::float64 tmp_qloop_45 =
                      jac_affine_inv_0_0_GRAY * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_30;
                  const walberla::float64 tmp_qloop_46 = jac_affine_inv_0_0_GRAY * tmp_qloop_44 - tmp_qloop_34;
                  const walberla::float64 tmp_qloop_47 = jac_affine_inv_0_1_GRAY * tmp_qloop_44 - tmp_qloop_37;
                  const walberla::float64 tmp_qloop_48 = tmp_qloop_7 * 0.66666666666666663 + tmp_qloop_8 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_49 = abs_det_jac_affine_GRAY * tmp_qloop_22;
                  const walberla::float64 tmp_qloop_50 = abs_det_jac_affine_GRAY * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_51 = jac_affine_inv_1_1_GRAY * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_52 = jac_affine_inv_0_1_GRAY * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_53 = tmp_qloop_51 + tmp_qloop_52;
                  const walberla::float64 tmp_qloop_54 = abs_det_jac_affine_GRAY * tmp_qloop_36;
                  const walberla::float64 tmp_qloop_55 =
                      jac_affine_inv_1_1_GRAY * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_52;
                  const walberla::float64 tmp_qloop_56 = abs_det_jac_affine_GRAY * tmp_qloop_42;
                  const walberla::float64 tmp_qloop_57 =
                      jac_affine_inv_0_1_GRAY * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_51;
                  const walberla::float64 tmp_qloop_58  = abs_det_jac_affine_GRAY * tmp_qloop_46;
                  const walberla::float64 tmp_qloop_59  = ( jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY );
                  const walberla::float64 tmp_qloop_60  = ( tmp_qloop_21 * tmp_qloop_21 );
                  const walberla::float64 tmp_qloop_61  = abs_det_jac_affine_GRAY * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_62  = tmp_qloop_60 * tmp_qloop_61;
                  const walberla::float64 tmp_qloop_63  = tmp_qloop_60 * 2.0;
                  const walberla::float64 tmp_qloop_64  = ( jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY );
                  const walberla::float64 tmp_qloop_65  = tmp_qloop_60 * 1.0;
                  const walberla::float64 tmp_qloop_66  = tmp_qloop_22 * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_67  = tmp_qloop_25 * 1.0;
                  const walberla::float64 tmp_qloop_68  = tmp_qloop_36 * 2.0;
                  const walberla::float64 tmp_qloop_69  = tmp_qloop_42 * 2.0;
                  const walberla::float64 tmp_qloop_70  = tmp_qloop_46 * 2.0;
                  const walberla::float64 tmp_qloop_71  = abs_det_jac_affine_GRAY * tmp_qloop_19 * 0.66666666666666674;
                  const walberla::float64 tmp_qloop_72  = abs_det_jac_affine_GRAY * tmp_qloop_57;
                  const walberla::float64 tmp_qloop_73  = ( jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY );
                  const walberla::float64 tmp_qloop_74  = ( tmp_qloop_26 * tmp_qloop_26 );
                  const walberla::float64 tmp_qloop_75  = tmp_qloop_61 * tmp_qloop_74;
                  const walberla::float64 tmp_qloop_76  = tmp_qloop_74 * 2.0;
                  const walberla::float64 tmp_qloop_77  = ( jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY );
                  const walberla::float64 tmp_qloop_78  = tmp_qloop_74 * 1.0;
                  const walberla::float64 tmp_qloop_79  = tmp_qloop_28 * 1.0;
                  const walberla::float64 tmp_qloop_80  = 2.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_81  = jac_affine_inv_1_1_GRAY * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_82  = 2.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_83  = jac_affine_inv_0_1_GRAY * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_84  = tmp_qloop_81 + tmp_qloop_83;
                  const walberla::float64 tmp_qloop_85  = tmp_qloop_84 * 2.0;
                  const walberla::float64 tmp_qloop_86  = tmp_qloop_54 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_87  = jac_affine_inv_1_1_GRAY * tmp_qloop_40 * 0.5 - tmp_qloop_83;
                  const walberla::float64 tmp_qloop_88  = tmp_qloop_87 * 2.0;
                  const walberla::float64 tmp_qloop_89  = tmp_qloop_56 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_90  = jac_affine_inv_0_1_GRAY * tmp_qloop_44 * 0.5 - tmp_qloop_81;
                  const walberla::float64 tmp_qloop_91  = tmp_qloop_90 * 2.0;
                  const walberla::float64 tmp_qloop_92  = tmp_qloop_58 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_93  = abs_det_jac_affine_GRAY * tmp_qloop_9;
                  const walberla::float64 tmp_qloop_94  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_95  = tmp_qloop_93 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_96  = tmp_qloop_9 * 2.0;
                  const walberla::float64 tmp_qloop_97  = tmp_qloop_25 * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_98  = tmp_qloop_22 * 1.0;
                  const walberla::float64 tmp_qloop_99  = abs_det_jac_affine_GRAY * tmp_qloop_25;
                  const walberla::float64 tmp_qloop_100 = tmp_qloop_39 * 2.0;
                  const walberla::float64 tmp_qloop_101 = tmp_qloop_43 * 2.0;
                  const walberla::float64 tmp_qloop_102 = tmp_qloop_47 * 2.0;
                  const walberla::float64 tmp_qloop_103 = abs_det_jac_affine_GRAY * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_104 = tmp_qloop_27 * 1.0;
                  const walberla::float64 tmp_qloop_105 = abs_det_jac_affine_GRAY * tmp_qloop_39;
                  const walberla::float64 tmp_qloop_106 = jac_affine_inv_1_0_GRAY * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_107 = jac_affine_inv_0_0_GRAY * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_108 = tmp_qloop_106 * 2.0 + tmp_qloop_107 * 2.0;
                  const walberla::float64 tmp_qloop_109 = abs_det_jac_affine_GRAY * tmp_qloop_43;
                  const walberla::float64 tmp_qloop_110 = jac_affine_inv_1_0_GRAY * tmp_qloop_40 + tmp_qloop_107 * -2.0;
                  const walberla::float64 q_tmp_0_0 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_9 + ( tmp_qloop_5 * tmp_qloop_5 ) * 2.0 ) -
                        tmp_qloop_6 * ( tmp_qloop_3 * 0.66666666666666663 + tmp_qloop_4 * 0.66666666666666663 ) );
                  const walberla::float64 q_tmp_0_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_25 + tmp_qloop_22 * tmp_qloop_24 ) -
                                       tmp_qloop_22 * tmp_qloop_23 );
                  const walberla::float64 q_tmp_0_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_28 + tmp_qloop_24 * tmp_qloop_27 ) -
                                       tmp_qloop_23 * tmp_qloop_27 );
                  const walberla::float64 q_tmp_0_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_39 + tmp_qloop_24 * tmp_qloop_36 ) -
                                       tmp_qloop_33 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_43 + tmp_qloop_24 * tmp_qloop_42 ) -
                                       tmp_qloop_41 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_10 * tmp_qloop_47 + tmp_qloop_24 * tmp_qloop_46 ) -
                                       tmp_qloop_45 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_6 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_6 - tmp_qloop_48 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_7 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_49 - tmp_qloop_23 * tmp_qloop_25 );
                  const walberla::float64 q_tmp_0_8 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_50 - tmp_qloop_23 * tmp_qloop_28 );
                  const walberla::float64 q_tmp_0_9 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_54 - tmp_qloop_53 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_10 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_56 - tmp_qloop_55 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_11 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_58 - tmp_qloop_57 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_1_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_59 * tmp_qloop_63 + tmp_qloop_64 * tmp_qloop_65 ) -
                                       tmp_qloop_59 * tmp_qloop_62 );
                  const walberla::float64 q_tmp_1_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_28 * tmp_qloop_67 + tmp_qloop_66 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_66 );
                  const walberla::float64 q_tmp_1_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_67 ) -
                                       tmp_qloop_33 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_67 ) -
                                       tmp_qloop_41 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_67 ) -
                                       tmp_qloop_45 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_7 =
                      jac_affine_inv_0_0_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_60 * tmp_qloop_71;
                  const walberla::float64 q_tmp_1_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * jac_affine_inv_1_0_GRAY *
                                           tmp_qloop_21 * tmp_qloop_26 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_28 * tmp_qloop_61 );
                  const walberla::float64 q_tmp_1_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_1_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_1_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_0_1_GRAY * tmp_qloop_21 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_2_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_73 * tmp_qloop_76 + tmp_qloop_77 * tmp_qloop_78 ) -
                                       tmp_qloop_73 * tmp_qloop_75 );
                  const walberla::float64 q_tmp_2_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_79 ) -
                                       tmp_qloop_33 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_79 ) -
                                       tmp_qloop_41 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_79 ) -
                                       tmp_qloop_45 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_27 * tmp_qloop_61 + tmp_qloop_49 * tmp_qloop_79 );
                  const walberla::float64 q_tmp_2_8 =
                      jac_affine_inv_1_0_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_71 * tmp_qloop_74;
                  const walberla::float64 q_tmp_2_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_2_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_2_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * jac_affine_inv_1_1_GRAY * tmp_qloop_26 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_27 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_3_3 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( ( tmp_qloop_36 * tmp_qloop_36 ) * 2.0 + tmp_qloop_39 * tmp_qloop_85 ) -
                        tmp_qloop_33 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_42 * tmp_qloop_68 + tmp_qloop_43 * tmp_qloop_85 ) -
                                       tmp_qloop_41 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_46 * tmp_qloop_68 + tmp_qloop_47 * tmp_qloop_85 ) -
                                       tmp_qloop_45 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * tmp_qloop_5 * tmp_qloop_84 * 2.0 - tmp_qloop_48 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_86 + tmp_qloop_49 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_86 + tmp_qloop_50 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_54 + tmp_qloop_54 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_42 * tmp_qloop_84 * 2.0 - tmp_qloop_54 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_3_11 =
                      tmp_qloop_20 * ( -tmp_qloop_54 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_4_4 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( ( tmp_qloop_42 * tmp_qloop_42 ) * 2.0 + tmp_qloop_43 * tmp_qloop_88 ) -
                        tmp_qloop_41 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_46 * tmp_qloop_69 + tmp_qloop_47 * tmp_qloop_88 ) -
                                       tmp_qloop_45 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * tmp_qloop_5 * tmp_qloop_87 * 2.0 - tmp_qloop_48 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_89 + tmp_qloop_49 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_89 + tmp_qloop_50 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_56 + tmp_qloop_54 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_42 * tmp_qloop_87 * 2.0 - tmp_qloop_55 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_11 =
                      tmp_qloop_20 * ( -tmp_qloop_56 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_5_5 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( ( tmp_qloop_46 * tmp_qloop_46 ) * 2.0 + tmp_qloop_47 * tmp_qloop_91 ) -
                        tmp_qloop_45 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * tmp_qloop_5 * tmp_qloop_90 * 2.0 - tmp_qloop_48 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_92 + tmp_qloop_49 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_92 + tmp_qloop_50 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_36 * tmp_qloop_90 * 2.0 - tmp_qloop_53 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_42 * tmp_qloop_90 * 2.0 - tmp_qloop_55 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * tmp_qloop_46 * tmp_qloop_90 * 2.0 - tmp_qloop_57 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_6_6 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_5 * tmp_qloop_94 + ( tmp_qloop_9 * tmp_qloop_9 ) * 2.0 ) -
                        tmp_qloop_48 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_22 * tmp_qloop_94 + tmp_qloop_25 * tmp_qloop_96 ) -
                                       tmp_qloop_25 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_94 + tmp_qloop_28 * tmp_qloop_96 ) -
                                       tmp_qloop_28 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_36 * tmp_qloop_94 + tmp_qloop_39 * tmp_qloop_96 ) -
                                       tmp_qloop_53 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_42 * tmp_qloop_94 + tmp_qloop_43 * tmp_qloop_96 ) -
                                       tmp_qloop_55 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_46 * tmp_qloop_94 + tmp_qloop_47 * tmp_qloop_96 ) -
                                       tmp_qloop_57 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_7_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_59 * tmp_qloop_65 + tmp_qloop_63 * tmp_qloop_64 ) -
                                       tmp_qloop_62 * tmp_qloop_64 );
                  const walberla::float64 q_tmp_7_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_27 * tmp_qloop_98 + tmp_qloop_97 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_97 );
                  const walberla::float64 q_tmp_7_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_25 + tmp_qloop_36 * tmp_qloop_98 ) -
                                       tmp_qloop_53 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_101 * tmp_qloop_25 + tmp_qloop_42 * tmp_qloop_98 ) -
                                       tmp_qloop_55 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_102 * tmp_qloop_25 + tmp_qloop_46 * tmp_qloop_98 ) -
                                       tmp_qloop_25 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_8_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_73 * tmp_qloop_78 + tmp_qloop_76 * tmp_qloop_77 ) -
                                       tmp_qloop_75 * tmp_qloop_77 );
                  const walberla::float64 q_tmp_8_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_36 ) -
                                       tmp_qloop_103 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_8_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_101 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_42 ) -
                                       tmp_qloop_103 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_8_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_102 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_46 ) -
                                       tmp_qloop_28 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_9_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_108 * tmp_qloop_36 + ( tmp_qloop_39 * tmp_qloop_39 ) * 2.0 ) -
                        tmp_qloop_105 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_9_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_43 + tmp_qloop_108 * tmp_qloop_42 ) -
                                       tmp_qloop_105 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_9_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_100 * tmp_qloop_47 + tmp_qloop_108 * tmp_qloop_46 ) -
                                       tmp_qloop_105 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_10_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY * ( tmp_qloop_110 * tmp_qloop_42 + ( tmp_qloop_43 * tmp_qloop_43 ) * 2.0 ) -
                        tmp_qloop_109 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_10_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_101 * tmp_qloop_47 + tmp_qloop_110 * tmp_qloop_46 ) -
                                       tmp_qloop_109 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_11_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_GRAY *
                            ( tmp_qloop_46 * ( jac_affine_inv_0_0_GRAY * tmp_qloop_44 * 0.5 - tmp_qloop_106 ) * 2.0 +
                              ( tmp_qloop_47 * tmp_qloop_47 ) * 2.0 ) -
                        tmp_qloop_47 * tmp_qloop_72 );
                  q_acc_0_0   = q_acc_0_0 + q_tmp_0_0;
                  q_acc_0_1   = q_acc_0_1 + q_tmp_0_1;
                  q_acc_0_2   = q_acc_0_2 + q_tmp_0_2;
                  q_acc_0_3   = q_acc_0_3 + q_tmp_0_3;
                  q_acc_0_4   = q_acc_0_4 + q_tmp_0_4;
                  q_acc_0_5   = q_acc_0_5 + q_tmp_0_5;
                  q_acc_0_6   = q_acc_0_6 + q_tmp_0_6;
                  q_acc_0_7   = q_acc_0_7 + q_tmp_0_7;
                  q_acc_0_8   = q_acc_0_8 + q_tmp_0_8;
                  q_acc_0_9   = q_acc_0_9 + q_tmp_0_9;
                  q_acc_0_10  = q_acc_0_10 + q_tmp_0_10;
                  q_acc_0_11  = q_acc_0_11 + q_tmp_0_11;
                  q_acc_1_1   = q_acc_1_1 + q_tmp_1_1;
                  q_acc_1_2   = q_acc_1_2 + q_tmp_1_2;
                  q_acc_1_3   = q_acc_1_3 + q_tmp_1_3;
                  q_acc_1_4   = q_acc_1_4 + q_tmp_1_4;
                  q_acc_1_5   = q_acc_1_5 + q_tmp_1_5;
                  q_acc_1_6   = q_acc_1_6 + q_tmp_1_6;
                  q_acc_1_7   = q_acc_1_7 + q_tmp_1_7;
                  q_acc_1_8   = q_acc_1_8 + q_tmp_1_8;
                  q_acc_1_9   = q_acc_1_9 + q_tmp_1_9;
                  q_acc_1_10  = q_acc_1_10 + q_tmp_1_10;
                  q_acc_1_11  = q_acc_1_11 + q_tmp_1_11;
                  q_acc_2_2   = q_acc_2_2 + q_tmp_2_2;
                  q_acc_2_3   = q_acc_2_3 + q_tmp_2_3;
                  q_acc_2_4   = q_acc_2_4 + q_tmp_2_4;
                  q_acc_2_5   = q_acc_2_5 + q_tmp_2_5;
                  q_acc_2_6   = q_acc_2_6 + q_tmp_2_6;
                  q_acc_2_7   = q_acc_2_7 + q_tmp_2_7;
                  q_acc_2_8   = q_acc_2_8 + q_tmp_2_8;
                  q_acc_2_9   = q_acc_2_9 + q_tmp_2_9;
                  q_acc_2_10  = q_acc_2_10 + q_tmp_2_10;
                  q_acc_2_11  = q_acc_2_11 + q_tmp_2_11;
                  q_acc_3_3   = q_acc_3_3 + q_tmp_3_3;
                  q_acc_3_4   = q_acc_3_4 + q_tmp_3_4;
                  q_acc_3_5   = q_acc_3_5 + q_tmp_3_5;
                  q_acc_3_6   = q_acc_3_6 + q_tmp_3_6;
                  q_acc_3_7   = q_acc_3_7 + q_tmp_3_7;
                  q_acc_3_8   = q_acc_3_8 + q_tmp_3_8;
                  q_acc_3_9   = q_acc_3_9 + q_tmp_3_9;
                  q_acc_3_10  = q_acc_3_10 + q_tmp_3_10;
                  q_acc_3_11  = q_acc_3_11 + q_tmp_3_11;
                  q_acc_4_4   = q_acc_4_4 + q_tmp_4_4;
                  q_acc_4_5   = q_acc_4_5 + q_tmp_4_5;
                  q_acc_4_6   = q_acc_4_6 + q_tmp_4_6;
                  q_acc_4_7   = q_acc_4_7 + q_tmp_4_7;
                  q_acc_4_8   = q_acc_4_8 + q_tmp_4_8;
                  q_acc_4_9   = q_acc_4_9 + q_tmp_4_9;
                  q_acc_4_10  = q_acc_4_10 + q_tmp_4_10;
                  q_acc_4_11  = q_acc_4_11 + q_tmp_4_11;
                  q_acc_5_5   = q_acc_5_5 + q_tmp_5_5;
                  q_acc_5_6   = q_acc_5_6 + q_tmp_5_6;
                  q_acc_5_7   = q_acc_5_7 + q_tmp_5_7;
                  q_acc_5_8   = q_acc_5_8 + q_tmp_5_8;
                  q_acc_5_9   = q_acc_5_9 + q_tmp_5_9;
                  q_acc_5_10  = q_acc_5_10 + q_tmp_5_10;
                  q_acc_5_11  = q_acc_5_11 + q_tmp_5_11;
                  q_acc_6_6   = q_acc_6_6 + q_tmp_6_6;
                  q_acc_6_7   = q_acc_6_7 + q_tmp_6_7;
                  q_acc_6_8   = q_acc_6_8 + q_tmp_6_8;
                  q_acc_6_9   = q_acc_6_9 + q_tmp_6_9;
                  q_acc_6_10  = q_acc_6_10 + q_tmp_6_10;
                  q_acc_6_11  = q_acc_6_11 + q_tmp_6_11;
                  q_acc_7_7   = q_acc_7_7 + q_tmp_7_7;
                  q_acc_7_8   = q_acc_7_8 + q_tmp_7_8;
                  q_acc_7_9   = q_acc_7_9 + q_tmp_7_9;
                  q_acc_7_10  = q_acc_7_10 + q_tmp_7_10;
                  q_acc_7_11  = q_acc_7_11 + q_tmp_7_11;
                  q_acc_8_8   = q_acc_8_8 + q_tmp_8_8;
                  q_acc_8_9   = q_acc_8_9 + q_tmp_8_9;
                  q_acc_8_10  = q_acc_8_10 + q_tmp_8_10;
                  q_acc_8_11  = q_acc_8_11 + q_tmp_8_11;
                  q_acc_9_9   = q_acc_9_9 + q_tmp_9_9;
                  q_acc_9_10  = q_acc_9_10 + q_tmp_9_10;
                  q_acc_9_11  = q_acc_9_11 + q_tmp_9_11;
                  q_acc_10_10 = q_acc_10_10 + q_tmp_10_10;
                  q_acc_10_11 = q_acc_10_11 + q_tmp_10_11;
                  q_acc_11_11 = q_acc_11_11 + q_tmp_11_11;
               }
               const walberla::float64 elMat_0_0   = q_acc_0_0;
               const walberla::float64 elMat_0_1   = q_acc_0_1;
               const walberla::float64 elMat_0_2   = q_acc_0_2;
               const walberla::float64 elMat_0_3   = q_acc_0_3;
               const walberla::float64 elMat_0_4   = q_acc_0_4;
               const walberla::float64 elMat_0_5   = q_acc_0_5;
               const walberla::float64 elMat_0_6   = q_acc_0_6;
               const walberla::float64 elMat_0_7   = q_acc_0_7;
               const walberla::float64 elMat_0_8   = q_acc_0_8;
               const walberla::float64 elMat_0_9   = q_acc_0_9;
               const walberla::float64 elMat_0_10  = q_acc_0_10;
               const walberla::float64 elMat_0_11  = q_acc_0_11;
               const walberla::float64 elMat_1_0   = q_acc_0_1;
               const walberla::float64 elMat_1_1   = q_acc_1_1;
               const walberla::float64 elMat_1_2   = q_acc_1_2;
               const walberla::float64 elMat_1_3   = q_acc_1_3;
               const walberla::float64 elMat_1_4   = q_acc_1_4;
               const walberla::float64 elMat_1_5   = q_acc_1_5;
               const walberla::float64 elMat_1_6   = q_acc_1_6;
               const walberla::float64 elMat_1_7   = q_acc_1_7;
               const walberla::float64 elMat_1_8   = q_acc_1_8;
               const walberla::float64 elMat_1_9   = q_acc_1_9;
               const walberla::float64 elMat_1_10  = q_acc_1_10;
               const walberla::float64 elMat_1_11  = q_acc_1_11;
               const walberla::float64 elMat_2_0   = q_acc_0_2;
               const walberla::float64 elMat_2_1   = q_acc_1_2;
               const walberla::float64 elMat_2_2   = q_acc_2_2;
               const walberla::float64 elMat_2_3   = q_acc_2_3;
               const walberla::float64 elMat_2_4   = q_acc_2_4;
               const walberla::float64 elMat_2_5   = q_acc_2_5;
               const walberla::float64 elMat_2_6   = q_acc_2_6;
               const walberla::float64 elMat_2_7   = q_acc_2_7;
               const walberla::float64 elMat_2_8   = q_acc_2_8;
               const walberla::float64 elMat_2_9   = q_acc_2_9;
               const walberla::float64 elMat_2_10  = q_acc_2_10;
               const walberla::float64 elMat_2_11  = q_acc_2_11;
               const walberla::float64 elMat_3_0   = q_acc_0_3;
               const walberla::float64 elMat_3_1   = q_acc_1_3;
               const walberla::float64 elMat_3_2   = q_acc_2_3;
               const walberla::float64 elMat_3_3   = q_acc_3_3;
               const walberla::float64 elMat_3_4   = q_acc_3_4;
               const walberla::float64 elMat_3_5   = q_acc_3_5;
               const walberla::float64 elMat_3_6   = q_acc_3_6;
               const walberla::float64 elMat_3_7   = q_acc_3_7;
               const walberla::float64 elMat_3_8   = q_acc_3_8;
               const walberla::float64 elMat_3_9   = q_acc_3_9;
               const walberla::float64 elMat_3_10  = q_acc_3_10;
               const walberla::float64 elMat_3_11  = q_acc_3_11;
               const walberla::float64 elMat_4_0   = q_acc_0_4;
               const walberla::float64 elMat_4_1   = q_acc_1_4;
               const walberla::float64 elMat_4_2   = q_acc_2_4;
               const walberla::float64 elMat_4_3   = q_acc_3_4;
               const walberla::float64 elMat_4_4   = q_acc_4_4;
               const walberla::float64 elMat_4_5   = q_acc_4_5;
               const walberla::float64 elMat_4_6   = q_acc_4_6;
               const walberla::float64 elMat_4_7   = q_acc_4_7;
               const walberla::float64 elMat_4_8   = q_acc_4_8;
               const walberla::float64 elMat_4_9   = q_acc_4_9;
               const walberla::float64 elMat_4_10  = q_acc_4_10;
               const walberla::float64 elMat_4_11  = q_acc_4_11;
               const walberla::float64 elMat_5_0   = q_acc_0_5;
               const walberla::float64 elMat_5_1   = q_acc_1_5;
               const walberla::float64 elMat_5_2   = q_acc_2_5;
               const walberla::float64 elMat_5_3   = q_acc_3_5;
               const walberla::float64 elMat_5_4   = q_acc_4_5;
               const walberla::float64 elMat_5_5   = q_acc_5_5;
               const walberla::float64 elMat_5_6   = q_acc_5_6;
               const walberla::float64 elMat_5_7   = q_acc_5_7;
               const walberla::float64 elMat_5_8   = q_acc_5_8;
               const walberla::float64 elMat_5_9   = q_acc_5_9;
               const walberla::float64 elMat_5_10  = q_acc_5_10;
               const walberla::float64 elMat_5_11  = q_acc_5_11;
               const walberla::float64 elMat_6_0   = q_acc_0_6;
               const walberla::float64 elMat_6_1   = q_acc_1_6;
               const walberla::float64 elMat_6_2   = q_acc_2_6;
               const walberla::float64 elMat_6_3   = q_acc_3_6;
               const walberla::float64 elMat_6_4   = q_acc_4_6;
               const walberla::float64 elMat_6_5   = q_acc_5_6;
               const walberla::float64 elMat_6_6   = q_acc_6_6;
               const walberla::float64 elMat_6_7   = q_acc_6_7;
               const walberla::float64 elMat_6_8   = q_acc_6_8;
               const walberla::float64 elMat_6_9   = q_acc_6_9;
               const walberla::float64 elMat_6_10  = q_acc_6_10;
               const walberla::float64 elMat_6_11  = q_acc_6_11;
               const walberla::float64 elMat_7_0   = q_acc_0_7;
               const walberla::float64 elMat_7_1   = q_acc_1_7;
               const walberla::float64 elMat_7_2   = q_acc_2_7;
               const walberla::float64 elMat_7_3   = q_acc_3_7;
               const walberla::float64 elMat_7_4   = q_acc_4_7;
               const walberla::float64 elMat_7_5   = q_acc_5_7;
               const walberla::float64 elMat_7_6   = q_acc_6_7;
               const walberla::float64 elMat_7_7   = q_acc_7_7;
               const walberla::float64 elMat_7_8   = q_acc_7_8;
               const walberla::float64 elMat_7_9   = q_acc_7_9;
               const walberla::float64 elMat_7_10  = q_acc_7_10;
               const walberla::float64 elMat_7_11  = q_acc_7_11;
               const walberla::float64 elMat_8_0   = q_acc_0_8;
               const walberla::float64 elMat_8_1   = q_acc_1_8;
               const walberla::float64 elMat_8_2   = q_acc_2_8;
               const walberla::float64 elMat_8_3   = q_acc_3_8;
               const walberla::float64 elMat_8_4   = q_acc_4_8;
               const walberla::float64 elMat_8_5   = q_acc_5_8;
               const walberla::float64 elMat_8_6   = q_acc_6_8;
               const walberla::float64 elMat_8_7   = q_acc_7_8;
               const walberla::float64 elMat_8_8   = q_acc_8_8;
               const walberla::float64 elMat_8_9   = q_acc_8_9;
               const walberla::float64 elMat_8_10  = q_acc_8_10;
               const walberla::float64 elMat_8_11  = q_acc_8_11;
               const walberla::float64 elMat_9_0   = q_acc_0_9;
               const walberla::float64 elMat_9_1   = q_acc_1_9;
               const walberla::float64 elMat_9_2   = q_acc_2_9;
               const walberla::float64 elMat_9_3   = q_acc_3_9;
               const walberla::float64 elMat_9_4   = q_acc_4_9;
               const walberla::float64 elMat_9_5   = q_acc_5_9;
               const walberla::float64 elMat_9_6   = q_acc_6_9;
               const walberla::float64 elMat_9_7   = q_acc_7_9;
               const walberla::float64 elMat_9_8   = q_acc_8_9;
               const walberla::float64 elMat_9_9   = q_acc_9_9;
               const walberla::float64 elMat_9_10  = q_acc_9_10;
               const walberla::float64 elMat_9_11  = q_acc_9_11;
               const walberla::float64 elMat_10_0  = q_acc_0_10;
               const walberla::float64 elMat_10_1  = q_acc_1_10;
               const walberla::float64 elMat_10_2  = q_acc_2_10;
               const walberla::float64 elMat_10_3  = q_acc_3_10;
               const walberla::float64 elMat_10_4  = q_acc_4_10;
               const walberla::float64 elMat_10_5  = q_acc_5_10;
               const walberla::float64 elMat_10_6  = q_acc_6_10;
               const walberla::float64 elMat_10_7  = q_acc_7_10;
               const walberla::float64 elMat_10_8  = q_acc_8_10;
               const walberla::float64 elMat_10_9  = q_acc_9_10;
               const walberla::float64 elMat_10_10 = q_acc_10_10;
               const walberla::float64 elMat_10_11 = q_acc_10_11;
               const walberla::float64 elMat_11_0  = q_acc_0_11;
               const walberla::float64 elMat_11_1  = q_acc_1_11;
               const walberla::float64 elMat_11_2  = q_acc_2_11;
               const walberla::float64 elMat_11_3  = q_acc_3_11;
               const walberla::float64 elMat_11_4  = q_acc_4_11;
               const walberla::float64 elMat_11_5  = q_acc_5_11;
               const walberla::float64 elMat_11_6  = q_acc_6_11;
               const walberla::float64 elMat_11_7  = q_acc_7_11;
               const walberla::float64 elMat_11_8  = q_acc_8_11;
               const walberla::float64 elMat_11_9  = q_acc_9_11;
               const walberla::float64 elMat_11_10 = q_acc_10_11;
               const walberla::float64 elMat_11_11 = q_acc_11_11;

               std::vector< uint_t > _data_rowIdx( 12 );
               std::vector< uint_t > _data_colIdx( 12 );
               std::vector< real_t > _data_mat( 144 );

               _data_rowIdx[0] = ( (uint64_t) ( _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_rowIdx[1] = ( (uint64_t) ( _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[2] = ( (uint64_t) ( _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_rowIdx[3] =
                   ( (uint64_t) ( _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );
               _data_rowIdx[4] =
                   ( (uint64_t) ( _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                         ( 2 ) )] ) );
               _data_rowIdx[5] = ( (uint64_t) ( _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                                 ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_rowIdx[6] = ( (uint64_t) ( _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_rowIdx[7] = ( (uint64_t) ( _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[8] = ( (uint64_t) ( _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_rowIdx[9] =
                   ( (uint64_t) ( _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );
               _data_rowIdx[10] =
                   ( (uint64_t) ( _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                         ( 2 ) )] ) );
               _data_rowIdx[11] = ( (uint64_t) ( _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                                  ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_colIdx[0]  = ( (uint64_t) ( _data_src_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_colIdx[1]  = ( (uint64_t) ( _data_src_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[2]  = ( (uint64_t) ( _data_src_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[3] =
                   ( (uint64_t) ( _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );
               _data_colIdx[4] =
                   ( (uint64_t) ( _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                         ( 2 ) )] ) );
               _data_colIdx[5] = ( (uint64_t) ( _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                                 ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_colIdx[6] = ( (uint64_t) ( _data_src_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );
               _data_colIdx[7] = ( (uint64_t) ( _data_src_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[8] = ( (uint64_t) ( _data_src_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[9] =
                   ( (uint64_t) ( _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );
               _data_colIdx[10] =
                   ( (uint64_t) ( _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                         ( 2 ) )] ) );
               _data_colIdx[11] = ( (uint64_t) ( _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                                  ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] ) );

               /* Apply basis transformation */

               _data_mat[0]   = ( (real_t) ( elMat_0_0 ) );
               _data_mat[1]   = ( (real_t) ( elMat_0_1 ) );
               _data_mat[2]   = ( (real_t) ( elMat_0_2 ) );
               _data_mat[3]   = ( (real_t) ( elMat_0_3 ) );
               _data_mat[4]   = ( (real_t) ( elMat_0_4 ) );
               _data_mat[5]   = ( (real_t) ( elMat_0_5 ) );
               _data_mat[6]   = ( (real_t) ( elMat_0_6 ) );
               _data_mat[7]   = ( (real_t) ( elMat_0_7 ) );
               _data_mat[8]   = ( (real_t) ( elMat_0_8 ) );
               _data_mat[9]   = ( (real_t) ( elMat_0_9 ) );
               _data_mat[10]  = ( (real_t) ( elMat_0_10 ) );
               _data_mat[11]  = ( (real_t) ( elMat_0_11 ) );
               _data_mat[12]  = ( (real_t) ( elMat_1_0 ) );
               _data_mat[13]  = ( (real_t) ( elMat_1_1 ) );
               _data_mat[14]  = ( (real_t) ( elMat_1_2 ) );
               _data_mat[15]  = ( (real_t) ( elMat_1_3 ) );
               _data_mat[16]  = ( (real_t) ( elMat_1_4 ) );
               _data_mat[17]  = ( (real_t) ( elMat_1_5 ) );
               _data_mat[18]  = ( (real_t) ( elMat_1_6 ) );
               _data_mat[19]  = ( (real_t) ( elMat_1_7 ) );
               _data_mat[20]  = ( (real_t) ( elMat_1_8 ) );
               _data_mat[21]  = ( (real_t) ( elMat_1_9 ) );
               _data_mat[22]  = ( (real_t) ( elMat_1_10 ) );
               _data_mat[23]  = ( (real_t) ( elMat_1_11 ) );
               _data_mat[24]  = ( (real_t) ( elMat_2_0 ) );
               _data_mat[25]  = ( (real_t) ( elMat_2_1 ) );
               _data_mat[26]  = ( (real_t) ( elMat_2_2 ) );
               _data_mat[27]  = ( (real_t) ( elMat_2_3 ) );
               _data_mat[28]  = ( (real_t) ( elMat_2_4 ) );
               _data_mat[29]  = ( (real_t) ( elMat_2_5 ) );
               _data_mat[30]  = ( (real_t) ( elMat_2_6 ) );
               _data_mat[31]  = ( (real_t) ( elMat_2_7 ) );
               _data_mat[32]  = ( (real_t) ( elMat_2_8 ) );
               _data_mat[33]  = ( (real_t) ( elMat_2_9 ) );
               _data_mat[34]  = ( (real_t) ( elMat_2_10 ) );
               _data_mat[35]  = ( (real_t) ( elMat_2_11 ) );
               _data_mat[36]  = ( (real_t) ( elMat_3_0 ) );
               _data_mat[37]  = ( (real_t) ( elMat_3_1 ) );
               _data_mat[38]  = ( (real_t) ( elMat_3_2 ) );
               _data_mat[39]  = ( (real_t) ( elMat_3_3 ) );
               _data_mat[40]  = ( (real_t) ( elMat_3_4 ) );
               _data_mat[41]  = ( (real_t) ( elMat_3_5 ) );
               _data_mat[42]  = ( (real_t) ( elMat_3_6 ) );
               _data_mat[43]  = ( (real_t) ( elMat_3_7 ) );
               _data_mat[44]  = ( (real_t) ( elMat_3_8 ) );
               _data_mat[45]  = ( (real_t) ( elMat_3_9 ) );
               _data_mat[46]  = ( (real_t) ( elMat_3_10 ) );
               _data_mat[47]  = ( (real_t) ( elMat_3_11 ) );
               _data_mat[48]  = ( (real_t) ( elMat_4_0 ) );
               _data_mat[49]  = ( (real_t) ( elMat_4_1 ) );
               _data_mat[50]  = ( (real_t) ( elMat_4_2 ) );
               _data_mat[51]  = ( (real_t) ( elMat_4_3 ) );
               _data_mat[52]  = ( (real_t) ( elMat_4_4 ) );
               _data_mat[53]  = ( (real_t) ( elMat_4_5 ) );
               _data_mat[54]  = ( (real_t) ( elMat_4_6 ) );
               _data_mat[55]  = ( (real_t) ( elMat_4_7 ) );
               _data_mat[56]  = ( (real_t) ( elMat_4_8 ) );
               _data_mat[57]  = ( (real_t) ( elMat_4_9 ) );
               _data_mat[58]  = ( (real_t) ( elMat_4_10 ) );
               _data_mat[59]  = ( (real_t) ( elMat_4_11 ) );
               _data_mat[60]  = ( (real_t) ( elMat_5_0 ) );
               _data_mat[61]  = ( (real_t) ( elMat_5_1 ) );
               _data_mat[62]  = ( (real_t) ( elMat_5_2 ) );
               _data_mat[63]  = ( (real_t) ( elMat_5_3 ) );
               _data_mat[64]  = ( (real_t) ( elMat_5_4 ) );
               _data_mat[65]  = ( (real_t) ( elMat_5_5 ) );
               _data_mat[66]  = ( (real_t) ( elMat_5_6 ) );
               _data_mat[67]  = ( (real_t) ( elMat_5_7 ) );
               _data_mat[68]  = ( (real_t) ( elMat_5_8 ) );
               _data_mat[69]  = ( (real_t) ( elMat_5_9 ) );
               _data_mat[70]  = ( (real_t) ( elMat_5_10 ) );
               _data_mat[71]  = ( (real_t) ( elMat_5_11 ) );
               _data_mat[72]  = ( (real_t) ( elMat_6_0 ) );
               _data_mat[73]  = ( (real_t) ( elMat_6_1 ) );
               _data_mat[74]  = ( (real_t) ( elMat_6_2 ) );
               _data_mat[75]  = ( (real_t) ( elMat_6_3 ) );
               _data_mat[76]  = ( (real_t) ( elMat_6_4 ) );
               _data_mat[77]  = ( (real_t) ( elMat_6_5 ) );
               _data_mat[78]  = ( (real_t) ( elMat_6_6 ) );
               _data_mat[79]  = ( (real_t) ( elMat_6_7 ) );
               _data_mat[80]  = ( (real_t) ( elMat_6_8 ) );
               _data_mat[81]  = ( (real_t) ( elMat_6_9 ) );
               _data_mat[82]  = ( (real_t) ( elMat_6_10 ) );
               _data_mat[83]  = ( (real_t) ( elMat_6_11 ) );
               _data_mat[84]  = ( (real_t) ( elMat_7_0 ) );
               _data_mat[85]  = ( (real_t) ( elMat_7_1 ) );
               _data_mat[86]  = ( (real_t) ( elMat_7_2 ) );
               _data_mat[87]  = ( (real_t) ( elMat_7_3 ) );
               _data_mat[88]  = ( (real_t) ( elMat_7_4 ) );
               _data_mat[89]  = ( (real_t) ( elMat_7_5 ) );
               _data_mat[90]  = ( (real_t) ( elMat_7_6 ) );
               _data_mat[91]  = ( (real_t) ( elMat_7_7 ) );
               _data_mat[92]  = ( (real_t) ( elMat_7_8 ) );
               _data_mat[93]  = ( (real_t) ( elMat_7_9 ) );
               _data_mat[94]  = ( (real_t) ( elMat_7_10 ) );
               _data_mat[95]  = ( (real_t) ( elMat_7_11 ) );
               _data_mat[96]  = ( (real_t) ( elMat_8_0 ) );
               _data_mat[97]  = ( (real_t) ( elMat_8_1 ) );
               _data_mat[98]  = ( (real_t) ( elMat_8_2 ) );
               _data_mat[99]  = ( (real_t) ( elMat_8_3 ) );
               _data_mat[100] = ( (real_t) ( elMat_8_4 ) );
               _data_mat[101] = ( (real_t) ( elMat_8_5 ) );
               _data_mat[102] = ( (real_t) ( elMat_8_6 ) );
               _data_mat[103] = ( (real_t) ( elMat_8_7 ) );
               _data_mat[104] = ( (real_t) ( elMat_8_8 ) );
               _data_mat[105] = ( (real_t) ( elMat_8_9 ) );
               _data_mat[106] = ( (real_t) ( elMat_8_10 ) );
               _data_mat[107] = ( (real_t) ( elMat_8_11 ) );
               _data_mat[108] = ( (real_t) ( elMat_9_0 ) );
               _data_mat[109] = ( (real_t) ( elMat_9_1 ) );
               _data_mat[110] = ( (real_t) ( elMat_9_2 ) );
               _data_mat[111] = ( (real_t) ( elMat_9_3 ) );
               _data_mat[112] = ( (real_t) ( elMat_9_4 ) );
               _data_mat[113] = ( (real_t) ( elMat_9_5 ) );
               _data_mat[114] = ( (real_t) ( elMat_9_6 ) );
               _data_mat[115] = ( (real_t) ( elMat_9_7 ) );
               _data_mat[116] = ( (real_t) ( elMat_9_8 ) );
               _data_mat[117] = ( (real_t) ( elMat_9_9 ) );
               _data_mat[118] = ( (real_t) ( elMat_9_10 ) );
               _data_mat[119] = ( (real_t) ( elMat_9_11 ) );
               _data_mat[120] = ( (real_t) ( elMat_10_0 ) );
               _data_mat[121] = ( (real_t) ( elMat_10_1 ) );
               _data_mat[122] = ( (real_t) ( elMat_10_2 ) );
               _data_mat[123] = ( (real_t) ( elMat_10_3 ) );
               _data_mat[124] = ( (real_t) ( elMat_10_4 ) );
               _data_mat[125] = ( (real_t) ( elMat_10_5 ) );
               _data_mat[126] = ( (real_t) ( elMat_10_6 ) );
               _data_mat[127] = ( (real_t) ( elMat_10_7 ) );
               _data_mat[128] = ( (real_t) ( elMat_10_8 ) );
               _data_mat[129] = ( (real_t) ( elMat_10_9 ) );
               _data_mat[130] = ( (real_t) ( elMat_10_10 ) );
               _data_mat[131] = ( (real_t) ( elMat_10_11 ) );
               _data_mat[132] = ( (real_t) ( elMat_11_0 ) );
               _data_mat[133] = ( (real_t) ( elMat_11_1 ) );
               _data_mat[134] = ( (real_t) ( elMat_11_2 ) );
               _data_mat[135] = ( (real_t) ( elMat_11_3 ) );
               _data_mat[136] = ( (real_t) ( elMat_11_4 ) );
               _data_mat[137] = ( (real_t) ( elMat_11_5 ) );
               _data_mat[138] = ( (real_t) ( elMat_11_6 ) );
               _data_mat[139] = ( (real_t) ( elMat_11_7 ) );
               _data_mat[140] = ( (real_t) ( elMat_11_8 ) );
               _data_mat[141] = ( (real_t) ( elMat_11_9 ) );
               _data_mat[142] = ( (real_t) ( elMat_11_10 ) );
               _data_mat[143] = ( (real_t) ( elMat_11_11 ) );

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
             jac_affine_0_0_BLUE * jac_affine_1_1_BLUE - jac_affine_0_1_BLUE * jac_affine_1_0_BLUE;
         const walberla::float64 tmp_coords_jac_6_BLUE   = 1.0 / ( tmp_coords_jac_5_BLUE );
         const walberla::float64 jac_affine_inv_0_0_BLUE = jac_affine_1_1_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_0_1_BLUE = -jac_affine_0_1_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_1_0_BLUE = -jac_affine_1_0_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_1_1_BLUE = jac_affine_0_0_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 abs_det_jac_affine_BLUE = abs( tmp_coords_jac_5_BLUE );
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
               const walberla::float64 uy_dof_1    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_2    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               walberla::float64       q_acc_0_0   = 0.0;
               walberla::float64       q_acc_0_1   = 0.0;
               walberla::float64       q_acc_0_2   = 0.0;
               walberla::float64       q_acc_0_3   = 0.0;
               walberla::float64       q_acc_0_4   = 0.0;
               walberla::float64       q_acc_0_5   = 0.0;
               walberla::float64       q_acc_0_6   = 0.0;
               walberla::float64       q_acc_0_7   = 0.0;
               walberla::float64       q_acc_0_8   = 0.0;
               walberla::float64       q_acc_0_9   = 0.0;
               walberla::float64       q_acc_0_10  = 0.0;
               walberla::float64       q_acc_0_11  = 0.0;
               walberla::float64       q_acc_1_1   = 0.0;
               walberla::float64       q_acc_1_2   = 0.0;
               walberla::float64       q_acc_1_3   = 0.0;
               walberla::float64       q_acc_1_4   = 0.0;
               walberla::float64       q_acc_1_5   = 0.0;
               walberla::float64       q_acc_1_6   = 0.0;
               walberla::float64       q_acc_1_7   = 0.0;
               walberla::float64       q_acc_1_8   = 0.0;
               walberla::float64       q_acc_1_9   = 0.0;
               walberla::float64       q_acc_1_10  = 0.0;
               walberla::float64       q_acc_1_11  = 0.0;
               walberla::float64       q_acc_2_2   = 0.0;
               walberla::float64       q_acc_2_3   = 0.0;
               walberla::float64       q_acc_2_4   = 0.0;
               walberla::float64       q_acc_2_5   = 0.0;
               walberla::float64       q_acc_2_6   = 0.0;
               walberla::float64       q_acc_2_7   = 0.0;
               walberla::float64       q_acc_2_8   = 0.0;
               walberla::float64       q_acc_2_9   = 0.0;
               walberla::float64       q_acc_2_10  = 0.0;
               walberla::float64       q_acc_2_11  = 0.0;
               walberla::float64       q_acc_3_3   = 0.0;
               walberla::float64       q_acc_3_4   = 0.0;
               walberla::float64       q_acc_3_5   = 0.0;
               walberla::float64       q_acc_3_6   = 0.0;
               walberla::float64       q_acc_3_7   = 0.0;
               walberla::float64       q_acc_3_8   = 0.0;
               walberla::float64       q_acc_3_9   = 0.0;
               walberla::float64       q_acc_3_10  = 0.0;
               walberla::float64       q_acc_3_11  = 0.0;
               walberla::float64       q_acc_4_4   = 0.0;
               walberla::float64       q_acc_4_5   = 0.0;
               walberla::float64       q_acc_4_6   = 0.0;
               walberla::float64       q_acc_4_7   = 0.0;
               walberla::float64       q_acc_4_8   = 0.0;
               walberla::float64       q_acc_4_9   = 0.0;
               walberla::float64       q_acc_4_10  = 0.0;
               walberla::float64       q_acc_4_11  = 0.0;
               walberla::float64       q_acc_5_5   = 0.0;
               walberla::float64       q_acc_5_6   = 0.0;
               walberla::float64       q_acc_5_7   = 0.0;
               walberla::float64       q_acc_5_8   = 0.0;
               walberla::float64       q_acc_5_9   = 0.0;
               walberla::float64       q_acc_5_10  = 0.0;
               walberla::float64       q_acc_5_11  = 0.0;
               walberla::float64       q_acc_6_6   = 0.0;
               walberla::float64       q_acc_6_7   = 0.0;
               walberla::float64       q_acc_6_8   = 0.0;
               walberla::float64       q_acc_6_9   = 0.0;
               walberla::float64       q_acc_6_10  = 0.0;
               walberla::float64       q_acc_6_11  = 0.0;
               walberla::float64       q_acc_7_7   = 0.0;
               walberla::float64       q_acc_7_8   = 0.0;
               walberla::float64       q_acc_7_9   = 0.0;
               walberla::float64       q_acc_7_10  = 0.0;
               walberla::float64       q_acc_7_11  = 0.0;
               walberla::float64       q_acc_8_8   = 0.0;
               walberla::float64       q_acc_8_9   = 0.0;
               walberla::float64       q_acc_8_10  = 0.0;
               walberla::float64       q_acc_8_11  = 0.0;
               walberla::float64       q_acc_9_9   = 0.0;
               walberla::float64       q_acc_9_10  = 0.0;
               walberla::float64       q_acc_9_11  = 0.0;
               walberla::float64       q_acc_10_10 = 0.0;
               walberla::float64       q_acc_10_11 = 0.0;
               walberla::float64       q_acc_11_11 = 0.0;
               for ( int64_t q = 0; q < 3; q += 1 )
               {
                  const walberla::float64 tmp_qloop_0  = 4.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_1  = 4.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_2  = tmp_qloop_0 + tmp_qloop_1 - 3.0;
                  const walberla::float64 tmp_qloop_3  = jac_affine_inv_0_0_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_4  = jac_affine_inv_1_0_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_5  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_6  = abs_det_jac_affine_BLUE * tmp_qloop_5;
                  const walberla::float64 tmp_qloop_7  = jac_affine_inv_0_1_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_8  = jac_affine_inv_1_1_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_9  = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_10 = tmp_qloop_7 + tmp_qloop_8;
                  const walberla::float64 tmp_qloop_11 = -ux_dof_0 + ux_dof_1;
                  const walberla::float64 tmp_qloop_12 = -ux_dof_0 + ux_dof_2;
                  const walberla::float64 tmp_qloop_13 = -uy_dof_0 + uy_dof_1;
                  const walberla::float64 tmp_qloop_14 = -uy_dof_0 + uy_dof_2;
                  const walberla::float64 tmp_qloop_15 = jac_affine_inv_0_0_BLUE * tmp_qloop_13;
                  const walberla::float64 tmp_qloop_16 = jac_affine_inv_0_1_BLUE * tmp_qloop_11;
                  const walberla::float64 tmp_qloop_17 = jac_affine_inv_1_0_BLUE * tmp_qloop_14;
                  const walberla::float64 tmp_qloop_18 = jac_affine_inv_1_1_BLUE * tmp_qloop_12;
                  const walberla::float64 tmp_qloop_19 =
                      1.0 /
                      ( 1.0 /
                            ( mu_star +
                              sigma_y * 1.0 /
                                  ( pow( ( ( jac_affine_inv_0_0_BLUE * tmp_qloop_11 + jac_affine_inv_1_0_BLUE * tmp_qloop_12 ) *
                                           ( jac_affine_inv_0_0_BLUE * tmp_qloop_11 + jac_affine_inv_1_0_BLUE * tmp_qloop_12 ) ) +
                                             ( ( jac_affine_inv_0_1_BLUE * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_BLUE * tmp_qloop_14 ) *
                                               ( jac_affine_inv_0_1_BLUE * tmp_qloop_13 +
                                                 jac_affine_inv_1_1_BLUE * tmp_qloop_14 ) ) +
                                             ( tmp_qloop_15 + tmp_qloop_16 + tmp_qloop_17 + tmp_qloop_18 ) *
                                                 ( tmp_qloop_15 * 0.5 + tmp_qloop_16 * 0.5 + tmp_qloop_17 * 0.5 +
                                                   tmp_qloop_18 * 0.5 ),
                                         0.50000000000000000 ) +
                                    9.9999999999999995e-21 ) ) *
                            1.0 +
                        1.0 /
                            ( mu_lin_dof_0 * ( 1.0 - _data_q_p_0[q] - _data_q_p_1[q] ) + mu_lin_dof_1 * _data_q_p_0[q] +
                              mu_lin_dof_2 * _data_q_p_1[q] ) *
                            1.0 ) *
                      _data_q_w[q];
                  const walberla::float64 tmp_qloop_20 = tmp_qloop_19 * 2.0;
                  const walberla::float64 tmp_qloop_21 = tmp_qloop_0 - 1.0;
                  const walberla::float64 tmp_qloop_22 = jac_affine_inv_0_0_BLUE * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_23 = tmp_qloop_6 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_24 = tmp_qloop_5 * 2.0;
                  const walberla::float64 tmp_qloop_25 = jac_affine_inv_0_1_BLUE * tmp_qloop_21;
                  const walberla::float64 tmp_qloop_26 = tmp_qloop_1 - 1.0;
                  const walberla::float64 tmp_qloop_27 = jac_affine_inv_1_0_BLUE * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_28 = jac_affine_inv_1_1_BLUE * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_29 = 2.6666666666666665 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_30 = jac_affine_inv_1_0_BLUE * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_31 = 2.6666666666666665 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_32 = jac_affine_inv_0_0_BLUE * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_33 = tmp_qloop_30 + tmp_qloop_32;
                  const walberla::float64 tmp_qloop_34 = jac_affine_inv_1_0_BLUE * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_35 = jac_affine_inv_0_0_BLUE * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_36 = tmp_qloop_34 + tmp_qloop_35;
                  const walberla::float64 tmp_qloop_37 = jac_affine_inv_1_1_BLUE * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_38 = jac_affine_inv_0_1_BLUE * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_39 = tmp_qloop_37 + tmp_qloop_38;
                  const walberla::float64 tmp_qloop_40 = -tmp_qloop_0 - 8.0 * _data_q_p_1[q] + 4.0;
                  const walberla::float64 tmp_qloop_41 =
                      jac_affine_inv_1_0_BLUE * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_32;
                  const walberla::float64 tmp_qloop_42 = jac_affine_inv_1_0_BLUE * tmp_qloop_40 - tmp_qloop_35;
                  const walberla::float64 tmp_qloop_43 = jac_affine_inv_1_1_BLUE * tmp_qloop_40 - tmp_qloop_38;
                  const walberla::float64 tmp_qloop_44 = -tmp_qloop_1 - 8.0 * _data_q_p_0[q] + 4.0;
                  const walberla::float64 tmp_qloop_45 =
                      jac_affine_inv_0_0_BLUE * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_30;
                  const walberla::float64 tmp_qloop_46 = jac_affine_inv_0_0_BLUE * tmp_qloop_44 - tmp_qloop_34;
                  const walberla::float64 tmp_qloop_47 = jac_affine_inv_0_1_BLUE * tmp_qloop_44 - tmp_qloop_37;
                  const walberla::float64 tmp_qloop_48 = tmp_qloop_7 * 0.66666666666666663 + tmp_qloop_8 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_49 = abs_det_jac_affine_BLUE * tmp_qloop_22;
                  const walberla::float64 tmp_qloop_50 = abs_det_jac_affine_BLUE * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_51 = jac_affine_inv_1_1_BLUE * tmp_qloop_29;
                  const walberla::float64 tmp_qloop_52 = jac_affine_inv_0_1_BLUE * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_53 = tmp_qloop_51 + tmp_qloop_52;
                  const walberla::float64 tmp_qloop_54 = abs_det_jac_affine_BLUE * tmp_qloop_36;
                  const walberla::float64 tmp_qloop_55 =
                      jac_affine_inv_1_1_BLUE * tmp_qloop_40 * 0.66666666666666663 - tmp_qloop_52;
                  const walberla::float64 tmp_qloop_56 = abs_det_jac_affine_BLUE * tmp_qloop_42;
                  const walberla::float64 tmp_qloop_57 =
                      jac_affine_inv_0_1_BLUE * tmp_qloop_44 * 0.66666666666666663 - tmp_qloop_51;
                  const walberla::float64 tmp_qloop_58  = abs_det_jac_affine_BLUE * tmp_qloop_46;
                  const walberla::float64 tmp_qloop_59  = ( jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE );
                  const walberla::float64 tmp_qloop_60  = ( tmp_qloop_21 * tmp_qloop_21 );
                  const walberla::float64 tmp_qloop_61  = abs_det_jac_affine_BLUE * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_62  = tmp_qloop_60 * tmp_qloop_61;
                  const walberla::float64 tmp_qloop_63  = tmp_qloop_60 * 2.0;
                  const walberla::float64 tmp_qloop_64  = ( jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE );
                  const walberla::float64 tmp_qloop_65  = tmp_qloop_60 * 1.0;
                  const walberla::float64 tmp_qloop_66  = tmp_qloop_22 * tmp_qloop_27;
                  const walberla::float64 tmp_qloop_67  = tmp_qloop_25 * 1.0;
                  const walberla::float64 tmp_qloop_68  = tmp_qloop_36 * 2.0;
                  const walberla::float64 tmp_qloop_69  = tmp_qloop_42 * 2.0;
                  const walberla::float64 tmp_qloop_70  = tmp_qloop_46 * 2.0;
                  const walberla::float64 tmp_qloop_71  = abs_det_jac_affine_BLUE * tmp_qloop_19 * 0.66666666666666674;
                  const walberla::float64 tmp_qloop_72  = abs_det_jac_affine_BLUE * tmp_qloop_57;
                  const walberla::float64 tmp_qloop_73  = ( jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE );
                  const walberla::float64 tmp_qloop_74  = ( tmp_qloop_26 * tmp_qloop_26 );
                  const walberla::float64 tmp_qloop_75  = tmp_qloop_61 * tmp_qloop_74;
                  const walberla::float64 tmp_qloop_76  = tmp_qloop_74 * 2.0;
                  const walberla::float64 tmp_qloop_77  = ( jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE );
                  const walberla::float64 tmp_qloop_78  = tmp_qloop_74 * 1.0;
                  const walberla::float64 tmp_qloop_79  = tmp_qloop_28 * 1.0;
                  const walberla::float64 tmp_qloop_80  = 2.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_81  = jac_affine_inv_1_1_BLUE * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_82  = 2.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_83  = jac_affine_inv_0_1_BLUE * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_84  = tmp_qloop_81 + tmp_qloop_83;
                  const walberla::float64 tmp_qloop_85  = tmp_qloop_84 * 2.0;
                  const walberla::float64 tmp_qloop_86  = tmp_qloop_54 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_87  = jac_affine_inv_1_1_BLUE * tmp_qloop_40 * 0.5 - tmp_qloop_83;
                  const walberla::float64 tmp_qloop_88  = tmp_qloop_87 * 2.0;
                  const walberla::float64 tmp_qloop_89  = tmp_qloop_56 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_90  = jac_affine_inv_0_1_BLUE * tmp_qloop_44 * 0.5 - tmp_qloop_81;
                  const walberla::float64 tmp_qloop_91  = tmp_qloop_90 * 2.0;
                  const walberla::float64 tmp_qloop_92  = tmp_qloop_58 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_93  = abs_det_jac_affine_BLUE * tmp_qloop_9;
                  const walberla::float64 tmp_qloop_94  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_95  = tmp_qloop_93 * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_96  = tmp_qloop_9 * 2.0;
                  const walberla::float64 tmp_qloop_97  = tmp_qloop_25 * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_98  = tmp_qloop_22 * 1.0;
                  const walberla::float64 tmp_qloop_99  = abs_det_jac_affine_BLUE * tmp_qloop_25;
                  const walberla::float64 tmp_qloop_100 = tmp_qloop_39 * 2.0;
                  const walberla::float64 tmp_qloop_101 = tmp_qloop_43 * 2.0;
                  const walberla::float64 tmp_qloop_102 = tmp_qloop_47 * 2.0;
                  const walberla::float64 tmp_qloop_103 = abs_det_jac_affine_BLUE * tmp_qloop_28;
                  const walberla::float64 tmp_qloop_104 = tmp_qloop_27 * 1.0;
                  const walberla::float64 tmp_qloop_105 = abs_det_jac_affine_BLUE * tmp_qloop_39;
                  const walberla::float64 tmp_qloop_106 = jac_affine_inv_1_0_BLUE * tmp_qloop_80;
                  const walberla::float64 tmp_qloop_107 = jac_affine_inv_0_0_BLUE * tmp_qloop_82;
                  const walberla::float64 tmp_qloop_108 = tmp_qloop_106 * 2.0 + tmp_qloop_107 * 2.0;
                  const walberla::float64 tmp_qloop_109 = abs_det_jac_affine_BLUE * tmp_qloop_43;
                  const walberla::float64 tmp_qloop_110 = jac_affine_inv_1_0_BLUE * tmp_qloop_40 + tmp_qloop_107 * -2.0;
                  const walberla::float64 q_tmp_0_0 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_9 + ( tmp_qloop_5 * tmp_qloop_5 ) * 2.0 ) -
                        tmp_qloop_6 * ( tmp_qloop_3 * 0.66666666666666663 + tmp_qloop_4 * 0.66666666666666663 ) );
                  const walberla::float64 q_tmp_0_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_25 + tmp_qloop_22 * tmp_qloop_24 ) -
                                       tmp_qloop_22 * tmp_qloop_23 );
                  const walberla::float64 q_tmp_0_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_28 + tmp_qloop_24 * tmp_qloop_27 ) -
                                       tmp_qloop_23 * tmp_qloop_27 );
                  const walberla::float64 q_tmp_0_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_39 + tmp_qloop_24 * tmp_qloop_36 ) -
                                       tmp_qloop_33 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_43 + tmp_qloop_24 * tmp_qloop_42 ) -
                                       tmp_qloop_41 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_10 * tmp_qloop_47 + tmp_qloop_24 * tmp_qloop_46 ) -
                                       tmp_qloop_45 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_6 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_6 - tmp_qloop_48 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_7 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_49 - tmp_qloop_23 * tmp_qloop_25 );
                  const walberla::float64 q_tmp_0_8 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_50 - tmp_qloop_23 * tmp_qloop_28 );
                  const walberla::float64 q_tmp_0_9 = tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_54 - tmp_qloop_53 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_10 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_56 - tmp_qloop_55 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_0_11 =
                      tmp_qloop_20 * ( tmp_qloop_10 * tmp_qloop_58 - tmp_qloop_57 * tmp_qloop_6 );
                  const walberla::float64 q_tmp_1_1 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_59 * tmp_qloop_63 + tmp_qloop_64 * tmp_qloop_65 ) -
                                       tmp_qloop_59 * tmp_qloop_62 );
                  const walberla::float64 q_tmp_1_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_28 * tmp_qloop_67 + tmp_qloop_66 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_66 );
                  const walberla::float64 q_tmp_1_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_67 ) -
                                       tmp_qloop_33 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_67 ) -
                                       tmp_qloop_41 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_67 ) -
                                       tmp_qloop_45 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_49 );
                  const walberla::float64 q_tmp_1_7 =
                      jac_affine_inv_0_0_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_60 * tmp_qloop_71;
                  const walberla::float64 q_tmp_1_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * jac_affine_inv_1_0_BLUE *
                                           tmp_qloop_21 * tmp_qloop_26 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_28 * tmp_qloop_61 );
                  const walberla::float64 q_tmp_1_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_1_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_49 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_1_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_0_1_BLUE * tmp_qloop_21 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_22 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_2_2 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_73 * tmp_qloop_76 + tmp_qloop_77 * tmp_qloop_78 ) -
                                       tmp_qloop_73 * tmp_qloop_75 );
                  const walberla::float64 q_tmp_2_3 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_68 + tmp_qloop_39 * tmp_qloop_79 ) -
                                       tmp_qloop_33 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_69 + tmp_qloop_43 * tmp_qloop_79 ) -
                                       tmp_qloop_41 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_70 + tmp_qloop_47 * tmp_qloop_79 ) -
                                       tmp_qloop_45 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_5 * 1.0 -
                                       tmp_qloop_48 * tmp_qloop_50 );
                  const walberla::float64 q_tmp_2_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_27 * tmp_qloop_61 + tmp_qloop_49 * tmp_qloop_79 );
                  const walberla::float64 q_tmp_2_8 =
                      jac_affine_inv_1_0_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_71 * tmp_qloop_74;
                  const walberla::float64 q_tmp_2_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_36 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_2_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_42 * 1.0 -
                                       tmp_qloop_50 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_2_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * jac_affine_inv_1_1_BLUE * tmp_qloop_26 * tmp_qloop_46 * 1.0 -
                                       tmp_qloop_27 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_3_3 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( ( tmp_qloop_36 * tmp_qloop_36 ) * 2.0 + tmp_qloop_39 * tmp_qloop_85 ) -
                        tmp_qloop_33 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_4 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_42 * tmp_qloop_68 + tmp_qloop_43 * tmp_qloop_85 ) -
                                       tmp_qloop_41 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_46 * tmp_qloop_68 + tmp_qloop_47 * tmp_qloop_85 ) -
                                       tmp_qloop_45 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * tmp_qloop_5 * tmp_qloop_84 * 2.0 - tmp_qloop_48 * tmp_qloop_54 );
                  const walberla::float64 q_tmp_3_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_86 + tmp_qloop_49 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_86 + tmp_qloop_50 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_54 + tmp_qloop_54 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_3_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_42 * tmp_qloop_84 * 2.0 - tmp_qloop_54 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_3_11 =
                      tmp_qloop_20 * ( -tmp_qloop_54 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_85 );
                  const walberla::float64 q_tmp_4_4 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( ( tmp_qloop_42 * tmp_qloop_42 ) * 2.0 + tmp_qloop_43 * tmp_qloop_88 ) -
                        tmp_qloop_41 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_5 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_46 * tmp_qloop_69 + tmp_qloop_47 * tmp_qloop_88 ) -
                                       tmp_qloop_45 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * tmp_qloop_5 * tmp_qloop_87 * 2.0 - tmp_qloop_48 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_89 + tmp_qloop_49 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_89 + tmp_qloop_50 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_9 =
                      tmp_qloop_20 * ( -tmp_qloop_53 * tmp_qloop_56 + tmp_qloop_54 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_4_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_42 * tmp_qloop_87 * 2.0 - tmp_qloop_55 * tmp_qloop_56 );
                  const walberla::float64 q_tmp_4_11 =
                      tmp_qloop_20 * ( -tmp_qloop_56 * tmp_qloop_57 + tmp_qloop_58 * tmp_qloop_88 );
                  const walberla::float64 q_tmp_5_5 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( ( tmp_qloop_46 * tmp_qloop_46 ) * 2.0 + tmp_qloop_47 * tmp_qloop_91 ) -
                        tmp_qloop_45 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_6 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * tmp_qloop_5 * tmp_qloop_90 * 2.0 - tmp_qloop_48 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_7 =
                      tmp_qloop_20 * ( -tmp_qloop_25 * tmp_qloop_92 + tmp_qloop_49 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_8 =
                      tmp_qloop_20 * ( -tmp_qloop_28 * tmp_qloop_92 + tmp_qloop_50 * tmp_qloop_91 );
                  const walberla::float64 q_tmp_5_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_36 * tmp_qloop_90 * 2.0 - tmp_qloop_53 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_42 * tmp_qloop_90 * 2.0 - tmp_qloop_55 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_5_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * tmp_qloop_46 * tmp_qloop_90 * 2.0 - tmp_qloop_57 * tmp_qloop_58 );
                  const walberla::float64 q_tmp_6_6 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_5 * tmp_qloop_94 + ( tmp_qloop_9 * tmp_qloop_9 ) * 2.0 ) -
                        tmp_qloop_48 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_22 * tmp_qloop_94 + tmp_qloop_25 * tmp_qloop_96 ) -
                                       tmp_qloop_25 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_94 + tmp_qloop_28 * tmp_qloop_96 ) -
                                       tmp_qloop_28 * tmp_qloop_95 );
                  const walberla::float64 q_tmp_6_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_36 * tmp_qloop_94 + tmp_qloop_39 * tmp_qloop_96 ) -
                                       tmp_qloop_53 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_42 * tmp_qloop_94 + tmp_qloop_43 * tmp_qloop_96 ) -
                                       tmp_qloop_55 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_6_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_46 * tmp_qloop_94 + tmp_qloop_47 * tmp_qloop_96 ) -
                                       tmp_qloop_57 * tmp_qloop_93 );
                  const walberla::float64 q_tmp_7_7 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_59 * tmp_qloop_65 + tmp_qloop_63 * tmp_qloop_64 ) -
                                       tmp_qloop_62 * tmp_qloop_64 );
                  const walberla::float64 q_tmp_7_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_27 * tmp_qloop_98 + tmp_qloop_97 * 2.0 ) -
                                       tmp_qloop_61 * tmp_qloop_97 );
                  const walberla::float64 q_tmp_7_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_25 + tmp_qloop_36 * tmp_qloop_98 ) -
                                       tmp_qloop_53 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_101 * tmp_qloop_25 + tmp_qloop_42 * tmp_qloop_98 ) -
                                       tmp_qloop_55 * tmp_qloop_99 );
                  const walberla::float64 q_tmp_7_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_102 * tmp_qloop_25 + tmp_qloop_46 * tmp_qloop_98 ) -
                                       tmp_qloop_25 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_8_8 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_73 * tmp_qloop_78 + tmp_qloop_76 * tmp_qloop_77 ) -
                                       tmp_qloop_75 * tmp_qloop_77 );
                  const walberla::float64 q_tmp_8_9 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_36 ) -
                                       tmp_qloop_103 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_8_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_101 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_42 ) -
                                       tmp_qloop_103 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_8_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_102 * tmp_qloop_28 + tmp_qloop_104 * tmp_qloop_46 ) -
                                       tmp_qloop_28 * tmp_qloop_72 );
                  const walberla::float64 q_tmp_9_9 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_108 * tmp_qloop_36 + ( tmp_qloop_39 * tmp_qloop_39 ) * 2.0 ) -
                        tmp_qloop_105 * tmp_qloop_53 );
                  const walberla::float64 q_tmp_9_10 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_43 + tmp_qloop_108 * tmp_qloop_42 ) -
                                       tmp_qloop_105 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_9_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_100 * tmp_qloop_47 + tmp_qloop_108 * tmp_qloop_46 ) -
                                       tmp_qloop_105 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_10_10 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE * ( tmp_qloop_110 * tmp_qloop_42 + ( tmp_qloop_43 * tmp_qloop_43 ) * 2.0 ) -
                        tmp_qloop_109 * tmp_qloop_55 );
                  const walberla::float64 q_tmp_10_11 =
                      tmp_qloop_20 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_101 * tmp_qloop_47 + tmp_qloop_110 * tmp_qloop_46 ) -
                                       tmp_qloop_109 * tmp_qloop_57 );
                  const walberla::float64 q_tmp_11_11 =
                      tmp_qloop_20 *
                      ( abs_det_jac_affine_BLUE *
                            ( tmp_qloop_46 * ( jac_affine_inv_0_0_BLUE * tmp_qloop_44 * 0.5 - tmp_qloop_106 ) * 2.0 +
                              ( tmp_qloop_47 * tmp_qloop_47 ) * 2.0 ) -
                        tmp_qloop_47 * tmp_qloop_72 );
                  q_acc_0_0   = q_acc_0_0 + q_tmp_0_0;
                  q_acc_0_1   = q_acc_0_1 + q_tmp_0_1;
                  q_acc_0_2   = q_acc_0_2 + q_tmp_0_2;
                  q_acc_0_3   = q_acc_0_3 + q_tmp_0_3;
                  q_acc_0_4   = q_acc_0_4 + q_tmp_0_4;
                  q_acc_0_5   = q_acc_0_5 + q_tmp_0_5;
                  q_acc_0_6   = q_acc_0_6 + q_tmp_0_6;
                  q_acc_0_7   = q_acc_0_7 + q_tmp_0_7;
                  q_acc_0_8   = q_acc_0_8 + q_tmp_0_8;
                  q_acc_0_9   = q_acc_0_9 + q_tmp_0_9;
                  q_acc_0_10  = q_acc_0_10 + q_tmp_0_10;
                  q_acc_0_11  = q_acc_0_11 + q_tmp_0_11;
                  q_acc_1_1   = q_acc_1_1 + q_tmp_1_1;
                  q_acc_1_2   = q_acc_1_2 + q_tmp_1_2;
                  q_acc_1_3   = q_acc_1_3 + q_tmp_1_3;
                  q_acc_1_4   = q_acc_1_4 + q_tmp_1_4;
                  q_acc_1_5   = q_acc_1_5 + q_tmp_1_5;
                  q_acc_1_6   = q_acc_1_6 + q_tmp_1_6;
                  q_acc_1_7   = q_acc_1_7 + q_tmp_1_7;
                  q_acc_1_8   = q_acc_1_8 + q_tmp_1_8;
                  q_acc_1_9   = q_acc_1_9 + q_tmp_1_9;
                  q_acc_1_10  = q_acc_1_10 + q_tmp_1_10;
                  q_acc_1_11  = q_acc_1_11 + q_tmp_1_11;
                  q_acc_2_2   = q_acc_2_2 + q_tmp_2_2;
                  q_acc_2_3   = q_acc_2_3 + q_tmp_2_3;
                  q_acc_2_4   = q_acc_2_4 + q_tmp_2_4;
                  q_acc_2_5   = q_acc_2_5 + q_tmp_2_5;
                  q_acc_2_6   = q_acc_2_6 + q_tmp_2_6;
                  q_acc_2_7   = q_acc_2_7 + q_tmp_2_7;
                  q_acc_2_8   = q_acc_2_8 + q_tmp_2_8;
                  q_acc_2_9   = q_acc_2_9 + q_tmp_2_9;
                  q_acc_2_10  = q_acc_2_10 + q_tmp_2_10;
                  q_acc_2_11  = q_acc_2_11 + q_tmp_2_11;
                  q_acc_3_3   = q_acc_3_3 + q_tmp_3_3;
                  q_acc_3_4   = q_acc_3_4 + q_tmp_3_4;
                  q_acc_3_5   = q_acc_3_5 + q_tmp_3_5;
                  q_acc_3_6   = q_acc_3_6 + q_tmp_3_6;
                  q_acc_3_7   = q_acc_3_7 + q_tmp_3_7;
                  q_acc_3_8   = q_acc_3_8 + q_tmp_3_8;
                  q_acc_3_9   = q_acc_3_9 + q_tmp_3_9;
                  q_acc_3_10  = q_acc_3_10 + q_tmp_3_10;
                  q_acc_3_11  = q_acc_3_11 + q_tmp_3_11;
                  q_acc_4_4   = q_acc_4_4 + q_tmp_4_4;
                  q_acc_4_5   = q_acc_4_5 + q_tmp_4_5;
                  q_acc_4_6   = q_acc_4_6 + q_tmp_4_6;
                  q_acc_4_7   = q_acc_4_7 + q_tmp_4_7;
                  q_acc_4_8   = q_acc_4_8 + q_tmp_4_8;
                  q_acc_4_9   = q_acc_4_9 + q_tmp_4_9;
                  q_acc_4_10  = q_acc_4_10 + q_tmp_4_10;
                  q_acc_4_11  = q_acc_4_11 + q_tmp_4_11;
                  q_acc_5_5   = q_acc_5_5 + q_tmp_5_5;
                  q_acc_5_6   = q_acc_5_6 + q_tmp_5_6;
                  q_acc_5_7   = q_acc_5_7 + q_tmp_5_7;
                  q_acc_5_8   = q_acc_5_8 + q_tmp_5_8;
                  q_acc_5_9   = q_acc_5_9 + q_tmp_5_9;
                  q_acc_5_10  = q_acc_5_10 + q_tmp_5_10;
                  q_acc_5_11  = q_acc_5_11 + q_tmp_5_11;
                  q_acc_6_6   = q_acc_6_6 + q_tmp_6_6;
                  q_acc_6_7   = q_acc_6_7 + q_tmp_6_7;
                  q_acc_6_8   = q_acc_6_8 + q_tmp_6_8;
                  q_acc_6_9   = q_acc_6_9 + q_tmp_6_9;
                  q_acc_6_10  = q_acc_6_10 + q_tmp_6_10;
                  q_acc_6_11  = q_acc_6_11 + q_tmp_6_11;
                  q_acc_7_7   = q_acc_7_7 + q_tmp_7_7;
                  q_acc_7_8   = q_acc_7_8 + q_tmp_7_8;
                  q_acc_7_9   = q_acc_7_9 + q_tmp_7_9;
                  q_acc_7_10  = q_acc_7_10 + q_tmp_7_10;
                  q_acc_7_11  = q_acc_7_11 + q_tmp_7_11;
                  q_acc_8_8   = q_acc_8_8 + q_tmp_8_8;
                  q_acc_8_9   = q_acc_8_9 + q_tmp_8_9;
                  q_acc_8_10  = q_acc_8_10 + q_tmp_8_10;
                  q_acc_8_11  = q_acc_8_11 + q_tmp_8_11;
                  q_acc_9_9   = q_acc_9_9 + q_tmp_9_9;
                  q_acc_9_10  = q_acc_9_10 + q_tmp_9_10;
                  q_acc_9_11  = q_acc_9_11 + q_tmp_9_11;
                  q_acc_10_10 = q_acc_10_10 + q_tmp_10_10;
                  q_acc_10_11 = q_acc_10_11 + q_tmp_10_11;
                  q_acc_11_11 = q_acc_11_11 + q_tmp_11_11;
               }
               const walberla::float64 elMat_0_0   = q_acc_0_0;
               const walberla::float64 elMat_0_1   = q_acc_0_1;
               const walberla::float64 elMat_0_2   = q_acc_0_2;
               const walberla::float64 elMat_0_3   = q_acc_0_3;
               const walberla::float64 elMat_0_4   = q_acc_0_4;
               const walberla::float64 elMat_0_5   = q_acc_0_5;
               const walberla::float64 elMat_0_6   = q_acc_0_6;
               const walberla::float64 elMat_0_7   = q_acc_0_7;
               const walberla::float64 elMat_0_8   = q_acc_0_8;
               const walberla::float64 elMat_0_9   = q_acc_0_9;
               const walberla::float64 elMat_0_10  = q_acc_0_10;
               const walberla::float64 elMat_0_11  = q_acc_0_11;
               const walberla::float64 elMat_1_0   = q_acc_0_1;
               const walberla::float64 elMat_1_1   = q_acc_1_1;
               const walberla::float64 elMat_1_2   = q_acc_1_2;
               const walberla::float64 elMat_1_3   = q_acc_1_3;
               const walberla::float64 elMat_1_4   = q_acc_1_4;
               const walberla::float64 elMat_1_5   = q_acc_1_5;
               const walberla::float64 elMat_1_6   = q_acc_1_6;
               const walberla::float64 elMat_1_7   = q_acc_1_7;
               const walberla::float64 elMat_1_8   = q_acc_1_8;
               const walberla::float64 elMat_1_9   = q_acc_1_9;
               const walberla::float64 elMat_1_10  = q_acc_1_10;
               const walberla::float64 elMat_1_11  = q_acc_1_11;
               const walberla::float64 elMat_2_0   = q_acc_0_2;
               const walberla::float64 elMat_2_1   = q_acc_1_2;
               const walberla::float64 elMat_2_2   = q_acc_2_2;
               const walberla::float64 elMat_2_3   = q_acc_2_3;
               const walberla::float64 elMat_2_4   = q_acc_2_4;
               const walberla::float64 elMat_2_5   = q_acc_2_5;
               const walberla::float64 elMat_2_6   = q_acc_2_6;
               const walberla::float64 elMat_2_7   = q_acc_2_7;
               const walberla::float64 elMat_2_8   = q_acc_2_8;
               const walberla::float64 elMat_2_9   = q_acc_2_9;
               const walberla::float64 elMat_2_10  = q_acc_2_10;
               const walberla::float64 elMat_2_11  = q_acc_2_11;
               const walberla::float64 elMat_3_0   = q_acc_0_3;
               const walberla::float64 elMat_3_1   = q_acc_1_3;
               const walberla::float64 elMat_3_2   = q_acc_2_3;
               const walberla::float64 elMat_3_3   = q_acc_3_3;
               const walberla::float64 elMat_3_4   = q_acc_3_4;
               const walberla::float64 elMat_3_5   = q_acc_3_5;
               const walberla::float64 elMat_3_6   = q_acc_3_6;
               const walberla::float64 elMat_3_7   = q_acc_3_7;
               const walberla::float64 elMat_3_8   = q_acc_3_8;
               const walberla::float64 elMat_3_9   = q_acc_3_9;
               const walberla::float64 elMat_3_10  = q_acc_3_10;
               const walberla::float64 elMat_3_11  = q_acc_3_11;
               const walberla::float64 elMat_4_0   = q_acc_0_4;
               const walberla::float64 elMat_4_1   = q_acc_1_4;
               const walberla::float64 elMat_4_2   = q_acc_2_4;
               const walberla::float64 elMat_4_3   = q_acc_3_4;
               const walberla::float64 elMat_4_4   = q_acc_4_4;
               const walberla::float64 elMat_4_5   = q_acc_4_5;
               const walberla::float64 elMat_4_6   = q_acc_4_6;
               const walberla::float64 elMat_4_7   = q_acc_4_7;
               const walberla::float64 elMat_4_8   = q_acc_4_8;
               const walberla::float64 elMat_4_9   = q_acc_4_9;
               const walberla::float64 elMat_4_10  = q_acc_4_10;
               const walberla::float64 elMat_4_11  = q_acc_4_11;
               const walberla::float64 elMat_5_0   = q_acc_0_5;
               const walberla::float64 elMat_5_1   = q_acc_1_5;
               const walberla::float64 elMat_5_2   = q_acc_2_5;
               const walberla::float64 elMat_5_3   = q_acc_3_5;
               const walberla::float64 elMat_5_4   = q_acc_4_5;
               const walberla::float64 elMat_5_5   = q_acc_5_5;
               const walberla::float64 elMat_5_6   = q_acc_5_6;
               const walberla::float64 elMat_5_7   = q_acc_5_7;
               const walberla::float64 elMat_5_8   = q_acc_5_8;
               const walberla::float64 elMat_5_9   = q_acc_5_9;
               const walberla::float64 elMat_5_10  = q_acc_5_10;
               const walberla::float64 elMat_5_11  = q_acc_5_11;
               const walberla::float64 elMat_6_0   = q_acc_0_6;
               const walberla::float64 elMat_6_1   = q_acc_1_6;
               const walberla::float64 elMat_6_2   = q_acc_2_6;
               const walberla::float64 elMat_6_3   = q_acc_3_6;
               const walberla::float64 elMat_6_4   = q_acc_4_6;
               const walberla::float64 elMat_6_5   = q_acc_5_6;
               const walberla::float64 elMat_6_6   = q_acc_6_6;
               const walberla::float64 elMat_6_7   = q_acc_6_7;
               const walberla::float64 elMat_6_8   = q_acc_6_8;
               const walberla::float64 elMat_6_9   = q_acc_6_9;
               const walberla::float64 elMat_6_10  = q_acc_6_10;
               const walberla::float64 elMat_6_11  = q_acc_6_11;
               const walberla::float64 elMat_7_0   = q_acc_0_7;
               const walberla::float64 elMat_7_1   = q_acc_1_7;
               const walberla::float64 elMat_7_2   = q_acc_2_7;
               const walberla::float64 elMat_7_3   = q_acc_3_7;
               const walberla::float64 elMat_7_4   = q_acc_4_7;
               const walberla::float64 elMat_7_5   = q_acc_5_7;
               const walberla::float64 elMat_7_6   = q_acc_6_7;
               const walberla::float64 elMat_7_7   = q_acc_7_7;
               const walberla::float64 elMat_7_8   = q_acc_7_8;
               const walberla::float64 elMat_7_9   = q_acc_7_9;
               const walberla::float64 elMat_7_10  = q_acc_7_10;
               const walberla::float64 elMat_7_11  = q_acc_7_11;
               const walberla::float64 elMat_8_0   = q_acc_0_8;
               const walberla::float64 elMat_8_1   = q_acc_1_8;
               const walberla::float64 elMat_8_2   = q_acc_2_8;
               const walberla::float64 elMat_8_3   = q_acc_3_8;
               const walberla::float64 elMat_8_4   = q_acc_4_8;
               const walberla::float64 elMat_8_5   = q_acc_5_8;
               const walberla::float64 elMat_8_6   = q_acc_6_8;
               const walberla::float64 elMat_8_7   = q_acc_7_8;
               const walberla::float64 elMat_8_8   = q_acc_8_8;
               const walberla::float64 elMat_8_9   = q_acc_8_9;
               const walberla::float64 elMat_8_10  = q_acc_8_10;
               const walberla::float64 elMat_8_11  = q_acc_8_11;
               const walberla::float64 elMat_9_0   = q_acc_0_9;
               const walberla::float64 elMat_9_1   = q_acc_1_9;
               const walberla::float64 elMat_9_2   = q_acc_2_9;
               const walberla::float64 elMat_9_3   = q_acc_3_9;
               const walberla::float64 elMat_9_4   = q_acc_4_9;
               const walberla::float64 elMat_9_5   = q_acc_5_9;
               const walberla::float64 elMat_9_6   = q_acc_6_9;
               const walberla::float64 elMat_9_7   = q_acc_7_9;
               const walberla::float64 elMat_9_8   = q_acc_8_9;
               const walberla::float64 elMat_9_9   = q_acc_9_9;
               const walberla::float64 elMat_9_10  = q_acc_9_10;
               const walberla::float64 elMat_9_11  = q_acc_9_11;
               const walberla::float64 elMat_10_0  = q_acc_0_10;
               const walberla::float64 elMat_10_1  = q_acc_1_10;
               const walberla::float64 elMat_10_2  = q_acc_2_10;
               const walberla::float64 elMat_10_3  = q_acc_3_10;
               const walberla::float64 elMat_10_4  = q_acc_4_10;
               const walberla::float64 elMat_10_5  = q_acc_5_10;
               const walberla::float64 elMat_10_6  = q_acc_6_10;
               const walberla::float64 elMat_10_7  = q_acc_7_10;
               const walberla::float64 elMat_10_8  = q_acc_8_10;
               const walberla::float64 elMat_10_9  = q_acc_9_10;
               const walberla::float64 elMat_10_10 = q_acc_10_10;
               const walberla::float64 elMat_10_11 = q_acc_10_11;
               const walberla::float64 elMat_11_0  = q_acc_0_11;
               const walberla::float64 elMat_11_1  = q_acc_1_11;
               const walberla::float64 elMat_11_2  = q_acc_2_11;
               const walberla::float64 elMat_11_3  = q_acc_3_11;
               const walberla::float64 elMat_11_4  = q_acc_4_11;
               const walberla::float64 elMat_11_5  = q_acc_5_11;
               const walberla::float64 elMat_11_6  = q_acc_6_11;
               const walberla::float64 elMat_11_7  = q_acc_7_11;
               const walberla::float64 elMat_11_8  = q_acc_8_11;
               const walberla::float64 elMat_11_9  = q_acc_9_11;
               const walberla::float64 elMat_11_10 = q_acc_10_11;
               const walberla::float64 elMat_11_11 = q_acc_11_11;

               std::vector< uint_t > _data_rowIdx( 12 );
               std::vector< uint_t > _data_colIdx( 12 );
               std::vector< real_t > _data_mat( 144 );

               _data_rowIdx[0] = ( (uint64_t) ( _data_dst_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[1] = ( (uint64_t) ( _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_rowIdx[2] = ( (uint64_t) ( _data_dst_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[3] = ( (uint64_t) ( _data_dst_edge_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                                 ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_rowIdx[4] = ( (
                   uint64_t) ( _data_dst_edge_0
                                   [ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[5] =
                   ( (uint64_t) ( _data_dst_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );
               _data_rowIdx[6]  = ( (uint64_t) ( _data_dst_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[7]  = ( (uint64_t) ( _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_rowIdx[8]  = ( (uint64_t) ( _data_dst_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[9]  = ( (uint64_t) ( _data_dst_edge_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                                 ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_rowIdx[10] = ( (
                   uint64_t) ( _data_dst_edge_1
                                   [ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_rowIdx[11] =
                   ( (uint64_t) ( _data_dst_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );
               _data_colIdx[0] = ( (uint64_t) ( _data_src_vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[1] = ( (uint64_t) ( _data_src_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[2] = ( (uint64_t) ( _data_src_vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[3] = ( (uint64_t) ( _data_src_edge_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                                 ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[4] = ( (
                   uint64_t) ( _data_src_edge_0
                                   [ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[5] =
                   ( (uint64_t) ( _data_src_edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );
               _data_colIdx[6]  = ( (uint64_t) ( _data_src_vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[7]  = ( (uint64_t) ( _data_src_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[8]  = ( (uint64_t) ( _data_src_vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                                   ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[9]  = ( (uint64_t) ( _data_src_edge_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                                 ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] ) );
               _data_colIdx[10] = ( (
                   uint64_t) ( _data_src_edge_1
                                   [ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                    2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] ) );
               _data_colIdx[11] =
                   ( (uint64_t) ( _data_src_edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                   ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                                   ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) /
                                                     ( 2 ) )] ) );

               /* Apply basis transformation */

               _data_mat[0]   = ( (real_t) ( elMat_0_0 ) );
               _data_mat[1]   = ( (real_t) ( elMat_0_1 ) );
               _data_mat[2]   = ( (real_t) ( elMat_0_2 ) );
               _data_mat[3]   = ( (real_t) ( elMat_0_3 ) );
               _data_mat[4]   = ( (real_t) ( elMat_0_4 ) );
               _data_mat[5]   = ( (real_t) ( elMat_0_5 ) );
               _data_mat[6]   = ( (real_t) ( elMat_0_6 ) );
               _data_mat[7]   = ( (real_t) ( elMat_0_7 ) );
               _data_mat[8]   = ( (real_t) ( elMat_0_8 ) );
               _data_mat[9]   = ( (real_t) ( elMat_0_9 ) );
               _data_mat[10]  = ( (real_t) ( elMat_0_10 ) );
               _data_mat[11]  = ( (real_t) ( elMat_0_11 ) );
               _data_mat[12]  = ( (real_t) ( elMat_1_0 ) );
               _data_mat[13]  = ( (real_t) ( elMat_1_1 ) );
               _data_mat[14]  = ( (real_t) ( elMat_1_2 ) );
               _data_mat[15]  = ( (real_t) ( elMat_1_3 ) );
               _data_mat[16]  = ( (real_t) ( elMat_1_4 ) );
               _data_mat[17]  = ( (real_t) ( elMat_1_5 ) );
               _data_mat[18]  = ( (real_t) ( elMat_1_6 ) );
               _data_mat[19]  = ( (real_t) ( elMat_1_7 ) );
               _data_mat[20]  = ( (real_t) ( elMat_1_8 ) );
               _data_mat[21]  = ( (real_t) ( elMat_1_9 ) );
               _data_mat[22]  = ( (real_t) ( elMat_1_10 ) );
               _data_mat[23]  = ( (real_t) ( elMat_1_11 ) );
               _data_mat[24]  = ( (real_t) ( elMat_2_0 ) );
               _data_mat[25]  = ( (real_t) ( elMat_2_1 ) );
               _data_mat[26]  = ( (real_t) ( elMat_2_2 ) );
               _data_mat[27]  = ( (real_t) ( elMat_2_3 ) );
               _data_mat[28]  = ( (real_t) ( elMat_2_4 ) );
               _data_mat[29]  = ( (real_t) ( elMat_2_5 ) );
               _data_mat[30]  = ( (real_t) ( elMat_2_6 ) );
               _data_mat[31]  = ( (real_t) ( elMat_2_7 ) );
               _data_mat[32]  = ( (real_t) ( elMat_2_8 ) );
               _data_mat[33]  = ( (real_t) ( elMat_2_9 ) );
               _data_mat[34]  = ( (real_t) ( elMat_2_10 ) );
               _data_mat[35]  = ( (real_t) ( elMat_2_11 ) );
               _data_mat[36]  = ( (real_t) ( elMat_3_0 ) );
               _data_mat[37]  = ( (real_t) ( elMat_3_1 ) );
               _data_mat[38]  = ( (real_t) ( elMat_3_2 ) );
               _data_mat[39]  = ( (real_t) ( elMat_3_3 ) );
               _data_mat[40]  = ( (real_t) ( elMat_3_4 ) );
               _data_mat[41]  = ( (real_t) ( elMat_3_5 ) );
               _data_mat[42]  = ( (real_t) ( elMat_3_6 ) );
               _data_mat[43]  = ( (real_t) ( elMat_3_7 ) );
               _data_mat[44]  = ( (real_t) ( elMat_3_8 ) );
               _data_mat[45]  = ( (real_t) ( elMat_3_9 ) );
               _data_mat[46]  = ( (real_t) ( elMat_3_10 ) );
               _data_mat[47]  = ( (real_t) ( elMat_3_11 ) );
               _data_mat[48]  = ( (real_t) ( elMat_4_0 ) );
               _data_mat[49]  = ( (real_t) ( elMat_4_1 ) );
               _data_mat[50]  = ( (real_t) ( elMat_4_2 ) );
               _data_mat[51]  = ( (real_t) ( elMat_4_3 ) );
               _data_mat[52]  = ( (real_t) ( elMat_4_4 ) );
               _data_mat[53]  = ( (real_t) ( elMat_4_5 ) );
               _data_mat[54]  = ( (real_t) ( elMat_4_6 ) );
               _data_mat[55]  = ( (real_t) ( elMat_4_7 ) );
               _data_mat[56]  = ( (real_t) ( elMat_4_8 ) );
               _data_mat[57]  = ( (real_t) ( elMat_4_9 ) );
               _data_mat[58]  = ( (real_t) ( elMat_4_10 ) );
               _data_mat[59]  = ( (real_t) ( elMat_4_11 ) );
               _data_mat[60]  = ( (real_t) ( elMat_5_0 ) );
               _data_mat[61]  = ( (real_t) ( elMat_5_1 ) );
               _data_mat[62]  = ( (real_t) ( elMat_5_2 ) );
               _data_mat[63]  = ( (real_t) ( elMat_5_3 ) );
               _data_mat[64]  = ( (real_t) ( elMat_5_4 ) );
               _data_mat[65]  = ( (real_t) ( elMat_5_5 ) );
               _data_mat[66]  = ( (real_t) ( elMat_5_6 ) );
               _data_mat[67]  = ( (real_t) ( elMat_5_7 ) );
               _data_mat[68]  = ( (real_t) ( elMat_5_8 ) );
               _data_mat[69]  = ( (real_t) ( elMat_5_9 ) );
               _data_mat[70]  = ( (real_t) ( elMat_5_10 ) );
               _data_mat[71]  = ( (real_t) ( elMat_5_11 ) );
               _data_mat[72]  = ( (real_t) ( elMat_6_0 ) );
               _data_mat[73]  = ( (real_t) ( elMat_6_1 ) );
               _data_mat[74]  = ( (real_t) ( elMat_6_2 ) );
               _data_mat[75]  = ( (real_t) ( elMat_6_3 ) );
               _data_mat[76]  = ( (real_t) ( elMat_6_4 ) );
               _data_mat[77]  = ( (real_t) ( elMat_6_5 ) );
               _data_mat[78]  = ( (real_t) ( elMat_6_6 ) );
               _data_mat[79]  = ( (real_t) ( elMat_6_7 ) );
               _data_mat[80]  = ( (real_t) ( elMat_6_8 ) );
               _data_mat[81]  = ( (real_t) ( elMat_6_9 ) );
               _data_mat[82]  = ( (real_t) ( elMat_6_10 ) );
               _data_mat[83]  = ( (real_t) ( elMat_6_11 ) );
               _data_mat[84]  = ( (real_t) ( elMat_7_0 ) );
               _data_mat[85]  = ( (real_t) ( elMat_7_1 ) );
               _data_mat[86]  = ( (real_t) ( elMat_7_2 ) );
               _data_mat[87]  = ( (real_t) ( elMat_7_3 ) );
               _data_mat[88]  = ( (real_t) ( elMat_7_4 ) );
               _data_mat[89]  = ( (real_t) ( elMat_7_5 ) );
               _data_mat[90]  = ( (real_t) ( elMat_7_6 ) );
               _data_mat[91]  = ( (real_t) ( elMat_7_7 ) );
               _data_mat[92]  = ( (real_t) ( elMat_7_8 ) );
               _data_mat[93]  = ( (real_t) ( elMat_7_9 ) );
               _data_mat[94]  = ( (real_t) ( elMat_7_10 ) );
               _data_mat[95]  = ( (real_t) ( elMat_7_11 ) );
               _data_mat[96]  = ( (real_t) ( elMat_8_0 ) );
               _data_mat[97]  = ( (real_t) ( elMat_8_1 ) );
               _data_mat[98]  = ( (real_t) ( elMat_8_2 ) );
               _data_mat[99]  = ( (real_t) ( elMat_8_3 ) );
               _data_mat[100] = ( (real_t) ( elMat_8_4 ) );
               _data_mat[101] = ( (real_t) ( elMat_8_5 ) );
               _data_mat[102] = ( (real_t) ( elMat_8_6 ) );
               _data_mat[103] = ( (real_t) ( elMat_8_7 ) );
               _data_mat[104] = ( (real_t) ( elMat_8_8 ) );
               _data_mat[105] = ( (real_t) ( elMat_8_9 ) );
               _data_mat[106] = ( (real_t) ( elMat_8_10 ) );
               _data_mat[107] = ( (real_t) ( elMat_8_11 ) );
               _data_mat[108] = ( (real_t) ( elMat_9_0 ) );
               _data_mat[109] = ( (real_t) ( elMat_9_1 ) );
               _data_mat[110] = ( (real_t) ( elMat_9_2 ) );
               _data_mat[111] = ( (real_t) ( elMat_9_3 ) );
               _data_mat[112] = ( (real_t) ( elMat_9_4 ) );
               _data_mat[113] = ( (real_t) ( elMat_9_5 ) );
               _data_mat[114] = ( (real_t) ( elMat_9_6 ) );
               _data_mat[115] = ( (real_t) ( elMat_9_7 ) );
               _data_mat[116] = ( (real_t) ( elMat_9_8 ) );
               _data_mat[117] = ( (real_t) ( elMat_9_9 ) );
               _data_mat[118] = ( (real_t) ( elMat_9_10 ) );
               _data_mat[119] = ( (real_t) ( elMat_9_11 ) );
               _data_mat[120] = ( (real_t) ( elMat_10_0 ) );
               _data_mat[121] = ( (real_t) ( elMat_10_1 ) );
               _data_mat[122] = ( (real_t) ( elMat_10_2 ) );
               _data_mat[123] = ( (real_t) ( elMat_10_3 ) );
               _data_mat[124] = ( (real_t) ( elMat_10_4 ) );
               _data_mat[125] = ( (real_t) ( elMat_10_5 ) );
               _data_mat[126] = ( (real_t) ( elMat_10_6 ) );
               _data_mat[127] = ( (real_t) ( elMat_10_7 ) );
               _data_mat[128] = ( (real_t) ( elMat_10_8 ) );
               _data_mat[129] = ( (real_t) ( elMat_10_9 ) );
               _data_mat[130] = ( (real_t) ( elMat_10_10 ) );
               _data_mat[131] = ( (real_t) ( elMat_10_11 ) );
               _data_mat[132] = ( (real_t) ( elMat_11_0 ) );
               _data_mat[133] = ( (real_t) ( elMat_11_1 ) );
               _data_mat[134] = ( (real_t) ( elMat_11_2 ) );
               _data_mat[135] = ( (real_t) ( elMat_11_3 ) );
               _data_mat[136] = ( (real_t) ( elMat_11_4 ) );
               _data_mat[137] = ( (real_t) ( elMat_11_5 ) );
               _data_mat[138] = ( (real_t) ( elMat_11_6 ) );
               _data_mat[139] = ( (real_t) ( elMat_11_7 ) );
               _data_mat[140] = ( (real_t) ( elMat_11_8 ) );
               _data_mat[141] = ( (real_t) ( elMat_11_9 ) );
               _data_mat[142] = ( (real_t) ( elMat_11_10 ) );
               _data_mat[143] = ( (real_t) ( elMat_11_11 ) );

               mat->addValues( _data_rowIdx, _data_colIdx, _data_mat );
            }
      }
   }
}
void P2VectorElementwiseFullStokesViscoplastic::
    computeInverseDiagonalOperatorValues_P2VectorElementwiseFullStokesViscoplastic_macro_2D(
        walberla::float64* RESTRICT _data_invDiag__edge_0,
        walberla::float64* RESTRICT _data_invDiag__edge_1,
        walberla::float64* RESTRICT _data_invDiag__vertex_0,
        walberla::float64* RESTRICT _data_invDiag__vertex_1,
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
      const walberla::float64 _data_q_w[] = { 0.16666666666666666, 0.16666666666666666, 0.16666666666666666 };

      const walberla::float64 _data_q_p_0[] = { 0.16666666666666666, 0.66666666666666663, 0.16666666666666666 };

      const walberla::float64 _data_q_p_1[] = { 0.66666666666666663, 0.16666666666666666, 0.16666666666666666 };

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
             jac_affine_0_0_GRAY * jac_affine_1_1_GRAY - jac_affine_0_1_GRAY * jac_affine_1_0_GRAY;
         const walberla::float64 tmp_coords_jac_2_GRAY   = 1.0 / ( tmp_coords_jac_1_GRAY );
         const walberla::float64 jac_affine_inv_0_0_GRAY = jac_affine_1_1_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_0_1_GRAY = -jac_affine_0_1_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_1_0_GRAY = -jac_affine_1_0_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 jac_affine_inv_1_1_GRAY = jac_affine_0_0_GRAY * tmp_coords_jac_2_GRAY;
         const walberla::float64 abs_det_jac_affine_GRAY = abs( tmp_coords_jac_1_GRAY );
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
               const walberla::float64 uy_dof_2    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               walberla::float64       q_acc_0_0   = 0.0;
               walberla::float64       q_acc_1_1   = 0.0;
               walberla::float64       q_acc_2_2   = 0.0;
               walberla::float64       q_acc_3_3   = 0.0;
               walberla::float64       q_acc_4_4   = 0.0;
               walberla::float64       q_acc_5_5   = 0.0;
               walberla::float64       q_acc_6_6   = 0.0;
               walberla::float64       q_acc_7_7   = 0.0;
               walberla::float64       q_acc_8_8   = 0.0;
               walberla::float64       q_acc_9_9   = 0.0;
               walberla::float64       q_acc_10_10 = 0.0;
               walberla::float64       q_acc_11_11 = 0.0;
               for ( int64_t q = 0; q < 3; q += 1 )
               {
                  const walberla::float64 tmp_qloop_0  = 4.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_1  = 4.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_2  = tmp_qloop_0 + tmp_qloop_1 - 3.0;
                  const walberla::float64 tmp_qloop_3  = jac_affine_inv_0_0_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_4  = jac_affine_inv_1_0_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_5  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_6  = jac_affine_inv_0_1_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_7  = jac_affine_inv_1_1_GRAY * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_8  = tmp_qloop_6 + tmp_qloop_7;
                  const walberla::float64 tmp_qloop_9  = -ux_dof_0 + ux_dof_1;
                  const walberla::float64 tmp_qloop_10 = -ux_dof_0 + ux_dof_2;
                  const walberla::float64 tmp_qloop_11 = -uy_dof_0 + uy_dof_1;
                  const walberla::float64 tmp_qloop_12 = -uy_dof_0 + uy_dof_2;
                  const walberla::float64 tmp_qloop_13 = jac_affine_inv_0_0_GRAY * tmp_qloop_11;
                  const walberla::float64 tmp_qloop_14 = jac_affine_inv_0_1_GRAY * tmp_qloop_9;
                  const walberla::float64 tmp_qloop_15 = jac_affine_inv_1_0_GRAY * tmp_qloop_12;
                  const walberla::float64 tmp_qloop_16 = jac_affine_inv_1_1_GRAY * tmp_qloop_10;
                  const walberla::float64 tmp_qloop_17 =
                      1.0 /
                      ( 1.0 /
                            ( mu_star +
                              sigma_y * 1.0 /
                                  ( pow( ( ( jac_affine_inv_0_0_GRAY * tmp_qloop_9 + jac_affine_inv_1_0_GRAY * tmp_qloop_10 ) *
                                           ( jac_affine_inv_0_0_GRAY * tmp_qloop_9 + jac_affine_inv_1_0_GRAY * tmp_qloop_10 ) ) +
                                             ( ( jac_affine_inv_0_1_GRAY * tmp_qloop_11 +
                                                 jac_affine_inv_1_1_GRAY * tmp_qloop_12 ) *
                                               ( jac_affine_inv_0_1_GRAY * tmp_qloop_11 +
                                                 jac_affine_inv_1_1_GRAY * tmp_qloop_12 ) ) +
                                             ( tmp_qloop_13 + tmp_qloop_14 + tmp_qloop_15 + tmp_qloop_16 ) *
                                                 ( tmp_qloop_13 * 0.5 + tmp_qloop_14 * 0.5 + tmp_qloop_15 * 0.5 +
                                                   tmp_qloop_16 * 0.5 ),
                                         0.50000000000000000 ) +
                                    9.9999999999999995e-21 ) ) *
                            1.0 +
                        1.0 /
                            ( mu_lin_dof_0 * ( 1.0 - _data_q_p_0[q] - _data_q_p_1[q] ) + mu_lin_dof_1 * _data_q_p_0[q] +
                              mu_lin_dof_2 * _data_q_p_1[q] ) *
                            1.0 ) *
                      2.0 * _data_q_w[q];
                  const walberla::float64 tmp_qloop_18 = ( jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY );
                  const walberla::float64 tmp_qloop_19 = ( ( tmp_qloop_0 - 1.0 ) * ( tmp_qloop_0 - 1.0 ) );
                  const walberla::float64 tmp_qloop_20 = abs_det_jac_affine_GRAY * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_21 = tmp_qloop_19 * tmp_qloop_20;
                  const walberla::float64 tmp_qloop_22 = tmp_qloop_19 * 2.0;
                  const walberla::float64 tmp_qloop_23 = ( jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY );
                  const walberla::float64 tmp_qloop_24 = tmp_qloop_19 * 1.0;
                  const walberla::float64 tmp_qloop_25 = ( jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY );
                  const walberla::float64 tmp_qloop_26 = ( ( tmp_qloop_1 - 1.0 ) * ( tmp_qloop_1 - 1.0 ) );
                  const walberla::float64 tmp_qloop_27 = tmp_qloop_20 * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_28 = tmp_qloop_26 * 2.0;
                  const walberla::float64 tmp_qloop_29 = ( jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY );
                  const walberla::float64 tmp_qloop_30 = tmp_qloop_26 * 1.0;
                  const walberla::float64 tmp_qloop_31 = 2.6666666666666665 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_32 = jac_affine_inv_1_0_GRAY * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_33 = 2.6666666666666665 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_34 = jac_affine_inv_0_0_GRAY * tmp_qloop_33;
                  const walberla::float64 tmp_qloop_35 = jac_affine_inv_1_0_GRAY * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_36 = jac_affine_inv_0_0_GRAY * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_37 = tmp_qloop_35 + tmp_qloop_36;
                  const walberla::float64 tmp_qloop_38 = 2.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_39 = jac_affine_inv_1_1_GRAY * tmp_qloop_38;
                  const walberla::float64 tmp_qloop_40 = 2.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_41 = jac_affine_inv_0_1_GRAY * tmp_qloop_40;
                  const walberla::float64 tmp_qloop_42 = jac_affine_inv_1_1_GRAY * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_43 = jac_affine_inv_0_1_GRAY * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_44 = tmp_qloop_42 + tmp_qloop_43;
                  const walberla::float64 tmp_qloop_45 = -tmp_qloop_0 - 8.0 * _data_q_p_1[q] + 4.0;
                  const walberla::float64 tmp_qloop_46 = jac_affine_inv_1_0_GRAY * tmp_qloop_45 - tmp_qloop_36;
                  const walberla::float64 tmp_qloop_47 = jac_affine_inv_1_1_GRAY * tmp_qloop_45 - tmp_qloop_43;
                  const walberla::float64 tmp_qloop_48 = -tmp_qloop_1 - 8.0 * _data_q_p_0[q] + 4.0;
                  const walberla::float64 tmp_qloop_49 = jac_affine_inv_0_0_GRAY * tmp_qloop_48 - tmp_qloop_35;
                  const walberla::float64 tmp_qloop_50 = jac_affine_inv_0_1_GRAY * tmp_qloop_48 - tmp_qloop_42;
                  const walberla::float64 tmp_qloop_51 = jac_affine_inv_1_1_GRAY * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_52 = jac_affine_inv_0_1_GRAY * tmp_qloop_33;
                  const walberla::float64 tmp_qloop_53 = jac_affine_inv_1_0_GRAY * tmp_qloop_38;
                  const walberla::float64 tmp_qloop_54 = jac_affine_inv_0_0_GRAY * tmp_qloop_40;
                  const walberla::float64 q_tmp_0_0 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_GRAY * tmp_qloop_5 *
                            ( tmp_qloop_3 * 0.66666666666666663 + tmp_qloop_4 * 0.66666666666666663 ) +
                        abs_det_jac_affine_GRAY * ( ( tmp_qloop_5 * tmp_qloop_5 ) * 2.0 +
                                                    tmp_qloop_8 * ( tmp_qloop_6 * 0.5 + tmp_qloop_7 * 0.5 ) * 2.0 ) );
                  const walberla::float64 q_tmp_1_1 =
                      tmp_qloop_17 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_18 * tmp_qloop_22 + tmp_qloop_23 * tmp_qloop_24 ) -
                                       tmp_qloop_18 * tmp_qloop_21 );
                  const walberla::float64 q_tmp_2_2 =
                      tmp_qloop_17 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_25 * tmp_qloop_28 + tmp_qloop_29 * tmp_qloop_30 ) -
                                       tmp_qloop_25 * tmp_qloop_27 );
                  const walberla::float64 q_tmp_3_3 =
                      tmp_qloop_17 * ( -abs_det_jac_affine_GRAY * tmp_qloop_37 * ( tmp_qloop_32 + tmp_qloop_34 ) +
                                       abs_det_jac_affine_GRAY * ( ( tmp_qloop_37 * tmp_qloop_37 ) * 2.0 +
                                                                   tmp_qloop_44 * ( tmp_qloop_39 + tmp_qloop_41 ) * 2.0 ) );
                  const walberla::float64 q_tmp_4_4 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_GRAY * tmp_qloop_46 *
                            ( jac_affine_inv_1_0_GRAY * tmp_qloop_45 * 0.66666666666666663 - tmp_qloop_34 ) +
                        abs_det_jac_affine_GRAY *
                            ( ( tmp_qloop_46 * tmp_qloop_46 ) * 2.0 +
                              tmp_qloop_47 * ( jac_affine_inv_1_1_GRAY * tmp_qloop_45 * 0.5 - tmp_qloop_41 ) * 2.0 ) );
                  const walberla::float64 q_tmp_5_5 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_GRAY * tmp_qloop_49 *
                            ( jac_affine_inv_0_0_GRAY * tmp_qloop_48 * 0.66666666666666663 - tmp_qloop_32 ) +
                        abs_det_jac_affine_GRAY *
                            ( ( tmp_qloop_49 * tmp_qloop_49 ) * 2.0 +
                              tmp_qloop_50 * ( jac_affine_inv_0_1_GRAY * tmp_qloop_48 * 0.5 - tmp_qloop_39 ) * 2.0 ) );
                  const walberla::float64 q_tmp_6_6 =
                      tmp_qloop_17 * ( -abs_det_jac_affine_GRAY * tmp_qloop_8 *
                                           ( tmp_qloop_6 * 0.66666666666666663 + tmp_qloop_7 * 0.66666666666666663 ) +
                                       abs_det_jac_affine_GRAY * ( tmp_qloop_5 * ( tmp_qloop_3 * 0.5 + tmp_qloop_4 * 0.5 ) * 2.0 +
                                                                   ( tmp_qloop_8 * tmp_qloop_8 ) * 2.0 ) );
                  const walberla::float64 q_tmp_7_7 =
                      tmp_qloop_17 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_18 * tmp_qloop_24 + tmp_qloop_22 * tmp_qloop_23 ) -
                                       tmp_qloop_21 * tmp_qloop_23 );
                  const walberla::float64 q_tmp_8_8 =
                      tmp_qloop_17 * ( abs_det_jac_affine_GRAY * ( tmp_qloop_25 * tmp_qloop_30 + tmp_qloop_28 * tmp_qloop_29 ) -
                                       tmp_qloop_27 * tmp_qloop_29 );
                  const walberla::float64 q_tmp_9_9 =
                      tmp_qloop_17 * ( -abs_det_jac_affine_GRAY * tmp_qloop_44 * ( tmp_qloop_51 + tmp_qloop_52 ) +
                                       abs_det_jac_affine_GRAY * ( tmp_qloop_37 * ( tmp_qloop_53 + tmp_qloop_54 ) * 2.0 +
                                                                   ( tmp_qloop_44 * tmp_qloop_44 ) * 2.0 ) );
                  const walberla::float64 q_tmp_10_10 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_GRAY * tmp_qloop_47 *
                            ( jac_affine_inv_1_1_GRAY * tmp_qloop_45 * 0.66666666666666663 - tmp_qloop_52 ) +
                        abs_det_jac_affine_GRAY *
                            ( tmp_qloop_46 * ( jac_affine_inv_1_0_GRAY * tmp_qloop_45 * 0.5 - tmp_qloop_54 ) * 2.0 +
                              ( tmp_qloop_47 * tmp_qloop_47 ) * 2.0 ) );
                  const walberla::float64 q_tmp_11_11 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_GRAY * tmp_qloop_50 *
                            ( jac_affine_inv_0_1_GRAY * tmp_qloop_48 * 0.66666666666666663 - tmp_qloop_51 ) +
                        abs_det_jac_affine_GRAY *
                            ( tmp_qloop_49 * ( jac_affine_inv_0_0_GRAY * tmp_qloop_48 * 0.5 - tmp_qloop_53 ) * 2.0 +
                              ( tmp_qloop_50 * tmp_qloop_50 ) * 2.0 ) );
                  q_acc_0_0   = q_acc_0_0 + q_tmp_0_0;
                  q_acc_1_1   = q_acc_1_1 + q_tmp_1_1;
                  q_acc_2_2   = q_acc_2_2 + q_tmp_2_2;
                  q_acc_3_3   = q_acc_3_3 + q_tmp_3_3;
                  q_acc_4_4   = q_acc_4_4 + q_tmp_4_4;
                  q_acc_5_5   = q_acc_5_5 + q_tmp_5_5;
                  q_acc_6_6   = q_acc_6_6 + q_tmp_6_6;
                  q_acc_7_7   = q_acc_7_7 + q_tmp_7_7;
                  q_acc_8_8   = q_acc_8_8 + q_tmp_8_8;
                  q_acc_9_9   = q_acc_9_9 + q_tmp_9_9;
                  q_acc_10_10 = q_acc_10_10 + q_tmp_10_10;
                  q_acc_11_11 = q_acc_11_11 + q_tmp_11_11;
               }
               const walberla::float64 elMatDiag_0  = q_acc_0_0;
               const walberla::float64 elMatDiag_1  = q_acc_1_1;
               const walberla::float64 elMatDiag_2  = q_acc_2_2;
               const walberla::float64 elMatDiag_3  = q_acc_3_3;
               const walberla::float64 elMatDiag_4  = q_acc_4_4;
               const walberla::float64 elMatDiag_5  = q_acc_5_5;
               const walberla::float64 elMatDiag_6  = q_acc_6_6;
               const walberla::float64 elMatDiag_7  = q_acc_7_7;
               const walberla::float64 elMatDiag_8  = q_acc_8_8;
               const walberla::float64 elMatDiag_9  = q_acc_9_9;
               const walberla::float64 elMatDiag_10 = q_acc_10_10;
               const walberla::float64 elMatDiag_11 = q_acc_11_11;
               _data_invDiag__vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatDiag_0 + _data_invDiag__vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_invDiag__vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_1 + _data_invDiag__vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag__vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_2 + _data_invDiag__vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatDiag_3 +
                   _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatDiag_4 +
                   _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatDiag_5 + _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_invDiag__vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatDiag_6 + _data_invDiag__vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
               _data_invDiag__vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_7 + _data_invDiag__vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag__vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_8 + _data_invDiag__vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatDiag_9 +
                   _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatDiag_10 +
                   _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )] =
                   elMatDiag_11 + _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                                        ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) )];
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
             jac_affine_0_0_BLUE * jac_affine_1_1_BLUE - jac_affine_0_1_BLUE * jac_affine_1_0_BLUE;
         const walberla::float64 tmp_coords_jac_6_BLUE   = 1.0 / ( tmp_coords_jac_5_BLUE );
         const walberla::float64 jac_affine_inv_0_0_BLUE = jac_affine_1_1_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_0_1_BLUE = -jac_affine_0_1_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_1_0_BLUE = -jac_affine_1_0_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 jac_affine_inv_1_1_BLUE = jac_affine_0_0_BLUE * tmp_coords_jac_6_BLUE;
         const walberla::float64 abs_det_jac_affine_BLUE = abs( tmp_coords_jac_5_BLUE );
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
               const walberla::float64 uy_dof_1    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               const walberla::float64 uy_dof_2    = _data_uy[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                           ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               walberla::float64       q_acc_0_0   = 0.0;
               walberla::float64       q_acc_1_1   = 0.0;
               walberla::float64       q_acc_2_2   = 0.0;
               walberla::float64       q_acc_3_3   = 0.0;
               walberla::float64       q_acc_4_4   = 0.0;
               walberla::float64       q_acc_5_5   = 0.0;
               walberla::float64       q_acc_6_6   = 0.0;
               walberla::float64       q_acc_7_7   = 0.0;
               walberla::float64       q_acc_8_8   = 0.0;
               walberla::float64       q_acc_9_9   = 0.0;
               walberla::float64       q_acc_10_10 = 0.0;
               walberla::float64       q_acc_11_11 = 0.0;
               for ( int64_t q = 0; q < 3; q += 1 )
               {
                  const walberla::float64 tmp_qloop_0  = 4.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_1  = 4.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_2  = tmp_qloop_0 + tmp_qloop_1 - 3.0;
                  const walberla::float64 tmp_qloop_3  = jac_affine_inv_0_0_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_4  = jac_affine_inv_1_0_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_5  = tmp_qloop_3 + tmp_qloop_4;
                  const walberla::float64 tmp_qloop_6  = jac_affine_inv_0_1_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_7  = jac_affine_inv_1_1_BLUE * tmp_qloop_2;
                  const walberla::float64 tmp_qloop_8  = tmp_qloop_6 + tmp_qloop_7;
                  const walberla::float64 tmp_qloop_9  = -ux_dof_0 + ux_dof_1;
                  const walberla::float64 tmp_qloop_10 = -ux_dof_0 + ux_dof_2;
                  const walberla::float64 tmp_qloop_11 = -uy_dof_0 + uy_dof_1;
                  const walberla::float64 tmp_qloop_12 = -uy_dof_0 + uy_dof_2;
                  const walberla::float64 tmp_qloop_13 = jac_affine_inv_0_0_BLUE * tmp_qloop_11;
                  const walberla::float64 tmp_qloop_14 = jac_affine_inv_0_1_BLUE * tmp_qloop_9;
                  const walberla::float64 tmp_qloop_15 = jac_affine_inv_1_0_BLUE * tmp_qloop_12;
                  const walberla::float64 tmp_qloop_16 = jac_affine_inv_1_1_BLUE * tmp_qloop_10;
                  const walberla::float64 tmp_qloop_17 =
                      1.0 /
                      ( 1.0 /
                            ( mu_star +
                              sigma_y * 1.0 /
                                  ( pow( ( ( jac_affine_inv_0_0_BLUE * tmp_qloop_9 + jac_affine_inv_1_0_BLUE * tmp_qloop_10 ) *
                                           ( jac_affine_inv_0_0_BLUE * tmp_qloop_9 + jac_affine_inv_1_0_BLUE * tmp_qloop_10 ) ) +
                                             ( ( jac_affine_inv_0_1_BLUE * tmp_qloop_11 +
                                                 jac_affine_inv_1_1_BLUE * tmp_qloop_12 ) *
                                               ( jac_affine_inv_0_1_BLUE * tmp_qloop_11 +
                                                 jac_affine_inv_1_1_BLUE * tmp_qloop_12 ) ) +
                                             ( tmp_qloop_13 + tmp_qloop_14 + tmp_qloop_15 + tmp_qloop_16 ) *
                                                 ( tmp_qloop_13 * 0.5 + tmp_qloop_14 * 0.5 + tmp_qloop_15 * 0.5 +
                                                   tmp_qloop_16 * 0.5 ),
                                         0.50000000000000000 ) +
                                    9.9999999999999995e-21 ) ) *
                            1.0 +
                        1.0 /
                            ( mu_lin_dof_0 * ( 1.0 - _data_q_p_0[q] - _data_q_p_1[q] ) + mu_lin_dof_1 * _data_q_p_0[q] +
                              mu_lin_dof_2 * _data_q_p_1[q] ) *
                            1.0 ) *
                      2.0 * _data_q_w[q];
                  const walberla::float64 tmp_qloop_18 = ( jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE );
                  const walberla::float64 tmp_qloop_19 = ( ( tmp_qloop_0 - 1.0 ) * ( tmp_qloop_0 - 1.0 ) );
                  const walberla::float64 tmp_qloop_20 = abs_det_jac_affine_BLUE * 0.66666666666666663;
                  const walberla::float64 tmp_qloop_21 = tmp_qloop_19 * tmp_qloop_20;
                  const walberla::float64 tmp_qloop_22 = tmp_qloop_19 * 2.0;
                  const walberla::float64 tmp_qloop_23 = ( jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE );
                  const walberla::float64 tmp_qloop_24 = tmp_qloop_19 * 1.0;
                  const walberla::float64 tmp_qloop_25 = ( jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE );
                  const walberla::float64 tmp_qloop_26 = ( ( tmp_qloop_1 - 1.0 ) * ( tmp_qloop_1 - 1.0 ) );
                  const walberla::float64 tmp_qloop_27 = tmp_qloop_20 * tmp_qloop_26;
                  const walberla::float64 tmp_qloop_28 = tmp_qloop_26 * 2.0;
                  const walberla::float64 tmp_qloop_29 = ( jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE );
                  const walberla::float64 tmp_qloop_30 = tmp_qloop_26 * 1.0;
                  const walberla::float64 tmp_qloop_31 = 2.6666666666666665 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_32 = jac_affine_inv_1_0_BLUE * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_33 = 2.6666666666666665 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_34 = jac_affine_inv_0_0_BLUE * tmp_qloop_33;
                  const walberla::float64 tmp_qloop_35 = jac_affine_inv_1_0_BLUE * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_36 = jac_affine_inv_0_0_BLUE * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_37 = tmp_qloop_35 + tmp_qloop_36;
                  const walberla::float64 tmp_qloop_38 = 2.0 * _data_q_p_0[q];
                  const walberla::float64 tmp_qloop_39 = jac_affine_inv_1_1_BLUE * tmp_qloop_38;
                  const walberla::float64 tmp_qloop_40 = 2.0 * _data_q_p_1[q];
                  const walberla::float64 tmp_qloop_41 = jac_affine_inv_0_1_BLUE * tmp_qloop_40;
                  const walberla::float64 tmp_qloop_42 = jac_affine_inv_1_1_BLUE * tmp_qloop_0;
                  const walberla::float64 tmp_qloop_43 = jac_affine_inv_0_1_BLUE * tmp_qloop_1;
                  const walberla::float64 tmp_qloop_44 = tmp_qloop_42 + tmp_qloop_43;
                  const walberla::float64 tmp_qloop_45 = -tmp_qloop_0 - 8.0 * _data_q_p_1[q] + 4.0;
                  const walberla::float64 tmp_qloop_46 = jac_affine_inv_1_0_BLUE * tmp_qloop_45 - tmp_qloop_36;
                  const walberla::float64 tmp_qloop_47 = jac_affine_inv_1_1_BLUE * tmp_qloop_45 - tmp_qloop_43;
                  const walberla::float64 tmp_qloop_48 = -tmp_qloop_1 - 8.0 * _data_q_p_0[q] + 4.0;
                  const walberla::float64 tmp_qloop_49 = jac_affine_inv_0_0_BLUE * tmp_qloop_48 - tmp_qloop_35;
                  const walberla::float64 tmp_qloop_50 = jac_affine_inv_0_1_BLUE * tmp_qloop_48 - tmp_qloop_42;
                  const walberla::float64 tmp_qloop_51 = jac_affine_inv_1_1_BLUE * tmp_qloop_31;
                  const walberla::float64 tmp_qloop_52 = jac_affine_inv_0_1_BLUE * tmp_qloop_33;
                  const walberla::float64 tmp_qloop_53 = jac_affine_inv_1_0_BLUE * tmp_qloop_38;
                  const walberla::float64 tmp_qloop_54 = jac_affine_inv_0_0_BLUE * tmp_qloop_40;
                  const walberla::float64 q_tmp_0_0 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_BLUE * tmp_qloop_5 *
                            ( tmp_qloop_3 * 0.66666666666666663 + tmp_qloop_4 * 0.66666666666666663 ) +
                        abs_det_jac_affine_BLUE * ( ( tmp_qloop_5 * tmp_qloop_5 ) * 2.0 +
                                                    tmp_qloop_8 * ( tmp_qloop_6 * 0.5 + tmp_qloop_7 * 0.5 ) * 2.0 ) );
                  const walberla::float64 q_tmp_1_1 =
                      tmp_qloop_17 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_18 * tmp_qloop_22 + tmp_qloop_23 * tmp_qloop_24 ) -
                                       tmp_qloop_18 * tmp_qloop_21 );
                  const walberla::float64 q_tmp_2_2 =
                      tmp_qloop_17 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_25 * tmp_qloop_28 + tmp_qloop_29 * tmp_qloop_30 ) -
                                       tmp_qloop_25 * tmp_qloop_27 );
                  const walberla::float64 q_tmp_3_3 =
                      tmp_qloop_17 * ( -abs_det_jac_affine_BLUE * tmp_qloop_37 * ( tmp_qloop_32 + tmp_qloop_34 ) +
                                       abs_det_jac_affine_BLUE * ( ( tmp_qloop_37 * tmp_qloop_37 ) * 2.0 +
                                                                   tmp_qloop_44 * ( tmp_qloop_39 + tmp_qloop_41 ) * 2.0 ) );
                  const walberla::float64 q_tmp_4_4 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_BLUE * tmp_qloop_46 *
                            ( jac_affine_inv_1_0_BLUE * tmp_qloop_45 * 0.66666666666666663 - tmp_qloop_34 ) +
                        abs_det_jac_affine_BLUE *
                            ( ( tmp_qloop_46 * tmp_qloop_46 ) * 2.0 +
                              tmp_qloop_47 * ( jac_affine_inv_1_1_BLUE * tmp_qloop_45 * 0.5 - tmp_qloop_41 ) * 2.0 ) );
                  const walberla::float64 q_tmp_5_5 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_BLUE * tmp_qloop_49 *
                            ( jac_affine_inv_0_0_BLUE * tmp_qloop_48 * 0.66666666666666663 - tmp_qloop_32 ) +
                        abs_det_jac_affine_BLUE *
                            ( ( tmp_qloop_49 * tmp_qloop_49 ) * 2.0 +
                              tmp_qloop_50 * ( jac_affine_inv_0_1_BLUE * tmp_qloop_48 * 0.5 - tmp_qloop_39 ) * 2.0 ) );
                  const walberla::float64 q_tmp_6_6 =
                      tmp_qloop_17 * ( -abs_det_jac_affine_BLUE * tmp_qloop_8 *
                                           ( tmp_qloop_6 * 0.66666666666666663 + tmp_qloop_7 * 0.66666666666666663 ) +
                                       abs_det_jac_affine_BLUE * ( tmp_qloop_5 * ( tmp_qloop_3 * 0.5 + tmp_qloop_4 * 0.5 ) * 2.0 +
                                                                   ( tmp_qloop_8 * tmp_qloop_8 ) * 2.0 ) );
                  const walberla::float64 q_tmp_7_7 =
                      tmp_qloop_17 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_18 * tmp_qloop_24 + tmp_qloop_22 * tmp_qloop_23 ) -
                                       tmp_qloop_21 * tmp_qloop_23 );
                  const walberla::float64 q_tmp_8_8 =
                      tmp_qloop_17 * ( abs_det_jac_affine_BLUE * ( tmp_qloop_25 * tmp_qloop_30 + tmp_qloop_28 * tmp_qloop_29 ) -
                                       tmp_qloop_27 * tmp_qloop_29 );
                  const walberla::float64 q_tmp_9_9 =
                      tmp_qloop_17 * ( -abs_det_jac_affine_BLUE * tmp_qloop_44 * ( tmp_qloop_51 + tmp_qloop_52 ) +
                                       abs_det_jac_affine_BLUE * ( tmp_qloop_37 * ( tmp_qloop_53 + tmp_qloop_54 ) * 2.0 +
                                                                   ( tmp_qloop_44 * tmp_qloop_44 ) * 2.0 ) );
                  const walberla::float64 q_tmp_10_10 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_BLUE * tmp_qloop_47 *
                            ( jac_affine_inv_1_1_BLUE * tmp_qloop_45 * 0.66666666666666663 - tmp_qloop_52 ) +
                        abs_det_jac_affine_BLUE *
                            ( tmp_qloop_46 * ( jac_affine_inv_1_0_BLUE * tmp_qloop_45 * 0.5 - tmp_qloop_54 ) * 2.0 +
                              ( tmp_qloop_47 * tmp_qloop_47 ) * 2.0 ) );
                  const walberla::float64 q_tmp_11_11 =
                      tmp_qloop_17 *
                      ( -abs_det_jac_affine_BLUE * tmp_qloop_50 *
                            ( jac_affine_inv_0_1_BLUE * tmp_qloop_48 * 0.66666666666666663 - tmp_qloop_51 ) +
                        abs_det_jac_affine_BLUE *
                            ( tmp_qloop_49 * ( jac_affine_inv_0_0_BLUE * tmp_qloop_48 * 0.5 - tmp_qloop_53 ) * 2.0 +
                              ( tmp_qloop_50 * tmp_qloop_50 ) * 2.0 ) );
                  q_acc_0_0   = q_acc_0_0 + q_tmp_0_0;
                  q_acc_1_1   = q_acc_1_1 + q_tmp_1_1;
                  q_acc_2_2   = q_acc_2_2 + q_tmp_2_2;
                  q_acc_3_3   = q_acc_3_3 + q_tmp_3_3;
                  q_acc_4_4   = q_acc_4_4 + q_tmp_4_4;
                  q_acc_5_5   = q_acc_5_5 + q_tmp_5_5;
                  q_acc_6_6   = q_acc_6_6 + q_tmp_6_6;
                  q_acc_7_7   = q_acc_7_7 + q_tmp_7_7;
                  q_acc_8_8   = q_acc_8_8 + q_tmp_8_8;
                  q_acc_9_9   = q_acc_9_9 + q_tmp_9_9;
                  q_acc_10_10 = q_acc_10_10 + q_tmp_10_10;
                  q_acc_11_11 = q_acc_11_11 + q_tmp_11_11;
               }
               const walberla::float64 elMatDiag_0  = q_acc_0_0;
               const walberla::float64 elMatDiag_1  = q_acc_1_1;
               const walberla::float64 elMatDiag_2  = q_acc_2_2;
               const walberla::float64 elMatDiag_3  = q_acc_3_3;
               const walberla::float64 elMatDiag_4  = q_acc_4_4;
               const walberla::float64 elMatDiag_5  = q_acc_5_5;
               const walberla::float64 elMatDiag_6  = q_acc_6_6;
               const walberla::float64 elMatDiag_7  = q_acc_7_7;
               const walberla::float64 elMatDiag_8  = q_acc_8_8;
               const walberla::float64 elMatDiag_9  = q_acc_9_9;
               const walberla::float64 elMatDiag_10 = q_acc_10_10;
               const walberla::float64 elMatDiag_11 = q_acc_11_11;
               _data_invDiag__vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_0 + _data_invDiag__vertex_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag__vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_1 + _data_invDiag__vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_invDiag__vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_2 + _data_invDiag__vertex_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               _data_invDiag__edge_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                     ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_3 + _data_invDiag__edge_0[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_4 +
                   _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatDiag_5 +
                   _data_invDiag__edge_0[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
               _data_invDiag__vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_6 + _data_invDiag__vertex_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag__vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_7 + _data_invDiag__vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_invDiag__vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_8 + _data_invDiag__vertex_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 2 ) -
                                                         ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) ) + 1];
               _data_invDiag__edge_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                     ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )] =
                   elMatDiag_9 + _data_invDiag__edge_1[ctr_0 + ( ctr_1 + 1 ) * ( micro_edges_per_macro_edge + 1 ) -
                                                       ( ( ( ctr_1 + 1 ) * ( ctr_1 + 2 ) ) / ( 2 ) )];
               _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1] =
                   elMatDiag_10 +
                   _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         2 * ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) ) + 1];
               _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) - ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                     ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )] =
                   elMatDiag_11 +
                   _data_invDiag__edge_1[ctr_0 + ctr_1 * ( micro_edges_per_macro_edge + 1 ) -
                                         ( ( ctr_1 * ( ctr_1 + 1 ) ) / ( 2 ) ) +
                                         ( ( micro_edges_per_macro_edge * ( micro_edges_per_macro_edge + 1 ) ) / ( 2 ) )];
            }
      }
   }
}

} // namespace operatorgeneration

} // namespace hyteg
