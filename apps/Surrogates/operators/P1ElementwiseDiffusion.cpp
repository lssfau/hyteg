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

#include "P1ElementwiseDiffusion.hpp"

#define FUNC_PREFIX

using std::max;

namespace hyteg {

namespace operatorgeneration {

P1ElementwiseDiffusion::P1ElementwiseDiffusion( const std::shared_ptr< PrimitiveStorage >& storage,
                                                size_t                                     minLevel,
                                                size_t                                     maxLevel )
: Operator( storage, minLevel, maxLevel )
{}

void P1ElementwiseDiffusion::apply( const P1Function< double >& src,
                                    const P1Function< double >& dst,
                                    uint_t                      level,
                                    DoFType                     flag,
                                    UpdateType                  updateType ) const
{
   this->startTiming( "apply" );

   // Make sure that halos are up-to-date
   this->timingTree_->start( "pre-communication" );
   if ( this->storage_->hasGlobalCells() )
   {
      // Note that the order of communication is important, since the face -> cell communication may overwrite
      // parts of the halos that carry the macro-vertex and macro-edge unknowns.
      src.communicate< Face, Cell >( level );
      src.communicate< Edge, Cell >( level );
      src.communicate< Vertex, Cell >( level );
   }
   else
   {
      communication::syncFunctionBetweenPrimitives( src, level, communication::syncDirection_t::LOW2HIGH );
   }
   this->timingTree_->stop( "pre-communication" );

   if ( updateType == Replace )
   {
      // We need to zero the destination array (including halos).
      // However, we must not zero out anything that is not flagged with the specified BCs.
      // Therefore, we first zero out everything that flagged, and then, later,
      // the halos of the highest dim primitives.
      dst.interpolate( walberla::numeric_cast< double >( 0 ), level, flag );
   }

   if ( storage_->hasGlobalCells() )
   {
      for ( auto& it : storage_->getCells() )
      {
         Cell& cell = *it.second;

         // get hold of the actual numerical data in the functions
         double* _data_src = cell.getData( src.getCellDataID() )->getPointer( level );
         double* _data_dst = cell.getData( dst.getCellDataID() )->getPointer( level );

         // Zero out dst halos only
         //
         // This is also necessary when using update type == Add.
         // During additive comm we then skip zeroing the data on the lower-dim primitives.
         for ( const auto& idx : vertexdof::macrocell::Iterator( level ) )
         {
            if ( !vertexdof::macrocell::isOnCellFace( idx, level ).empty() )
            {
               auto arrayIdx       = vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() );
               _data_dst[arrayIdx] = double( 0 );
            }
         }

         const auto   micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
         const auto   num_microfaces_per_face          = (int64_t) levelinfo::num_microfaces_per_face( level );
         const auto   micro_edges_per_macro_edge_float = (double) levelinfo::num_microedges_per_edge( level );
         const double macro_vertex_coord_id_0comp0     = (double) cell.getCoordinates()[0][0];
         const double macro_vertex_coord_id_0comp1     = (double) cell.getCoordinates()[0][1];
         const double macro_vertex_coord_id_0comp2     = (double) cell.getCoordinates()[0][2];
         const double macro_vertex_coord_id_1comp0     = (double) cell.getCoordinates()[1][0];
         const double macro_vertex_coord_id_1comp1     = (double) cell.getCoordinates()[1][1];
         const double macro_vertex_coord_id_1comp2     = (double) cell.getCoordinates()[1][2];
         const double macro_vertex_coord_id_2comp0     = (double) cell.getCoordinates()[2][0];
         const double macro_vertex_coord_id_2comp1     = (double) cell.getCoordinates()[2][1];
         const double macro_vertex_coord_id_2comp2     = (double) cell.getCoordinates()[2][2];
         const double macro_vertex_coord_id_3comp0     = (double) cell.getCoordinates()[3][0];
         const double macro_vertex_coord_id_3comp1     = (double) cell.getCoordinates()[3][1];
         const double macro_vertex_coord_id_3comp2     = (double) cell.getCoordinates()[3][2];

         this->timingTree_->start( "kernel" );

         apply_P1ElementwiseDiffusion_macro_3D(

             _data_dst,
             _data_src,
             macro_vertex_coord_id_0comp0,
             macro_vertex_coord_id_0comp1,
             macro_vertex_coord_id_0comp2,
             macro_vertex_coord_id_1comp0,
             macro_vertex_coord_id_1comp1,
             macro_vertex_coord_id_1comp2,
             macro_vertex_coord_id_2comp0,
             macro_vertex_coord_id_2comp1,
             macro_vertex_coord_id_2comp2,
             macro_vertex_coord_id_3comp0,
             macro_vertex_coord_id_3comp1,
             macro_vertex_coord_id_3comp2,
             micro_edges_per_macro_edge,
             micro_edges_per_macro_edge_float );

         this->timingTree_->stop( "kernel" );
      }

      // Push result to lower-dimensional primitives
      //
      this->timingTree_->start( "post-communication" );
      // Note: We could avoid communication here by implementing the apply() also for the respective
      //       lower dimensional primitives!
      dst.communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      dst.communicateAdditively< Cell, Vertex >( level, DoFType::All ^ flag, *storage_, updateType == Replace );
      this->timingTree_->stop( "post-communication" );
   }
   else
   {
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         // get hold of the actual numerical data in the functions
         double* _data_src = face.getData( src.getFaceDataID() )->getPointer( level );
         double* _data_dst = face.getData( dst.getFaceDataID() )->getPointer( level );

         // Zero out dst halos only
         //
         // This is also necessary when using update type == Add.
         // During additive comm we then skip zeroing the data on the lower-dim primitives.
         for ( const auto& idx : vertexdof::macroface::Iterator( level ) )
         {
            if ( vertexdof::macroface::isVertexOnBoundary( level, idx ) )
            {
               auto arrayIdx       = vertexdof::macroface::index( level, idx.x(), idx.y() );
               _data_dst[arrayIdx] = double( 0 );
            }
         }

         const auto   micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
         const auto   num_microfaces_per_face          = (int64_t) levelinfo::num_microfaces_per_face( level );
         const auto   micro_edges_per_macro_edge_float = (double) levelinfo::num_microedges_per_edge( level );
         const double macro_vertex_coord_id_0comp0     = (double) face.getCoordinates()[0][0];
         const double macro_vertex_coord_id_0comp1     = (double) face.getCoordinates()[0][1];
         const double macro_vertex_coord_id_1comp0     = (double) face.getCoordinates()[1][0];
         const double macro_vertex_coord_id_1comp1     = (double) face.getCoordinates()[1][1];
         const double macro_vertex_coord_id_2comp0     = (double) face.getCoordinates()[2][0];
         const double macro_vertex_coord_id_2comp1     = (double) face.getCoordinates()[2][1];

         this->timingTree_->start( "kernel" );

         apply_P1ElementwiseDiffusion_macro_2D(

             _data_dst,
             _data_src,
             macro_vertex_coord_id_0comp0,
             macro_vertex_coord_id_0comp1,
             macro_vertex_coord_id_1comp0,
             macro_vertex_coord_id_1comp1,
             macro_vertex_coord_id_2comp0,
             macro_vertex_coord_id_2comp1,
             micro_edges_per_macro_edge,
             micro_edges_per_macro_edge_float );

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
void P1ElementwiseDiffusion::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                       const P1Function< int64_t >&                src,
                                       const P1Function< int64_t >&                dst,
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

      this->timingTree_->stop( "pre-communication" );

      for ( auto& it : storage_->getCells() )
      {
         Cell& cell = *it.second;

         // get hold of the actual numerical data
         int64_t* _data_src = cell.getData( src.getCellDataID() )->getPointer( level );
         int64_t* _data_dst = cell.getData( dst.getCellDataID() )->getPointer( level );

         const auto   micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
         const auto   num_microfaces_per_face          = (int64_t) levelinfo::num_microfaces_per_face( level );
         const auto   micro_edges_per_macro_edge_float = (double) levelinfo::num_microedges_per_edge( level );
         const double macro_vertex_coord_id_0comp0     = (double) cell.getCoordinates()[0][0];
         const double macro_vertex_coord_id_0comp1     = (double) cell.getCoordinates()[0][1];
         const double macro_vertex_coord_id_0comp2     = (double) cell.getCoordinates()[0][2];
         const double macro_vertex_coord_id_1comp0     = (double) cell.getCoordinates()[1][0];
         const double macro_vertex_coord_id_1comp1     = (double) cell.getCoordinates()[1][1];
         const double macro_vertex_coord_id_1comp2     = (double) cell.getCoordinates()[1][2];
         const double macro_vertex_coord_id_2comp0     = (double) cell.getCoordinates()[2][0];
         const double macro_vertex_coord_id_2comp1     = (double) cell.getCoordinates()[2][1];
         const double macro_vertex_coord_id_2comp2     = (double) cell.getCoordinates()[2][2];
         const double macro_vertex_coord_id_3comp0     = (double) cell.getCoordinates()[3][0];
         const double macro_vertex_coord_id_3comp1     = (double) cell.getCoordinates()[3][1];
         const double macro_vertex_coord_id_3comp2     = (double) cell.getCoordinates()[3][2];

         this->timingTree_->start( "kernel" );

         toMatrix_P1ElementwiseDiffusion_macro_3D(

             _data_dst,
             _data_src,
             macro_vertex_coord_id_0comp0,
             macro_vertex_coord_id_0comp1,
             macro_vertex_coord_id_0comp2,
             macro_vertex_coord_id_1comp0,
             macro_vertex_coord_id_1comp1,
             macro_vertex_coord_id_1comp2,
             macro_vertex_coord_id_2comp0,
             macro_vertex_coord_id_2comp1,
             macro_vertex_coord_id_2comp2,
             macro_vertex_coord_id_3comp0,
             macro_vertex_coord_id_3comp1,
             macro_vertex_coord_id_3comp2,
             mat,
             micro_edges_per_macro_edge,
             micro_edges_per_macro_edge_float );

         this->timingTree_->stop( "kernel" );
      }
   }
   else
   {
      this->timingTree_->start( "pre-communication" );

      this->timingTree_->stop( "pre-communication" );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         // get hold of the actual numerical data
         int64_t* _data_src = face.getData( src.getFaceDataID() )->getPointer( level );
         int64_t* _data_dst = face.getData( dst.getFaceDataID() )->getPointer( level );

         const auto   micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
         const auto   num_microfaces_per_face          = (int64_t) levelinfo::num_microfaces_per_face( level );
         const auto   micro_edges_per_macro_edge_float = (double) levelinfo::num_microedges_per_edge( level );
         const double macro_vertex_coord_id_0comp0     = (double) face.getCoordinates()[0][0];
         const double macro_vertex_coord_id_0comp1     = (double) face.getCoordinates()[0][1];
         const double macro_vertex_coord_id_1comp0     = (double) face.getCoordinates()[1][0];
         const double macro_vertex_coord_id_1comp1     = (double) face.getCoordinates()[1][1];
         const double macro_vertex_coord_id_2comp0     = (double) face.getCoordinates()[2][0];
         const double macro_vertex_coord_id_2comp1     = (double) face.getCoordinates()[2][1];

         this->timingTree_->start( "kernel" );

         toMatrix_P1ElementwiseDiffusion_macro_2D(

             _data_dst,
             _data_src,
             macro_vertex_coord_id_0comp0,
             macro_vertex_coord_id_0comp1,
             macro_vertex_coord_id_1comp0,
             macro_vertex_coord_id_1comp1,
             macro_vertex_coord_id_2comp0,
             macro_vertex_coord_id_2comp1,
             mat,
             micro_edges_per_macro_edge,
             micro_edges_per_macro_edge_float );

         this->timingTree_->stop( "kernel" );
      }
   }
   this->stopTiming( "toMatrix" );
}
void P1ElementwiseDiffusion::computeInverseDiagonalOperatorValues()
{
   this->startTiming( "computeInverseDiagonalOperatorValues" );

   if ( invDiag_ == nullptr )
   {
      invDiag_ = std::make_shared< P1Function< double > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
   }

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      invDiag_->setToZero( level );

      if ( storage_->hasGlobalCells() )
      {
         this->timingTree_->start( "pre-communication" );

         this->timingTree_->stop( "pre-communication" );

         for ( auto& it : storage_->getCells() )
         {
            Cell& cell = *it.second;

            // get hold of the actual numerical data
            double* _data_invDiag_ = cell.getData( ( *invDiag_ ).getCellDataID() )->getPointer( level );

            const auto   micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
            const auto   num_microfaces_per_face          = (int64_t) levelinfo::num_microfaces_per_face( level );
            const auto   micro_edges_per_macro_edge_float = (double) levelinfo::num_microedges_per_edge( level );
            const double macro_vertex_coord_id_0comp0     = (double) cell.getCoordinates()[0][0];
            const double macro_vertex_coord_id_0comp1     = (double) cell.getCoordinates()[0][1];
            const double macro_vertex_coord_id_0comp2     = (double) cell.getCoordinates()[0][2];
            const double macro_vertex_coord_id_1comp0     = (double) cell.getCoordinates()[1][0];
            const double macro_vertex_coord_id_1comp1     = (double) cell.getCoordinates()[1][1];
            const double macro_vertex_coord_id_1comp2     = (double) cell.getCoordinates()[1][2];
            const double macro_vertex_coord_id_2comp0     = (double) cell.getCoordinates()[2][0];
            const double macro_vertex_coord_id_2comp1     = (double) cell.getCoordinates()[2][1];
            const double macro_vertex_coord_id_2comp2     = (double) cell.getCoordinates()[2][2];
            const double macro_vertex_coord_id_3comp0     = (double) cell.getCoordinates()[3][0];
            const double macro_vertex_coord_id_3comp1     = (double) cell.getCoordinates()[3][1];
            const double macro_vertex_coord_id_3comp2     = (double) cell.getCoordinates()[3][2];

            this->timingTree_->start( "kernel" );

            computeInverseDiagonalOperatorValues_P1ElementwiseDiffusion_macro_3D(

                _data_invDiag_,
                macro_vertex_coord_id_0comp0,
                macro_vertex_coord_id_0comp1,
                macro_vertex_coord_id_0comp2,
                macro_vertex_coord_id_1comp0,
                macro_vertex_coord_id_1comp1,
                macro_vertex_coord_id_1comp2,
                macro_vertex_coord_id_2comp0,
                macro_vertex_coord_id_2comp1,
                macro_vertex_coord_id_2comp2,
                macro_vertex_coord_id_3comp0,
                macro_vertex_coord_id_3comp1,
                macro_vertex_coord_id_3comp2,
                micro_edges_per_macro_edge,
                micro_edges_per_macro_edge_float );

            this->timingTree_->stop( "kernel" );
         }

         // Push result to lower-dimensional primitives
         //
         this->timingTree_->start( "post-communication" );
         // Note: We could avoid communication here by implementing the apply() also for the respective
         //       lower dimensional primitives!
         ( *invDiag_ ).communicateAdditively< Cell, Face >( level );
         ( *invDiag_ ).communicateAdditively< Cell, Edge >( level );
         ( *invDiag_ ).communicateAdditively< Cell, Vertex >( level );
         this->timingTree_->stop( "post-communication" );
         ( *invDiag_ ).invertElementwise( level );
      }
      else
      {
         this->timingTree_->start( "pre-communication" );

         this->timingTree_->stop( "pre-communication" );

         for ( auto& it : storage_->getFaces() )
         {
            Face& face = *it.second;

            // get hold of the actual numerical data
            double* _data_invDiag_ = face.getData( ( *invDiag_ ).getFaceDataID() )->getPointer( level );

            const auto   micro_edges_per_macro_edge       = (int64_t) levelinfo::num_microedges_per_edge( level );
            const auto   num_microfaces_per_face          = (int64_t) levelinfo::num_microfaces_per_face( level );
            const auto   micro_edges_per_macro_edge_float = (double) levelinfo::num_microedges_per_edge( level );
            const double macro_vertex_coord_id_0comp0     = (double) face.getCoordinates()[0][0];
            const double macro_vertex_coord_id_0comp1     = (double) face.getCoordinates()[0][1];
            const double macro_vertex_coord_id_1comp0     = (double) face.getCoordinates()[1][0];
            const double macro_vertex_coord_id_1comp1     = (double) face.getCoordinates()[1][1];
            const double macro_vertex_coord_id_2comp0     = (double) face.getCoordinates()[2][0];
            const double macro_vertex_coord_id_2comp1     = (double) face.getCoordinates()[2][1];

            this->timingTree_->start( "kernel" );

            computeInverseDiagonalOperatorValues_P1ElementwiseDiffusion_macro_2D(

                _data_invDiag_,
                macro_vertex_coord_id_0comp0,
                macro_vertex_coord_id_0comp1,
                macro_vertex_coord_id_1comp0,
                macro_vertex_coord_id_1comp1,
                macro_vertex_coord_id_2comp0,
                macro_vertex_coord_id_2comp1,
                micro_edges_per_macro_edge,
                micro_edges_per_macro_edge_float );

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
std::shared_ptr< P1Function< double > > P1ElementwiseDiffusion::getInverseDiagonalValues() const
{
   return invDiag_;
}
void P1ElementwiseDiffusion::apply_P1ElementwiseDiffusion_macro_3D( double* RESTRICT const _data_dst,
                                                                    double* RESTRICT const _data_src,
                                                                    const double           macro_vertex_coord_id_0comp0,
                                                                    const double           macro_vertex_coord_id_0comp1,
                                                                    const double           macro_vertex_coord_id_0comp2,
                                                                    const double           macro_vertex_coord_id_1comp0,
                                                                    const double           macro_vertex_coord_id_1comp1,
                                                                    const double           macro_vertex_coord_id_1comp2,
                                                                    const double           macro_vertex_coord_id_2comp0,
                                                                    const double           macro_vertex_coord_id_2comp1,
                                                                    const double           macro_vertex_coord_id_2comp2,
                                                                    const double           macro_vertex_coord_id_3comp0,
                                                                    const double           macro_vertex_coord_id_3comp1,
                                                                    const double           macro_vertex_coord_id_3comp2,
                                                                    const int64_t          micro_edges_per_macro_edge,
                                                                    const double micro_edges_per_macro_edge_float ) const
{
   {
      {
         /* CellType.WHITE_UP */
         const double tmp_coords_jac_0__4   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__4   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__4   = tmp_coords_jac_1__4 * 0.0;
         const double tmp_coords_jac_3__4   = tmp_coords_jac_0__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_4__4   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__4   = tmp_coords_jac_2__4 * tmp_coords_jac_4__4;
         const double tmp_coords_jac_6__4   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__4   = tmp_coords_jac_2__4 * tmp_coords_jac_6__4;
         const double tmp_coords_jac_8__4   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_5__4 + tmp_coords_jac_7__4;
         const double tmp_coords_jac_9__4   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__4  = tmp_coords_jac_2__4 * tmp_coords_jac_9__4;
         const double tmp_coords_jac_11__4  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12__4  = tmp_coords_jac_11__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_13__4  = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_14__4  = tmp_coords_jac_13__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_15__4  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12__4 + tmp_coords_jac_14__4;
         const double tmp_coords_jac_16__4  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_17__4  = tmp_coords_jac_16__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_18__4  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__4  = tmp_coords_jac_18__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_20__4  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_21__4  = tmp_coords_jac_2__4 * tmp_coords_jac_20__4;
         const double tmp_coords_jac_22__4  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_19__4 + tmp_coords_jac_21__4;
         const double tmp_coords_jac_23__4  = tmp_coords_jac_1__4 * 1.0;
         const double tmp_coords_jac_24__4  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__4;
         const double tmp_coords_jac_25__4  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10__4;
         const double tmp_coords_jac_26__4  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_17__4;
         const double p_affine_const_0_0__4 = tmp_coords_jac_3__4 + tmp_coords_jac_8__4;
         const double p_affine_const_0_1__4 = tmp_coords_jac_10__4 + tmp_coords_jac_15__4;
         const double p_affine_const_0_2__4 = tmp_coords_jac_17__4 + tmp_coords_jac_22__4;
         const double p_affine_const_1_0__4 = tmp_coords_jac_8__4 + tmp_coords_jac_0__4 * tmp_coords_jac_23__4;
         const double p_affine_const_1_1__4 = tmp_coords_jac_15__4 + tmp_coords_jac_23__4 * tmp_coords_jac_9__4;
         const double p_affine_const_1_2__4 = tmp_coords_jac_22__4 + tmp_coords_jac_16__4 * tmp_coords_jac_23__4;
         const double p_affine_const_2_0__4 =
             tmp_coords_jac_24__4 + tmp_coords_jac_7__4 + tmp_coords_jac_23__4 * tmp_coords_jac_4__4;
         const double p_affine_const_2_1__4 =
             tmp_coords_jac_14__4 + tmp_coords_jac_25__4 + tmp_coords_jac_11__4 * tmp_coords_jac_23__4;
         const double p_affine_const_2_2__4 =
             tmp_coords_jac_21__4 + tmp_coords_jac_26__4 + tmp_coords_jac_18__4 * tmp_coords_jac_23__4;
         const double p_affine_const_3_0__4 =
             tmp_coords_jac_24__4 + tmp_coords_jac_5__4 + tmp_coords_jac_23__4 * tmp_coords_jac_6__4;
         const double p_affine_const_3_1__4 =
             tmp_coords_jac_12__4 + tmp_coords_jac_25__4 + tmp_coords_jac_13__4 * tmp_coords_jac_23__4;
         const double p_affine_const_3_2__4 =
             tmp_coords_jac_19__4 + tmp_coords_jac_26__4 + tmp_coords_jac_20__4 * tmp_coords_jac_23__4;
         const double jac_affine_0_0__4    = p_affine_const_1_0__4 - p_affine_const_0_0__4;
         const double jac_affine_0_1__4    = p_affine_const_2_0__4 - p_affine_const_0_0__4;
         const double jac_affine_0_2__4    = p_affine_const_3_0__4 - p_affine_const_0_0__4;
         const double jac_affine_1_0__4    = p_affine_const_1_1__4 - p_affine_const_0_1__4;
         const double jac_affine_1_1__4    = p_affine_const_2_1__4 - p_affine_const_0_1__4;
         const double tmp_coords_jac_31__2 = jac_affine_0_2__4 * jac_affine_1_1__4;
         const double jac_affine_1_2__4    = p_affine_const_3_1__4 - p_affine_const_0_1__4;
         const double tmp_coords_jac_29__4 = jac_affine_0_1__4 * jac_affine_1_2__4;
         const double jac_affine_2_0__4    = p_affine_const_1_2__4 - p_affine_const_0_2__4;
         const double jac_affine_2_1__4    = p_affine_const_2_2__4 - p_affine_const_0_2__4;
         const double tmp_coords_jac_28__4 = jac_affine_1_2__4 * jac_affine_2_1__4;
         const double jac_affine_2_2__4    = p_affine_const_3_2__4 - p_affine_const_0_2__4;
         const double tmp_coords_jac_27__4 = jac_affine_1_1__4 * jac_affine_2_2__4;
         const double tmp_coords_jac_30__4 = jac_affine_0_1__4 * jac_affine_2_2__4;
         const double tmp_coords_jac_32__2 = jac_affine_0_0__4 * tmp_coords_jac_27__4 + jac_affine_2_0__4 * tmp_coords_jac_29__4 -
                                             jac_affine_0_0__4 * tmp_coords_jac_28__4 - jac_affine_1_0__4 * tmp_coords_jac_30__4 -
                                             jac_affine_2_0__4 * tmp_coords_jac_31__2 +
                                             jac_affine_0_2__4 * jac_affine_1_0__4 * jac_affine_2_1__4;
         const double tmp_coords_jac_33__2  = 1.0 / tmp_coords_jac_32__2;
         const double jac_affine_inv_0_0__4 = tmp_coords_jac_33__2 * ( tmp_coords_jac_27__4 - tmp_coords_jac_28__4 );
         const double jac_affine_inv_0_1__4 =
             tmp_coords_jac_33__2 * ( -1.0 * tmp_coords_jac_30__4 + jac_affine_0_2__4 * jac_affine_2_1__4 );
         const double jac_affine_inv_0_2__4 = tmp_coords_jac_33__2 * ( tmp_coords_jac_29__4 - tmp_coords_jac_31__2 );
         const double jac_affine_inv_1_0__4 =
             tmp_coords_jac_33__2 * ( jac_affine_1_2__4 * jac_affine_2_0__4 - jac_affine_1_0__4 * jac_affine_2_2__4 );
         const double jac_affine_inv_1_1__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_0__4 * jac_affine_2_2__4 - jac_affine_0_2__4 * jac_affine_2_0__4 );
         const double jac_affine_inv_1_2__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_2__4 * jac_affine_1_0__4 - jac_affine_0_0__4 * jac_affine_1_2__4 );
         const double jac_affine_inv_2_0__4 =
             tmp_coords_jac_33__2 * ( jac_affine_1_0__4 * jac_affine_2_1__4 - jac_affine_1_1__4 * jac_affine_2_0__4 );
         const double jac_affine_inv_2_1__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_1__4 * jac_affine_2_0__4 - jac_affine_0_0__4 * jac_affine_2_1__4 );
         const double jac_affine_inv_2_2__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_0__4 * jac_affine_1_1__4 - jac_affine_0_1__4 * jac_affine_1_0__4 );
         const double abs_det_jac_affine__4 = abs( tmp_coords_jac_32__2 );
         for ( int64_t ctr_2__4 = 0LL; ctr_2__4 < micro_edges_per_macro_edge; ctr_2__4 += 1LL )
         {
            for ( int64_t ctr_1__4 = 0LL; ctr_1__4 < -1LL * ctr_2__4 + micro_edges_per_macro_edge; ctr_1__4 += 1LL )
            {
               for ( int64_t ctr_0__4 = 0LL; ctr_0__4 < -1LL * ctr_1__4 - ctr_2__4 + micro_edges_per_macro_edge; ctr_0__4 += 1LL )
               {
                  const double p_affine_0_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_0_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_0_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_1_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_1_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_1_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_2_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_2_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_2_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_3_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__4 );
                  const double p_affine_3_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__4 );
                  const double p_affine_3_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__4 );
                  const double src_dof_0__4 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_1__4 =
                      _data_src[1LL - ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL -
                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_2__4 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__4 ) * ( 2LL + ctr_1__4 ) / 2LL ) -
                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__4 ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) + ctr_0__4 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_3__4 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                ( -1LL * ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double tmp_kernel_op_0__4 = -1.0 * jac_affine_inv_0_0__4 - jac_affine_inv_1_0__4 - jac_affine_inv_2_0__4;
                  const double tmp_kernel_op_1__4 = -1.0 * jac_affine_inv_0_1__4 - jac_affine_inv_1_1__4 - jac_affine_inv_2_1__4;
                  const double tmp_kernel_op_2__4 = -1.0 * jac_affine_inv_0_2__4 - jac_affine_inv_1_2__4 - jac_affine_inv_2_2__4;
                  const double tmp_kernel_op_3__4 = 0.16666666666666663 * abs_det_jac_affine__4;
                  const double tmp_kernel_op_4__4 = src_dof_0__4 * tmp_kernel_op_3__4;
                  const double tmp_kernel_op_5__4 = jac_affine_inv_0_0__4 * tmp_kernel_op_0__4 +
                                                    jac_affine_inv_0_1__4 * tmp_kernel_op_1__4 +
                                                    jac_affine_inv_0_2__4 * tmp_kernel_op_2__4;
                  const double tmp_kernel_op_6__4 = src_dof_1__4 * tmp_kernel_op_3__4;
                  const double tmp_kernel_op_7__4 = jac_affine_inv_1_0__4 * tmp_kernel_op_0__4 +
                                                    jac_affine_inv_1_1__4 * tmp_kernel_op_1__4 +
                                                    jac_affine_inv_1_2__4 * tmp_kernel_op_2__4;
                  const double tmp_kernel_op_8__4 = src_dof_2__4 * tmp_kernel_op_3__4;
                  const double tmp_kernel_op_9__4 = jac_affine_inv_2_0__4 * tmp_kernel_op_0__4 +
                                                    jac_affine_inv_2_1__4 * tmp_kernel_op_1__4 +
                                                    jac_affine_inv_2_2__4 * tmp_kernel_op_2__4;
                  const double tmp_kernel_op_10__4 = src_dof_3__4 * tmp_kernel_op_3__4;
                  const double tmp_kernel_op_11__4 = jac_affine_inv_0_0__4 * jac_affine_inv_1_0__4 +
                                                     jac_affine_inv_0_1__4 * jac_affine_inv_1_1__4 +
                                                     jac_affine_inv_0_2__4 * jac_affine_inv_1_2__4;
                  const double tmp_kernel_op_12__4 = jac_affine_inv_0_0__4 * jac_affine_inv_2_0__4 +
                                                     jac_affine_inv_0_1__4 * jac_affine_inv_2_1__4 +
                                                     jac_affine_inv_0_2__4 * jac_affine_inv_2_2__4;
                  const double tmp_kernel_op_13__4 = jac_affine_inv_1_0__4 * jac_affine_inv_2_0__4 +
                                                     jac_affine_inv_1_1__4 * jac_affine_inv_2_1__4 +
                                                     jac_affine_inv_1_2__4 * jac_affine_inv_2_2__4;
                  const double elMatVec_0__4 =
                      tmp_kernel_op_10__4 * tmp_kernel_op_9__4 +
                      tmp_kernel_op_4__4 * ( tmp_kernel_op_0__4 * tmp_kernel_op_0__4 + tmp_kernel_op_1__4 * tmp_kernel_op_1__4 +
                                             tmp_kernel_op_2__4 * tmp_kernel_op_2__4 ) +
                      tmp_kernel_op_5__4 * tmp_kernel_op_6__4 + tmp_kernel_op_7__4 * tmp_kernel_op_8__4;
                  const double elMatVec_1__4 =
                      tmp_kernel_op_10__4 * tmp_kernel_op_12__4 + tmp_kernel_op_11__4 * tmp_kernel_op_8__4 +
                      tmp_kernel_op_4__4 * tmp_kernel_op_5__4 +
                      tmp_kernel_op_6__4 *
                          ( jac_affine_inv_0_0__4 * jac_affine_inv_0_0__4 + jac_affine_inv_0_1__4 * jac_affine_inv_0_1__4 +
                            jac_affine_inv_0_2__4 * jac_affine_inv_0_2__4 );
                  const double elMatVec_2__4 =
                      tmp_kernel_op_10__4 * tmp_kernel_op_13__4 + tmp_kernel_op_11__4 * tmp_kernel_op_6__4 +
                      tmp_kernel_op_4__4 * tmp_kernel_op_7__4 +
                      tmp_kernel_op_8__4 *
                          ( jac_affine_inv_1_0__4 * jac_affine_inv_1_0__4 + jac_affine_inv_1_1__4 * jac_affine_inv_1_1__4 +
                            jac_affine_inv_1_2__4 * jac_affine_inv_1_2__4 );
                  const double elMatVec_3__4 = tmp_kernel_op_10__4 * ( jac_affine_inv_2_0__4 * jac_affine_inv_2_0__4 +
                                                                       jac_affine_inv_2_1__4 * jac_affine_inv_2_1__4 +
                                                                       jac_affine_inv_2_2__4 * jac_affine_inv_2_2__4 ) +
                                               tmp_kernel_op_12__4 * tmp_kernel_op_6__4 +
                                               tmp_kernel_op_13__4 * tmp_kernel_op_8__4 + tmp_kernel_op_4__4 * tmp_kernel_op_9__4;
                  _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                            ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                            ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_0__4 + _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL -
                            ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                            ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_1__4 + _data_dst[1LL - ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL -
                                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ( 2LL + ctr_1__4 ) / 2LL ) -
                            ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__4 ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) + ctr_0__4 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_2__4 + _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ( 2LL + ctr_1__4 ) / 2LL ) -
                                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__4 ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) + ctr_0__4 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                            ( -1LL * ctr_2__4 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_3__4 + _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                                ( -1LL * ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.WHITE_DOWN */
         const double tmp_coords_jac_0__3   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__3   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__3   = tmp_coords_jac_1__3 * 0.0;
         const double tmp_coords_jac_3__3   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4__3   = tmp_coords_jac_1__3 * 1.0;
         const double tmp_coords_jac_5__3   = tmp_coords_jac_3__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_6__3   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__3   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_4__3 * tmp_coords_jac_6__3;
         const double tmp_coords_jac_8__3   = tmp_coords_jac_5__3 + tmp_coords_jac_7__3;
         const double tmp_coords_jac_9__3   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__3  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__3  = tmp_coords_jac_10__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_12__3  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__3  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_14__3  = tmp_coords_jac_11__3 + tmp_coords_jac_13__3;
         const double tmp_coords_jac_15__3  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__3  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_17__3  = tmp_coords_jac_16__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_18__3  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__3  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_20__3  = tmp_coords_jac_17__3 + tmp_coords_jac_19__3;
         const double tmp_coords_jac_21__3  = tmp_coords_jac_0__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_22__3  = tmp_coords_jac_4__3 * tmp_coords_jac_9__3;
         const double tmp_coords_jac_23__3  = tmp_coords_jac_15__3 * tmp_coords_jac_4__3;
         const double p_affine_const_0_0__3 = tmp_coords_jac_8__3 + tmp_coords_jac_0__3 * tmp_coords_jac_2__3;
         const double p_affine_const_0_1__3 = tmp_coords_jac_14__3 + tmp_coords_jac_2__3 * tmp_coords_jac_9__3;
         const double p_affine_const_0_2__3 = tmp_coords_jac_20__3 + tmp_coords_jac_15__3 * tmp_coords_jac_2__3;
         const double p_affine_const_1_0__3 =
             tmp_coords_jac_21__3 + tmp_coords_jac_7__3 + tmp_coords_jac_2__3 * tmp_coords_jac_3__3;
         const double p_affine_const_1_1__3 =
             tmp_coords_jac_13__3 + tmp_coords_jac_22__3 + tmp_coords_jac_10__3 * tmp_coords_jac_2__3;
         const double p_affine_const_1_2__3 =
             tmp_coords_jac_19__3 + tmp_coords_jac_23__3 + tmp_coords_jac_16__3 * tmp_coords_jac_2__3;
         const double p_affine_const_2_0__3 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_21__3 + tmp_coords_jac_5__3 +
                                              tmp_coords_jac_2__3 * tmp_coords_jac_6__3;
         const double p_affine_const_2_1__3 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_11__3 + tmp_coords_jac_22__3 +
                                              tmp_coords_jac_12__3 * tmp_coords_jac_2__3;
         const double p_affine_const_2_2__3 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_17__3 + tmp_coords_jac_23__3 +
                                              tmp_coords_jac_18__3 * tmp_coords_jac_2__3;
         const double p_affine_const_3_0__3 = tmp_coords_jac_21__3 + tmp_coords_jac_8__3;
         const double p_affine_const_3_1__3 = tmp_coords_jac_14__3 + tmp_coords_jac_22__3;
         const double p_affine_const_3_2__3 = tmp_coords_jac_20__3 + tmp_coords_jac_23__3;
         const double jac_affine_0_0__3     = p_affine_const_1_0__3 - p_affine_const_0_0__3;
         const double jac_affine_0_1__3     = p_affine_const_2_0__3 - p_affine_const_0_0__3;
         const double jac_affine_0_2__3     = p_affine_const_3_0__3 - p_affine_const_0_0__3;
         const double jac_affine_1_0__3     = p_affine_const_1_1__3 - p_affine_const_0_1__3;
         const double jac_affine_1_1__3     = p_affine_const_2_1__3 - p_affine_const_0_1__3;
         const double tmp_coords_jac_28__3  = jac_affine_0_2__3 * jac_affine_1_1__3;
         const double jac_affine_1_2__3     = p_affine_const_3_1__3 - p_affine_const_0_1__3;
         const double tmp_coords_jac_26__3  = jac_affine_0_1__3 * jac_affine_1_2__3;
         const double jac_affine_2_0__3     = p_affine_const_1_2__3 - p_affine_const_0_2__3;
         const double jac_affine_2_1__3     = p_affine_const_2_2__3 - p_affine_const_0_2__3;
         const double tmp_coords_jac_25__3  = jac_affine_1_2__3 * jac_affine_2_1__3;
         const double jac_affine_2_2__3     = p_affine_const_3_2__3 - p_affine_const_0_2__3;
         const double tmp_coords_jac_24__3  = jac_affine_1_1__3 * jac_affine_2_2__3;
         const double tmp_coords_jac_27__3  = jac_affine_0_1__3 * jac_affine_2_2__3;
         const double tmp_coords_jac_29__3 = jac_affine_0_0__3 * tmp_coords_jac_24__3 + jac_affine_2_0__3 * tmp_coords_jac_26__3 -
                                             jac_affine_0_0__3 * tmp_coords_jac_25__3 - jac_affine_1_0__3 * tmp_coords_jac_27__3 -
                                             jac_affine_2_0__3 * tmp_coords_jac_28__3 +
                                             jac_affine_0_2__3 * jac_affine_1_0__3 * jac_affine_2_1__3;
         const double tmp_coords_jac_30__3  = 1.0 / tmp_coords_jac_29__3;
         const double jac_affine_inv_0_0__3 = tmp_coords_jac_30__3 * ( tmp_coords_jac_24__3 - tmp_coords_jac_25__3 );
         const double jac_affine_inv_0_1__3 =
             tmp_coords_jac_30__3 * ( -1.0 * tmp_coords_jac_27__3 + jac_affine_0_2__3 * jac_affine_2_1__3 );
         const double jac_affine_inv_0_2__3 = tmp_coords_jac_30__3 * ( tmp_coords_jac_26__3 - tmp_coords_jac_28__3 );
         const double jac_affine_inv_1_0__3 =
             tmp_coords_jac_30__3 * ( jac_affine_1_2__3 * jac_affine_2_0__3 - jac_affine_1_0__3 * jac_affine_2_2__3 );
         const double jac_affine_inv_1_1__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_0__3 * jac_affine_2_2__3 - jac_affine_0_2__3 * jac_affine_2_0__3 );
         const double jac_affine_inv_1_2__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_2__3 * jac_affine_1_0__3 - jac_affine_0_0__3 * jac_affine_1_2__3 );
         const double jac_affine_inv_2_0__3 =
             tmp_coords_jac_30__3 * ( jac_affine_1_0__3 * jac_affine_2_1__3 - jac_affine_1_1__3 * jac_affine_2_0__3 );
         const double jac_affine_inv_2_1__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_1__3 * jac_affine_2_0__3 - jac_affine_0_0__3 * jac_affine_2_1__3 );
         const double jac_affine_inv_2_2__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_0__3 * jac_affine_1_1__3 - jac_affine_0_1__3 * jac_affine_1_0__3 );
         const double abs_det_jac_affine__3 = abs( tmp_coords_jac_29__3 );
         for ( int64_t ctr_2__3 = 0LL; ctr_2__3 < micro_edges_per_macro_edge; ctr_2__3 += 1LL )
         {
            for ( int64_t ctr_1__3 = 0LL; ctr_1__3 < -1LL * ctr_2__3 + micro_edges_per_macro_edge; ctr_1__3 += 1LL )
            {
               for ( int64_t ctr_0__3 = 0LL; ctr_0__3 < -2LL - ctr_1__3 - ctr_2__3 + micro_edges_per_macro_edge; ctr_0__3 += 1LL )
               {
                  const double p_affine_0_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__3;
                  const double p_affine_0_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__3;
                  const double p_affine_0_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__3;
                  const double p_affine_1_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_1_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_1_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double src_dof_0__3 =
                      _data_src[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__3 ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_1__3 =
                      _data_src[1LL - ( 1LL + ctr_1__3 ) * ctr_1__3 / 2LL -
                                ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ctr_1__3 + ctr_0__3 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_2__3 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL ) -
                                ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_3__3 =
                      _data_src[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double tmp_kernel_op_0__3 = -1.0 * jac_affine_inv_0_0__3 - jac_affine_inv_1_0__3 - jac_affine_inv_2_0__3;
                  const double tmp_kernel_op_1__3 = -1.0 * jac_affine_inv_0_1__3 - jac_affine_inv_1_1__3 - jac_affine_inv_2_1__3;
                  const double tmp_kernel_op_2__3 = -1.0 * jac_affine_inv_0_2__3 - jac_affine_inv_1_2__3 - jac_affine_inv_2_2__3;
                  const double tmp_kernel_op_3__3 = 0.16666666666666663 * abs_det_jac_affine__3;
                  const double tmp_kernel_op_4__3 = src_dof_0__3 * tmp_kernel_op_3__3;
                  const double tmp_kernel_op_5__3 = jac_affine_inv_0_0__3 * tmp_kernel_op_0__3 +
                                                    jac_affine_inv_0_1__3 * tmp_kernel_op_1__3 +
                                                    jac_affine_inv_0_2__3 * tmp_kernel_op_2__3;
                  const double tmp_kernel_op_6__3 = src_dof_1__3 * tmp_kernel_op_3__3;
                  const double tmp_kernel_op_7__3 = jac_affine_inv_1_0__3 * tmp_kernel_op_0__3 +
                                                    jac_affine_inv_1_1__3 * tmp_kernel_op_1__3 +
                                                    jac_affine_inv_1_2__3 * tmp_kernel_op_2__3;
                  const double tmp_kernel_op_8__3 = src_dof_2__3 * tmp_kernel_op_3__3;
                  const double tmp_kernel_op_9__3 = jac_affine_inv_2_0__3 * tmp_kernel_op_0__3 +
                                                    jac_affine_inv_2_1__3 * tmp_kernel_op_1__3 +
                                                    jac_affine_inv_2_2__3 * tmp_kernel_op_2__3;
                  const double tmp_kernel_op_10__3 = src_dof_3__3 * tmp_kernel_op_3__3;
                  const double tmp_kernel_op_11__3 = jac_affine_inv_0_0__3 * jac_affine_inv_1_0__3 +
                                                     jac_affine_inv_0_1__3 * jac_affine_inv_1_1__3 +
                                                     jac_affine_inv_0_2__3 * jac_affine_inv_1_2__3;
                  const double tmp_kernel_op_12__3 = jac_affine_inv_0_0__3 * jac_affine_inv_2_0__3 +
                                                     jac_affine_inv_0_1__3 * jac_affine_inv_2_1__3 +
                                                     jac_affine_inv_0_2__3 * jac_affine_inv_2_2__3;
                  const double tmp_kernel_op_13__3 = jac_affine_inv_1_0__3 * jac_affine_inv_2_0__3 +
                                                     jac_affine_inv_1_1__3 * jac_affine_inv_2_1__3 +
                                                     jac_affine_inv_1_2__3 * jac_affine_inv_2_2__3;
                  const double elMatVec_0__3 =
                      tmp_kernel_op_10__3 * tmp_kernel_op_9__3 +
                      tmp_kernel_op_4__3 * ( tmp_kernel_op_0__3 * tmp_kernel_op_0__3 + tmp_kernel_op_1__3 * tmp_kernel_op_1__3 +
                                             tmp_kernel_op_2__3 * tmp_kernel_op_2__3 ) +
                      tmp_kernel_op_5__3 * tmp_kernel_op_6__3 + tmp_kernel_op_7__3 * tmp_kernel_op_8__3;
                  const double elMatVec_1__3 =
                      tmp_kernel_op_10__3 * tmp_kernel_op_12__3 + tmp_kernel_op_11__3 * tmp_kernel_op_8__3 +
                      tmp_kernel_op_4__3 * tmp_kernel_op_5__3 +
                      tmp_kernel_op_6__3 *
                          ( jac_affine_inv_0_0__3 * jac_affine_inv_0_0__3 + jac_affine_inv_0_1__3 * jac_affine_inv_0_1__3 +
                            jac_affine_inv_0_2__3 * jac_affine_inv_0_2__3 );
                  const double elMatVec_2__3 =
                      tmp_kernel_op_10__3 * tmp_kernel_op_13__3 + tmp_kernel_op_11__3 * tmp_kernel_op_6__3 +
                      tmp_kernel_op_4__3 * tmp_kernel_op_7__3 +
                      tmp_kernel_op_8__3 *
                          ( jac_affine_inv_1_0__3 * jac_affine_inv_1_0__3 + jac_affine_inv_1_1__3 * jac_affine_inv_1_1__3 +
                            jac_affine_inv_1_2__3 * jac_affine_inv_1_2__3 );
                  const double elMatVec_3__3 = tmp_kernel_op_10__3 * ( jac_affine_inv_2_0__3 * jac_affine_inv_2_0__3 +
                                                                       jac_affine_inv_2_1__3 * jac_affine_inv_2_1__3 +
                                                                       jac_affine_inv_2_2__3 * jac_affine_inv_2_2__3 ) +
                                               tmp_kernel_op_12__3 * tmp_kernel_op_6__3 +
                                               tmp_kernel_op_13__3 * tmp_kernel_op_8__3 + tmp_kernel_op_4__3 * tmp_kernel_op_9__3;
                  _data_dst[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                            ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__3 ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_0__3 + _data_dst[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                                ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__3 ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1__3 ) * ctr_1__3 / 2LL -
                            ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ctr_1__3 + ctr_0__3 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_1__3 + _data_dst[1LL - ( 1LL + ctr_1__3 ) * ctr_1__3 / 2LL -
                                                ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ctr_1__3 + ctr_0__3 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL ) -
                            ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_2__3 + _data_dst[-1LL * ( ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL ) -
                                                ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                            ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_3__3 + _data_dst[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                                ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.BLUE_UP */
         const double tmp_coords_jac_0__2   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__2   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__2   = tmp_coords_jac_1__2 * 0.0;
         const double tmp_coords_jac_3__2   = tmp_coords_jac_0__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_4__2   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__2   = tmp_coords_jac_1__2 * 1.0;
         const double tmp_coords_jac_6__2   = tmp_coords_jac_4__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_7__2   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_8__2   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__2 * tmp_coords_jac_7__2;
         const double tmp_coords_jac_9__2   = tmp_coords_jac_6__2 + tmp_coords_jac_8__2;
         const double tmp_coords_jac_10__2  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__2  = tmp_coords_jac_10__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_12__2  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__2  = tmp_coords_jac_12__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_14__2  = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_15__2  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_14__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_16__2  = tmp_coords_jac_13__2 + tmp_coords_jac_15__2;
         const double tmp_coords_jac_17__2  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18__2  = tmp_coords_jac_17__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_19__2  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_20__2  = tmp_coords_jac_19__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_21__2  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_22__2  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_2__2 * tmp_coords_jac_21__2;
         const double tmp_coords_jac_23__2  = tmp_coords_jac_20__2 + tmp_coords_jac_22__2;
         const double tmp_coords_jac_24__2  = tmp_coords_jac_0__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_25__2  = tmp_coords_jac_10__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_26__2  = tmp_coords_jac_17__2 * tmp_coords_jac_5__2;
         const double p_affine_const_0_0__2 = tmp_coords_jac_3__2 + tmp_coords_jac_9__2;
         const double p_affine_const_0_1__2 = tmp_coords_jac_11__2 + tmp_coords_jac_16__2;
         const double p_affine_const_0_2__2 = tmp_coords_jac_18__2 + tmp_coords_jac_23__2;
         const double p_affine_const_1_0__2 =
             tmp_coords_jac_24__2 + tmp_coords_jac_8__2 + tmp_coords_jac_2__2 * tmp_coords_jac_4__2;
         const double p_affine_const_1_1__2 =
             tmp_coords_jac_15__2 + tmp_coords_jac_25__2 + tmp_coords_jac_12__2 * tmp_coords_jac_2__2;
         const double p_affine_const_1_2__2 =
             tmp_coords_jac_22__2 + tmp_coords_jac_26__2 + tmp_coords_jac_19__2 * tmp_coords_jac_2__2;
         const double p_affine_const_2_0__2 = tmp_coords_jac_24__2 + tmp_coords_jac_9__2;
         const double p_affine_const_2_1__2 = tmp_coords_jac_16__2 + tmp_coords_jac_25__2;
         const double p_affine_const_2_2__2 = tmp_coords_jac_23__2 + tmp_coords_jac_26__2;
         const double p_affine_const_3_0__2 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__2 + tmp_coords_jac_6__2 + tmp_coords_jac_5__2 * tmp_coords_jac_7__2;
         const double p_affine_const_3_1__2 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_11__2 + tmp_coords_jac_13__2 +
                                              tmp_coords_jac_14__2 * tmp_coords_jac_5__2;
         const double p_affine_const_3_2__2 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18__2 + tmp_coords_jac_20__2 +
                                              tmp_coords_jac_21__2 * tmp_coords_jac_5__2;
         const double jac_affine_0_0__2    = p_affine_const_1_0__2 - p_affine_const_0_0__2;
         const double jac_affine_0_1__2    = p_affine_const_2_0__2 - p_affine_const_0_0__2;
         const double jac_affine_0_2__2    = p_affine_const_3_0__2 - p_affine_const_0_0__2;
         const double jac_affine_1_0__2    = p_affine_const_1_1__2 - p_affine_const_0_1__2;
         const double jac_affine_1_1__2    = p_affine_const_2_1__2 - p_affine_const_0_1__2;
         const double tmp_coords_jac_31__1 = jac_affine_0_2__2 * jac_affine_1_1__2;
         const double jac_affine_1_2__2    = p_affine_const_3_1__2 - p_affine_const_0_1__2;
         const double tmp_coords_jac_29__2 = jac_affine_0_1__2 * jac_affine_1_2__2;
         const double jac_affine_2_0__2    = p_affine_const_1_2__2 - p_affine_const_0_2__2;
         const double jac_affine_2_1__2    = p_affine_const_2_2__2 - p_affine_const_0_2__2;
         const double tmp_coords_jac_28__2 = jac_affine_1_2__2 * jac_affine_2_1__2;
         const double jac_affine_2_2__2    = p_affine_const_3_2__2 - p_affine_const_0_2__2;
         const double tmp_coords_jac_27__2 = jac_affine_1_1__2 * jac_affine_2_2__2;
         const double tmp_coords_jac_30__2 = jac_affine_0_1__2 * jac_affine_2_2__2;
         const double tmp_coords_jac_32__1 = jac_affine_0_0__2 * tmp_coords_jac_27__2 + jac_affine_2_0__2 * tmp_coords_jac_29__2 -
                                             jac_affine_0_0__2 * tmp_coords_jac_28__2 - jac_affine_1_0__2 * tmp_coords_jac_30__2 -
                                             jac_affine_2_0__2 * tmp_coords_jac_31__1 +
                                             jac_affine_0_2__2 * jac_affine_1_0__2 * jac_affine_2_1__2;
         const double tmp_coords_jac_33__1  = 1.0 / tmp_coords_jac_32__1;
         const double jac_affine_inv_0_0__2 = tmp_coords_jac_33__1 * ( tmp_coords_jac_27__2 - tmp_coords_jac_28__2 );
         const double jac_affine_inv_0_1__2 =
             tmp_coords_jac_33__1 * ( -1.0 * tmp_coords_jac_30__2 + jac_affine_0_2__2 * jac_affine_2_1__2 );
         const double jac_affine_inv_0_2__2 = tmp_coords_jac_33__1 * ( tmp_coords_jac_29__2 - tmp_coords_jac_31__1 );
         const double jac_affine_inv_1_0__2 =
             tmp_coords_jac_33__1 * ( jac_affine_1_2__2 * jac_affine_2_0__2 - jac_affine_1_0__2 * jac_affine_2_2__2 );
         const double jac_affine_inv_1_1__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_0__2 * jac_affine_2_2__2 - jac_affine_0_2__2 * jac_affine_2_0__2 );
         const double jac_affine_inv_1_2__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_2__2 * jac_affine_1_0__2 - jac_affine_0_0__2 * jac_affine_1_2__2 );
         const double jac_affine_inv_2_0__2 =
             tmp_coords_jac_33__1 * ( jac_affine_1_0__2 * jac_affine_2_1__2 - jac_affine_1_1__2 * jac_affine_2_0__2 );
         const double jac_affine_inv_2_1__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_1__2 * jac_affine_2_0__2 - jac_affine_0_0__2 * jac_affine_2_1__2 );
         const double jac_affine_inv_2_2__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_0__2 * jac_affine_1_1__2 - jac_affine_0_1__2 * jac_affine_1_0__2 );
         const double abs_det_jac_affine__2 = abs( tmp_coords_jac_32__1 );
         for ( int64_t ctr_2__2 = 0LL; ctr_2__2 < micro_edges_per_macro_edge; ctr_2__2 += 1LL )
         {
            for ( int64_t ctr_1__2 = 0LL; ctr_1__2 < -1LL * ctr_2__2 + micro_edges_per_macro_edge; ctr_1__2 += 1LL )
            {
               for ( int64_t ctr_0__2 = 0LL; ctr_0__2 < -1LL - ctr_1__2 - ctr_2__2 + micro_edges_per_macro_edge; ctr_0__2 += 1LL )
               {
                  const double p_affine_0_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_0_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_0_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_1_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_1_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_1_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_2_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_2_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_2_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_3_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__2 );
                  const double p_affine_3_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__2 );
                  const double p_affine_3_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__2 );
                  const double src_dof_0__2 =
                      _data_src[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_1__2 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL ) -
                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_2__2 =
                      _data_src[1LL - ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL -
                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_3__2 =
                      _data_src[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                ( -1LL * ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double tmp_kernel_op_0__2 = -1.0 * jac_affine_inv_0_0__2 - jac_affine_inv_1_0__2 - jac_affine_inv_2_0__2;
                  const double tmp_kernel_op_1__2 = -1.0 * jac_affine_inv_0_1__2 - jac_affine_inv_1_1__2 - jac_affine_inv_2_1__2;
                  const double tmp_kernel_op_2__2 = -1.0 * jac_affine_inv_0_2__2 - jac_affine_inv_1_2__2 - jac_affine_inv_2_2__2;
                  const double tmp_kernel_op_3__2 = 0.16666666666666663 * abs_det_jac_affine__2;
                  const double tmp_kernel_op_4__2 = src_dof_0__2 * tmp_kernel_op_3__2;
                  const double tmp_kernel_op_5__2 = jac_affine_inv_0_0__2 * tmp_kernel_op_0__2 +
                                                    jac_affine_inv_0_1__2 * tmp_kernel_op_1__2 +
                                                    jac_affine_inv_0_2__2 * tmp_kernel_op_2__2;
                  const double tmp_kernel_op_6__2 = src_dof_1__2 * tmp_kernel_op_3__2;
                  const double tmp_kernel_op_7__2 = jac_affine_inv_1_0__2 * tmp_kernel_op_0__2 +
                                                    jac_affine_inv_1_1__2 * tmp_kernel_op_1__2 +
                                                    jac_affine_inv_1_2__2 * tmp_kernel_op_2__2;
                  const double tmp_kernel_op_8__2 = src_dof_2__2 * tmp_kernel_op_3__2;
                  const double tmp_kernel_op_9__2 = jac_affine_inv_2_0__2 * tmp_kernel_op_0__2 +
                                                    jac_affine_inv_2_1__2 * tmp_kernel_op_1__2 +
                                                    jac_affine_inv_2_2__2 * tmp_kernel_op_2__2;
                  const double tmp_kernel_op_10__2 = src_dof_3__2 * tmp_kernel_op_3__2;
                  const double tmp_kernel_op_11__2 = jac_affine_inv_0_0__2 * jac_affine_inv_1_0__2 +
                                                     jac_affine_inv_0_1__2 * jac_affine_inv_1_1__2 +
                                                     jac_affine_inv_0_2__2 * jac_affine_inv_1_2__2;
                  const double tmp_kernel_op_12__2 = jac_affine_inv_0_0__2 * jac_affine_inv_2_0__2 +
                                                     jac_affine_inv_0_1__2 * jac_affine_inv_2_1__2 +
                                                     jac_affine_inv_0_2__2 * jac_affine_inv_2_2__2;
                  const double tmp_kernel_op_13__2 = jac_affine_inv_1_0__2 * jac_affine_inv_2_0__2 +
                                                     jac_affine_inv_1_1__2 * jac_affine_inv_2_1__2 +
                                                     jac_affine_inv_1_2__2 * jac_affine_inv_2_2__2;
                  const double elMatVec_0__2 =
                      tmp_kernel_op_10__2 * tmp_kernel_op_9__2 +
                      tmp_kernel_op_4__2 * ( tmp_kernel_op_0__2 * tmp_kernel_op_0__2 + tmp_kernel_op_1__2 * tmp_kernel_op_1__2 +
                                             tmp_kernel_op_2__2 * tmp_kernel_op_2__2 ) +
                      tmp_kernel_op_5__2 * tmp_kernel_op_6__2 + tmp_kernel_op_7__2 * tmp_kernel_op_8__2;
                  const double elMatVec_1__2 =
                      tmp_kernel_op_10__2 * tmp_kernel_op_12__2 + tmp_kernel_op_11__2 * tmp_kernel_op_8__2 +
                      tmp_kernel_op_4__2 * tmp_kernel_op_5__2 +
                      tmp_kernel_op_6__2 *
                          ( jac_affine_inv_0_0__2 * jac_affine_inv_0_0__2 + jac_affine_inv_0_1__2 * jac_affine_inv_0_1__2 +
                            jac_affine_inv_0_2__2 * jac_affine_inv_0_2__2 );
                  const double elMatVec_2__2 =
                      tmp_kernel_op_10__2 * tmp_kernel_op_13__2 + tmp_kernel_op_11__2 * tmp_kernel_op_6__2 +
                      tmp_kernel_op_4__2 * tmp_kernel_op_7__2 +
                      tmp_kernel_op_8__2 *
                          ( jac_affine_inv_1_0__2 * jac_affine_inv_1_0__2 + jac_affine_inv_1_1__2 * jac_affine_inv_1_1__2 +
                            jac_affine_inv_1_2__2 * jac_affine_inv_1_2__2 );
                  const double elMatVec_3__2 = tmp_kernel_op_10__2 * ( jac_affine_inv_2_0__2 * jac_affine_inv_2_0__2 +
                                                                       jac_affine_inv_2_1__2 * jac_affine_inv_2_1__2 +
                                                                       jac_affine_inv_2_2__2 * jac_affine_inv_2_2__2 ) +
                                               tmp_kernel_op_12__2 * tmp_kernel_op_6__2 +
                                               tmp_kernel_op_13__2 * tmp_kernel_op_8__2 + tmp_kernel_op_4__2 * tmp_kernel_op_9__2;
                  _data_dst[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                            ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_0__2 + _data_dst[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL ) -
                            ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_1__2 + _data_dst[-1LL * ( ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL ) -
                                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL -
                            ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_2__2 + _data_dst[1LL - ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL -
                                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                            ( -1LL * ctr_2__2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_3__2 + _data_dst[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                                ( -1LL * ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.BLUE_DOWN */
         const double tmp_coords_jac_0__1   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__1   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__1   = tmp_coords_jac_1__1 * 0.0;
         const double tmp_coords_jac_3__1   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4__1   = tmp_coords_jac_1__1 * 1.0;
         const double tmp_coords_jac_5__1   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_6__1   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__1 * tmp_coords_jac_5__1;
         const double tmp_coords_jac_7__1   = tmp_coords_jac_6__1 + tmp_coords_jac_3__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_8__1   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_9__1   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__1  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__1  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10__1 * tmp_coords_jac_2__1;
         const double tmp_coords_jac_12__1  = tmp_coords_jac_11__1 + tmp_coords_jac_4__1 * tmp_coords_jac_9__1;
         const double tmp_coords_jac_13__1  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_14__1  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_15__1  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__1  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__1 * tmp_coords_jac_2__1;
         const double tmp_coords_jac_17__1  = tmp_coords_jac_16__1 + tmp_coords_jac_14__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_18__1  = tmp_coords_jac_0__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_19__1  = tmp_coords_jac_18__1 + tmp_coords_jac_2__1 * tmp_coords_jac_3__1;
         const double tmp_coords_jac_20__1  = tmp_coords_jac_4__1 * tmp_coords_jac_8__1;
         const double tmp_coords_jac_21__1  = tmp_coords_jac_20__1 + tmp_coords_jac_2__1 * tmp_coords_jac_9__1;
         const double tmp_coords_jac_22__1  = tmp_coords_jac_13__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_23__1  = tmp_coords_jac_22__1 + tmp_coords_jac_14__1 * tmp_coords_jac_2__1;
         const double p_affine_const_0_0__1 = tmp_coords_jac_7__1 + tmp_coords_jac_0__1 * tmp_coords_jac_2__1;
         const double p_affine_const_0_1__1 = tmp_coords_jac_12__1 + tmp_coords_jac_2__1 * tmp_coords_jac_8__1;
         const double p_affine_const_0_2__1 = tmp_coords_jac_17__1 + tmp_coords_jac_13__1 * tmp_coords_jac_2__1;
         const double p_affine_const_1_0__1 = tmp_coords_jac_19__1 + tmp_coords_jac_6__1;
         const double p_affine_const_1_1__1 = tmp_coords_jac_11__1 + tmp_coords_jac_21__1;
         const double p_affine_const_1_2__1 = tmp_coords_jac_16__1 + tmp_coords_jac_23__1;
         const double p_affine_const_2_0__1 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_19__1 + tmp_coords_jac_4__1 * tmp_coords_jac_5__1;
         const double p_affine_const_2_1__1 =
             macro_vertex_coord_id_0comp1 + tmp_coords_jac_21__1 + tmp_coords_jac_10__1 * tmp_coords_jac_4__1;
         const double p_affine_const_2_2__1 =
             macro_vertex_coord_id_0comp2 + tmp_coords_jac_23__1 + tmp_coords_jac_15__1 * tmp_coords_jac_4__1;
         const double p_affine_const_3_0__1 = tmp_coords_jac_18__1 + tmp_coords_jac_7__1;
         const double p_affine_const_3_1__1 = tmp_coords_jac_12__1 + tmp_coords_jac_20__1;
         const double p_affine_const_3_2__1 = tmp_coords_jac_17__1 + tmp_coords_jac_22__1;
         const double jac_affine_0_0__1     = p_affine_const_1_0__1 - p_affine_const_0_0__1;
         const double jac_affine_0_1__1     = p_affine_const_2_0__1 - p_affine_const_0_0__1;
         const double jac_affine_0_2__1     = p_affine_const_3_0__1 - p_affine_const_0_0__1;
         const double jac_affine_1_0__1     = p_affine_const_1_1__1 - p_affine_const_0_1__1;
         const double jac_affine_1_1__1     = p_affine_const_2_1__1 - p_affine_const_0_1__1;
         const double tmp_coords_jac_28__1  = jac_affine_0_2__1 * jac_affine_1_1__1;
         const double jac_affine_1_2__1     = p_affine_const_3_1__1 - p_affine_const_0_1__1;
         const double tmp_coords_jac_26__1  = jac_affine_0_1__1 * jac_affine_1_2__1;
         const double jac_affine_2_0__1     = p_affine_const_1_2__1 - p_affine_const_0_2__1;
         const double jac_affine_2_1__1     = p_affine_const_2_2__1 - p_affine_const_0_2__1;
         const double tmp_coords_jac_25__1  = jac_affine_1_2__1 * jac_affine_2_1__1;
         const double jac_affine_2_2__1     = p_affine_const_3_2__1 - p_affine_const_0_2__1;
         const double tmp_coords_jac_24__1  = jac_affine_1_1__1 * jac_affine_2_2__1;
         const double tmp_coords_jac_27__1  = jac_affine_0_1__1 * jac_affine_2_2__1;
         const double tmp_coords_jac_29__1 = jac_affine_0_0__1 * tmp_coords_jac_24__1 + jac_affine_2_0__1 * tmp_coords_jac_26__1 -
                                             jac_affine_0_0__1 * tmp_coords_jac_25__1 - jac_affine_1_0__1 * tmp_coords_jac_27__1 -
                                             jac_affine_2_0__1 * tmp_coords_jac_28__1 +
                                             jac_affine_0_2__1 * jac_affine_1_0__1 * jac_affine_2_1__1;
         const double tmp_coords_jac_30__1  = 1.0 / tmp_coords_jac_29__1;
         const double jac_affine_inv_0_0__1 = tmp_coords_jac_30__1 * ( tmp_coords_jac_24__1 - tmp_coords_jac_25__1 );
         const double jac_affine_inv_0_1__1 =
             tmp_coords_jac_30__1 * ( -1.0 * tmp_coords_jac_27__1 + jac_affine_0_2__1 * jac_affine_2_1__1 );
         const double jac_affine_inv_0_2__1 = tmp_coords_jac_30__1 * ( tmp_coords_jac_26__1 - tmp_coords_jac_28__1 );
         const double jac_affine_inv_1_0__1 =
             tmp_coords_jac_30__1 * ( jac_affine_1_2__1 * jac_affine_2_0__1 - jac_affine_1_0__1 * jac_affine_2_2__1 );
         const double jac_affine_inv_1_1__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_0__1 * jac_affine_2_2__1 - jac_affine_0_2__1 * jac_affine_2_0__1 );
         const double jac_affine_inv_1_2__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_2__1 * jac_affine_1_0__1 - jac_affine_0_0__1 * jac_affine_1_2__1 );
         const double jac_affine_inv_2_0__1 =
             tmp_coords_jac_30__1 * ( jac_affine_1_0__1 * jac_affine_2_1__1 - jac_affine_1_1__1 * jac_affine_2_0__1 );
         const double jac_affine_inv_2_1__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_1__1 * jac_affine_2_0__1 - jac_affine_0_0__1 * jac_affine_2_1__1 );
         const double jac_affine_inv_2_2__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_0__1 * jac_affine_1_1__1 - jac_affine_0_1__1 * jac_affine_1_0__1 );
         const double abs_det_jac_affine__1 = abs( tmp_coords_jac_29__1 );
         for ( int64_t ctr_2__1 = 0LL; ctr_2__1 < micro_edges_per_macro_edge; ctr_2__1 += 1LL )
         {
            for ( int64_t ctr_1__1 = 0LL; ctr_1__1 < -1LL * ctr_2__1 + micro_edges_per_macro_edge; ctr_1__1 += 1LL )
            {
               for ( int64_t ctr_0__1 = 0LL; ctr_0__1 < -1LL - ctr_1__1 - ctr_2__1 + micro_edges_per_macro_edge; ctr_0__1 += 1LL )
               {
                  const double p_affine_0_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__1;
                  const double p_affine_0_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__1;
                  const double p_affine_0_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__1;
                  const double p_affine_1_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_1_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_1_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double src_dof_0__1 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__1 ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_1__1 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL ) -
                                ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_2__1 =
                      _data_src[1LL - ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL -
                                ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_3__1 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__1 ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double tmp_kernel_op_0__1 = -1.0 * jac_affine_inv_0_0__1 - jac_affine_inv_1_0__1 - jac_affine_inv_2_0__1;
                  const double tmp_kernel_op_1__1 = -1.0 * jac_affine_inv_0_1__1 - jac_affine_inv_1_1__1 - jac_affine_inv_2_1__1;
                  const double tmp_kernel_op_2__1 = -1.0 * jac_affine_inv_0_2__1 - jac_affine_inv_1_2__1 - jac_affine_inv_2_2__1;
                  const double tmp_kernel_op_3__1 = 0.16666666666666663 * abs_det_jac_affine__1;
                  const double tmp_kernel_op_4__1 = src_dof_0__1 * tmp_kernel_op_3__1;
                  const double tmp_kernel_op_5__1 = jac_affine_inv_0_0__1 * tmp_kernel_op_0__1 +
                                                    jac_affine_inv_0_1__1 * tmp_kernel_op_1__1 +
                                                    jac_affine_inv_0_2__1 * tmp_kernel_op_2__1;
                  const double tmp_kernel_op_6__1 = src_dof_1__1 * tmp_kernel_op_3__1;
                  const double tmp_kernel_op_7__1 = jac_affine_inv_1_0__1 * tmp_kernel_op_0__1 +
                                                    jac_affine_inv_1_1__1 * tmp_kernel_op_1__1 +
                                                    jac_affine_inv_1_2__1 * tmp_kernel_op_2__1;
                  const double tmp_kernel_op_8__1 = src_dof_2__1 * tmp_kernel_op_3__1;
                  const double tmp_kernel_op_9__1 = jac_affine_inv_2_0__1 * tmp_kernel_op_0__1 +
                                                    jac_affine_inv_2_1__1 * tmp_kernel_op_1__1 +
                                                    jac_affine_inv_2_2__1 * tmp_kernel_op_2__1;
                  const double tmp_kernel_op_10__1 = src_dof_3__1 * tmp_kernel_op_3__1;
                  const double tmp_kernel_op_11__1 = jac_affine_inv_0_0__1 * jac_affine_inv_1_0__1 +
                                                     jac_affine_inv_0_1__1 * jac_affine_inv_1_1__1 +
                                                     jac_affine_inv_0_2__1 * jac_affine_inv_1_2__1;
                  const double tmp_kernel_op_12__1 = jac_affine_inv_0_0__1 * jac_affine_inv_2_0__1 +
                                                     jac_affine_inv_0_1__1 * jac_affine_inv_2_1__1 +
                                                     jac_affine_inv_0_2__1 * jac_affine_inv_2_2__1;
                  const double tmp_kernel_op_13__1 = jac_affine_inv_1_0__1 * jac_affine_inv_2_0__1 +
                                                     jac_affine_inv_1_1__1 * jac_affine_inv_2_1__1 +
                                                     jac_affine_inv_1_2__1 * jac_affine_inv_2_2__1;
                  const double elMatVec_0__1 =
                      tmp_kernel_op_10__1 * tmp_kernel_op_9__1 +
                      tmp_kernel_op_4__1 * ( tmp_kernel_op_0__1 * tmp_kernel_op_0__1 + tmp_kernel_op_1__1 * tmp_kernel_op_1__1 +
                                             tmp_kernel_op_2__1 * tmp_kernel_op_2__1 ) +
                      tmp_kernel_op_5__1 * tmp_kernel_op_6__1 + tmp_kernel_op_7__1 * tmp_kernel_op_8__1;
                  const double elMatVec_1__1 =
                      tmp_kernel_op_10__1 * tmp_kernel_op_12__1 + tmp_kernel_op_11__1 * tmp_kernel_op_8__1 +
                      tmp_kernel_op_4__1 * tmp_kernel_op_5__1 +
                      tmp_kernel_op_6__1 *
                          ( jac_affine_inv_0_0__1 * jac_affine_inv_0_0__1 + jac_affine_inv_0_1__1 * jac_affine_inv_0_1__1 +
                            jac_affine_inv_0_2__1 * jac_affine_inv_0_2__1 );
                  const double elMatVec_2__1 =
                      tmp_kernel_op_10__1 * tmp_kernel_op_13__1 + tmp_kernel_op_11__1 * tmp_kernel_op_6__1 +
                      tmp_kernel_op_4__1 * tmp_kernel_op_7__1 +
                      tmp_kernel_op_8__1 *
                          ( jac_affine_inv_1_0__1 * jac_affine_inv_1_0__1 + jac_affine_inv_1_1__1 * jac_affine_inv_1_1__1 +
                            jac_affine_inv_1_2__1 * jac_affine_inv_1_2__1 );
                  const double elMatVec_3__1 = tmp_kernel_op_10__1 * ( jac_affine_inv_2_0__1 * jac_affine_inv_2_0__1 +
                                                                       jac_affine_inv_2_1__1 * jac_affine_inv_2_1__1 +
                                                                       jac_affine_inv_2_2__1 * jac_affine_inv_2_2__1 ) +
                                               tmp_kernel_op_12__1 * tmp_kernel_op_6__1 +
                                               tmp_kernel_op_13__1 * tmp_kernel_op_8__1 + tmp_kernel_op_4__1 * tmp_kernel_op_9__1;
                  _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                            ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__1 ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_0__1 + _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                                ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__1 ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL ) -
                            ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_1__1 + _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL ) -
                                                ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL -
                            ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_2__1 + _data_dst[1LL - ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL -
                                                ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                            ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__1 ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_3__1 + _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                                ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__1 ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.GREEN_UP */
         const double tmp_coords_jac_0__0   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__0   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__0   = tmp_coords_jac_1__0 * 0.0;
         const double tmp_coords_jac_3__0   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_0__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_4__0   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__0   = tmp_coords_jac_1__0 * 1.0;
         const double tmp_coords_jac_6__0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__0   = tmp_coords_jac_2__0 * tmp_coords_jac_6__0;
         const double tmp_coords_jac_8__0   = tmp_coords_jac_7__0 + tmp_coords_jac_4__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_9__0   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__0  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_2__0 * tmp_coords_jac_9__0;
         const double tmp_coords_jac_11__0  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12__0  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__0  = tmp_coords_jac_12__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_14__0  = tmp_coords_jac_13__0 + tmp_coords_jac_11__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_15__0  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__0  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_17__0  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18__0  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__0  = tmp_coords_jac_18__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_20__0  = tmp_coords_jac_19__0 + tmp_coords_jac_17__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_21__0  = tmp_coords_jac_2__0 * tmp_coords_jac_4__0;
         const double tmp_coords_jac_22__0  = tmp_coords_jac_11__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_23__0  = tmp_coords_jac_17__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_24__0  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_0__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_25__0  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_5__0 * tmp_coords_jac_9__0;
         const double tmp_coords_jac_26__0  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__0 * tmp_coords_jac_5__0;
         const double p_affine_const_0_0__0 = tmp_coords_jac_3__0 + tmp_coords_jac_8__0;
         const double p_affine_const_0_1__0 = tmp_coords_jac_10__0 + tmp_coords_jac_14__0;
         const double p_affine_const_0_2__0 = tmp_coords_jac_16__0 + tmp_coords_jac_20__0;
         const double p_affine_const_1_0__0 =
             tmp_coords_jac_21__0 + tmp_coords_jac_3__0 + tmp_coords_jac_5__0 * tmp_coords_jac_6__0;
         const double p_affine_const_1_1__0 =
             tmp_coords_jac_10__0 + tmp_coords_jac_22__0 + tmp_coords_jac_12__0 * tmp_coords_jac_5__0;
         const double p_affine_const_1_2__0 =
             tmp_coords_jac_16__0 + tmp_coords_jac_23__0 + tmp_coords_jac_18__0 * tmp_coords_jac_5__0;
         const double p_affine_const_2_0__0 = tmp_coords_jac_21__0 + tmp_coords_jac_24__0 + tmp_coords_jac_7__0;
         const double p_affine_const_2_1__0 = tmp_coords_jac_13__0 + tmp_coords_jac_22__0 + tmp_coords_jac_25__0;
         const double p_affine_const_2_2__0 = tmp_coords_jac_19__0 + tmp_coords_jac_23__0 + tmp_coords_jac_26__0;
         const double p_affine_const_3_0__0 = tmp_coords_jac_24__0 + tmp_coords_jac_8__0;
         const double p_affine_const_3_1__0 = tmp_coords_jac_14__0 + tmp_coords_jac_25__0;
         const double p_affine_const_3_2__0 = tmp_coords_jac_20__0 + tmp_coords_jac_26__0;
         const double jac_affine_0_0__0     = p_affine_const_1_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_1__0     = p_affine_const_2_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_2__0     = p_affine_const_3_0__0 - p_affine_const_0_0__0;
         const double jac_affine_1_0__0     = p_affine_const_1_1__0 - p_affine_const_0_1__0;
         const double jac_affine_1_1__0     = p_affine_const_2_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_31__0  = jac_affine_0_2__0 * jac_affine_1_1__0;
         const double jac_affine_1_2__0     = p_affine_const_3_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_29__0  = jac_affine_0_1__0 * jac_affine_1_2__0;
         const double jac_affine_2_0__0     = p_affine_const_1_2__0 - p_affine_const_0_2__0;
         const double jac_affine_2_1__0     = p_affine_const_2_2__0 - p_affine_const_0_2__0;
         const double tmp_coords_jac_28__0  = jac_affine_1_2__0 * jac_affine_2_1__0;
         const double jac_affine_2_2__0     = p_affine_const_3_2__0 - p_affine_const_0_2__0;
         const double tmp_coords_jac_27__0  = jac_affine_1_1__0 * jac_affine_2_2__0;
         const double tmp_coords_jac_30__0  = jac_affine_0_1__0 * jac_affine_2_2__0;
         const double tmp_coords_jac_32__0 = jac_affine_0_0__0 * tmp_coords_jac_27__0 + jac_affine_2_0__0 * tmp_coords_jac_29__0 -
                                             jac_affine_0_0__0 * tmp_coords_jac_28__0 - jac_affine_1_0__0 * tmp_coords_jac_30__0 -
                                             jac_affine_2_0__0 * tmp_coords_jac_31__0 +
                                             jac_affine_0_2__0 * jac_affine_1_0__0 * jac_affine_2_1__0;
         const double tmp_coords_jac_33__0  = 1.0 / tmp_coords_jac_32__0;
         const double jac_affine_inv_0_0__0 = tmp_coords_jac_33__0 * ( tmp_coords_jac_27__0 - tmp_coords_jac_28__0 );
         const double jac_affine_inv_0_1__0 =
             tmp_coords_jac_33__0 * ( -1.0 * tmp_coords_jac_30__0 + jac_affine_0_2__0 * jac_affine_2_1__0 );
         const double jac_affine_inv_0_2__0 = tmp_coords_jac_33__0 * ( tmp_coords_jac_29__0 - tmp_coords_jac_31__0 );
         const double jac_affine_inv_1_0__0 =
             tmp_coords_jac_33__0 * ( jac_affine_1_2__0 * jac_affine_2_0__0 - jac_affine_1_0__0 * jac_affine_2_2__0 );
         const double jac_affine_inv_1_1__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_0__0 * jac_affine_2_2__0 - jac_affine_0_2__0 * jac_affine_2_0__0 );
         const double jac_affine_inv_1_2__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_2__0 * jac_affine_1_0__0 - jac_affine_0_0__0 * jac_affine_1_2__0 );
         const double jac_affine_inv_2_0__0 =
             tmp_coords_jac_33__0 * ( jac_affine_1_0__0 * jac_affine_2_1__0 - jac_affine_1_1__0 * jac_affine_2_0__0 );
         const double jac_affine_inv_2_1__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_1__0 * jac_affine_2_0__0 - jac_affine_0_0__0 * jac_affine_2_1__0 );
         const double jac_affine_inv_2_2__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_0__0 * jac_affine_1_1__0 - jac_affine_0_1__0 * jac_affine_1_0__0 );
         const double abs_det_jac_affine__0 = abs( tmp_coords_jac_32__0 );
         for ( int64_t ctr_2__0 = 0LL; ctr_2__0 < micro_edges_per_macro_edge; ctr_2__0 += 1LL )
         {
            for ( int64_t ctr_1__0 = 0LL; ctr_1__0 < -1LL * ctr_2__0 + micro_edges_per_macro_edge; ctr_1__0 += 1LL )
            {
               for ( int64_t ctr_0__0 = 0LL; ctr_0__0 < -1LL - ctr_1__0 - ctr_2__0 + micro_edges_per_macro_edge; ctr_0__0 += 1LL )
               {
                  const double p_affine_0_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__0;
                  const double p_affine_0_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__0;
                  const double p_affine_0_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__0;
                  const double p_affine_1_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__0;
                  const double p_affine_1_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__0;
                  const double p_affine_1_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__0;
                  const double p_affine_2_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_2_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_2_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__0 );
                  const double src_dof_0__0 =
                      _data_src[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_1__0 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) -
                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1__0 ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) + ctr_0__0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_2__0 =
                      _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) -
                                ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_3__0 =
                      _data_src[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double tmp_kernel_op_0__0 = -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0 - jac_affine_inv_2_0__0;
                  const double tmp_kernel_op_1__0 = -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0 - jac_affine_inv_2_1__0;
                  const double tmp_kernel_op_2__0 = -1.0 * jac_affine_inv_0_2__0 - jac_affine_inv_1_2__0 - jac_affine_inv_2_2__0;
                  const double tmp_kernel_op_3__0 = 0.16666666666666663 * abs_det_jac_affine__0;
                  const double tmp_kernel_op_4__0 = src_dof_0__0 * tmp_kernel_op_3__0;
                  const double tmp_kernel_op_5__0 = jac_affine_inv_0_0__0 * tmp_kernel_op_0__0 +
                                                    jac_affine_inv_0_1__0 * tmp_kernel_op_1__0 +
                                                    jac_affine_inv_0_2__0 * tmp_kernel_op_2__0;
                  const double tmp_kernel_op_6__0 = src_dof_1__0 * tmp_kernel_op_3__0;
                  const double tmp_kernel_op_7__0 = jac_affine_inv_1_0__0 * tmp_kernel_op_0__0 +
                                                    jac_affine_inv_1_1__0 * tmp_kernel_op_1__0 +
                                                    jac_affine_inv_1_2__0 * tmp_kernel_op_2__0;
                  const double tmp_kernel_op_8__0 = src_dof_2__0 * tmp_kernel_op_3__0;
                  const double tmp_kernel_op_9__0 = jac_affine_inv_2_0__0 * tmp_kernel_op_0__0 +
                                                    jac_affine_inv_2_1__0 * tmp_kernel_op_1__0 +
                                                    jac_affine_inv_2_2__0 * tmp_kernel_op_2__0;
                  const double tmp_kernel_op_10__0 = src_dof_3__0 * tmp_kernel_op_3__0;
                  const double tmp_kernel_op_11__0 = jac_affine_inv_0_0__0 * jac_affine_inv_1_0__0 +
                                                     jac_affine_inv_0_1__0 * jac_affine_inv_1_1__0 +
                                                     jac_affine_inv_0_2__0 * jac_affine_inv_1_2__0;
                  const double tmp_kernel_op_12__0 = jac_affine_inv_0_0__0 * jac_affine_inv_2_0__0 +
                                                     jac_affine_inv_0_1__0 * jac_affine_inv_2_1__0 +
                                                     jac_affine_inv_0_2__0 * jac_affine_inv_2_2__0;
                  const double tmp_kernel_op_13__0 = jac_affine_inv_1_0__0 * jac_affine_inv_2_0__0 +
                                                     jac_affine_inv_1_1__0 * jac_affine_inv_2_1__0 +
                                                     jac_affine_inv_1_2__0 * jac_affine_inv_2_2__0;
                  const double elMatVec_0__0 =
                      tmp_kernel_op_10__0 * tmp_kernel_op_9__0 +
                      tmp_kernel_op_4__0 * ( tmp_kernel_op_0__0 * tmp_kernel_op_0__0 + tmp_kernel_op_1__0 * tmp_kernel_op_1__0 +
                                             tmp_kernel_op_2__0 * tmp_kernel_op_2__0 ) +
                      tmp_kernel_op_5__0 * tmp_kernel_op_6__0 + tmp_kernel_op_7__0 * tmp_kernel_op_8__0;
                  const double elMatVec_1__0 =
                      tmp_kernel_op_10__0 * tmp_kernel_op_12__0 + tmp_kernel_op_11__0 * tmp_kernel_op_8__0 +
                      tmp_kernel_op_4__0 * tmp_kernel_op_5__0 +
                      tmp_kernel_op_6__0 *
                          ( jac_affine_inv_0_0__0 * jac_affine_inv_0_0__0 + jac_affine_inv_0_1__0 * jac_affine_inv_0_1__0 +
                            jac_affine_inv_0_2__0 * jac_affine_inv_0_2__0 );
                  const double elMatVec_2__0 =
                      tmp_kernel_op_10__0 * tmp_kernel_op_13__0 + tmp_kernel_op_11__0 * tmp_kernel_op_6__0 +
                      tmp_kernel_op_4__0 * tmp_kernel_op_7__0 +
                      tmp_kernel_op_8__0 *
                          ( jac_affine_inv_1_0__0 * jac_affine_inv_1_0__0 + jac_affine_inv_1_1__0 * jac_affine_inv_1_1__0 +
                            jac_affine_inv_1_2__0 * jac_affine_inv_1_2__0 );
                  const double elMatVec_3__0 = tmp_kernel_op_10__0 * ( jac_affine_inv_2_0__0 * jac_affine_inv_2_0__0 +
                                                                       jac_affine_inv_2_1__0 * jac_affine_inv_2_1__0 +
                                                                       jac_affine_inv_2_2__0 * jac_affine_inv_2_2__0 ) +
                                               tmp_kernel_op_12__0 * tmp_kernel_op_6__0 +
                                               tmp_kernel_op_13__0 * tmp_kernel_op_8__0 + tmp_kernel_op_4__0 * tmp_kernel_op_9__0;
                  _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                            ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                            ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_0__0 + _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) -
                            ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1__0 ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) + ctr_0__0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_1__0 + _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) -
                                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL + ctr_1__0 ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) + ctr_0__0 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) -
                            ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_2__0 + _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) -
                                                ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                            ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_3__0 + _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                                ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                                    ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                                ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.GREEN_DOWN */
         const double tmp_coords_jac_0  = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1  = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2  = tmp_coords_jac_1 * 0.0;
         const double tmp_coords_jac_3  = tmp_coords_jac_0 * tmp_coords_jac_2;
         const double tmp_coords_jac_4  = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5  = tmp_coords_jac_1 * 1.0;
         const double tmp_coords_jac_6  = tmp_coords_jac_4 * tmp_coords_jac_5;
         const double tmp_coords_jac_7  = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_8  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_6 + tmp_coords_jac_2 * tmp_coords_jac_7;
         const double tmp_coords_jac_9  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10 = tmp_coords_jac_2 * tmp_coords_jac_9;
         const double tmp_coords_jac_11 = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12 = tmp_coords_jac_11 * tmp_coords_jac_5;
         const double tmp_coords_jac_13 = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_14 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12 + tmp_coords_jac_13 * tmp_coords_jac_2;
         const double tmp_coords_jac_15 = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16 = tmp_coords_jac_15 * tmp_coords_jac_2;
         const double tmp_coords_jac_17 = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18 = tmp_coords_jac_17 * tmp_coords_jac_5;
         const double tmp_coords_jac_19 = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_20 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18 + tmp_coords_jac_19 * tmp_coords_jac_2;
         const double tmp_coords_jac_21 = tmp_coords_jac_0 * tmp_coords_jac_5;
         const double tmp_coords_jac_22 = tmp_coords_jac_5 * tmp_coords_jac_9;
         const double tmp_coords_jac_23 = tmp_coords_jac_15 * tmp_coords_jac_5;
         const double tmp_coords_jac_24 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_5 * tmp_coords_jac_7;
         const double tmp_coords_jac_25 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_13 * tmp_coords_jac_5;
         const double tmp_coords_jac_26 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_19 * tmp_coords_jac_5;
         const double p_affine_const_0_0 = tmp_coords_jac_3 + tmp_coords_jac_8;
         const double p_affine_const_0_1 = tmp_coords_jac_10 + tmp_coords_jac_14;
         const double p_affine_const_0_2 = tmp_coords_jac_16 + tmp_coords_jac_20;
         const double p_affine_const_1_0 = tmp_coords_jac_21 + tmp_coords_jac_8;
         const double p_affine_const_1_1 = tmp_coords_jac_14 + tmp_coords_jac_22;
         const double p_affine_const_1_2 = tmp_coords_jac_20 + tmp_coords_jac_23;
         const double p_affine_const_2_0 = tmp_coords_jac_21 + tmp_coords_jac_24 + tmp_coords_jac_2 * tmp_coords_jac_4;
         const double p_affine_const_2_1 = tmp_coords_jac_22 + tmp_coords_jac_25 + tmp_coords_jac_11 * tmp_coords_jac_2;
         const double p_affine_const_2_2 = tmp_coords_jac_23 + tmp_coords_jac_26 + tmp_coords_jac_17 * tmp_coords_jac_2;
         const double p_affine_const_3_0 = tmp_coords_jac_24 + tmp_coords_jac_3 + tmp_coords_jac_6;
         const double p_affine_const_3_1 = tmp_coords_jac_10 + tmp_coords_jac_12 + tmp_coords_jac_25;
         const double p_affine_const_3_2 = tmp_coords_jac_16 + tmp_coords_jac_18 + tmp_coords_jac_26;
         const double jac_affine_0_0     = p_affine_const_1_0 - p_affine_const_0_0;
         const double jac_affine_0_1     = p_affine_const_2_0 - p_affine_const_0_0;
         const double jac_affine_0_2     = p_affine_const_3_0 - p_affine_const_0_0;
         const double jac_affine_1_0     = p_affine_const_1_1 - p_affine_const_0_1;
         const double jac_affine_1_1     = p_affine_const_2_1 - p_affine_const_0_1;
         const double tmp_coords_jac_31  = jac_affine_0_2 * jac_affine_1_1;
         const double jac_affine_1_2     = p_affine_const_3_1 - p_affine_const_0_1;
         const double tmp_coords_jac_29  = jac_affine_0_1 * jac_affine_1_2;
         const double jac_affine_2_0     = p_affine_const_1_2 - p_affine_const_0_2;
         const double jac_affine_2_1     = p_affine_const_2_2 - p_affine_const_0_2;
         const double tmp_coords_jac_28  = jac_affine_1_2 * jac_affine_2_1;
         const double jac_affine_2_2     = p_affine_const_3_2 - p_affine_const_0_2;
         const double tmp_coords_jac_27  = jac_affine_1_1 * jac_affine_2_2;
         const double tmp_coords_jac_30  = jac_affine_0_1 * jac_affine_2_2;
         const double tmp_coords_jac_32  = jac_affine_0_0 * tmp_coords_jac_27 + jac_affine_2_0 * tmp_coords_jac_29 -
                                          jac_affine_0_0 * tmp_coords_jac_28 - jac_affine_1_0 * tmp_coords_jac_30 -
                                          jac_affine_2_0 * tmp_coords_jac_31 + jac_affine_0_2 * jac_affine_1_0 * jac_affine_2_1;
         const double tmp_coords_jac_33  = 1.0 / tmp_coords_jac_32;
         const double jac_affine_inv_0_0 = tmp_coords_jac_33 * ( tmp_coords_jac_27 - tmp_coords_jac_28 );
         const double jac_affine_inv_0_1 = tmp_coords_jac_33 * ( -1.0 * tmp_coords_jac_30 + jac_affine_0_2 * jac_affine_2_1 );
         const double jac_affine_inv_0_2 = tmp_coords_jac_33 * ( tmp_coords_jac_29 - tmp_coords_jac_31 );
         const double jac_affine_inv_1_0 =
             tmp_coords_jac_33 * ( jac_affine_1_2 * jac_affine_2_0 - jac_affine_1_0 * jac_affine_2_2 );
         const double jac_affine_inv_1_1 =
             tmp_coords_jac_33 * ( jac_affine_0_0 * jac_affine_2_2 - jac_affine_0_2 * jac_affine_2_0 );
         const double jac_affine_inv_1_2 =
             tmp_coords_jac_33 * ( jac_affine_0_2 * jac_affine_1_0 - jac_affine_0_0 * jac_affine_1_2 );
         const double jac_affine_inv_2_0 =
             tmp_coords_jac_33 * ( jac_affine_1_0 * jac_affine_2_1 - jac_affine_1_1 * jac_affine_2_0 );
         const double jac_affine_inv_2_1 =
             tmp_coords_jac_33 * ( jac_affine_0_1 * jac_affine_2_0 - jac_affine_0_0 * jac_affine_2_1 );
         const double jac_affine_inv_2_2 =
             tmp_coords_jac_33 * ( jac_affine_0_0 * jac_affine_1_1 - jac_affine_0_1 * jac_affine_1_0 );
         const double abs_det_jac_affine = abs( tmp_coords_jac_32 );
         for ( int64_t ctr_2 = 0LL; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1LL )
         {
            for ( int64_t ctr_1 = 0LL; ctr_1 < -1LL * ctr_2 + micro_edges_per_macro_edge; ctr_1 += 1LL )
            {
               for ( int64_t ctr_0 = 0LL; ctr_0 < -1LL - ctr_1 - ctr_2 + micro_edges_per_macro_edge; ctr_0 += 1LL )
               {
                  const double p_affine_0_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2;
                  const double p_affine_0_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2;
                  const double p_affine_0_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2;
                  const double p_affine_1_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2;
                  const double p_affine_1_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2;
                  const double p_affine_1_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2;
                  const double p_affine_2_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_2_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_2_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2 );
                  const double src_dof_0 =
                      _data_src[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_1 =
                      _data_src[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_2 =
                      _data_src[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ctr_1 + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double src_dof_3 =
                      _data_src[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  const double tmp_kernel_op_0 = -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0 - jac_affine_inv_2_0;
                  const double tmp_kernel_op_1 = -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1 - jac_affine_inv_2_1;
                  const double tmp_kernel_op_2 = -1.0 * jac_affine_inv_0_2 - jac_affine_inv_1_2 - jac_affine_inv_2_2;
                  const double tmp_kernel_op_3 = 0.16666666666666663 * abs_det_jac_affine;
                  const double tmp_kernel_op_4 = src_dof_0 * tmp_kernel_op_3;
                  const double tmp_kernel_op_5 = jac_affine_inv_0_0 * tmp_kernel_op_0 + jac_affine_inv_0_1 * tmp_kernel_op_1 +
                                                 jac_affine_inv_0_2 * tmp_kernel_op_2;
                  const double tmp_kernel_op_6 = src_dof_1 * tmp_kernel_op_3;
                  const double tmp_kernel_op_7 = jac_affine_inv_1_0 * tmp_kernel_op_0 + jac_affine_inv_1_1 * tmp_kernel_op_1 +
                                                 jac_affine_inv_1_2 * tmp_kernel_op_2;
                  const double tmp_kernel_op_8 = src_dof_2 * tmp_kernel_op_3;
                  const double tmp_kernel_op_9 = jac_affine_inv_2_0 * tmp_kernel_op_0 + jac_affine_inv_2_1 * tmp_kernel_op_1 +
                                                 jac_affine_inv_2_2 * tmp_kernel_op_2;
                  const double tmp_kernel_op_10 = src_dof_3 * tmp_kernel_op_3;
                  const double tmp_kernel_op_11 = jac_affine_inv_0_0 * jac_affine_inv_1_0 +
                                                  jac_affine_inv_0_1 * jac_affine_inv_1_1 +
                                                  jac_affine_inv_0_2 * jac_affine_inv_1_2;
                  const double tmp_kernel_op_12 = jac_affine_inv_0_0 * jac_affine_inv_2_0 +
                                                  jac_affine_inv_0_1 * jac_affine_inv_2_1 +
                                                  jac_affine_inv_0_2 * jac_affine_inv_2_2;
                  const double tmp_kernel_op_13 = jac_affine_inv_1_0 * jac_affine_inv_2_0 +
                                                  jac_affine_inv_1_1 * jac_affine_inv_2_1 +
                                                  jac_affine_inv_1_2 * jac_affine_inv_2_2;
                  const double elMatVec_0 =
                      tmp_kernel_op_10 * tmp_kernel_op_9 +
                      tmp_kernel_op_4 * ( tmp_kernel_op_0 * tmp_kernel_op_0 + tmp_kernel_op_1 * tmp_kernel_op_1 +
                                          tmp_kernel_op_2 * tmp_kernel_op_2 ) +
                      tmp_kernel_op_5 * tmp_kernel_op_6 + tmp_kernel_op_7 * tmp_kernel_op_8;
                  const double elMatVec_1 =
                      tmp_kernel_op_10 * tmp_kernel_op_12 + tmp_kernel_op_11 * tmp_kernel_op_8 +
                      tmp_kernel_op_4 * tmp_kernel_op_5 +
                      tmp_kernel_op_6 * ( jac_affine_inv_0_0 * jac_affine_inv_0_0 + jac_affine_inv_0_1 * jac_affine_inv_0_1 +
                                          jac_affine_inv_0_2 * jac_affine_inv_0_2 );
                  const double elMatVec_2 =
                      tmp_kernel_op_10 * tmp_kernel_op_13 + tmp_kernel_op_11 * tmp_kernel_op_6 +
                      tmp_kernel_op_4 * tmp_kernel_op_7 +
                      tmp_kernel_op_8 * ( jac_affine_inv_1_0 * jac_affine_inv_1_0 + jac_affine_inv_1_1 * jac_affine_inv_1_1 +
                                          jac_affine_inv_1_2 * jac_affine_inv_1_2 );
                  const double elMatVec_3 =
                      tmp_kernel_op_10 * ( jac_affine_inv_2_0 * jac_affine_inv_2_0 + jac_affine_inv_2_1 * jac_affine_inv_2_1 +
                                           jac_affine_inv_2_2 * jac_affine_inv_2_2 ) +
                      tmp_kernel_op_12 * tmp_kernel_op_6 + tmp_kernel_op_13 * tmp_kernel_op_8 + tmp_kernel_op_4 * tmp_kernel_op_9;
                  _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                            ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_0 +
                      _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL -
                            ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_1 +
                      _data_dst[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL -
                            ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ctr_1 + ctr_0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_2 +
                      _data_dst[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ctr_1 + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                            ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                            ( 1LL + ctr_1 ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                            ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatVec_3 +
                      _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
   }
}
void P1ElementwiseDiffusion::apply_P1ElementwiseDiffusion_macro_2D( double* RESTRICT const _data_dst,
                                                                    double* RESTRICT const _data_src,
                                                                    const double           macro_vertex_coord_id_0comp0,
                                                                    const double           macro_vertex_coord_id_0comp1,
                                                                    const double           macro_vertex_coord_id_1comp0,
                                                                    const double           macro_vertex_coord_id_1comp1,
                                                                    const double           macro_vertex_coord_id_2comp0,
                                                                    const double           macro_vertex_coord_id_2comp1,
                                                                    const int64_t          micro_edges_per_macro_edge,
                                                                    const double micro_edges_per_macro_edge_float ) const
{
   {
      {
         /* FaceType.GRAY */
         const double tmp_coords_jac_0__0   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__0   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__0   = tmp_coords_jac_1__0 * 0.0;
         const double tmp_coords_jac_3__0   = tmp_coords_jac_0__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_4__0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__0   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__0 * tmp_coords_jac_4__0;
         const double tmp_coords_jac_6__0   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_7__0   = tmp_coords_jac_2__0 * tmp_coords_jac_6__0;
         const double tmp_coords_jac_8__0   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_9__0   = macro_vertex_coord_id_0comp1 + tmp_coords_jac_2__0 * tmp_coords_jac_8__0;
         const double tmp_coords_jac_10__0  = tmp_coords_jac_1__0 * 1.0;
         const double p_affine_const_0_0__0 = tmp_coords_jac_3__0 + tmp_coords_jac_5__0;
         const double p_affine_const_0_1__0 = tmp_coords_jac_7__0 + tmp_coords_jac_9__0;
         const double p_affine_const_1_0__0 = tmp_coords_jac_5__0 + tmp_coords_jac_0__0 * tmp_coords_jac_10__0;
         const double p_affine_const_1_1__0 = tmp_coords_jac_9__0 + tmp_coords_jac_10__0 * tmp_coords_jac_6__0;
         const double p_affine_const_2_0__0 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__0 + tmp_coords_jac_10__0 * tmp_coords_jac_4__0;
         const double p_affine_const_2_1__0 =
             macro_vertex_coord_id_0comp1 + tmp_coords_jac_7__0 + tmp_coords_jac_10__0 * tmp_coords_jac_8__0;
         const double jac_affine_0_0__0     = p_affine_const_1_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_1__0     = p_affine_const_2_0__0 - p_affine_const_0_0__0;
         const double jac_affine_1_0__0     = p_affine_const_1_1__0 - p_affine_const_0_1__0;
         const double jac_affine_1_1__0     = p_affine_const_2_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_11__0  = jac_affine_0_0__0 * jac_affine_1_1__0 - jac_affine_0_1__0 * jac_affine_1_0__0;
         const double tmp_coords_jac_12__0  = 1.0 / tmp_coords_jac_11__0;
         const double jac_affine_inv_0_0__0 = jac_affine_1_1__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_0_1__0 = -1.0 * jac_affine_0_1__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_1_0__0 = -1.0 * jac_affine_1_0__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_1_1__0 = jac_affine_0_0__0 * tmp_coords_jac_12__0;
         const double abs_det_jac_affine__0 = abs( tmp_coords_jac_11__0 );
         for ( int64_t ctr_1__0 = 0LL; ctr_1__0 < micro_edges_per_macro_edge; ctr_1__0 += 1LL )
         {
            for ( int64_t ctr_0__0 = 0LL; ctr_0__0 < -1LL * ctr_1__0 + micro_edges_per_macro_edge; ctr_0__0 += 1LL )
            {
               const double p_affine_0_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0;
               const double p_affine_0_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0;
               const double p_affine_1_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0;
               const double p_affine_1_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0;
               const double p_affine_2_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__0 );
               const double p_affine_2_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__0 );
               const double src_dof_0__0       = _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) +
                                                     ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               const double src_dof_1__0       = _data_src[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL +
                                                     ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               const double src_dof_2__0       = _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) +
                                                     ( 1LL + ctr_1__0 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0__0];
               const double tmp_kernel_op_0__0 = -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0;
               const double tmp_kernel_op_1__0 = -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0;
               const double tmp_kernel_op_2__0 = 0.5 * abs_det_jac_affine__0;
               const double tmp_kernel_op_3__0 = src_dof_0__0 * tmp_kernel_op_2__0;
               const double tmp_kernel_op_4__0 =
                   jac_affine_inv_0_0__0 * tmp_kernel_op_0__0 + jac_affine_inv_0_1__0 * tmp_kernel_op_1__0;
               const double tmp_kernel_op_5__0 = src_dof_1__0 * tmp_kernel_op_2__0;
               const double tmp_kernel_op_6__0 =
                   jac_affine_inv_1_0__0 * tmp_kernel_op_0__0 + jac_affine_inv_1_1__0 * tmp_kernel_op_1__0;
               const double tmp_kernel_op_7__0 = src_dof_2__0 * tmp_kernel_op_2__0;
               const double tmp_kernel_op_8__0 =
                   jac_affine_inv_0_0__0 * jac_affine_inv_1_0__0 + jac_affine_inv_0_1__0 * jac_affine_inv_1_1__0;
               const double elMatVec_0__0 =
                   tmp_kernel_op_3__0 * ( tmp_kernel_op_0__0 * tmp_kernel_op_0__0 + tmp_kernel_op_1__0 * tmp_kernel_op_1__0 ) +
                   tmp_kernel_op_4__0 * tmp_kernel_op_5__0 + tmp_kernel_op_6__0 * tmp_kernel_op_7__0;
               const double elMatVec_1__0 = tmp_kernel_op_3__0 * tmp_kernel_op_4__0 +
                                            tmp_kernel_op_5__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_0_0__0 +
                                                                   jac_affine_inv_0_1__0 * jac_affine_inv_0_1__0 ) +
                                            tmp_kernel_op_7__0 * tmp_kernel_op_8__0;
               const double elMatVec_2__0 = tmp_kernel_op_3__0 * tmp_kernel_op_6__0 + tmp_kernel_op_5__0 * tmp_kernel_op_8__0 +
                                            tmp_kernel_op_7__0 * ( jac_affine_inv_1_0__0 * jac_affine_inv_1_0__0 +
                                                                   jac_affine_inv_1_1__0 * jac_affine_inv_1_1__0 );
               _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) + ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 +
                         ctr_0__0] = elMatVec_0__0 + _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) +
                                                               ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0] =
                   elMatVec_1__0 + _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL +
                                             ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) +
                         ( 1LL + ctr_1__0 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0__0] =
                   elMatVec_2__0 + _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) +
                                             ( 1LL + ctr_1__0 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0__0];
            }
         }
      }
      {
         /* FaceType.BLUE */
         const double tmp_coords_jac_0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2   = tmp_coords_jac_1 * 0.0;
         const double tmp_coords_jac_3   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4   = tmp_coords_jac_1 * 1.0;
         const double tmp_coords_jac_5   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3 * tmp_coords_jac_4;
         const double tmp_coords_jac_6   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_7   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_8   = macro_vertex_coord_id_0comp1 + tmp_coords_jac_4 * tmp_coords_jac_7;
         const double tmp_coords_jac_9   = tmp_coords_jac_0 * tmp_coords_jac_4;
         const double tmp_coords_jac_10  = tmp_coords_jac_4 * tmp_coords_jac_6;
         const double p_affine_const_0_0 = tmp_coords_jac_5 + tmp_coords_jac_0 * tmp_coords_jac_2;
         const double p_affine_const_0_1 = tmp_coords_jac_8 + tmp_coords_jac_2 * tmp_coords_jac_6;
         const double p_affine_const_1_0 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_9 + tmp_coords_jac_2 * tmp_coords_jac_3;
         const double p_affine_const_1_1 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10 + tmp_coords_jac_2 * tmp_coords_jac_7;
         const double p_affine_const_2_0 = tmp_coords_jac_5 + tmp_coords_jac_9;
         const double p_affine_const_2_1 = tmp_coords_jac_10 + tmp_coords_jac_8;
         const double jac_affine_0_0     = p_affine_const_1_0 - p_affine_const_0_0;
         const double jac_affine_0_1     = p_affine_const_2_0 - p_affine_const_0_0;
         const double jac_affine_1_0     = p_affine_const_1_1 - p_affine_const_0_1;
         const double jac_affine_1_1     = p_affine_const_2_1 - p_affine_const_0_1;
         const double tmp_coords_jac_11  = jac_affine_0_0 * jac_affine_1_1 - jac_affine_0_1 * jac_affine_1_0;
         const double tmp_coords_jac_12  = 1.0 / tmp_coords_jac_11;
         const double jac_affine_inv_0_0 = jac_affine_1_1 * tmp_coords_jac_12;
         const double jac_affine_inv_0_1 = -1.0 * jac_affine_0_1 * tmp_coords_jac_12;
         const double jac_affine_inv_1_0 = -1.0 * jac_affine_1_0 * tmp_coords_jac_12;
         const double jac_affine_inv_1_1 = jac_affine_0_0 * tmp_coords_jac_12;
         const double abs_det_jac_affine = abs( tmp_coords_jac_11 );
         for ( int64_t ctr_1 = 0LL; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1LL )
         {
            for ( int64_t ctr_0 = 0LL; ctr_0 < -1LL - ctr_1 + micro_edges_per_macro_edge; ctr_0 += 1LL )
            {
               const double p_affine_0_0 = macro_vertex_coord_id_0comp0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) *
                                               (double) ( 1LL + ctr_0 ) +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1;
               const double p_affine_0_1 = macro_vertex_coord_id_0comp1 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) *
                                               (double) ( 1LL + ctr_0 ) +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1;
               const double p_affine_1_0 = macro_vertex_coord_id_0comp0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) *
                                               (double) ( 1LL + ctr_1 );
               const double p_affine_1_1 = macro_vertex_coord_id_0comp1 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) *
                                               (double) ( 1LL + ctr_1 );
               const double p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 );
               const double p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 );
               const double src_dof_0 =
                   _data_src[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1 + ctr_0];
               const double src_dof_1       = _data_src[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) +
                                                  ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               const double src_dof_2       = _data_src[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL +
                                                  ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               const double tmp_kernel_op_0 = -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0;
               const double tmp_kernel_op_1 = -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1;
               const double tmp_kernel_op_2 = 0.5 * abs_det_jac_affine;
               const double tmp_kernel_op_3 = src_dof_0 * tmp_kernel_op_2;
               const double tmp_kernel_op_4 = jac_affine_inv_0_0 * tmp_kernel_op_0 + jac_affine_inv_0_1 * tmp_kernel_op_1;
               const double tmp_kernel_op_5 = src_dof_1 * tmp_kernel_op_2;
               const double tmp_kernel_op_6 = jac_affine_inv_1_0 * tmp_kernel_op_0 + jac_affine_inv_1_1 * tmp_kernel_op_1;
               const double tmp_kernel_op_7 = src_dof_2 * tmp_kernel_op_2;
               const double tmp_kernel_op_8 = jac_affine_inv_0_0 * jac_affine_inv_1_0 + jac_affine_inv_0_1 * jac_affine_inv_1_1;
               const double elMatVec_0 =
                   tmp_kernel_op_3 * ( tmp_kernel_op_0 * tmp_kernel_op_0 + tmp_kernel_op_1 * tmp_kernel_op_1 ) +
                   tmp_kernel_op_4 * tmp_kernel_op_5 + tmp_kernel_op_6 * tmp_kernel_op_7;
               const double elMatVec_1 =
                   tmp_kernel_op_3 * tmp_kernel_op_4 +
                   tmp_kernel_op_5 * ( jac_affine_inv_0_0 * jac_affine_inv_0_0 + jac_affine_inv_0_1 * jac_affine_inv_0_1 ) +
                   tmp_kernel_op_7 * tmp_kernel_op_8;
               const double elMatVec_2 =
                   tmp_kernel_op_3 * tmp_kernel_op_6 + tmp_kernel_op_5 * tmp_kernel_op_8 +
                   tmp_kernel_op_7 * ( jac_affine_inv_1_0 * jac_affine_inv_1_0 + jac_affine_inv_1_1 * jac_affine_inv_1_1 );
               _data_dst[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1 + ctr_0] =
                   elMatVec_0 +
                   _data_dst[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1 + ctr_0];
               _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) +
                         ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0] =
                   elMatVec_1 + _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) +
                                          ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               _data_dst[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL + ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) +
                         ctr_0] = elMatVec_2 + _data_dst[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL +
                                                         ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
            }
         }
      }
   }
}
void P1ElementwiseDiffusion::toMatrix_P1ElementwiseDiffusion_macro_3D( int64_t* RESTRICT const _data_dst,
                                                                       int64_t* RESTRICT const _data_src,
                                                                       const double            macro_vertex_coord_id_0comp0,
                                                                       const double            macro_vertex_coord_id_0comp1,
                                                                       const double            macro_vertex_coord_id_0comp2,
                                                                       const double            macro_vertex_coord_id_1comp0,
                                                                       const double            macro_vertex_coord_id_1comp1,
                                                                       const double            macro_vertex_coord_id_1comp2,
                                                                       const double            macro_vertex_coord_id_2comp0,
                                                                       const double            macro_vertex_coord_id_2comp1,
                                                                       const double            macro_vertex_coord_id_2comp2,
                                                                       const double            macro_vertex_coord_id_3comp0,
                                                                       const double            macro_vertex_coord_id_3comp1,
                                                                       const double            macro_vertex_coord_id_3comp2,
                                                                       const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                       const int64_t micro_edges_per_macro_edge,
                                                                       const double  micro_edges_per_macro_edge_float ) const
{
   {
      {
         /* CellType.WHITE_UP */
         const double tmp_coords_jac_0__4   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__4   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__4   = tmp_coords_jac_1__4 * 0.0;
         const double tmp_coords_jac_3__4   = tmp_coords_jac_0__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_4__4   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__4   = tmp_coords_jac_2__4 * tmp_coords_jac_4__4;
         const double tmp_coords_jac_6__4   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__4   = tmp_coords_jac_2__4 * tmp_coords_jac_6__4;
         const double tmp_coords_jac_8__4   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_5__4 + tmp_coords_jac_7__4;
         const double tmp_coords_jac_9__4   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__4  = tmp_coords_jac_2__4 * tmp_coords_jac_9__4;
         const double tmp_coords_jac_11__4  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12__4  = tmp_coords_jac_11__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_13__4  = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_14__4  = tmp_coords_jac_13__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_15__4  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12__4 + tmp_coords_jac_14__4;
         const double tmp_coords_jac_16__4  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_17__4  = tmp_coords_jac_16__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_18__4  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__4  = tmp_coords_jac_18__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_20__4  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_21__4  = tmp_coords_jac_2__4 * tmp_coords_jac_20__4;
         const double tmp_coords_jac_22__4  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_19__4 + tmp_coords_jac_21__4;
         const double tmp_coords_jac_23__4  = tmp_coords_jac_1__4 * 1.0;
         const double tmp_coords_jac_24__4  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__4;
         const double tmp_coords_jac_25__4  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10__4;
         const double tmp_coords_jac_26__4  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_17__4;
         const double p_affine_const_0_0__4 = tmp_coords_jac_3__4 + tmp_coords_jac_8__4;
         const double p_affine_const_0_1__4 = tmp_coords_jac_10__4 + tmp_coords_jac_15__4;
         const double p_affine_const_0_2__4 = tmp_coords_jac_17__4 + tmp_coords_jac_22__4;
         const double p_affine_const_1_0__4 = tmp_coords_jac_8__4 + tmp_coords_jac_0__4 * tmp_coords_jac_23__4;
         const double p_affine_const_1_1__4 = tmp_coords_jac_15__4 + tmp_coords_jac_23__4 * tmp_coords_jac_9__4;
         const double p_affine_const_1_2__4 = tmp_coords_jac_22__4 + tmp_coords_jac_16__4 * tmp_coords_jac_23__4;
         const double p_affine_const_2_0__4 =
             tmp_coords_jac_24__4 + tmp_coords_jac_7__4 + tmp_coords_jac_23__4 * tmp_coords_jac_4__4;
         const double p_affine_const_2_1__4 =
             tmp_coords_jac_14__4 + tmp_coords_jac_25__4 + tmp_coords_jac_11__4 * tmp_coords_jac_23__4;
         const double p_affine_const_2_2__4 =
             tmp_coords_jac_21__4 + tmp_coords_jac_26__4 + tmp_coords_jac_18__4 * tmp_coords_jac_23__4;
         const double p_affine_const_3_0__4 =
             tmp_coords_jac_24__4 + tmp_coords_jac_5__4 + tmp_coords_jac_23__4 * tmp_coords_jac_6__4;
         const double p_affine_const_3_1__4 =
             tmp_coords_jac_12__4 + tmp_coords_jac_25__4 + tmp_coords_jac_13__4 * tmp_coords_jac_23__4;
         const double p_affine_const_3_2__4 =
             tmp_coords_jac_19__4 + tmp_coords_jac_26__4 + tmp_coords_jac_20__4 * tmp_coords_jac_23__4;
         const double jac_affine_0_0__4    = p_affine_const_1_0__4 - p_affine_const_0_0__4;
         const double jac_affine_0_1__4    = p_affine_const_2_0__4 - p_affine_const_0_0__4;
         const double jac_affine_0_2__4    = p_affine_const_3_0__4 - p_affine_const_0_0__4;
         const double jac_affine_1_0__4    = p_affine_const_1_1__4 - p_affine_const_0_1__4;
         const double jac_affine_1_1__4    = p_affine_const_2_1__4 - p_affine_const_0_1__4;
         const double tmp_coords_jac_31__2 = jac_affine_0_2__4 * jac_affine_1_1__4;
         const double jac_affine_1_2__4    = p_affine_const_3_1__4 - p_affine_const_0_1__4;
         const double tmp_coords_jac_29__4 = jac_affine_0_1__4 * jac_affine_1_2__4;
         const double jac_affine_2_0__4    = p_affine_const_1_2__4 - p_affine_const_0_2__4;
         const double jac_affine_2_1__4    = p_affine_const_2_2__4 - p_affine_const_0_2__4;
         const double tmp_coords_jac_28__4 = jac_affine_1_2__4 * jac_affine_2_1__4;
         const double jac_affine_2_2__4    = p_affine_const_3_2__4 - p_affine_const_0_2__4;
         const double tmp_coords_jac_27__4 = jac_affine_1_1__4 * jac_affine_2_2__4;
         const double tmp_coords_jac_30__4 = jac_affine_0_1__4 * jac_affine_2_2__4;
         const double tmp_coords_jac_32__2 = jac_affine_0_0__4 * tmp_coords_jac_27__4 + jac_affine_2_0__4 * tmp_coords_jac_29__4 -
                                             jac_affine_0_0__4 * tmp_coords_jac_28__4 - jac_affine_1_0__4 * tmp_coords_jac_30__4 -
                                             jac_affine_2_0__4 * tmp_coords_jac_31__2 +
                                             jac_affine_0_2__4 * jac_affine_1_0__4 * jac_affine_2_1__4;
         const double tmp_coords_jac_33__2  = 1.0 / tmp_coords_jac_32__2;
         const double jac_affine_inv_0_0__4 = tmp_coords_jac_33__2 * ( tmp_coords_jac_27__4 - tmp_coords_jac_28__4 );
         const double jac_affine_inv_0_1__4 =
             tmp_coords_jac_33__2 * ( -1.0 * tmp_coords_jac_30__4 + jac_affine_0_2__4 * jac_affine_2_1__4 );
         const double jac_affine_inv_0_2__4 = tmp_coords_jac_33__2 * ( tmp_coords_jac_29__4 - tmp_coords_jac_31__2 );
         const double jac_affine_inv_1_0__4 =
             tmp_coords_jac_33__2 * ( jac_affine_1_2__4 * jac_affine_2_0__4 - jac_affine_1_0__4 * jac_affine_2_2__4 );
         const double jac_affine_inv_1_1__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_0__4 * jac_affine_2_2__4 - jac_affine_0_2__4 * jac_affine_2_0__4 );
         const double jac_affine_inv_1_2__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_2__4 * jac_affine_1_0__4 - jac_affine_0_0__4 * jac_affine_1_2__4 );
         const double jac_affine_inv_2_0__4 =
             tmp_coords_jac_33__2 * ( jac_affine_1_0__4 * jac_affine_2_1__4 - jac_affine_1_1__4 * jac_affine_2_0__4 );
         const double jac_affine_inv_2_1__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_1__4 * jac_affine_2_0__4 - jac_affine_0_0__4 * jac_affine_2_1__4 );
         const double jac_affine_inv_2_2__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_0__4 * jac_affine_1_1__4 - jac_affine_0_1__4 * jac_affine_1_0__4 );
         const double abs_det_jac_affine__4 = abs( tmp_coords_jac_32__2 );
         for ( int64_t ctr_2__4 = 0LL; ctr_2__4 < micro_edges_per_macro_edge; ctr_2__4 += 1LL )
         {
            for ( int64_t ctr_1__4 = 0LL; ctr_1__4 < -1LL * ctr_2__4 + micro_edges_per_macro_edge; ctr_1__4 += 1LL )
            {
               for ( int64_t ctr_0__4 = 0LL; ctr_0__4 < -1LL * ctr_1__4 - ctr_2__4 + micro_edges_per_macro_edge; ctr_0__4 += 1LL )
               {
                  const double p_affine_0_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_0_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_0_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_1_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_1_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_1_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_2_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_2_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_2_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_3_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__4 );
                  const double p_affine_3_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__4 );
                  const double p_affine_3_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__4 );
                  const double tmp_kernel_op_0__4 = -1.0 * jac_affine_inv_0_0__4 - jac_affine_inv_1_0__4 - jac_affine_inv_2_0__4;
                  const double tmp_kernel_op_1__4 = -1.0 * jac_affine_inv_0_1__4 - jac_affine_inv_1_1__4 - jac_affine_inv_2_1__4;
                  const double tmp_kernel_op_2__4 = -1.0 * jac_affine_inv_0_2__4 - jac_affine_inv_1_2__4 - jac_affine_inv_2_2__4;
                  const double tmp_kernel_op_3__4 = 0.16666666666666663 * abs_det_jac_affine__4;
                  const double tmp_kernel_op_4__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_0_0__4 * tmp_kernel_op_0__4 +
                                                                           jac_affine_inv_0_1__4 * tmp_kernel_op_1__4 +
                                                                           jac_affine_inv_0_2__4 * tmp_kernel_op_2__4 );
                  const double tmp_kernel_op_5__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_1_0__4 * tmp_kernel_op_0__4 +
                                                                           jac_affine_inv_1_1__4 * tmp_kernel_op_1__4 +
                                                                           jac_affine_inv_1_2__4 * tmp_kernel_op_2__4 );
                  const double tmp_kernel_op_6__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_2_0__4 * tmp_kernel_op_0__4 +
                                                                           jac_affine_inv_2_1__4 * tmp_kernel_op_1__4 +
                                                                           jac_affine_inv_2_2__4 * tmp_kernel_op_2__4 );
                  const double tmp_kernel_op_7__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_0_0__4 * jac_affine_inv_1_0__4 +
                                                                           jac_affine_inv_0_1__4 * jac_affine_inv_1_1__4 +
                                                                           jac_affine_inv_0_2__4 * jac_affine_inv_1_2__4 );
                  const double tmp_kernel_op_8__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_0_0__4 * jac_affine_inv_2_0__4 +
                                                                           jac_affine_inv_0_1__4 * jac_affine_inv_2_1__4 +
                                                                           jac_affine_inv_0_2__4 * jac_affine_inv_2_2__4 );
                  const double tmp_kernel_op_9__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_1_0__4 * jac_affine_inv_2_0__4 +
                                                                           jac_affine_inv_1_1__4 * jac_affine_inv_2_1__4 +
                                                                           jac_affine_inv_1_2__4 * jac_affine_inv_2_2__4 );
                  const double elMat_0_0__4 =
                      tmp_kernel_op_3__4 * ( tmp_kernel_op_0__4 * tmp_kernel_op_0__4 + tmp_kernel_op_1__4 * tmp_kernel_op_1__4 +
                                             tmp_kernel_op_2__4 * tmp_kernel_op_2__4 );
                  const double elMat_0_1__4 = tmp_kernel_op_4__4;
                  const double elMat_0_2__4 = tmp_kernel_op_5__4;
                  const double elMat_0_3__4 = tmp_kernel_op_6__4;
                  const double elMat_1_0__4 = tmp_kernel_op_4__4;
                  const double elMat_1_1__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_0_0__4 * jac_affine_inv_0_0__4 +
                                                                     jac_affine_inv_0_1__4 * jac_affine_inv_0_1__4 +
                                                                     jac_affine_inv_0_2__4 * jac_affine_inv_0_2__4 );
                  const double elMat_1_2__4 = tmp_kernel_op_7__4;
                  const double elMat_1_3__4 = tmp_kernel_op_8__4;
                  const double elMat_2_0__4 = tmp_kernel_op_5__4;
                  const double elMat_2_1__4 = tmp_kernel_op_7__4;
                  const double elMat_2_2__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_1_0__4 * jac_affine_inv_1_0__4 +
                                                                     jac_affine_inv_1_1__4 * jac_affine_inv_1_1__4 +
                                                                     jac_affine_inv_1_2__4 * jac_affine_inv_1_2__4 );
                  const double elMat_2_3__4 = tmp_kernel_op_9__4;
                  const double elMat_3_0__4 = tmp_kernel_op_6__4;
                  const double elMat_3_1__4 = tmp_kernel_op_8__4;
                  const double elMat_3_2__4 = tmp_kernel_op_9__4;
                  const double elMat_3_3__4 = tmp_kernel_op_3__4 * ( jac_affine_inv_2_0__4 * jac_affine_inv_2_0__4 +
                                                                     jac_affine_inv_2_1__4 * jac_affine_inv_2_1__4 +
                                                                     jac_affine_inv_2_2__4 * jac_affine_inv_2_2__4 );
                  /*  */
                  /* Apply basis transformation */
                  /*  */
                  std::vector< uint_t > rowIdx__9   = std::vector< uint_t >( 4LL );
                  std::vector< uint_t > colIdx__9   = std::vector< uint_t >( 4LL );
                  std::vector< real_t > matData__21 = std::vector< real_t >( 16LL );
                  rowIdx__9.operator[]( 0 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__9.operator[]( 1 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL -
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__9.operator[]( 2 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ( 2LL + ctr_1__4 ) / 2LL ) -
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__4 ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__9.operator[]( 3 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                         ( -1LL * ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__9.operator[]( 0 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__9.operator[]( 1 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL -
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__9.operator[]( 2 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__4 ) * ( 2LL + ctr_1__4 ) / 2LL ) -
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__4 ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__9.operator[]( 3 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                         ( -1LL * ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  matData__21.operator[]( 0 )  = (real_t) elMat_0_0__4;
                  matData__21.operator[]( 1 )  = (real_t) elMat_0_1__4;
                  matData__21.operator[]( 2 )  = (real_t) elMat_0_2__4;
                  matData__21.operator[]( 3 )  = (real_t) elMat_0_3__4;
                  matData__21.operator[]( 4 )  = (real_t) elMat_1_0__4;
                  matData__21.operator[]( 5 )  = (real_t) elMat_1_1__4;
                  matData__21.operator[]( 6 )  = (real_t) elMat_1_2__4;
                  matData__21.operator[]( 7 )  = (real_t) elMat_1_3__4;
                  matData__21.operator[]( 8 )  = (real_t) elMat_2_0__4;
                  matData__21.operator[]( 9 )  = (real_t) elMat_2_1__4;
                  matData__21.operator[]( 10 ) = (real_t) elMat_2_2__4;
                  matData__21.operator[]( 11 ) = (real_t) elMat_2_3__4;
                  matData__21.operator[]( 12 ) = (real_t) elMat_3_0__4;
                  matData__21.operator[]( 13 ) = (real_t) elMat_3_1__4;
                  matData__21.operator[]( 14 ) = (real_t) elMat_3_2__4;
                  matData__21.operator[]( 15 ) = (real_t) elMat_3_3__4;
                  /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
                  (void) mat;
                  mat->addValues( rowIdx__9, colIdx__9, matData__21 );
               }
            }
         }
      }
      {
         /* CellType.WHITE_DOWN */
         const double tmp_coords_jac_0__3   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__3   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__3   = tmp_coords_jac_1__3 * 0.0;
         const double tmp_coords_jac_3__3   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4__3   = tmp_coords_jac_1__3 * 1.0;
         const double tmp_coords_jac_5__3   = tmp_coords_jac_3__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_6__3   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__3   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_4__3 * tmp_coords_jac_6__3;
         const double tmp_coords_jac_8__3   = tmp_coords_jac_5__3 + tmp_coords_jac_7__3;
         const double tmp_coords_jac_9__3   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__3  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__3  = tmp_coords_jac_10__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_12__3  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__3  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_14__3  = tmp_coords_jac_11__3 + tmp_coords_jac_13__3;
         const double tmp_coords_jac_15__3  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__3  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_17__3  = tmp_coords_jac_16__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_18__3  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__3  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_20__3  = tmp_coords_jac_17__3 + tmp_coords_jac_19__3;
         const double tmp_coords_jac_21__3  = tmp_coords_jac_0__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_22__3  = tmp_coords_jac_4__3 * tmp_coords_jac_9__3;
         const double tmp_coords_jac_23__3  = tmp_coords_jac_15__3 * tmp_coords_jac_4__3;
         const double p_affine_const_0_0__3 = tmp_coords_jac_8__3 + tmp_coords_jac_0__3 * tmp_coords_jac_2__3;
         const double p_affine_const_0_1__3 = tmp_coords_jac_14__3 + tmp_coords_jac_2__3 * tmp_coords_jac_9__3;
         const double p_affine_const_0_2__3 = tmp_coords_jac_20__3 + tmp_coords_jac_15__3 * tmp_coords_jac_2__3;
         const double p_affine_const_1_0__3 =
             tmp_coords_jac_21__3 + tmp_coords_jac_7__3 + tmp_coords_jac_2__3 * tmp_coords_jac_3__3;
         const double p_affine_const_1_1__3 =
             tmp_coords_jac_13__3 + tmp_coords_jac_22__3 + tmp_coords_jac_10__3 * tmp_coords_jac_2__3;
         const double p_affine_const_1_2__3 =
             tmp_coords_jac_19__3 + tmp_coords_jac_23__3 + tmp_coords_jac_16__3 * tmp_coords_jac_2__3;
         const double p_affine_const_2_0__3 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_21__3 + tmp_coords_jac_5__3 +
                                              tmp_coords_jac_2__3 * tmp_coords_jac_6__3;
         const double p_affine_const_2_1__3 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_11__3 + tmp_coords_jac_22__3 +
                                              tmp_coords_jac_12__3 * tmp_coords_jac_2__3;
         const double p_affine_const_2_2__3 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_17__3 + tmp_coords_jac_23__3 +
                                              tmp_coords_jac_18__3 * tmp_coords_jac_2__3;
         const double p_affine_const_3_0__3 = tmp_coords_jac_21__3 + tmp_coords_jac_8__3;
         const double p_affine_const_3_1__3 = tmp_coords_jac_14__3 + tmp_coords_jac_22__3;
         const double p_affine_const_3_2__3 = tmp_coords_jac_20__3 + tmp_coords_jac_23__3;
         const double jac_affine_0_0__3     = p_affine_const_1_0__3 - p_affine_const_0_0__3;
         const double jac_affine_0_1__3     = p_affine_const_2_0__3 - p_affine_const_0_0__3;
         const double jac_affine_0_2__3     = p_affine_const_3_0__3 - p_affine_const_0_0__3;
         const double jac_affine_1_0__3     = p_affine_const_1_1__3 - p_affine_const_0_1__3;
         const double jac_affine_1_1__3     = p_affine_const_2_1__3 - p_affine_const_0_1__3;
         const double tmp_coords_jac_28__3  = jac_affine_0_2__3 * jac_affine_1_1__3;
         const double jac_affine_1_2__3     = p_affine_const_3_1__3 - p_affine_const_0_1__3;
         const double tmp_coords_jac_26__3  = jac_affine_0_1__3 * jac_affine_1_2__3;
         const double jac_affine_2_0__3     = p_affine_const_1_2__3 - p_affine_const_0_2__3;
         const double jac_affine_2_1__3     = p_affine_const_2_2__3 - p_affine_const_0_2__3;
         const double tmp_coords_jac_25__3  = jac_affine_1_2__3 * jac_affine_2_1__3;
         const double jac_affine_2_2__3     = p_affine_const_3_2__3 - p_affine_const_0_2__3;
         const double tmp_coords_jac_24__3  = jac_affine_1_1__3 * jac_affine_2_2__3;
         const double tmp_coords_jac_27__3  = jac_affine_0_1__3 * jac_affine_2_2__3;
         const double tmp_coords_jac_29__3 = jac_affine_0_0__3 * tmp_coords_jac_24__3 + jac_affine_2_0__3 * tmp_coords_jac_26__3 -
                                             jac_affine_0_0__3 * tmp_coords_jac_25__3 - jac_affine_1_0__3 * tmp_coords_jac_27__3 -
                                             jac_affine_2_0__3 * tmp_coords_jac_28__3 +
                                             jac_affine_0_2__3 * jac_affine_1_0__3 * jac_affine_2_1__3;
         const double tmp_coords_jac_30__3  = 1.0 / tmp_coords_jac_29__3;
         const double jac_affine_inv_0_0__3 = tmp_coords_jac_30__3 * ( tmp_coords_jac_24__3 - tmp_coords_jac_25__3 );
         const double jac_affine_inv_0_1__3 =
             tmp_coords_jac_30__3 * ( -1.0 * tmp_coords_jac_27__3 + jac_affine_0_2__3 * jac_affine_2_1__3 );
         const double jac_affine_inv_0_2__3 = tmp_coords_jac_30__3 * ( tmp_coords_jac_26__3 - tmp_coords_jac_28__3 );
         const double jac_affine_inv_1_0__3 =
             tmp_coords_jac_30__3 * ( jac_affine_1_2__3 * jac_affine_2_0__3 - jac_affine_1_0__3 * jac_affine_2_2__3 );
         const double jac_affine_inv_1_1__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_0__3 * jac_affine_2_2__3 - jac_affine_0_2__3 * jac_affine_2_0__3 );
         const double jac_affine_inv_1_2__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_2__3 * jac_affine_1_0__3 - jac_affine_0_0__3 * jac_affine_1_2__3 );
         const double jac_affine_inv_2_0__3 =
             tmp_coords_jac_30__3 * ( jac_affine_1_0__3 * jac_affine_2_1__3 - jac_affine_1_1__3 * jac_affine_2_0__3 );
         const double jac_affine_inv_2_1__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_1__3 * jac_affine_2_0__3 - jac_affine_0_0__3 * jac_affine_2_1__3 );
         const double jac_affine_inv_2_2__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_0__3 * jac_affine_1_1__3 - jac_affine_0_1__3 * jac_affine_1_0__3 );
         const double abs_det_jac_affine__3 = abs( tmp_coords_jac_29__3 );
         for ( int64_t ctr_2__3 = 0LL; ctr_2__3 < micro_edges_per_macro_edge; ctr_2__3 += 1LL )
         {
            for ( int64_t ctr_1__3 = 0LL; ctr_1__3 < -1LL * ctr_2__3 + micro_edges_per_macro_edge; ctr_1__3 += 1LL )
            {
               for ( int64_t ctr_0__3 = 0LL; ctr_0__3 < -2LL - ctr_1__3 - ctr_2__3 + micro_edges_per_macro_edge; ctr_0__3 += 1LL )
               {
                  const double p_affine_0_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__3;
                  const double p_affine_0_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__3;
                  const double p_affine_0_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__3;
                  const double p_affine_1_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_1_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_1_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double tmp_kernel_op_0__3 = -1.0 * jac_affine_inv_0_0__3 - jac_affine_inv_1_0__3 - jac_affine_inv_2_0__3;
                  const double tmp_kernel_op_1__3 = -1.0 * jac_affine_inv_0_1__3 - jac_affine_inv_1_1__3 - jac_affine_inv_2_1__3;
                  const double tmp_kernel_op_2__3 = -1.0 * jac_affine_inv_0_2__3 - jac_affine_inv_1_2__3 - jac_affine_inv_2_2__3;
                  const double tmp_kernel_op_3__3 = 0.16666666666666663 * abs_det_jac_affine__3;
                  const double tmp_kernel_op_4__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_0_0__3 * tmp_kernel_op_0__3 +
                                                                           jac_affine_inv_0_1__3 * tmp_kernel_op_1__3 +
                                                                           jac_affine_inv_0_2__3 * tmp_kernel_op_2__3 );
                  const double tmp_kernel_op_5__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_1_0__3 * tmp_kernel_op_0__3 +
                                                                           jac_affine_inv_1_1__3 * tmp_kernel_op_1__3 +
                                                                           jac_affine_inv_1_2__3 * tmp_kernel_op_2__3 );
                  const double tmp_kernel_op_6__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_2_0__3 * tmp_kernel_op_0__3 +
                                                                           jac_affine_inv_2_1__3 * tmp_kernel_op_1__3 +
                                                                           jac_affine_inv_2_2__3 * tmp_kernel_op_2__3 );
                  const double tmp_kernel_op_7__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_0_0__3 * jac_affine_inv_1_0__3 +
                                                                           jac_affine_inv_0_1__3 * jac_affine_inv_1_1__3 +
                                                                           jac_affine_inv_0_2__3 * jac_affine_inv_1_2__3 );
                  const double tmp_kernel_op_8__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_0_0__3 * jac_affine_inv_2_0__3 +
                                                                           jac_affine_inv_0_1__3 * jac_affine_inv_2_1__3 +
                                                                           jac_affine_inv_0_2__3 * jac_affine_inv_2_2__3 );
                  const double tmp_kernel_op_9__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_1_0__3 * jac_affine_inv_2_0__3 +
                                                                           jac_affine_inv_1_1__3 * jac_affine_inv_2_1__3 +
                                                                           jac_affine_inv_1_2__3 * jac_affine_inv_2_2__3 );
                  const double elMat_0_0__3 =
                      tmp_kernel_op_3__3 * ( tmp_kernel_op_0__3 * tmp_kernel_op_0__3 + tmp_kernel_op_1__3 * tmp_kernel_op_1__3 +
                                             tmp_kernel_op_2__3 * tmp_kernel_op_2__3 );
                  const double elMat_0_1__3 = tmp_kernel_op_4__3;
                  const double elMat_0_2__3 = tmp_kernel_op_5__3;
                  const double elMat_0_3__3 = tmp_kernel_op_6__3;
                  const double elMat_1_0__3 = tmp_kernel_op_4__3;
                  const double elMat_1_1__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_0_0__3 * jac_affine_inv_0_0__3 +
                                                                     jac_affine_inv_0_1__3 * jac_affine_inv_0_1__3 +
                                                                     jac_affine_inv_0_2__3 * jac_affine_inv_0_2__3 );
                  const double elMat_1_2__3 = tmp_kernel_op_7__3;
                  const double elMat_1_3__3 = tmp_kernel_op_8__3;
                  const double elMat_2_0__3 = tmp_kernel_op_5__3;
                  const double elMat_2_1__3 = tmp_kernel_op_7__3;
                  const double elMat_2_2__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_1_0__3 * jac_affine_inv_1_0__3 +
                                                                     jac_affine_inv_1_1__3 * jac_affine_inv_1_1__3 +
                                                                     jac_affine_inv_1_2__3 * jac_affine_inv_1_2__3 );
                  const double elMat_2_3__3 = tmp_kernel_op_9__3;
                  const double elMat_3_0__3 = tmp_kernel_op_6__3;
                  const double elMat_3_1__3 = tmp_kernel_op_8__3;
                  const double elMat_3_2__3 = tmp_kernel_op_9__3;
                  const double elMat_3_3__3 = tmp_kernel_op_3__3 * ( jac_affine_inv_2_0__3 * jac_affine_inv_2_0__3 +
                                                                     jac_affine_inv_2_1__3 * jac_affine_inv_2_1__3 +
                                                                     jac_affine_inv_2_2__3 * jac_affine_inv_2_2__3 );
                  /*  */
                  /* Apply basis transformation */
                  /*  */
                  std::vector< uint_t > rowIdx__8   = std::vector< uint_t >( 4LL );
                  std::vector< uint_t > colIdx__8   = std::vector< uint_t >( 4LL );
                  std::vector< real_t > matData__20 = std::vector< real_t >( 16LL );
                  rowIdx__8.operator[]( 0 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                         ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__3 ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__8.operator[]( 1 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__3 ) * ctr_1__3 / 2LL -
                                         ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ctr_1__3 + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__8.operator[]( 2 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL ) -
                                         ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__8.operator[]( 3 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                         ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__8.operator[]( 0 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                         ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__3 ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__8.operator[]( 1 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__3 ) * ctr_1__3 / 2LL -
                                         ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ctr_1__3 + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__8.operator[]( 2 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL ) -
                                         ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__8.operator[]( 3 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                         ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  matData__20.operator[]( 0 )  = (real_t) elMat_0_0__3;
                  matData__20.operator[]( 1 )  = (real_t) elMat_0_1__3;
                  matData__20.operator[]( 2 )  = (real_t) elMat_0_2__3;
                  matData__20.operator[]( 3 )  = (real_t) elMat_0_3__3;
                  matData__20.operator[]( 4 )  = (real_t) elMat_1_0__3;
                  matData__20.operator[]( 5 )  = (real_t) elMat_1_1__3;
                  matData__20.operator[]( 6 )  = (real_t) elMat_1_2__3;
                  matData__20.operator[]( 7 )  = (real_t) elMat_1_3__3;
                  matData__20.operator[]( 8 )  = (real_t) elMat_2_0__3;
                  matData__20.operator[]( 9 )  = (real_t) elMat_2_1__3;
                  matData__20.operator[]( 10 ) = (real_t) elMat_2_2__3;
                  matData__20.operator[]( 11 ) = (real_t) elMat_2_3__3;
                  matData__20.operator[]( 12 ) = (real_t) elMat_3_0__3;
                  matData__20.operator[]( 13 ) = (real_t) elMat_3_1__3;
                  matData__20.operator[]( 14 ) = (real_t) elMat_3_2__3;
                  matData__20.operator[]( 15 ) = (real_t) elMat_3_3__3;
                  /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
                  (void) mat;
                  mat->addValues( rowIdx__8, colIdx__8, matData__20 );
               }
            }
         }
      }
      {
         /* CellType.BLUE_UP */
         const double tmp_coords_jac_0__2   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__2   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__2   = tmp_coords_jac_1__2 * 0.0;
         const double tmp_coords_jac_3__2   = tmp_coords_jac_0__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_4__2   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__2   = tmp_coords_jac_1__2 * 1.0;
         const double tmp_coords_jac_6__2   = tmp_coords_jac_4__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_7__2   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_8__2   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__2 * tmp_coords_jac_7__2;
         const double tmp_coords_jac_9__2   = tmp_coords_jac_6__2 + tmp_coords_jac_8__2;
         const double tmp_coords_jac_10__2  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__2  = tmp_coords_jac_10__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_12__2  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__2  = tmp_coords_jac_12__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_14__2  = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_15__2  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_14__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_16__2  = tmp_coords_jac_13__2 + tmp_coords_jac_15__2;
         const double tmp_coords_jac_17__2  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18__2  = tmp_coords_jac_17__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_19__2  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_20__2  = tmp_coords_jac_19__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_21__2  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_22__2  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_2__2 * tmp_coords_jac_21__2;
         const double tmp_coords_jac_23__2  = tmp_coords_jac_20__2 + tmp_coords_jac_22__2;
         const double tmp_coords_jac_24__2  = tmp_coords_jac_0__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_25__2  = tmp_coords_jac_10__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_26__2  = tmp_coords_jac_17__2 * tmp_coords_jac_5__2;
         const double p_affine_const_0_0__2 = tmp_coords_jac_3__2 + tmp_coords_jac_9__2;
         const double p_affine_const_0_1__2 = tmp_coords_jac_11__2 + tmp_coords_jac_16__2;
         const double p_affine_const_0_2__2 = tmp_coords_jac_18__2 + tmp_coords_jac_23__2;
         const double p_affine_const_1_0__2 =
             tmp_coords_jac_24__2 + tmp_coords_jac_8__2 + tmp_coords_jac_2__2 * tmp_coords_jac_4__2;
         const double p_affine_const_1_1__2 =
             tmp_coords_jac_15__2 + tmp_coords_jac_25__2 + tmp_coords_jac_12__2 * tmp_coords_jac_2__2;
         const double p_affine_const_1_2__2 =
             tmp_coords_jac_22__2 + tmp_coords_jac_26__2 + tmp_coords_jac_19__2 * tmp_coords_jac_2__2;
         const double p_affine_const_2_0__2 = tmp_coords_jac_24__2 + tmp_coords_jac_9__2;
         const double p_affine_const_2_1__2 = tmp_coords_jac_16__2 + tmp_coords_jac_25__2;
         const double p_affine_const_2_2__2 = tmp_coords_jac_23__2 + tmp_coords_jac_26__2;
         const double p_affine_const_3_0__2 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__2 + tmp_coords_jac_6__2 + tmp_coords_jac_5__2 * tmp_coords_jac_7__2;
         const double p_affine_const_3_1__2 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_11__2 + tmp_coords_jac_13__2 +
                                              tmp_coords_jac_14__2 * tmp_coords_jac_5__2;
         const double p_affine_const_3_2__2 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18__2 + tmp_coords_jac_20__2 +
                                              tmp_coords_jac_21__2 * tmp_coords_jac_5__2;
         const double jac_affine_0_0__2    = p_affine_const_1_0__2 - p_affine_const_0_0__2;
         const double jac_affine_0_1__2    = p_affine_const_2_0__2 - p_affine_const_0_0__2;
         const double jac_affine_0_2__2    = p_affine_const_3_0__2 - p_affine_const_0_0__2;
         const double jac_affine_1_0__2    = p_affine_const_1_1__2 - p_affine_const_0_1__2;
         const double jac_affine_1_1__2    = p_affine_const_2_1__2 - p_affine_const_0_1__2;
         const double tmp_coords_jac_31__1 = jac_affine_0_2__2 * jac_affine_1_1__2;
         const double jac_affine_1_2__2    = p_affine_const_3_1__2 - p_affine_const_0_1__2;
         const double tmp_coords_jac_29__2 = jac_affine_0_1__2 * jac_affine_1_2__2;
         const double jac_affine_2_0__2    = p_affine_const_1_2__2 - p_affine_const_0_2__2;
         const double jac_affine_2_1__2    = p_affine_const_2_2__2 - p_affine_const_0_2__2;
         const double tmp_coords_jac_28__2 = jac_affine_1_2__2 * jac_affine_2_1__2;
         const double jac_affine_2_2__2    = p_affine_const_3_2__2 - p_affine_const_0_2__2;
         const double tmp_coords_jac_27__2 = jac_affine_1_1__2 * jac_affine_2_2__2;
         const double tmp_coords_jac_30__2 = jac_affine_0_1__2 * jac_affine_2_2__2;
         const double tmp_coords_jac_32__1 = jac_affine_0_0__2 * tmp_coords_jac_27__2 + jac_affine_2_0__2 * tmp_coords_jac_29__2 -
                                             jac_affine_0_0__2 * tmp_coords_jac_28__2 - jac_affine_1_0__2 * tmp_coords_jac_30__2 -
                                             jac_affine_2_0__2 * tmp_coords_jac_31__1 +
                                             jac_affine_0_2__2 * jac_affine_1_0__2 * jac_affine_2_1__2;
         const double tmp_coords_jac_33__1  = 1.0 / tmp_coords_jac_32__1;
         const double jac_affine_inv_0_0__2 = tmp_coords_jac_33__1 * ( tmp_coords_jac_27__2 - tmp_coords_jac_28__2 );
         const double jac_affine_inv_0_1__2 =
             tmp_coords_jac_33__1 * ( -1.0 * tmp_coords_jac_30__2 + jac_affine_0_2__2 * jac_affine_2_1__2 );
         const double jac_affine_inv_0_2__2 = tmp_coords_jac_33__1 * ( tmp_coords_jac_29__2 - tmp_coords_jac_31__1 );
         const double jac_affine_inv_1_0__2 =
             tmp_coords_jac_33__1 * ( jac_affine_1_2__2 * jac_affine_2_0__2 - jac_affine_1_0__2 * jac_affine_2_2__2 );
         const double jac_affine_inv_1_1__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_0__2 * jac_affine_2_2__2 - jac_affine_0_2__2 * jac_affine_2_0__2 );
         const double jac_affine_inv_1_2__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_2__2 * jac_affine_1_0__2 - jac_affine_0_0__2 * jac_affine_1_2__2 );
         const double jac_affine_inv_2_0__2 =
             tmp_coords_jac_33__1 * ( jac_affine_1_0__2 * jac_affine_2_1__2 - jac_affine_1_1__2 * jac_affine_2_0__2 );
         const double jac_affine_inv_2_1__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_1__2 * jac_affine_2_0__2 - jac_affine_0_0__2 * jac_affine_2_1__2 );
         const double jac_affine_inv_2_2__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_0__2 * jac_affine_1_1__2 - jac_affine_0_1__2 * jac_affine_1_0__2 );
         const double abs_det_jac_affine__2 = abs( tmp_coords_jac_32__1 );
         for ( int64_t ctr_2__2 = 0LL; ctr_2__2 < micro_edges_per_macro_edge; ctr_2__2 += 1LL )
         {
            for ( int64_t ctr_1__2 = 0LL; ctr_1__2 < -1LL * ctr_2__2 + micro_edges_per_macro_edge; ctr_1__2 += 1LL )
            {
               for ( int64_t ctr_0__2 = 0LL; ctr_0__2 < -1LL - ctr_1__2 - ctr_2__2 + micro_edges_per_macro_edge; ctr_0__2 += 1LL )
               {
                  const double p_affine_0_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_0_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_0_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_1_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_1_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_1_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_2_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_2_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_2_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_3_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__2 );
                  const double p_affine_3_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__2 );
                  const double p_affine_3_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__2 );
                  const double tmp_kernel_op_0__2 = -1.0 * jac_affine_inv_0_0__2 - jac_affine_inv_1_0__2 - jac_affine_inv_2_0__2;
                  const double tmp_kernel_op_1__2 = -1.0 * jac_affine_inv_0_1__2 - jac_affine_inv_1_1__2 - jac_affine_inv_2_1__2;
                  const double tmp_kernel_op_2__2 = -1.0 * jac_affine_inv_0_2__2 - jac_affine_inv_1_2__2 - jac_affine_inv_2_2__2;
                  const double tmp_kernel_op_3__2 = 0.16666666666666663 * abs_det_jac_affine__2;
                  const double tmp_kernel_op_4__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_0_0__2 * tmp_kernel_op_0__2 +
                                                                           jac_affine_inv_0_1__2 * tmp_kernel_op_1__2 +
                                                                           jac_affine_inv_0_2__2 * tmp_kernel_op_2__2 );
                  const double tmp_kernel_op_5__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_1_0__2 * tmp_kernel_op_0__2 +
                                                                           jac_affine_inv_1_1__2 * tmp_kernel_op_1__2 +
                                                                           jac_affine_inv_1_2__2 * tmp_kernel_op_2__2 );
                  const double tmp_kernel_op_6__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_2_0__2 * tmp_kernel_op_0__2 +
                                                                           jac_affine_inv_2_1__2 * tmp_kernel_op_1__2 +
                                                                           jac_affine_inv_2_2__2 * tmp_kernel_op_2__2 );
                  const double tmp_kernel_op_7__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_0_0__2 * jac_affine_inv_1_0__2 +
                                                                           jac_affine_inv_0_1__2 * jac_affine_inv_1_1__2 +
                                                                           jac_affine_inv_0_2__2 * jac_affine_inv_1_2__2 );
                  const double tmp_kernel_op_8__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_0_0__2 * jac_affine_inv_2_0__2 +
                                                                           jac_affine_inv_0_1__2 * jac_affine_inv_2_1__2 +
                                                                           jac_affine_inv_0_2__2 * jac_affine_inv_2_2__2 );
                  const double tmp_kernel_op_9__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_1_0__2 * jac_affine_inv_2_0__2 +
                                                                           jac_affine_inv_1_1__2 * jac_affine_inv_2_1__2 +
                                                                           jac_affine_inv_1_2__2 * jac_affine_inv_2_2__2 );
                  const double elMat_0_0__2 =
                      tmp_kernel_op_3__2 * ( tmp_kernel_op_0__2 * tmp_kernel_op_0__2 + tmp_kernel_op_1__2 * tmp_kernel_op_1__2 +
                                             tmp_kernel_op_2__2 * tmp_kernel_op_2__2 );
                  const double elMat_0_1__2 = tmp_kernel_op_4__2;
                  const double elMat_0_2__2 = tmp_kernel_op_5__2;
                  const double elMat_0_3__2 = tmp_kernel_op_6__2;
                  const double elMat_1_0__2 = tmp_kernel_op_4__2;
                  const double elMat_1_1__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_0_0__2 * jac_affine_inv_0_0__2 +
                                                                     jac_affine_inv_0_1__2 * jac_affine_inv_0_1__2 +
                                                                     jac_affine_inv_0_2__2 * jac_affine_inv_0_2__2 );
                  const double elMat_1_2__2 = tmp_kernel_op_7__2;
                  const double elMat_1_3__2 = tmp_kernel_op_8__2;
                  const double elMat_2_0__2 = tmp_kernel_op_5__2;
                  const double elMat_2_1__2 = tmp_kernel_op_7__2;
                  const double elMat_2_2__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_1_0__2 * jac_affine_inv_1_0__2 +
                                                                     jac_affine_inv_1_1__2 * jac_affine_inv_1_1__2 +
                                                                     jac_affine_inv_1_2__2 * jac_affine_inv_1_2__2 );
                  const double elMat_2_3__2 = tmp_kernel_op_9__2;
                  const double elMat_3_0__2 = tmp_kernel_op_6__2;
                  const double elMat_3_1__2 = tmp_kernel_op_8__2;
                  const double elMat_3_2__2 = tmp_kernel_op_9__2;
                  const double elMat_3_3__2 = tmp_kernel_op_3__2 * ( jac_affine_inv_2_0__2 * jac_affine_inv_2_0__2 +
                                                                     jac_affine_inv_2_1__2 * jac_affine_inv_2_1__2 +
                                                                     jac_affine_inv_2_2__2 * jac_affine_inv_2_2__2 );
                  /*  */
                  /* Apply basis transformation */
                  /*  */
                  std::vector< uint_t > rowIdx__7   = std::vector< uint_t >( 4LL );
                  std::vector< uint_t > colIdx__7   = std::vector< uint_t >( 4LL );
                  std::vector< real_t > matData__19 = std::vector< real_t >( 16LL );
                  rowIdx__7.operator[]( 0 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__7.operator[]( 1 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL ) -
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__7.operator[]( 2 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL -
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__7.operator[]( 3 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                         ( -1LL * ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__7.operator[]( 0 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__7.operator[]( 1 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL ) -
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__7.operator[]( 2 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL -
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__7.operator[]( 3 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                         ( -1LL * ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  matData__19.operator[]( 0 )  = (real_t) elMat_0_0__2;
                  matData__19.operator[]( 1 )  = (real_t) elMat_0_1__2;
                  matData__19.operator[]( 2 )  = (real_t) elMat_0_2__2;
                  matData__19.operator[]( 3 )  = (real_t) elMat_0_3__2;
                  matData__19.operator[]( 4 )  = (real_t) elMat_1_0__2;
                  matData__19.operator[]( 5 )  = (real_t) elMat_1_1__2;
                  matData__19.operator[]( 6 )  = (real_t) elMat_1_2__2;
                  matData__19.operator[]( 7 )  = (real_t) elMat_1_3__2;
                  matData__19.operator[]( 8 )  = (real_t) elMat_2_0__2;
                  matData__19.operator[]( 9 )  = (real_t) elMat_2_1__2;
                  matData__19.operator[]( 10 ) = (real_t) elMat_2_2__2;
                  matData__19.operator[]( 11 ) = (real_t) elMat_2_3__2;
                  matData__19.operator[]( 12 ) = (real_t) elMat_3_0__2;
                  matData__19.operator[]( 13 ) = (real_t) elMat_3_1__2;
                  matData__19.operator[]( 14 ) = (real_t) elMat_3_2__2;
                  matData__19.operator[]( 15 ) = (real_t) elMat_3_3__2;
                  /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
                  (void) mat;
                  mat->addValues( rowIdx__7, colIdx__7, matData__19 );
               }
            }
         }
      }
      {
         /* CellType.BLUE_DOWN */
         const double tmp_coords_jac_0__1   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__1   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__1   = tmp_coords_jac_1__1 * 0.0;
         const double tmp_coords_jac_3__1   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4__1   = tmp_coords_jac_1__1 * 1.0;
         const double tmp_coords_jac_5__1   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_6__1   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__1 * tmp_coords_jac_5__1;
         const double tmp_coords_jac_7__1   = tmp_coords_jac_6__1 + tmp_coords_jac_3__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_8__1   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_9__1   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__1  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__1  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10__1 * tmp_coords_jac_2__1;
         const double tmp_coords_jac_12__1  = tmp_coords_jac_11__1 + tmp_coords_jac_4__1 * tmp_coords_jac_9__1;
         const double tmp_coords_jac_13__1  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_14__1  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_15__1  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__1  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__1 * tmp_coords_jac_2__1;
         const double tmp_coords_jac_17__1  = tmp_coords_jac_16__1 + tmp_coords_jac_14__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_18__1  = tmp_coords_jac_0__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_19__1  = tmp_coords_jac_18__1 + tmp_coords_jac_2__1 * tmp_coords_jac_3__1;
         const double tmp_coords_jac_20__1  = tmp_coords_jac_4__1 * tmp_coords_jac_8__1;
         const double tmp_coords_jac_21__1  = tmp_coords_jac_20__1 + tmp_coords_jac_2__1 * tmp_coords_jac_9__1;
         const double tmp_coords_jac_22__1  = tmp_coords_jac_13__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_23__1  = tmp_coords_jac_22__1 + tmp_coords_jac_14__1 * tmp_coords_jac_2__1;
         const double p_affine_const_0_0__1 = tmp_coords_jac_7__1 + tmp_coords_jac_0__1 * tmp_coords_jac_2__1;
         const double p_affine_const_0_1__1 = tmp_coords_jac_12__1 + tmp_coords_jac_2__1 * tmp_coords_jac_8__1;
         const double p_affine_const_0_2__1 = tmp_coords_jac_17__1 + tmp_coords_jac_13__1 * tmp_coords_jac_2__1;
         const double p_affine_const_1_0__1 = tmp_coords_jac_19__1 + tmp_coords_jac_6__1;
         const double p_affine_const_1_1__1 = tmp_coords_jac_11__1 + tmp_coords_jac_21__1;
         const double p_affine_const_1_2__1 = tmp_coords_jac_16__1 + tmp_coords_jac_23__1;
         const double p_affine_const_2_0__1 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_19__1 + tmp_coords_jac_4__1 * tmp_coords_jac_5__1;
         const double p_affine_const_2_1__1 =
             macro_vertex_coord_id_0comp1 + tmp_coords_jac_21__1 + tmp_coords_jac_10__1 * tmp_coords_jac_4__1;
         const double p_affine_const_2_2__1 =
             macro_vertex_coord_id_0comp2 + tmp_coords_jac_23__1 + tmp_coords_jac_15__1 * tmp_coords_jac_4__1;
         const double p_affine_const_3_0__1 = tmp_coords_jac_18__1 + tmp_coords_jac_7__1;
         const double p_affine_const_3_1__1 = tmp_coords_jac_12__1 + tmp_coords_jac_20__1;
         const double p_affine_const_3_2__1 = tmp_coords_jac_17__1 + tmp_coords_jac_22__1;
         const double jac_affine_0_0__1     = p_affine_const_1_0__1 - p_affine_const_0_0__1;
         const double jac_affine_0_1__1     = p_affine_const_2_0__1 - p_affine_const_0_0__1;
         const double jac_affine_0_2__1     = p_affine_const_3_0__1 - p_affine_const_0_0__1;
         const double jac_affine_1_0__1     = p_affine_const_1_1__1 - p_affine_const_0_1__1;
         const double jac_affine_1_1__1     = p_affine_const_2_1__1 - p_affine_const_0_1__1;
         const double tmp_coords_jac_28__1  = jac_affine_0_2__1 * jac_affine_1_1__1;
         const double jac_affine_1_2__1     = p_affine_const_3_1__1 - p_affine_const_0_1__1;
         const double tmp_coords_jac_26__1  = jac_affine_0_1__1 * jac_affine_1_2__1;
         const double jac_affine_2_0__1     = p_affine_const_1_2__1 - p_affine_const_0_2__1;
         const double jac_affine_2_1__1     = p_affine_const_2_2__1 - p_affine_const_0_2__1;
         const double tmp_coords_jac_25__1  = jac_affine_1_2__1 * jac_affine_2_1__1;
         const double jac_affine_2_2__1     = p_affine_const_3_2__1 - p_affine_const_0_2__1;
         const double tmp_coords_jac_24__1  = jac_affine_1_1__1 * jac_affine_2_2__1;
         const double tmp_coords_jac_27__1  = jac_affine_0_1__1 * jac_affine_2_2__1;
         const double tmp_coords_jac_29__1 = jac_affine_0_0__1 * tmp_coords_jac_24__1 + jac_affine_2_0__1 * tmp_coords_jac_26__1 -
                                             jac_affine_0_0__1 * tmp_coords_jac_25__1 - jac_affine_1_0__1 * tmp_coords_jac_27__1 -
                                             jac_affine_2_0__1 * tmp_coords_jac_28__1 +
                                             jac_affine_0_2__1 * jac_affine_1_0__1 * jac_affine_2_1__1;
         const double tmp_coords_jac_30__1  = 1.0 / tmp_coords_jac_29__1;
         const double jac_affine_inv_0_0__1 = tmp_coords_jac_30__1 * ( tmp_coords_jac_24__1 - tmp_coords_jac_25__1 );
         const double jac_affine_inv_0_1__1 =
             tmp_coords_jac_30__1 * ( -1.0 * tmp_coords_jac_27__1 + jac_affine_0_2__1 * jac_affine_2_1__1 );
         const double jac_affine_inv_0_2__1 = tmp_coords_jac_30__1 * ( tmp_coords_jac_26__1 - tmp_coords_jac_28__1 );
         const double jac_affine_inv_1_0__1 =
             tmp_coords_jac_30__1 * ( jac_affine_1_2__1 * jac_affine_2_0__1 - jac_affine_1_0__1 * jac_affine_2_2__1 );
         const double jac_affine_inv_1_1__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_0__1 * jac_affine_2_2__1 - jac_affine_0_2__1 * jac_affine_2_0__1 );
         const double jac_affine_inv_1_2__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_2__1 * jac_affine_1_0__1 - jac_affine_0_0__1 * jac_affine_1_2__1 );
         const double jac_affine_inv_2_0__1 =
             tmp_coords_jac_30__1 * ( jac_affine_1_0__1 * jac_affine_2_1__1 - jac_affine_1_1__1 * jac_affine_2_0__1 );
         const double jac_affine_inv_2_1__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_1__1 * jac_affine_2_0__1 - jac_affine_0_0__1 * jac_affine_2_1__1 );
         const double jac_affine_inv_2_2__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_0__1 * jac_affine_1_1__1 - jac_affine_0_1__1 * jac_affine_1_0__1 );
         const double abs_det_jac_affine__1 = abs( tmp_coords_jac_29__1 );
         for ( int64_t ctr_2__1 = 0LL; ctr_2__1 < micro_edges_per_macro_edge; ctr_2__1 += 1LL )
         {
            for ( int64_t ctr_1__1 = 0LL; ctr_1__1 < -1LL * ctr_2__1 + micro_edges_per_macro_edge; ctr_1__1 += 1LL )
            {
               for ( int64_t ctr_0__1 = 0LL; ctr_0__1 < -1LL - ctr_1__1 - ctr_2__1 + micro_edges_per_macro_edge; ctr_0__1 += 1LL )
               {
                  const double p_affine_0_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__1;
                  const double p_affine_0_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__1;
                  const double p_affine_0_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__1;
                  const double p_affine_1_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_1_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_1_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double tmp_kernel_op_0__1 = -1.0 * jac_affine_inv_0_0__1 - jac_affine_inv_1_0__1 - jac_affine_inv_2_0__1;
                  const double tmp_kernel_op_1__1 = -1.0 * jac_affine_inv_0_1__1 - jac_affine_inv_1_1__1 - jac_affine_inv_2_1__1;
                  const double tmp_kernel_op_2__1 = -1.0 * jac_affine_inv_0_2__1 - jac_affine_inv_1_2__1 - jac_affine_inv_2_2__1;
                  const double tmp_kernel_op_3__1 = 0.16666666666666663 * abs_det_jac_affine__1;
                  const double tmp_kernel_op_4__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_0_0__1 * tmp_kernel_op_0__1 +
                                                                           jac_affine_inv_0_1__1 * tmp_kernel_op_1__1 +
                                                                           jac_affine_inv_0_2__1 * tmp_kernel_op_2__1 );
                  const double tmp_kernel_op_5__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_1_0__1 * tmp_kernel_op_0__1 +
                                                                           jac_affine_inv_1_1__1 * tmp_kernel_op_1__1 +
                                                                           jac_affine_inv_1_2__1 * tmp_kernel_op_2__1 );
                  const double tmp_kernel_op_6__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_2_0__1 * tmp_kernel_op_0__1 +
                                                                           jac_affine_inv_2_1__1 * tmp_kernel_op_1__1 +
                                                                           jac_affine_inv_2_2__1 * tmp_kernel_op_2__1 );
                  const double tmp_kernel_op_7__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_0_0__1 * jac_affine_inv_1_0__1 +
                                                                           jac_affine_inv_0_1__1 * jac_affine_inv_1_1__1 +
                                                                           jac_affine_inv_0_2__1 * jac_affine_inv_1_2__1 );
                  const double tmp_kernel_op_8__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_0_0__1 * jac_affine_inv_2_0__1 +
                                                                           jac_affine_inv_0_1__1 * jac_affine_inv_2_1__1 +
                                                                           jac_affine_inv_0_2__1 * jac_affine_inv_2_2__1 );
                  const double tmp_kernel_op_9__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_1_0__1 * jac_affine_inv_2_0__1 +
                                                                           jac_affine_inv_1_1__1 * jac_affine_inv_2_1__1 +
                                                                           jac_affine_inv_1_2__1 * jac_affine_inv_2_2__1 );
                  const double elMat_0_0__1 =
                      tmp_kernel_op_3__1 * ( tmp_kernel_op_0__1 * tmp_kernel_op_0__1 + tmp_kernel_op_1__1 * tmp_kernel_op_1__1 +
                                             tmp_kernel_op_2__1 * tmp_kernel_op_2__1 );
                  const double elMat_0_1__1 = tmp_kernel_op_4__1;
                  const double elMat_0_2__1 = tmp_kernel_op_5__1;
                  const double elMat_0_3__1 = tmp_kernel_op_6__1;
                  const double elMat_1_0__1 = tmp_kernel_op_4__1;
                  const double elMat_1_1__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_0_0__1 * jac_affine_inv_0_0__1 +
                                                                     jac_affine_inv_0_1__1 * jac_affine_inv_0_1__1 +
                                                                     jac_affine_inv_0_2__1 * jac_affine_inv_0_2__1 );
                  const double elMat_1_2__1 = tmp_kernel_op_7__1;
                  const double elMat_1_3__1 = tmp_kernel_op_8__1;
                  const double elMat_2_0__1 = tmp_kernel_op_5__1;
                  const double elMat_2_1__1 = tmp_kernel_op_7__1;
                  const double elMat_2_2__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_1_0__1 * jac_affine_inv_1_0__1 +
                                                                     jac_affine_inv_1_1__1 * jac_affine_inv_1_1__1 +
                                                                     jac_affine_inv_1_2__1 * jac_affine_inv_1_2__1 );
                  const double elMat_2_3__1 = tmp_kernel_op_9__1;
                  const double elMat_3_0__1 = tmp_kernel_op_6__1;
                  const double elMat_3_1__1 = tmp_kernel_op_8__1;
                  const double elMat_3_2__1 = tmp_kernel_op_9__1;
                  const double elMat_3_3__1 = tmp_kernel_op_3__1 * ( jac_affine_inv_2_0__1 * jac_affine_inv_2_0__1 +
                                                                     jac_affine_inv_2_1__1 * jac_affine_inv_2_1__1 +
                                                                     jac_affine_inv_2_2__1 * jac_affine_inv_2_2__1 );
                  /*  */
                  /* Apply basis transformation */
                  /*  */
                  std::vector< uint_t > rowIdx__6   = std::vector< uint_t >( 4LL );
                  std::vector< uint_t > colIdx__6   = std::vector< uint_t >( 4LL );
                  std::vector< real_t > matData__18 = std::vector< real_t >( 16LL );
                  rowIdx__6.operator[]( 0 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__1 ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__6.operator[]( 1 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL ) -
                                         ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__6.operator[]( 2 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL -
                                         ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__6.operator[]( 3 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                         ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__1 ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__6.operator[]( 0 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__1 ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__6.operator[]( 1 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL ) -
                                         ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__6.operator[]( 2 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL -
                                         ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__6.operator[]( 3 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                         ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__1 ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  matData__18.operator[]( 0 )  = (real_t) elMat_0_0__1;
                  matData__18.operator[]( 1 )  = (real_t) elMat_0_1__1;
                  matData__18.operator[]( 2 )  = (real_t) elMat_0_2__1;
                  matData__18.operator[]( 3 )  = (real_t) elMat_0_3__1;
                  matData__18.operator[]( 4 )  = (real_t) elMat_1_0__1;
                  matData__18.operator[]( 5 )  = (real_t) elMat_1_1__1;
                  matData__18.operator[]( 6 )  = (real_t) elMat_1_2__1;
                  matData__18.operator[]( 7 )  = (real_t) elMat_1_3__1;
                  matData__18.operator[]( 8 )  = (real_t) elMat_2_0__1;
                  matData__18.operator[]( 9 )  = (real_t) elMat_2_1__1;
                  matData__18.operator[]( 10 ) = (real_t) elMat_2_2__1;
                  matData__18.operator[]( 11 ) = (real_t) elMat_2_3__1;
                  matData__18.operator[]( 12 ) = (real_t) elMat_3_0__1;
                  matData__18.operator[]( 13 ) = (real_t) elMat_3_1__1;
                  matData__18.operator[]( 14 ) = (real_t) elMat_3_2__1;
                  matData__18.operator[]( 15 ) = (real_t) elMat_3_3__1;
                  /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
                  (void) mat;
                  mat->addValues( rowIdx__6, colIdx__6, matData__18 );
               }
            }
         }
      }
      {
         /* CellType.GREEN_UP */
         const double tmp_coords_jac_0__0   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__0   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__0   = tmp_coords_jac_1__0 * 0.0;
         const double tmp_coords_jac_3__0   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_0__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_4__0   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__0   = tmp_coords_jac_1__0 * 1.0;
         const double tmp_coords_jac_6__0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__0   = tmp_coords_jac_2__0 * tmp_coords_jac_6__0;
         const double tmp_coords_jac_8__0   = tmp_coords_jac_7__0 + tmp_coords_jac_4__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_9__0   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__0  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_2__0 * tmp_coords_jac_9__0;
         const double tmp_coords_jac_11__0  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12__0  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__0  = tmp_coords_jac_12__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_14__0  = tmp_coords_jac_13__0 + tmp_coords_jac_11__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_15__0  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__0  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_17__0  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18__0  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__0  = tmp_coords_jac_18__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_20__0  = tmp_coords_jac_19__0 + tmp_coords_jac_17__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_21__0  = tmp_coords_jac_2__0 * tmp_coords_jac_4__0;
         const double tmp_coords_jac_22__0  = tmp_coords_jac_11__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_23__0  = tmp_coords_jac_17__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_24__0  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_0__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_25__0  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_5__0 * tmp_coords_jac_9__0;
         const double tmp_coords_jac_26__0  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__0 * tmp_coords_jac_5__0;
         const double p_affine_const_0_0__0 = tmp_coords_jac_3__0 + tmp_coords_jac_8__0;
         const double p_affine_const_0_1__0 = tmp_coords_jac_10__0 + tmp_coords_jac_14__0;
         const double p_affine_const_0_2__0 = tmp_coords_jac_16__0 + tmp_coords_jac_20__0;
         const double p_affine_const_1_0__0 =
             tmp_coords_jac_21__0 + tmp_coords_jac_3__0 + tmp_coords_jac_5__0 * tmp_coords_jac_6__0;
         const double p_affine_const_1_1__0 =
             tmp_coords_jac_10__0 + tmp_coords_jac_22__0 + tmp_coords_jac_12__0 * tmp_coords_jac_5__0;
         const double p_affine_const_1_2__0 =
             tmp_coords_jac_16__0 + tmp_coords_jac_23__0 + tmp_coords_jac_18__0 * tmp_coords_jac_5__0;
         const double p_affine_const_2_0__0 = tmp_coords_jac_21__0 + tmp_coords_jac_24__0 + tmp_coords_jac_7__0;
         const double p_affine_const_2_1__0 = tmp_coords_jac_13__0 + tmp_coords_jac_22__0 + tmp_coords_jac_25__0;
         const double p_affine_const_2_2__0 = tmp_coords_jac_19__0 + tmp_coords_jac_23__0 + tmp_coords_jac_26__0;
         const double p_affine_const_3_0__0 = tmp_coords_jac_24__0 + tmp_coords_jac_8__0;
         const double p_affine_const_3_1__0 = tmp_coords_jac_14__0 + tmp_coords_jac_25__0;
         const double p_affine_const_3_2__0 = tmp_coords_jac_20__0 + tmp_coords_jac_26__0;
         const double jac_affine_0_0__0     = p_affine_const_1_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_1__0     = p_affine_const_2_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_2__0     = p_affine_const_3_0__0 - p_affine_const_0_0__0;
         const double jac_affine_1_0__0     = p_affine_const_1_1__0 - p_affine_const_0_1__0;
         const double jac_affine_1_1__0     = p_affine_const_2_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_31__0  = jac_affine_0_2__0 * jac_affine_1_1__0;
         const double jac_affine_1_2__0     = p_affine_const_3_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_29__0  = jac_affine_0_1__0 * jac_affine_1_2__0;
         const double jac_affine_2_0__0     = p_affine_const_1_2__0 - p_affine_const_0_2__0;
         const double jac_affine_2_1__0     = p_affine_const_2_2__0 - p_affine_const_0_2__0;
         const double tmp_coords_jac_28__0  = jac_affine_1_2__0 * jac_affine_2_1__0;
         const double jac_affine_2_2__0     = p_affine_const_3_2__0 - p_affine_const_0_2__0;
         const double tmp_coords_jac_27__0  = jac_affine_1_1__0 * jac_affine_2_2__0;
         const double tmp_coords_jac_30__0  = jac_affine_0_1__0 * jac_affine_2_2__0;
         const double tmp_coords_jac_32__0 = jac_affine_0_0__0 * tmp_coords_jac_27__0 + jac_affine_2_0__0 * tmp_coords_jac_29__0 -
                                             jac_affine_0_0__0 * tmp_coords_jac_28__0 - jac_affine_1_0__0 * tmp_coords_jac_30__0 -
                                             jac_affine_2_0__0 * tmp_coords_jac_31__0 +
                                             jac_affine_0_2__0 * jac_affine_1_0__0 * jac_affine_2_1__0;
         const double tmp_coords_jac_33__0  = 1.0 / tmp_coords_jac_32__0;
         const double jac_affine_inv_0_0__0 = tmp_coords_jac_33__0 * ( tmp_coords_jac_27__0 - tmp_coords_jac_28__0 );
         const double jac_affine_inv_0_1__0 =
             tmp_coords_jac_33__0 * ( -1.0 * tmp_coords_jac_30__0 + jac_affine_0_2__0 * jac_affine_2_1__0 );
         const double jac_affine_inv_0_2__0 = tmp_coords_jac_33__0 * ( tmp_coords_jac_29__0 - tmp_coords_jac_31__0 );
         const double jac_affine_inv_1_0__0 =
             tmp_coords_jac_33__0 * ( jac_affine_1_2__0 * jac_affine_2_0__0 - jac_affine_1_0__0 * jac_affine_2_2__0 );
         const double jac_affine_inv_1_1__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_0__0 * jac_affine_2_2__0 - jac_affine_0_2__0 * jac_affine_2_0__0 );
         const double jac_affine_inv_1_2__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_2__0 * jac_affine_1_0__0 - jac_affine_0_0__0 * jac_affine_1_2__0 );
         const double jac_affine_inv_2_0__0 =
             tmp_coords_jac_33__0 * ( jac_affine_1_0__0 * jac_affine_2_1__0 - jac_affine_1_1__0 * jac_affine_2_0__0 );
         const double jac_affine_inv_2_1__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_1__0 * jac_affine_2_0__0 - jac_affine_0_0__0 * jac_affine_2_1__0 );
         const double jac_affine_inv_2_2__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_0__0 * jac_affine_1_1__0 - jac_affine_0_1__0 * jac_affine_1_0__0 );
         const double abs_det_jac_affine__0 = abs( tmp_coords_jac_32__0 );
         for ( int64_t ctr_2__0 = 0LL; ctr_2__0 < micro_edges_per_macro_edge; ctr_2__0 += 1LL )
         {
            for ( int64_t ctr_1__0 = 0LL; ctr_1__0 < -1LL * ctr_2__0 + micro_edges_per_macro_edge; ctr_1__0 += 1LL )
            {
               for ( int64_t ctr_0__0 = 0LL; ctr_0__0 < -1LL - ctr_1__0 - ctr_2__0 + micro_edges_per_macro_edge; ctr_0__0 += 1LL )
               {
                  const double p_affine_0_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__0;
                  const double p_affine_0_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__0;
                  const double p_affine_0_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__0;
                  const double p_affine_1_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__0;
                  const double p_affine_1_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__0;
                  const double p_affine_1_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__0;
                  const double p_affine_2_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_2_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_2_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__0 );
                  const double tmp_kernel_op_0__0 = -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0 - jac_affine_inv_2_0__0;
                  const double tmp_kernel_op_1__0 = -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0 - jac_affine_inv_2_1__0;
                  const double tmp_kernel_op_2__0 = -1.0 * jac_affine_inv_0_2__0 - jac_affine_inv_1_2__0 - jac_affine_inv_2_2__0;
                  const double tmp_kernel_op_3__0 = 0.16666666666666663 * abs_det_jac_affine__0;
                  const double tmp_kernel_op_4__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_0_0__0 * tmp_kernel_op_0__0 +
                                                                           jac_affine_inv_0_1__0 * tmp_kernel_op_1__0 +
                                                                           jac_affine_inv_0_2__0 * tmp_kernel_op_2__0 );
                  const double tmp_kernel_op_5__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_1_0__0 * tmp_kernel_op_0__0 +
                                                                           jac_affine_inv_1_1__0 * tmp_kernel_op_1__0 +
                                                                           jac_affine_inv_1_2__0 * tmp_kernel_op_2__0 );
                  const double tmp_kernel_op_6__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_2_0__0 * tmp_kernel_op_0__0 +
                                                                           jac_affine_inv_2_1__0 * tmp_kernel_op_1__0 +
                                                                           jac_affine_inv_2_2__0 * tmp_kernel_op_2__0 );
                  const double tmp_kernel_op_7__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_1_0__0 +
                                                                           jac_affine_inv_0_1__0 * jac_affine_inv_1_1__0 +
                                                                           jac_affine_inv_0_2__0 * jac_affine_inv_1_2__0 );
                  const double tmp_kernel_op_8__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_2_0__0 +
                                                                           jac_affine_inv_0_1__0 * jac_affine_inv_2_1__0 +
                                                                           jac_affine_inv_0_2__0 * jac_affine_inv_2_2__0 );
                  const double tmp_kernel_op_9__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_1_0__0 * jac_affine_inv_2_0__0 +
                                                                           jac_affine_inv_1_1__0 * jac_affine_inv_2_1__0 +
                                                                           jac_affine_inv_1_2__0 * jac_affine_inv_2_2__0 );
                  const double elMat_0_0__0 =
                      tmp_kernel_op_3__0 * ( tmp_kernel_op_0__0 * tmp_kernel_op_0__0 + tmp_kernel_op_1__0 * tmp_kernel_op_1__0 +
                                             tmp_kernel_op_2__0 * tmp_kernel_op_2__0 );
                  const double elMat_0_1__0 = tmp_kernel_op_4__0;
                  const double elMat_0_2__0 = tmp_kernel_op_5__0;
                  const double elMat_0_3__0 = tmp_kernel_op_6__0;
                  const double elMat_1_0__0 = tmp_kernel_op_4__0;
                  const double elMat_1_1__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_0_0__0 +
                                                                     jac_affine_inv_0_1__0 * jac_affine_inv_0_1__0 +
                                                                     jac_affine_inv_0_2__0 * jac_affine_inv_0_2__0 );
                  const double elMat_1_2__0 = tmp_kernel_op_7__0;
                  const double elMat_1_3__0 = tmp_kernel_op_8__0;
                  const double elMat_2_0__0 = tmp_kernel_op_5__0;
                  const double elMat_2_1__0 = tmp_kernel_op_7__0;
                  const double elMat_2_2__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_1_0__0 * jac_affine_inv_1_0__0 +
                                                                     jac_affine_inv_1_1__0 * jac_affine_inv_1_1__0 +
                                                                     jac_affine_inv_1_2__0 * jac_affine_inv_1_2__0 );
                  const double elMat_2_3__0 = tmp_kernel_op_9__0;
                  const double elMat_3_0__0 = tmp_kernel_op_6__0;
                  const double elMat_3_1__0 = tmp_kernel_op_8__0;
                  const double elMat_3_2__0 = tmp_kernel_op_9__0;
                  const double elMat_3_3__0 = tmp_kernel_op_3__0 * ( jac_affine_inv_2_0__0 * jac_affine_inv_2_0__0 +
                                                                     jac_affine_inv_2_1__0 * jac_affine_inv_2_1__0 +
                                                                     jac_affine_inv_2_2__0 * jac_affine_inv_2_2__0 );
                  /*  */
                  /* Apply basis transformation */
                  /*  */
                  std::vector< uint_t > rowIdx__5   = std::vector< uint_t >( 4LL );
                  std::vector< uint_t > colIdx__5   = std::vector< uint_t >( 4LL );
                  std::vector< real_t > matData__17 = std::vector< real_t >( 16LL );
                  rowIdx__5.operator[]( 0 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__5.operator[]( 1 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) -
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__0 ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__5.operator[]( 2 ) =
                      (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) -
                                         ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx__5.operator[]( 3 ) =
                      (uint_t) _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                         ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__5.operator[]( 0 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__5.operator[]( 1 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) -
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL + ctr_1__0 ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__5.operator[]( 2 ) =
                      (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) -
                                         ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx__5.operator[]( 3 ) =
                      (uint_t) _data_src[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                         ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                             ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                         ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                             ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  matData__17.operator[]( 0 )  = (real_t) elMat_0_0__0;
                  matData__17.operator[]( 1 )  = (real_t) elMat_0_1__0;
                  matData__17.operator[]( 2 )  = (real_t) elMat_0_2__0;
                  matData__17.operator[]( 3 )  = (real_t) elMat_0_3__0;
                  matData__17.operator[]( 4 )  = (real_t) elMat_1_0__0;
                  matData__17.operator[]( 5 )  = (real_t) elMat_1_1__0;
                  matData__17.operator[]( 6 )  = (real_t) elMat_1_2__0;
                  matData__17.operator[]( 7 )  = (real_t) elMat_1_3__0;
                  matData__17.operator[]( 8 )  = (real_t) elMat_2_0__0;
                  matData__17.operator[]( 9 )  = (real_t) elMat_2_1__0;
                  matData__17.operator[]( 10 ) = (real_t) elMat_2_2__0;
                  matData__17.operator[]( 11 ) = (real_t) elMat_2_3__0;
                  matData__17.operator[]( 12 ) = (real_t) elMat_3_0__0;
                  matData__17.operator[]( 13 ) = (real_t) elMat_3_1__0;
                  matData__17.operator[]( 14 ) = (real_t) elMat_3_2__0;
                  matData__17.operator[]( 15 ) = (real_t) elMat_3_3__0;
                  /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
                  (void) mat;
                  mat->addValues( rowIdx__5, colIdx__5, matData__17 );
               }
            }
         }
      }
      {
         /* CellType.GREEN_DOWN */
         const double tmp_coords_jac_0  = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1  = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2  = tmp_coords_jac_1 * 0.0;
         const double tmp_coords_jac_3  = tmp_coords_jac_0 * tmp_coords_jac_2;
         const double tmp_coords_jac_4  = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5  = tmp_coords_jac_1 * 1.0;
         const double tmp_coords_jac_6  = tmp_coords_jac_4 * tmp_coords_jac_5;
         const double tmp_coords_jac_7  = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_8  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_6 + tmp_coords_jac_2 * tmp_coords_jac_7;
         const double tmp_coords_jac_9  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10 = tmp_coords_jac_2 * tmp_coords_jac_9;
         const double tmp_coords_jac_11 = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12 = tmp_coords_jac_11 * tmp_coords_jac_5;
         const double tmp_coords_jac_13 = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_14 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12 + tmp_coords_jac_13 * tmp_coords_jac_2;
         const double tmp_coords_jac_15 = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16 = tmp_coords_jac_15 * tmp_coords_jac_2;
         const double tmp_coords_jac_17 = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18 = tmp_coords_jac_17 * tmp_coords_jac_5;
         const double tmp_coords_jac_19 = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_20 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18 + tmp_coords_jac_19 * tmp_coords_jac_2;
         const double tmp_coords_jac_21 = tmp_coords_jac_0 * tmp_coords_jac_5;
         const double tmp_coords_jac_22 = tmp_coords_jac_5 * tmp_coords_jac_9;
         const double tmp_coords_jac_23 = tmp_coords_jac_15 * tmp_coords_jac_5;
         const double tmp_coords_jac_24 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_5 * tmp_coords_jac_7;
         const double tmp_coords_jac_25 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_13 * tmp_coords_jac_5;
         const double tmp_coords_jac_26 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_19 * tmp_coords_jac_5;
         const double p_affine_const_0_0 = tmp_coords_jac_3 + tmp_coords_jac_8;
         const double p_affine_const_0_1 = tmp_coords_jac_10 + tmp_coords_jac_14;
         const double p_affine_const_0_2 = tmp_coords_jac_16 + tmp_coords_jac_20;
         const double p_affine_const_1_0 = tmp_coords_jac_21 + tmp_coords_jac_8;
         const double p_affine_const_1_1 = tmp_coords_jac_14 + tmp_coords_jac_22;
         const double p_affine_const_1_2 = tmp_coords_jac_20 + tmp_coords_jac_23;
         const double p_affine_const_2_0 = tmp_coords_jac_21 + tmp_coords_jac_24 + tmp_coords_jac_2 * tmp_coords_jac_4;
         const double p_affine_const_2_1 = tmp_coords_jac_22 + tmp_coords_jac_25 + tmp_coords_jac_11 * tmp_coords_jac_2;
         const double p_affine_const_2_2 = tmp_coords_jac_23 + tmp_coords_jac_26 + tmp_coords_jac_17 * tmp_coords_jac_2;
         const double p_affine_const_3_0 = tmp_coords_jac_24 + tmp_coords_jac_3 + tmp_coords_jac_6;
         const double p_affine_const_3_1 = tmp_coords_jac_10 + tmp_coords_jac_12 + tmp_coords_jac_25;
         const double p_affine_const_3_2 = tmp_coords_jac_16 + tmp_coords_jac_18 + tmp_coords_jac_26;
         const double jac_affine_0_0     = p_affine_const_1_0 - p_affine_const_0_0;
         const double jac_affine_0_1     = p_affine_const_2_0 - p_affine_const_0_0;
         const double jac_affine_0_2     = p_affine_const_3_0 - p_affine_const_0_0;
         const double jac_affine_1_0     = p_affine_const_1_1 - p_affine_const_0_1;
         const double jac_affine_1_1     = p_affine_const_2_1 - p_affine_const_0_1;
         const double tmp_coords_jac_31  = jac_affine_0_2 * jac_affine_1_1;
         const double jac_affine_1_2     = p_affine_const_3_1 - p_affine_const_0_1;
         const double tmp_coords_jac_29  = jac_affine_0_1 * jac_affine_1_2;
         const double jac_affine_2_0     = p_affine_const_1_2 - p_affine_const_0_2;
         const double jac_affine_2_1     = p_affine_const_2_2 - p_affine_const_0_2;
         const double tmp_coords_jac_28  = jac_affine_1_2 * jac_affine_2_1;
         const double jac_affine_2_2     = p_affine_const_3_2 - p_affine_const_0_2;
         const double tmp_coords_jac_27  = jac_affine_1_1 * jac_affine_2_2;
         const double tmp_coords_jac_30  = jac_affine_0_1 * jac_affine_2_2;
         const double tmp_coords_jac_32  = jac_affine_0_0 * tmp_coords_jac_27 + jac_affine_2_0 * tmp_coords_jac_29 -
                                          jac_affine_0_0 * tmp_coords_jac_28 - jac_affine_1_0 * tmp_coords_jac_30 -
                                          jac_affine_2_0 * tmp_coords_jac_31 + jac_affine_0_2 * jac_affine_1_0 * jac_affine_2_1;
         const double tmp_coords_jac_33  = 1.0 / tmp_coords_jac_32;
         const double jac_affine_inv_0_0 = tmp_coords_jac_33 * ( tmp_coords_jac_27 - tmp_coords_jac_28 );
         const double jac_affine_inv_0_1 = tmp_coords_jac_33 * ( -1.0 * tmp_coords_jac_30 + jac_affine_0_2 * jac_affine_2_1 );
         const double jac_affine_inv_0_2 = tmp_coords_jac_33 * ( tmp_coords_jac_29 - tmp_coords_jac_31 );
         const double jac_affine_inv_1_0 =
             tmp_coords_jac_33 * ( jac_affine_1_2 * jac_affine_2_0 - jac_affine_1_0 * jac_affine_2_2 );
         const double jac_affine_inv_1_1 =
             tmp_coords_jac_33 * ( jac_affine_0_0 * jac_affine_2_2 - jac_affine_0_2 * jac_affine_2_0 );
         const double jac_affine_inv_1_2 =
             tmp_coords_jac_33 * ( jac_affine_0_2 * jac_affine_1_0 - jac_affine_0_0 * jac_affine_1_2 );
         const double jac_affine_inv_2_0 =
             tmp_coords_jac_33 * ( jac_affine_1_0 * jac_affine_2_1 - jac_affine_1_1 * jac_affine_2_0 );
         const double jac_affine_inv_2_1 =
             tmp_coords_jac_33 * ( jac_affine_0_1 * jac_affine_2_0 - jac_affine_0_0 * jac_affine_2_1 );
         const double jac_affine_inv_2_2 =
             tmp_coords_jac_33 * ( jac_affine_0_0 * jac_affine_1_1 - jac_affine_0_1 * jac_affine_1_0 );
         const double abs_det_jac_affine = abs( tmp_coords_jac_32 );
         for ( int64_t ctr_2 = 0LL; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1LL )
         {
            for ( int64_t ctr_1 = 0LL; ctr_1 < -1LL * ctr_2 + micro_edges_per_macro_edge; ctr_1 += 1LL )
            {
               for ( int64_t ctr_0 = 0LL; ctr_0 < -1LL - ctr_1 - ctr_2 + micro_edges_per_macro_edge; ctr_0 += 1LL )
               {
                  const double p_affine_0_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2;
                  const double p_affine_0_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2;
                  const double p_affine_0_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2;
                  const double p_affine_1_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2;
                  const double p_affine_1_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2;
                  const double p_affine_1_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2;
                  const double p_affine_2_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_2_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_2_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2 );
                  const double tmp_kernel_op_0 = -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0 - jac_affine_inv_2_0;
                  const double tmp_kernel_op_1 = -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1 - jac_affine_inv_2_1;
                  const double tmp_kernel_op_2 = -1.0 * jac_affine_inv_0_2 - jac_affine_inv_1_2 - jac_affine_inv_2_2;
                  const double tmp_kernel_op_3 = 0.16666666666666663 * abs_det_jac_affine;
                  const double tmp_kernel_op_4 =
                      tmp_kernel_op_3 * ( jac_affine_inv_0_0 * tmp_kernel_op_0 + jac_affine_inv_0_1 * tmp_kernel_op_1 +
                                          jac_affine_inv_0_2 * tmp_kernel_op_2 );
                  const double tmp_kernel_op_5 =
                      tmp_kernel_op_3 * ( jac_affine_inv_1_0 * tmp_kernel_op_0 + jac_affine_inv_1_1 * tmp_kernel_op_1 +
                                          jac_affine_inv_1_2 * tmp_kernel_op_2 );
                  const double tmp_kernel_op_6 =
                      tmp_kernel_op_3 * ( jac_affine_inv_2_0 * tmp_kernel_op_0 + jac_affine_inv_2_1 * tmp_kernel_op_1 +
                                          jac_affine_inv_2_2 * tmp_kernel_op_2 );
                  const double tmp_kernel_op_7 =
                      tmp_kernel_op_3 * ( jac_affine_inv_0_0 * jac_affine_inv_1_0 + jac_affine_inv_0_1 * jac_affine_inv_1_1 +
                                          jac_affine_inv_0_2 * jac_affine_inv_1_2 );
                  const double tmp_kernel_op_8 =
                      tmp_kernel_op_3 * ( jac_affine_inv_0_0 * jac_affine_inv_2_0 + jac_affine_inv_0_1 * jac_affine_inv_2_1 +
                                          jac_affine_inv_0_2 * jac_affine_inv_2_2 );
                  const double tmp_kernel_op_9 =
                      tmp_kernel_op_3 * ( jac_affine_inv_1_0 * jac_affine_inv_2_0 + jac_affine_inv_1_1 * jac_affine_inv_2_1 +
                                          jac_affine_inv_1_2 * jac_affine_inv_2_2 );
                  const double elMat_0_0 =
                      tmp_kernel_op_3 * ( tmp_kernel_op_0 * tmp_kernel_op_0 + tmp_kernel_op_1 * tmp_kernel_op_1 +
                                          tmp_kernel_op_2 * tmp_kernel_op_2 );
                  const double elMat_0_1 = tmp_kernel_op_4;
                  const double elMat_0_2 = tmp_kernel_op_5;
                  const double elMat_0_3 = tmp_kernel_op_6;
                  const double elMat_1_0 = tmp_kernel_op_4;
                  const double elMat_1_1 =
                      tmp_kernel_op_3 * ( jac_affine_inv_0_0 * jac_affine_inv_0_0 + jac_affine_inv_0_1 * jac_affine_inv_0_1 +
                                          jac_affine_inv_0_2 * jac_affine_inv_0_2 );
                  const double elMat_1_2 = tmp_kernel_op_7;
                  const double elMat_1_3 = tmp_kernel_op_8;
                  const double elMat_2_0 = tmp_kernel_op_5;
                  const double elMat_2_1 = tmp_kernel_op_7;
                  const double elMat_2_2 =
                      tmp_kernel_op_3 * ( jac_affine_inv_1_0 * jac_affine_inv_1_0 + jac_affine_inv_1_1 * jac_affine_inv_1_1 +
                                          jac_affine_inv_1_2 * jac_affine_inv_1_2 );
                  const double elMat_2_3 = tmp_kernel_op_9;
                  const double elMat_3_0 = tmp_kernel_op_6;
                  const double elMat_3_1 = tmp_kernel_op_8;
                  const double elMat_3_2 = tmp_kernel_op_9;
                  const double elMat_3_3 =
                      tmp_kernel_op_3 * ( jac_affine_inv_2_0 * jac_affine_inv_2_0 + jac_affine_inv_2_1 * jac_affine_inv_2_1 +
                                          jac_affine_inv_2_2 * jac_affine_inv_2_2 );
                  /*  */
                  /* Apply basis transformation */
                  /*  */
                  std::vector< uint_t > rowIdx  = std::vector< uint_t >( 4LL );
                  std::vector< uint_t > colIdx  = std::vector< uint_t >( 4LL );
                  std::vector< real_t > matData = std::vector< real_t >( 16LL );
                  rowIdx.operator[]( 0 )        = (uint_t)
                      _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx.operator[]( 1 ) = (uint_t)
                      _data_dst[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx.operator[]( 2 ) = (uint_t)
                      _data_dst[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ctr_1 + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  rowIdx.operator[]( 3 ) = (uint_t)
                      _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx.operator[]( 0 ) = (uint_t)
                      _data_src[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx.operator[]( 1 ) = (uint_t)
                      _data_src[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL -
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx.operator[]( 2 ) = (uint_t)
                      _data_src[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ctr_1 + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  colIdx.operator[]( 3 ) = (uint_t)
                      _data_src[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                    ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                ( 1LL + ctr_1 ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                    ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  matData.operator[]( 0 )  = (real_t) elMat_0_0;
                  matData.operator[]( 1 )  = (real_t) elMat_0_1;
                  matData.operator[]( 2 )  = (real_t) elMat_0_2;
                  matData.operator[]( 3 )  = (real_t) elMat_0_3;
                  matData.operator[]( 4 )  = (real_t) elMat_1_0;
                  matData.operator[]( 5 )  = (real_t) elMat_1_1;
                  matData.operator[]( 6 )  = (real_t) elMat_1_2;
                  matData.operator[]( 7 )  = (real_t) elMat_1_3;
                  matData.operator[]( 8 )  = (real_t) elMat_2_0;
                  matData.operator[]( 9 )  = (real_t) elMat_2_1;
                  matData.operator[]( 10 ) = (real_t) elMat_2_2;
                  matData.operator[]( 11 ) = (real_t) elMat_2_3;
                  matData.operator[]( 12 ) = (real_t) elMat_3_0;
                  matData.operator[]( 13 ) = (real_t) elMat_3_1;
                  matData.operator[]( 14 ) = (real_t) elMat_3_2;
                  matData.operator[]( 15 ) = (real_t) elMat_3_3;
                  /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
                  (void) mat;
                  mat->addValues( rowIdx, colIdx, matData );
               }
            }
         }
      }
   }
}
void P1ElementwiseDiffusion::toMatrix_P1ElementwiseDiffusion_macro_2D( int64_t* RESTRICT const _data_dst,
                                                                       int64_t* RESTRICT const _data_src,
                                                                       const double            macro_vertex_coord_id_0comp0,
                                                                       const double            macro_vertex_coord_id_0comp1,
                                                                       const double            macro_vertex_coord_id_1comp0,
                                                                       const double            macro_vertex_coord_id_1comp1,
                                                                       const double            macro_vertex_coord_id_2comp0,
                                                                       const double            macro_vertex_coord_id_2comp1,
                                                                       const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                       const int64_t micro_edges_per_macro_edge,
                                                                       const double  micro_edges_per_macro_edge_float ) const
{
   {
      {
         /* FaceType.GRAY */
         const double tmp_coords_jac_0__0   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__0   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__0   = tmp_coords_jac_1__0 * 0.0;
         const double tmp_coords_jac_3__0   = tmp_coords_jac_0__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_4__0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__0   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__0 * tmp_coords_jac_4__0;
         const double tmp_coords_jac_6__0   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_7__0   = tmp_coords_jac_2__0 * tmp_coords_jac_6__0;
         const double tmp_coords_jac_8__0   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_9__0   = macro_vertex_coord_id_0comp1 + tmp_coords_jac_2__0 * tmp_coords_jac_8__0;
         const double tmp_coords_jac_10__0  = tmp_coords_jac_1__0 * 1.0;
         const double p_affine_const_0_0__0 = tmp_coords_jac_3__0 + tmp_coords_jac_5__0;
         const double p_affine_const_0_1__0 = tmp_coords_jac_7__0 + tmp_coords_jac_9__0;
         const double p_affine_const_1_0__0 = tmp_coords_jac_5__0 + tmp_coords_jac_0__0 * tmp_coords_jac_10__0;
         const double p_affine_const_1_1__0 = tmp_coords_jac_9__0 + tmp_coords_jac_10__0 * tmp_coords_jac_6__0;
         const double p_affine_const_2_0__0 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__0 + tmp_coords_jac_10__0 * tmp_coords_jac_4__0;
         const double p_affine_const_2_1__0 =
             macro_vertex_coord_id_0comp1 + tmp_coords_jac_7__0 + tmp_coords_jac_10__0 * tmp_coords_jac_8__0;
         const double jac_affine_0_0__0     = p_affine_const_1_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_1__0     = p_affine_const_2_0__0 - p_affine_const_0_0__0;
         const double jac_affine_1_0__0     = p_affine_const_1_1__0 - p_affine_const_0_1__0;
         const double jac_affine_1_1__0     = p_affine_const_2_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_11__0  = jac_affine_0_0__0 * jac_affine_1_1__0 - jac_affine_0_1__0 * jac_affine_1_0__0;
         const double tmp_coords_jac_12__0  = 1.0 / tmp_coords_jac_11__0;
         const double jac_affine_inv_0_0__0 = jac_affine_1_1__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_0_1__0 = -1.0 * jac_affine_0_1__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_1_0__0 = -1.0 * jac_affine_1_0__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_1_1__0 = jac_affine_0_0__0 * tmp_coords_jac_12__0;
         const double abs_det_jac_affine__0 = abs( tmp_coords_jac_11__0 );
         for ( int64_t ctr_1__0 = 0LL; ctr_1__0 < micro_edges_per_macro_edge; ctr_1__0 += 1LL )
         {
            for ( int64_t ctr_0__0 = 0LL; ctr_0__0 < -1LL * ctr_1__0 + micro_edges_per_macro_edge; ctr_0__0 += 1LL )
            {
               const double p_affine_0_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0;
               const double p_affine_0_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0;
               const double p_affine_1_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0;
               const double p_affine_1_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0;
               const double p_affine_2_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__0 );
               const double p_affine_2_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__0 );
               const double tmp_kernel_op_0__0 = -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0;
               const double tmp_kernel_op_1__0 = -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0;
               const double tmp_kernel_op_2__0 = 0.5 * abs_det_jac_affine__0;
               const double tmp_kernel_op_3__0 = tmp_kernel_op_2__0 * ( jac_affine_inv_0_0__0 * tmp_kernel_op_0__0 +
                                                                        jac_affine_inv_0_1__0 * tmp_kernel_op_1__0 );
               const double tmp_kernel_op_4__0 = tmp_kernel_op_2__0 * ( jac_affine_inv_1_0__0 * tmp_kernel_op_0__0 +
                                                                        jac_affine_inv_1_1__0 * tmp_kernel_op_1__0 );
               const double tmp_kernel_op_5__0 = tmp_kernel_op_2__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_1_0__0 +
                                                                        jac_affine_inv_0_1__0 * jac_affine_inv_1_1__0 );
               const double elMat_0_0__0 =
                   tmp_kernel_op_2__0 * ( tmp_kernel_op_0__0 * tmp_kernel_op_0__0 + tmp_kernel_op_1__0 * tmp_kernel_op_1__0 );
               const double elMat_0_1__0 = tmp_kernel_op_3__0;
               const double elMat_0_2__0 = tmp_kernel_op_4__0;
               const double elMat_1_0__0 = tmp_kernel_op_3__0;
               const double elMat_1_1__0 = tmp_kernel_op_2__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_0_0__0 +
                                                                  jac_affine_inv_0_1__0 * jac_affine_inv_0_1__0 );
               const double elMat_1_2__0 = tmp_kernel_op_5__0;
               const double elMat_2_0__0 = tmp_kernel_op_4__0;
               const double elMat_2_1__0 = tmp_kernel_op_5__0;
               const double elMat_2_2__0 = tmp_kernel_op_2__0 * ( jac_affine_inv_1_0__0 * jac_affine_inv_1_0__0 +
                                                                  jac_affine_inv_1_1__0 * jac_affine_inv_1_1__0 );
               /*  */
               /* Apply basis transformation */
               /*  */
               std::vector< uint_t > rowIdx__4   = std::vector< uint_t >( 3LL );
               std::vector< uint_t > colIdx__4   = std::vector< uint_t >( 3LL );
               std::vector< real_t > matData__10 = std::vector< real_t >( 9LL );
               rowIdx__4.operator[]( 0 )         = (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) +
                                                              ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               rowIdx__4.operator[]( 1 )         = (uint_t) _data_dst[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL +
                                                              ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               rowIdx__4.operator[]( 2 ) =
                   (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) +
                                      ( 1LL + ctr_1__0 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0__0];
               colIdx__4.operator[]( 0 ) = (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) +
                                                              ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               colIdx__4.operator[]( 1 ) = (uint_t) _data_src[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL +
                                                              ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               colIdx__4.operator[]( 2 ) =
                   (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) +
                                      ( 1LL + ctr_1__0 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0__0];
               matData__10.operator[]( 0 ) = (real_t) elMat_0_0__0;
               matData__10.operator[]( 1 ) = (real_t) elMat_0_1__0;
               matData__10.operator[]( 2 ) = (real_t) elMat_0_2__0;
               matData__10.operator[]( 3 ) = (real_t) elMat_1_0__0;
               matData__10.operator[]( 4 ) = (real_t) elMat_1_1__0;
               matData__10.operator[]( 5 ) = (real_t) elMat_1_2__0;
               matData__10.operator[]( 6 ) = (real_t) elMat_2_0__0;
               matData__10.operator[]( 7 ) = (real_t) elMat_2_1__0;
               matData__10.operator[]( 8 ) = (real_t) elMat_2_2__0;
               /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
               (void) mat;
               mat->addValues( rowIdx__4, colIdx__4, matData__10 );
            }
         }
      }
      {
         /* FaceType.BLUE */
         const double tmp_coords_jac_0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2   = tmp_coords_jac_1 * 0.0;
         const double tmp_coords_jac_3   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4   = tmp_coords_jac_1 * 1.0;
         const double tmp_coords_jac_5   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3 * tmp_coords_jac_4;
         const double tmp_coords_jac_6   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_7   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_8   = macro_vertex_coord_id_0comp1 + tmp_coords_jac_4 * tmp_coords_jac_7;
         const double tmp_coords_jac_9   = tmp_coords_jac_0 * tmp_coords_jac_4;
         const double tmp_coords_jac_10  = tmp_coords_jac_4 * tmp_coords_jac_6;
         const double p_affine_const_0_0 = tmp_coords_jac_5 + tmp_coords_jac_0 * tmp_coords_jac_2;
         const double p_affine_const_0_1 = tmp_coords_jac_8 + tmp_coords_jac_2 * tmp_coords_jac_6;
         const double p_affine_const_1_0 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_9 + tmp_coords_jac_2 * tmp_coords_jac_3;
         const double p_affine_const_1_1 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10 + tmp_coords_jac_2 * tmp_coords_jac_7;
         const double p_affine_const_2_0 = tmp_coords_jac_5 + tmp_coords_jac_9;
         const double p_affine_const_2_1 = tmp_coords_jac_10 + tmp_coords_jac_8;
         const double jac_affine_0_0     = p_affine_const_1_0 - p_affine_const_0_0;
         const double jac_affine_0_1     = p_affine_const_2_0 - p_affine_const_0_0;
         const double jac_affine_1_0     = p_affine_const_1_1 - p_affine_const_0_1;
         const double jac_affine_1_1     = p_affine_const_2_1 - p_affine_const_0_1;
         const double tmp_coords_jac_11  = jac_affine_0_0 * jac_affine_1_1 - jac_affine_0_1 * jac_affine_1_0;
         const double tmp_coords_jac_12  = 1.0 / tmp_coords_jac_11;
         const double jac_affine_inv_0_0 = jac_affine_1_1 * tmp_coords_jac_12;
         const double jac_affine_inv_0_1 = -1.0 * jac_affine_0_1 * tmp_coords_jac_12;
         const double jac_affine_inv_1_0 = -1.0 * jac_affine_1_0 * tmp_coords_jac_12;
         const double jac_affine_inv_1_1 = jac_affine_0_0 * tmp_coords_jac_12;
         const double abs_det_jac_affine = abs( tmp_coords_jac_11 );
         for ( int64_t ctr_1 = 0LL; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1LL )
         {
            for ( int64_t ctr_0 = 0LL; ctr_0 < -1LL - ctr_1 + micro_edges_per_macro_edge; ctr_0 += 1LL )
            {
               const double p_affine_0_0 = macro_vertex_coord_id_0comp0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) *
                                               (double) ( 1LL + ctr_0 ) +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1;
               const double p_affine_0_1 = macro_vertex_coord_id_0comp1 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) *
                                               (double) ( 1LL + ctr_0 ) +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1;
               const double p_affine_1_0 = macro_vertex_coord_id_0comp0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) *
                                               (double) ( 1LL + ctr_1 );
               const double p_affine_1_1 = macro_vertex_coord_id_0comp1 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) *
                                               (double) ( 1LL + ctr_1 );
               const double p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 );
               const double p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 );
               const double tmp_kernel_op_0 = -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0;
               const double tmp_kernel_op_1 = -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1;
               const double tmp_kernel_op_2 = 0.5 * abs_det_jac_affine;
               const double tmp_kernel_op_3 =
                   tmp_kernel_op_2 * ( jac_affine_inv_0_0 * tmp_kernel_op_0 + jac_affine_inv_0_1 * tmp_kernel_op_1 );
               const double tmp_kernel_op_4 =
                   tmp_kernel_op_2 * ( jac_affine_inv_1_0 * tmp_kernel_op_0 + jac_affine_inv_1_1 * tmp_kernel_op_1 );
               const double tmp_kernel_op_5 =
                   tmp_kernel_op_2 * ( jac_affine_inv_0_0 * jac_affine_inv_1_0 + jac_affine_inv_0_1 * jac_affine_inv_1_1 );
               const double elMat_0_0 =
                   tmp_kernel_op_2 * ( tmp_kernel_op_0 * tmp_kernel_op_0 + tmp_kernel_op_1 * tmp_kernel_op_1 );
               const double elMat_0_1 = tmp_kernel_op_3;
               const double elMat_0_2 = tmp_kernel_op_4;
               const double elMat_1_0 = tmp_kernel_op_3;
               const double elMat_1_1 =
                   tmp_kernel_op_2 * ( jac_affine_inv_0_0 * jac_affine_inv_0_0 + jac_affine_inv_0_1 * jac_affine_inv_0_1 );
               const double elMat_1_2 = tmp_kernel_op_5;
               const double elMat_2_0 = tmp_kernel_op_4;
               const double elMat_2_1 = tmp_kernel_op_5;
               const double elMat_2_2 =
                   tmp_kernel_op_2 * ( jac_affine_inv_1_0 * jac_affine_inv_1_0 + jac_affine_inv_1_1 * jac_affine_inv_1_1 );
               /*  */
               /* Apply basis transformation */
               /*  */
               std::vector< uint_t > rowIdx  = std::vector< uint_t >( 3LL );
               std::vector< uint_t > colIdx  = std::vector< uint_t >( 3LL );
               std::vector< real_t > matData = std::vector< real_t >( 9LL );
               rowIdx.operator[]( 0 ) =
                   (uint_t) _data_dst[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1 + ctr_0];
               rowIdx.operator[]( 1 ) = (uint_t) _data_dst[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) +
                                                           ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               rowIdx.operator[]( 2 ) = (uint_t) _data_dst[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL +
                                                           ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               colIdx.operator[]( 0 ) =
                   (uint_t) _data_src[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1 + ctr_0];
               colIdx.operator[]( 1 )  = (uint_t) _data_src[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) +
                                                           ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               colIdx.operator[]( 2 )  = (uint_t) _data_src[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL +
                                                           ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               matData.operator[]( 0 ) = (real_t) elMat_0_0;
               matData.operator[]( 1 ) = (real_t) elMat_0_1;
               matData.operator[]( 2 ) = (real_t) elMat_0_2;
               matData.operator[]( 3 ) = (real_t) elMat_1_0;
               matData.operator[]( 4 ) = (real_t) elMat_1_1;
               matData.operator[]( 5 ) = (real_t) elMat_1_2;
               matData.operator[]( 6 ) = (real_t) elMat_2_0;
               matData.operator[]( 7 ) = (real_t) elMat_2_1;
               matData.operator[]( 8 ) = (real_t) elMat_2_2;
               /* Artifact from code generation: workaround to add `mat` to the kernel parameter list. */
               (void) mat;
               mat->addValues( rowIdx, colIdx, matData );
            }
         }
      }
   }
}
void P1ElementwiseDiffusion::computeInverseDiagonalOperatorValues_P1ElementwiseDiffusion_macro_3D(
    double* RESTRICT const _data_invDiag_,
    const double           macro_vertex_coord_id_0comp0,
    const double           macro_vertex_coord_id_0comp1,
    const double           macro_vertex_coord_id_0comp2,
    const double           macro_vertex_coord_id_1comp0,
    const double           macro_vertex_coord_id_1comp1,
    const double           macro_vertex_coord_id_1comp2,
    const double           macro_vertex_coord_id_2comp0,
    const double           macro_vertex_coord_id_2comp1,
    const double           macro_vertex_coord_id_2comp2,
    const double           macro_vertex_coord_id_3comp0,
    const double           macro_vertex_coord_id_3comp1,
    const double           macro_vertex_coord_id_3comp2,
    const int64_t          micro_edges_per_macro_edge,
    const double           micro_edges_per_macro_edge_float ) const
{
   {
      {
         /* CellType.WHITE_UP */
         const double tmp_coords_jac_0__4   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__4   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__4   = tmp_coords_jac_1__4 * 0.0;
         const double tmp_coords_jac_3__4   = tmp_coords_jac_0__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_4__4   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__4   = tmp_coords_jac_2__4 * tmp_coords_jac_4__4;
         const double tmp_coords_jac_6__4   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__4   = tmp_coords_jac_2__4 * tmp_coords_jac_6__4;
         const double tmp_coords_jac_8__4   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_5__4 + tmp_coords_jac_7__4;
         const double tmp_coords_jac_9__4   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__4  = tmp_coords_jac_2__4 * tmp_coords_jac_9__4;
         const double tmp_coords_jac_11__4  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12__4  = tmp_coords_jac_11__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_13__4  = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_14__4  = tmp_coords_jac_13__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_15__4  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12__4 + tmp_coords_jac_14__4;
         const double tmp_coords_jac_16__4  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_17__4  = tmp_coords_jac_16__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_18__4  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__4  = tmp_coords_jac_18__4 * tmp_coords_jac_2__4;
         const double tmp_coords_jac_20__4  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_21__4  = tmp_coords_jac_2__4 * tmp_coords_jac_20__4;
         const double tmp_coords_jac_22__4  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_19__4 + tmp_coords_jac_21__4;
         const double tmp_coords_jac_23__4  = tmp_coords_jac_1__4 * 1.0;
         const double tmp_coords_jac_24__4  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__4;
         const double tmp_coords_jac_25__4  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10__4;
         const double tmp_coords_jac_26__4  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_17__4;
         const double p_affine_const_0_0__4 = tmp_coords_jac_3__4 + tmp_coords_jac_8__4;
         const double p_affine_const_0_1__4 = tmp_coords_jac_10__4 + tmp_coords_jac_15__4;
         const double p_affine_const_0_2__4 = tmp_coords_jac_17__4 + tmp_coords_jac_22__4;
         const double p_affine_const_1_0__4 = tmp_coords_jac_8__4 + tmp_coords_jac_0__4 * tmp_coords_jac_23__4;
         const double p_affine_const_1_1__4 = tmp_coords_jac_15__4 + tmp_coords_jac_23__4 * tmp_coords_jac_9__4;
         const double p_affine_const_1_2__4 = tmp_coords_jac_22__4 + tmp_coords_jac_16__4 * tmp_coords_jac_23__4;
         const double p_affine_const_2_0__4 =
             tmp_coords_jac_24__4 + tmp_coords_jac_7__4 + tmp_coords_jac_23__4 * tmp_coords_jac_4__4;
         const double p_affine_const_2_1__4 =
             tmp_coords_jac_14__4 + tmp_coords_jac_25__4 + tmp_coords_jac_11__4 * tmp_coords_jac_23__4;
         const double p_affine_const_2_2__4 =
             tmp_coords_jac_21__4 + tmp_coords_jac_26__4 + tmp_coords_jac_18__4 * tmp_coords_jac_23__4;
         const double p_affine_const_3_0__4 =
             tmp_coords_jac_24__4 + tmp_coords_jac_5__4 + tmp_coords_jac_23__4 * tmp_coords_jac_6__4;
         const double p_affine_const_3_1__4 =
             tmp_coords_jac_12__4 + tmp_coords_jac_25__4 + tmp_coords_jac_13__4 * tmp_coords_jac_23__4;
         const double p_affine_const_3_2__4 =
             tmp_coords_jac_19__4 + tmp_coords_jac_26__4 + tmp_coords_jac_20__4 * tmp_coords_jac_23__4;
         const double jac_affine_0_0__4    = p_affine_const_1_0__4 - p_affine_const_0_0__4;
         const double jac_affine_0_1__4    = p_affine_const_2_0__4 - p_affine_const_0_0__4;
         const double jac_affine_0_2__4    = p_affine_const_3_0__4 - p_affine_const_0_0__4;
         const double jac_affine_1_0__4    = p_affine_const_1_1__4 - p_affine_const_0_1__4;
         const double jac_affine_1_1__4    = p_affine_const_2_1__4 - p_affine_const_0_1__4;
         const double tmp_coords_jac_31__2 = jac_affine_0_2__4 * jac_affine_1_1__4;
         const double jac_affine_1_2__4    = p_affine_const_3_1__4 - p_affine_const_0_1__4;
         const double tmp_coords_jac_29__4 = jac_affine_0_1__4 * jac_affine_1_2__4;
         const double jac_affine_2_0__4    = p_affine_const_1_2__4 - p_affine_const_0_2__4;
         const double jac_affine_2_1__4    = p_affine_const_2_2__4 - p_affine_const_0_2__4;
         const double tmp_coords_jac_28__4 = jac_affine_1_2__4 * jac_affine_2_1__4;
         const double jac_affine_2_2__4    = p_affine_const_3_2__4 - p_affine_const_0_2__4;
         const double tmp_coords_jac_27__4 = jac_affine_1_1__4 * jac_affine_2_2__4;
         const double tmp_coords_jac_30__4 = jac_affine_0_1__4 * jac_affine_2_2__4;
         const double tmp_coords_jac_32__2 = jac_affine_0_0__4 * tmp_coords_jac_27__4 + jac_affine_2_0__4 * tmp_coords_jac_29__4 -
                                             jac_affine_0_0__4 * tmp_coords_jac_28__4 - jac_affine_1_0__4 * tmp_coords_jac_30__4 -
                                             jac_affine_2_0__4 * tmp_coords_jac_31__2 +
                                             jac_affine_0_2__4 * jac_affine_1_0__4 * jac_affine_2_1__4;
         const double tmp_coords_jac_33__2  = 1.0 / tmp_coords_jac_32__2;
         const double jac_affine_inv_0_0__4 = tmp_coords_jac_33__2 * ( tmp_coords_jac_27__4 - tmp_coords_jac_28__4 );
         const double jac_affine_inv_0_1__4 =
             tmp_coords_jac_33__2 * ( -1.0 * tmp_coords_jac_30__4 + jac_affine_0_2__4 * jac_affine_2_1__4 );
         const double jac_affine_inv_0_2__4 = tmp_coords_jac_33__2 * ( tmp_coords_jac_29__4 - tmp_coords_jac_31__2 );
         const double jac_affine_inv_1_0__4 =
             tmp_coords_jac_33__2 * ( jac_affine_1_2__4 * jac_affine_2_0__4 - jac_affine_1_0__4 * jac_affine_2_2__4 );
         const double jac_affine_inv_1_1__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_0__4 * jac_affine_2_2__4 - jac_affine_0_2__4 * jac_affine_2_0__4 );
         const double jac_affine_inv_1_2__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_2__4 * jac_affine_1_0__4 - jac_affine_0_0__4 * jac_affine_1_2__4 );
         const double jac_affine_inv_2_0__4 =
             tmp_coords_jac_33__2 * ( jac_affine_1_0__4 * jac_affine_2_1__4 - jac_affine_1_1__4 * jac_affine_2_0__4 );
         const double jac_affine_inv_2_1__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_1__4 * jac_affine_2_0__4 - jac_affine_0_0__4 * jac_affine_2_1__4 );
         const double jac_affine_inv_2_2__4 =
             tmp_coords_jac_33__2 * ( jac_affine_0_0__4 * jac_affine_1_1__4 - jac_affine_0_1__4 * jac_affine_1_0__4 );
         const double abs_det_jac_affine__4 = abs( tmp_coords_jac_32__2 );
         for ( int64_t ctr_2__4 = 0LL; ctr_2__4 < micro_edges_per_macro_edge; ctr_2__4 += 1LL )
         {
            for ( int64_t ctr_1__4 = 0LL; ctr_1__4 < -1LL * ctr_2__4 + micro_edges_per_macro_edge; ctr_1__4 += 1LL )
            {
               for ( int64_t ctr_0__4 = 0LL; ctr_0__4 < -1LL * ctr_1__4 - ctr_2__4 + micro_edges_per_macro_edge; ctr_0__4 += 1LL )
               {
                  const double p_affine_0_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_0_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_0_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_1_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_1_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_1_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_2_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__4;
                  const double p_affine_2_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__4;
                  const double p_affine_2_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__4 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__4;
                  const double p_affine_3_0__4 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__4 );
                  const double p_affine_3_1__4 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__4 );
                  const double p_affine_3_2__4 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__4 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__4 );
                  const double tmp_kernel_op_0__4 = 0.16666666666666663 * abs_det_jac_affine__4;
                  const double elMatDiag_0__4 =
                      tmp_kernel_op_0__4 *
                      ( ( -1.0 * jac_affine_inv_0_0__4 - jac_affine_inv_1_0__4 - jac_affine_inv_2_0__4 ) *
                            ( -1.0 * jac_affine_inv_0_0__4 - jac_affine_inv_1_0__4 - jac_affine_inv_2_0__4 ) +
                        ( -1.0 * jac_affine_inv_0_1__4 - jac_affine_inv_1_1__4 - jac_affine_inv_2_1__4 ) *
                            ( -1.0 * jac_affine_inv_0_1__4 - jac_affine_inv_1_1__4 - jac_affine_inv_2_1__4 ) +
                        ( -1.0 * jac_affine_inv_0_2__4 - jac_affine_inv_1_2__4 - jac_affine_inv_2_2__4 ) *
                            ( -1.0 * jac_affine_inv_0_2__4 - jac_affine_inv_1_2__4 - jac_affine_inv_2_2__4 ) );
                  const double elMatDiag_1__4 = tmp_kernel_op_0__4 * ( jac_affine_inv_0_0__4 * jac_affine_inv_0_0__4 +
                                                                       jac_affine_inv_0_1__4 * jac_affine_inv_0_1__4 +
                                                                       jac_affine_inv_0_2__4 * jac_affine_inv_0_2__4 );
                  const double elMatDiag_2__4 = tmp_kernel_op_0__4 * ( jac_affine_inv_1_0__4 * jac_affine_inv_1_0__4 +
                                                                       jac_affine_inv_1_1__4 * jac_affine_inv_1_1__4 +
                                                                       jac_affine_inv_1_2__4 * jac_affine_inv_1_2__4 );
                  const double elMatDiag_3__4 = tmp_kernel_op_0__4 * ( jac_affine_inv_2_0__4 * jac_affine_inv_2_0__4 +
                                                                       jac_affine_inv_2_1__4 * jac_affine_inv_2_1__4 +
                                                                       jac_affine_inv_2_2__4 * jac_affine_inv_2_2__4 );
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                 ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_0__4 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                     ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL -
                                 ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_1__4 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL -
                                     ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__4 ) * ( 2LL + ctr_1__4 ) / 2LL ) -
                                 ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__4 ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) + ctr_0__4 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_2__4 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__4 ) * ( 2LL + ctr_1__4 ) / 2LL ) -
                                     ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__4 ) * ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) + ctr_0__4 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                 ( -1LL * ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_3__4 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__4 ) * ctr_1__4 / 2LL ) -
                                     ( -1LL * ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__4 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL - ctr_2__4 + micro_edges_per_macro_edge ) * ctr_1__4 + ctr_0__4 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.WHITE_DOWN */
         const double tmp_coords_jac_0__3   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__3   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__3   = tmp_coords_jac_1__3 * 0.0;
         const double tmp_coords_jac_3__3   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4__3   = tmp_coords_jac_1__3 * 1.0;
         const double tmp_coords_jac_5__3   = tmp_coords_jac_3__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_6__3   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__3   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_4__3 * tmp_coords_jac_6__3;
         const double tmp_coords_jac_8__3   = tmp_coords_jac_5__3 + tmp_coords_jac_7__3;
         const double tmp_coords_jac_9__3   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__3  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__3  = tmp_coords_jac_10__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_12__3  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__3  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_14__3  = tmp_coords_jac_11__3 + tmp_coords_jac_13__3;
         const double tmp_coords_jac_15__3  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__3  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_17__3  = tmp_coords_jac_16__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_18__3  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__3  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_20__3  = tmp_coords_jac_17__3 + tmp_coords_jac_19__3;
         const double tmp_coords_jac_21__3  = tmp_coords_jac_0__3 * tmp_coords_jac_4__3;
         const double tmp_coords_jac_22__3  = tmp_coords_jac_4__3 * tmp_coords_jac_9__3;
         const double tmp_coords_jac_23__3  = tmp_coords_jac_15__3 * tmp_coords_jac_4__3;
         const double p_affine_const_0_0__3 = tmp_coords_jac_8__3 + tmp_coords_jac_0__3 * tmp_coords_jac_2__3;
         const double p_affine_const_0_1__3 = tmp_coords_jac_14__3 + tmp_coords_jac_2__3 * tmp_coords_jac_9__3;
         const double p_affine_const_0_2__3 = tmp_coords_jac_20__3 + tmp_coords_jac_15__3 * tmp_coords_jac_2__3;
         const double p_affine_const_1_0__3 =
             tmp_coords_jac_21__3 + tmp_coords_jac_7__3 + tmp_coords_jac_2__3 * tmp_coords_jac_3__3;
         const double p_affine_const_1_1__3 =
             tmp_coords_jac_13__3 + tmp_coords_jac_22__3 + tmp_coords_jac_10__3 * tmp_coords_jac_2__3;
         const double p_affine_const_1_2__3 =
             tmp_coords_jac_19__3 + tmp_coords_jac_23__3 + tmp_coords_jac_16__3 * tmp_coords_jac_2__3;
         const double p_affine_const_2_0__3 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_21__3 + tmp_coords_jac_5__3 +
                                              tmp_coords_jac_2__3 * tmp_coords_jac_6__3;
         const double p_affine_const_2_1__3 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_11__3 + tmp_coords_jac_22__3 +
                                              tmp_coords_jac_12__3 * tmp_coords_jac_2__3;
         const double p_affine_const_2_2__3 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_17__3 + tmp_coords_jac_23__3 +
                                              tmp_coords_jac_18__3 * tmp_coords_jac_2__3;
         const double p_affine_const_3_0__3 = tmp_coords_jac_21__3 + tmp_coords_jac_8__3;
         const double p_affine_const_3_1__3 = tmp_coords_jac_14__3 + tmp_coords_jac_22__3;
         const double p_affine_const_3_2__3 = tmp_coords_jac_20__3 + tmp_coords_jac_23__3;
         const double jac_affine_0_0__3     = p_affine_const_1_0__3 - p_affine_const_0_0__3;
         const double jac_affine_0_1__3     = p_affine_const_2_0__3 - p_affine_const_0_0__3;
         const double jac_affine_0_2__3     = p_affine_const_3_0__3 - p_affine_const_0_0__3;
         const double jac_affine_1_0__3     = p_affine_const_1_1__3 - p_affine_const_0_1__3;
         const double jac_affine_1_1__3     = p_affine_const_2_1__3 - p_affine_const_0_1__3;
         const double tmp_coords_jac_28__3  = jac_affine_0_2__3 * jac_affine_1_1__3;
         const double jac_affine_1_2__3     = p_affine_const_3_1__3 - p_affine_const_0_1__3;
         const double tmp_coords_jac_26__3  = jac_affine_0_1__3 * jac_affine_1_2__3;
         const double jac_affine_2_0__3     = p_affine_const_1_2__3 - p_affine_const_0_2__3;
         const double jac_affine_2_1__3     = p_affine_const_2_2__3 - p_affine_const_0_2__3;
         const double tmp_coords_jac_25__3  = jac_affine_1_2__3 * jac_affine_2_1__3;
         const double jac_affine_2_2__3     = p_affine_const_3_2__3 - p_affine_const_0_2__3;
         const double tmp_coords_jac_24__3  = jac_affine_1_1__3 * jac_affine_2_2__3;
         const double tmp_coords_jac_27__3  = jac_affine_0_1__3 * jac_affine_2_2__3;
         const double tmp_coords_jac_29__3 = jac_affine_0_0__3 * tmp_coords_jac_24__3 + jac_affine_2_0__3 * tmp_coords_jac_26__3 -
                                             jac_affine_0_0__3 * tmp_coords_jac_25__3 - jac_affine_1_0__3 * tmp_coords_jac_27__3 -
                                             jac_affine_2_0__3 * tmp_coords_jac_28__3 +
                                             jac_affine_0_2__3 * jac_affine_1_0__3 * jac_affine_2_1__3;
         const double tmp_coords_jac_30__3  = 1.0 / tmp_coords_jac_29__3;
         const double jac_affine_inv_0_0__3 = tmp_coords_jac_30__3 * ( tmp_coords_jac_24__3 - tmp_coords_jac_25__3 );
         const double jac_affine_inv_0_1__3 =
             tmp_coords_jac_30__3 * ( -1.0 * tmp_coords_jac_27__3 + jac_affine_0_2__3 * jac_affine_2_1__3 );
         const double jac_affine_inv_0_2__3 = tmp_coords_jac_30__3 * ( tmp_coords_jac_26__3 - tmp_coords_jac_28__3 );
         const double jac_affine_inv_1_0__3 =
             tmp_coords_jac_30__3 * ( jac_affine_1_2__3 * jac_affine_2_0__3 - jac_affine_1_0__3 * jac_affine_2_2__3 );
         const double jac_affine_inv_1_1__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_0__3 * jac_affine_2_2__3 - jac_affine_0_2__3 * jac_affine_2_0__3 );
         const double jac_affine_inv_1_2__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_2__3 * jac_affine_1_0__3 - jac_affine_0_0__3 * jac_affine_1_2__3 );
         const double jac_affine_inv_2_0__3 =
             tmp_coords_jac_30__3 * ( jac_affine_1_0__3 * jac_affine_2_1__3 - jac_affine_1_1__3 * jac_affine_2_0__3 );
         const double jac_affine_inv_2_1__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_1__3 * jac_affine_2_0__3 - jac_affine_0_0__3 * jac_affine_2_1__3 );
         const double jac_affine_inv_2_2__3 =
             tmp_coords_jac_30__3 * ( jac_affine_0_0__3 * jac_affine_1_1__3 - jac_affine_0_1__3 * jac_affine_1_0__3 );
         const double abs_det_jac_affine__3 = abs( tmp_coords_jac_29__3 );
         for ( int64_t ctr_2__3 = 0LL; ctr_2__3 < micro_edges_per_macro_edge; ctr_2__3 += 1LL )
         {
            for ( int64_t ctr_1__3 = 0LL; ctr_1__3 < -1LL * ctr_2__3 + micro_edges_per_macro_edge; ctr_1__3 += 1LL )
            {
               for ( int64_t ctr_0__3 = 0LL; ctr_0__3 < -2LL - ctr_1__3 - ctr_2__3 + micro_edges_per_macro_edge; ctr_0__3 += 1LL )
               {
                  const double p_affine_0_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__3;
                  const double p_affine_0_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__3;
                  const double p_affine_0_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__3;
                  const double p_affine_1_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_1_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_1_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_2_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__3 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_0__3 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_1__3 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__3 );
                  const double p_affine_3_2__3 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__3 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__3 );
                  const double tmp_kernel_op_0__3 = 0.16666666666666663 * abs_det_jac_affine__3;
                  const double elMatDiag_0__3 =
                      tmp_kernel_op_0__3 *
                      ( ( -1.0 * jac_affine_inv_0_0__3 - jac_affine_inv_1_0__3 - jac_affine_inv_2_0__3 ) *
                            ( -1.0 * jac_affine_inv_0_0__3 - jac_affine_inv_1_0__3 - jac_affine_inv_2_0__3 ) +
                        ( -1.0 * jac_affine_inv_0_1__3 - jac_affine_inv_1_1__3 - jac_affine_inv_2_1__3 ) *
                            ( -1.0 * jac_affine_inv_0_1__3 - jac_affine_inv_1_1__3 - jac_affine_inv_2_1__3 ) +
                        ( -1.0 * jac_affine_inv_0_2__3 - jac_affine_inv_1_2__3 - jac_affine_inv_2_2__3 ) *
                            ( -1.0 * jac_affine_inv_0_2__3 - jac_affine_inv_1_2__3 - jac_affine_inv_2_2__3 ) );
                  const double elMatDiag_1__3 = tmp_kernel_op_0__3 * ( jac_affine_inv_0_0__3 * jac_affine_inv_0_0__3 +
                                                                       jac_affine_inv_0_1__3 * jac_affine_inv_0_1__3 +
                                                                       jac_affine_inv_0_2__3 * jac_affine_inv_0_2__3 );
                  const double elMatDiag_2__3 = tmp_kernel_op_0__3 * ( jac_affine_inv_1_0__3 * jac_affine_inv_1_0__3 +
                                                                       jac_affine_inv_1_1__3 * jac_affine_inv_1_1__3 +
                                                                       jac_affine_inv_1_2__3 * jac_affine_inv_1_2__3 );
                  const double elMatDiag_3__3 = tmp_kernel_op_0__3 * ( jac_affine_inv_2_0__3 * jac_affine_inv_2_0__3 +
                                                                       jac_affine_inv_2_1__3 * jac_affine_inv_2_1__3 +
                                                                       jac_affine_inv_2_2__3 * jac_affine_inv_2_2__3 );
                  _data_invDiag_[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                 ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__3 ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_0__3 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                     ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__3 ) * ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1__3 ) * ctr_1__3 / 2LL -
                                 ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ctr_1__3 + ctr_0__3 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_1__3 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__3 ) * ctr_1__3 / 2LL -
                                     ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) * ctr_1__3 + ctr_0__3 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL ) -
                                 ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_2__3 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL ) -
                                     ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                 ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_3__3 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__3 ) * ( 2LL + ctr_1__3 ) / 2LL -
                                     ( -1LL * ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__3 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__3 ) * ( 1LL - ctr_2__3 + micro_edges_per_macro_edge ) + ctr_0__3 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.BLUE_UP */
         const double tmp_coords_jac_0__2   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__2   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__2   = tmp_coords_jac_1__2 * 0.0;
         const double tmp_coords_jac_3__2   = tmp_coords_jac_0__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_4__2   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__2   = tmp_coords_jac_1__2 * 1.0;
         const double tmp_coords_jac_6__2   = tmp_coords_jac_4__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_7__2   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_8__2   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__2 * tmp_coords_jac_7__2;
         const double tmp_coords_jac_9__2   = tmp_coords_jac_6__2 + tmp_coords_jac_8__2;
         const double tmp_coords_jac_10__2  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__2  = tmp_coords_jac_10__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_12__2  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__2  = tmp_coords_jac_12__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_14__2  = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_15__2  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_14__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_16__2  = tmp_coords_jac_13__2 + tmp_coords_jac_15__2;
         const double tmp_coords_jac_17__2  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18__2  = tmp_coords_jac_17__2 * tmp_coords_jac_2__2;
         const double tmp_coords_jac_19__2  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_20__2  = tmp_coords_jac_19__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_21__2  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_22__2  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_2__2 * tmp_coords_jac_21__2;
         const double tmp_coords_jac_23__2  = tmp_coords_jac_20__2 + tmp_coords_jac_22__2;
         const double tmp_coords_jac_24__2  = tmp_coords_jac_0__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_25__2  = tmp_coords_jac_10__2 * tmp_coords_jac_5__2;
         const double tmp_coords_jac_26__2  = tmp_coords_jac_17__2 * tmp_coords_jac_5__2;
         const double p_affine_const_0_0__2 = tmp_coords_jac_3__2 + tmp_coords_jac_9__2;
         const double p_affine_const_0_1__2 = tmp_coords_jac_11__2 + tmp_coords_jac_16__2;
         const double p_affine_const_0_2__2 = tmp_coords_jac_18__2 + tmp_coords_jac_23__2;
         const double p_affine_const_1_0__2 =
             tmp_coords_jac_24__2 + tmp_coords_jac_8__2 + tmp_coords_jac_2__2 * tmp_coords_jac_4__2;
         const double p_affine_const_1_1__2 =
             tmp_coords_jac_15__2 + tmp_coords_jac_25__2 + tmp_coords_jac_12__2 * tmp_coords_jac_2__2;
         const double p_affine_const_1_2__2 =
             tmp_coords_jac_22__2 + tmp_coords_jac_26__2 + tmp_coords_jac_19__2 * tmp_coords_jac_2__2;
         const double p_affine_const_2_0__2 = tmp_coords_jac_24__2 + tmp_coords_jac_9__2;
         const double p_affine_const_2_1__2 = tmp_coords_jac_16__2 + tmp_coords_jac_25__2;
         const double p_affine_const_2_2__2 = tmp_coords_jac_23__2 + tmp_coords_jac_26__2;
         const double p_affine_const_3_0__2 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__2 + tmp_coords_jac_6__2 + tmp_coords_jac_5__2 * tmp_coords_jac_7__2;
         const double p_affine_const_3_1__2 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_11__2 + tmp_coords_jac_13__2 +
                                              tmp_coords_jac_14__2 * tmp_coords_jac_5__2;
         const double p_affine_const_3_2__2 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18__2 + tmp_coords_jac_20__2 +
                                              tmp_coords_jac_21__2 * tmp_coords_jac_5__2;
         const double jac_affine_0_0__2    = p_affine_const_1_0__2 - p_affine_const_0_0__2;
         const double jac_affine_0_1__2    = p_affine_const_2_0__2 - p_affine_const_0_0__2;
         const double jac_affine_0_2__2    = p_affine_const_3_0__2 - p_affine_const_0_0__2;
         const double jac_affine_1_0__2    = p_affine_const_1_1__2 - p_affine_const_0_1__2;
         const double jac_affine_1_1__2    = p_affine_const_2_1__2 - p_affine_const_0_1__2;
         const double tmp_coords_jac_31__1 = jac_affine_0_2__2 * jac_affine_1_1__2;
         const double jac_affine_1_2__2    = p_affine_const_3_1__2 - p_affine_const_0_1__2;
         const double tmp_coords_jac_29__2 = jac_affine_0_1__2 * jac_affine_1_2__2;
         const double jac_affine_2_0__2    = p_affine_const_1_2__2 - p_affine_const_0_2__2;
         const double jac_affine_2_1__2    = p_affine_const_2_2__2 - p_affine_const_0_2__2;
         const double tmp_coords_jac_28__2 = jac_affine_1_2__2 * jac_affine_2_1__2;
         const double jac_affine_2_2__2    = p_affine_const_3_2__2 - p_affine_const_0_2__2;
         const double tmp_coords_jac_27__2 = jac_affine_1_1__2 * jac_affine_2_2__2;
         const double tmp_coords_jac_30__2 = jac_affine_0_1__2 * jac_affine_2_2__2;
         const double tmp_coords_jac_32__1 = jac_affine_0_0__2 * tmp_coords_jac_27__2 + jac_affine_2_0__2 * tmp_coords_jac_29__2 -
                                             jac_affine_0_0__2 * tmp_coords_jac_28__2 - jac_affine_1_0__2 * tmp_coords_jac_30__2 -
                                             jac_affine_2_0__2 * tmp_coords_jac_31__1 +
                                             jac_affine_0_2__2 * jac_affine_1_0__2 * jac_affine_2_1__2;
         const double tmp_coords_jac_33__1  = 1.0 / tmp_coords_jac_32__1;
         const double jac_affine_inv_0_0__2 = tmp_coords_jac_33__1 * ( tmp_coords_jac_27__2 - tmp_coords_jac_28__2 );
         const double jac_affine_inv_0_1__2 =
             tmp_coords_jac_33__1 * ( -1.0 * tmp_coords_jac_30__2 + jac_affine_0_2__2 * jac_affine_2_1__2 );
         const double jac_affine_inv_0_2__2 = tmp_coords_jac_33__1 * ( tmp_coords_jac_29__2 - tmp_coords_jac_31__1 );
         const double jac_affine_inv_1_0__2 =
             tmp_coords_jac_33__1 * ( jac_affine_1_2__2 * jac_affine_2_0__2 - jac_affine_1_0__2 * jac_affine_2_2__2 );
         const double jac_affine_inv_1_1__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_0__2 * jac_affine_2_2__2 - jac_affine_0_2__2 * jac_affine_2_0__2 );
         const double jac_affine_inv_1_2__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_2__2 * jac_affine_1_0__2 - jac_affine_0_0__2 * jac_affine_1_2__2 );
         const double jac_affine_inv_2_0__2 =
             tmp_coords_jac_33__1 * ( jac_affine_1_0__2 * jac_affine_2_1__2 - jac_affine_1_1__2 * jac_affine_2_0__2 );
         const double jac_affine_inv_2_1__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_1__2 * jac_affine_2_0__2 - jac_affine_0_0__2 * jac_affine_2_1__2 );
         const double jac_affine_inv_2_2__2 =
             tmp_coords_jac_33__1 * ( jac_affine_0_0__2 * jac_affine_1_1__2 - jac_affine_0_1__2 * jac_affine_1_0__2 );
         const double abs_det_jac_affine__2 = abs( tmp_coords_jac_32__1 );
         for ( int64_t ctr_2__2 = 0LL; ctr_2__2 < micro_edges_per_macro_edge; ctr_2__2 += 1LL )
         {
            for ( int64_t ctr_1__2 = 0LL; ctr_1__2 < -1LL * ctr_2__2 + micro_edges_per_macro_edge; ctr_1__2 += 1LL )
            {
               for ( int64_t ctr_0__2 = 0LL; ctr_0__2 < -1LL - ctr_1__2 - ctr_2__2 + micro_edges_per_macro_edge; ctr_0__2 += 1LL )
               {
                  const double p_affine_0_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_0_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_0_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_1_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_1_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_1_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_2_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__2;
                  const double p_affine_2_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__2;
                  const double p_affine_2_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__2;
                  const double p_affine_3_0__2 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__2 );
                  const double p_affine_3_1__2 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__2 );
                  const double p_affine_3_2__2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__2 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__2 );
                  const double tmp_kernel_op_0__2 = 0.16666666666666663 * abs_det_jac_affine__2;
                  const double elMatDiag_0__2 =
                      tmp_kernel_op_0__2 *
                      ( ( -1.0 * jac_affine_inv_0_0__2 - jac_affine_inv_1_0__2 - jac_affine_inv_2_0__2 ) *
                            ( -1.0 * jac_affine_inv_0_0__2 - jac_affine_inv_1_0__2 - jac_affine_inv_2_0__2 ) +
                        ( -1.0 * jac_affine_inv_0_1__2 - jac_affine_inv_1_1__2 - jac_affine_inv_2_1__2 ) *
                            ( -1.0 * jac_affine_inv_0_1__2 - jac_affine_inv_1_1__2 - jac_affine_inv_2_1__2 ) +
                        ( -1.0 * jac_affine_inv_0_2__2 - jac_affine_inv_1_2__2 - jac_affine_inv_2_2__2 ) *
                            ( -1.0 * jac_affine_inv_0_2__2 - jac_affine_inv_1_2__2 - jac_affine_inv_2_2__2 ) );
                  const double elMatDiag_1__2 = tmp_kernel_op_0__2 * ( jac_affine_inv_0_0__2 * jac_affine_inv_0_0__2 +
                                                                       jac_affine_inv_0_1__2 * jac_affine_inv_0_1__2 +
                                                                       jac_affine_inv_0_2__2 * jac_affine_inv_0_2__2 );
                  const double elMatDiag_2__2 = tmp_kernel_op_0__2 * ( jac_affine_inv_1_0__2 * jac_affine_inv_1_0__2 +
                                                                       jac_affine_inv_1_1__2 * jac_affine_inv_1_1__2 +
                                                                       jac_affine_inv_1_2__2 * jac_affine_inv_1_2__2 );
                  const double elMatDiag_3__2 = tmp_kernel_op_0__2 * ( jac_affine_inv_2_0__2 * jac_affine_inv_2_0__2 +
                                                                       jac_affine_inv_2_1__2 * jac_affine_inv_2_1__2 +
                                                                       jac_affine_inv_2_2__2 * jac_affine_inv_2_2__2 );
                  _data_invDiag_[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                 ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_0__2 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                     ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL ) -
                                 ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_1__2 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL ) -
                                     ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL -
                                 ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_2__2 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__2 ) * ( 2LL + ctr_1__2 ) / 2LL -
                                     ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__2 ) * ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) + ctr_0__2 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                 ( -1LL * ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_3__2 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__2 ) * ctr_1__2 / 2LL -
                                     ( -1LL * ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__2 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL - ctr_2__2 + micro_edges_per_macro_edge ) * ctr_1__2 + ctr_0__2 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.BLUE_DOWN */
         const double tmp_coords_jac_0__1   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__1   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__1   = tmp_coords_jac_1__1 * 0.0;
         const double tmp_coords_jac_3__1   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4__1   = tmp_coords_jac_1__1 * 1.0;
         const double tmp_coords_jac_5__1   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_6__1   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__1 * tmp_coords_jac_5__1;
         const double tmp_coords_jac_7__1   = tmp_coords_jac_6__1 + tmp_coords_jac_3__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_8__1   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_9__1   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__1  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_11__1  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10__1 * tmp_coords_jac_2__1;
         const double tmp_coords_jac_12__1  = tmp_coords_jac_11__1 + tmp_coords_jac_4__1 * tmp_coords_jac_9__1;
         const double tmp_coords_jac_13__1  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_14__1  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_15__1  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__1  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__1 * tmp_coords_jac_2__1;
         const double tmp_coords_jac_17__1  = tmp_coords_jac_16__1 + tmp_coords_jac_14__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_18__1  = tmp_coords_jac_0__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_19__1  = tmp_coords_jac_18__1 + tmp_coords_jac_2__1 * tmp_coords_jac_3__1;
         const double tmp_coords_jac_20__1  = tmp_coords_jac_4__1 * tmp_coords_jac_8__1;
         const double tmp_coords_jac_21__1  = tmp_coords_jac_20__1 + tmp_coords_jac_2__1 * tmp_coords_jac_9__1;
         const double tmp_coords_jac_22__1  = tmp_coords_jac_13__1 * tmp_coords_jac_4__1;
         const double tmp_coords_jac_23__1  = tmp_coords_jac_22__1 + tmp_coords_jac_14__1 * tmp_coords_jac_2__1;
         const double p_affine_const_0_0__1 = tmp_coords_jac_7__1 + tmp_coords_jac_0__1 * tmp_coords_jac_2__1;
         const double p_affine_const_0_1__1 = tmp_coords_jac_12__1 + tmp_coords_jac_2__1 * tmp_coords_jac_8__1;
         const double p_affine_const_0_2__1 = tmp_coords_jac_17__1 + tmp_coords_jac_13__1 * tmp_coords_jac_2__1;
         const double p_affine_const_1_0__1 = tmp_coords_jac_19__1 + tmp_coords_jac_6__1;
         const double p_affine_const_1_1__1 = tmp_coords_jac_11__1 + tmp_coords_jac_21__1;
         const double p_affine_const_1_2__1 = tmp_coords_jac_16__1 + tmp_coords_jac_23__1;
         const double p_affine_const_2_0__1 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_19__1 + tmp_coords_jac_4__1 * tmp_coords_jac_5__1;
         const double p_affine_const_2_1__1 =
             macro_vertex_coord_id_0comp1 + tmp_coords_jac_21__1 + tmp_coords_jac_10__1 * tmp_coords_jac_4__1;
         const double p_affine_const_2_2__1 =
             macro_vertex_coord_id_0comp2 + tmp_coords_jac_23__1 + tmp_coords_jac_15__1 * tmp_coords_jac_4__1;
         const double p_affine_const_3_0__1 = tmp_coords_jac_18__1 + tmp_coords_jac_7__1;
         const double p_affine_const_3_1__1 = tmp_coords_jac_12__1 + tmp_coords_jac_20__1;
         const double p_affine_const_3_2__1 = tmp_coords_jac_17__1 + tmp_coords_jac_22__1;
         const double jac_affine_0_0__1     = p_affine_const_1_0__1 - p_affine_const_0_0__1;
         const double jac_affine_0_1__1     = p_affine_const_2_0__1 - p_affine_const_0_0__1;
         const double jac_affine_0_2__1     = p_affine_const_3_0__1 - p_affine_const_0_0__1;
         const double jac_affine_1_0__1     = p_affine_const_1_1__1 - p_affine_const_0_1__1;
         const double jac_affine_1_1__1     = p_affine_const_2_1__1 - p_affine_const_0_1__1;
         const double tmp_coords_jac_28__1  = jac_affine_0_2__1 * jac_affine_1_1__1;
         const double jac_affine_1_2__1     = p_affine_const_3_1__1 - p_affine_const_0_1__1;
         const double tmp_coords_jac_26__1  = jac_affine_0_1__1 * jac_affine_1_2__1;
         const double jac_affine_2_0__1     = p_affine_const_1_2__1 - p_affine_const_0_2__1;
         const double jac_affine_2_1__1     = p_affine_const_2_2__1 - p_affine_const_0_2__1;
         const double tmp_coords_jac_25__1  = jac_affine_1_2__1 * jac_affine_2_1__1;
         const double jac_affine_2_2__1     = p_affine_const_3_2__1 - p_affine_const_0_2__1;
         const double tmp_coords_jac_24__1  = jac_affine_1_1__1 * jac_affine_2_2__1;
         const double tmp_coords_jac_27__1  = jac_affine_0_1__1 * jac_affine_2_2__1;
         const double tmp_coords_jac_29__1 = jac_affine_0_0__1 * tmp_coords_jac_24__1 + jac_affine_2_0__1 * tmp_coords_jac_26__1 -
                                             jac_affine_0_0__1 * tmp_coords_jac_25__1 - jac_affine_1_0__1 * tmp_coords_jac_27__1 -
                                             jac_affine_2_0__1 * tmp_coords_jac_28__1 +
                                             jac_affine_0_2__1 * jac_affine_1_0__1 * jac_affine_2_1__1;
         const double tmp_coords_jac_30__1  = 1.0 / tmp_coords_jac_29__1;
         const double jac_affine_inv_0_0__1 = tmp_coords_jac_30__1 * ( tmp_coords_jac_24__1 - tmp_coords_jac_25__1 );
         const double jac_affine_inv_0_1__1 =
             tmp_coords_jac_30__1 * ( -1.0 * tmp_coords_jac_27__1 + jac_affine_0_2__1 * jac_affine_2_1__1 );
         const double jac_affine_inv_0_2__1 = tmp_coords_jac_30__1 * ( tmp_coords_jac_26__1 - tmp_coords_jac_28__1 );
         const double jac_affine_inv_1_0__1 =
             tmp_coords_jac_30__1 * ( jac_affine_1_2__1 * jac_affine_2_0__1 - jac_affine_1_0__1 * jac_affine_2_2__1 );
         const double jac_affine_inv_1_1__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_0__1 * jac_affine_2_2__1 - jac_affine_0_2__1 * jac_affine_2_0__1 );
         const double jac_affine_inv_1_2__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_2__1 * jac_affine_1_0__1 - jac_affine_0_0__1 * jac_affine_1_2__1 );
         const double jac_affine_inv_2_0__1 =
             tmp_coords_jac_30__1 * ( jac_affine_1_0__1 * jac_affine_2_1__1 - jac_affine_1_1__1 * jac_affine_2_0__1 );
         const double jac_affine_inv_2_1__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_1__1 * jac_affine_2_0__1 - jac_affine_0_0__1 * jac_affine_2_1__1 );
         const double jac_affine_inv_2_2__1 =
             tmp_coords_jac_30__1 * ( jac_affine_0_0__1 * jac_affine_1_1__1 - jac_affine_0_1__1 * jac_affine_1_0__1 );
         const double abs_det_jac_affine__1 = abs( tmp_coords_jac_29__1 );
         for ( int64_t ctr_2__1 = 0LL; ctr_2__1 < micro_edges_per_macro_edge; ctr_2__1 += 1LL )
         {
            for ( int64_t ctr_1__1 = 0LL; ctr_1__1 < -1LL * ctr_2__1 + micro_edges_per_macro_edge; ctr_1__1 += 1LL )
            {
               for ( int64_t ctr_0__1 = 0LL; ctr_0__1 < -1LL - ctr_1__1 - ctr_2__1 + micro_edges_per_macro_edge; ctr_0__1 += 1LL )
               {
                  const double p_affine_0_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__1;
                  const double p_affine_0_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__1;
                  const double p_affine_0_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__1;
                  const double p_affine_1_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_1_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_1_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_2_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_0__1 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_1__1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__1 );
                  const double p_affine_3_2__1 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__1 );
                  const double tmp_kernel_op_0__1 = 0.16666666666666663 * abs_det_jac_affine__1;
                  const double elMatDiag_0__1 =
                      tmp_kernel_op_0__1 *
                      ( ( -1.0 * jac_affine_inv_0_0__1 - jac_affine_inv_1_0__1 - jac_affine_inv_2_0__1 ) *
                            ( -1.0 * jac_affine_inv_0_0__1 - jac_affine_inv_1_0__1 - jac_affine_inv_2_0__1 ) +
                        ( -1.0 * jac_affine_inv_0_1__1 - jac_affine_inv_1_1__1 - jac_affine_inv_2_1__1 ) *
                            ( -1.0 * jac_affine_inv_0_1__1 - jac_affine_inv_1_1__1 - jac_affine_inv_2_1__1 ) +
                        ( -1.0 * jac_affine_inv_0_2__1 - jac_affine_inv_1_2__1 - jac_affine_inv_2_2__1 ) *
                            ( -1.0 * jac_affine_inv_0_2__1 - jac_affine_inv_1_2__1 - jac_affine_inv_2_2__1 ) );
                  const double elMatDiag_1__1 = tmp_kernel_op_0__1 * ( jac_affine_inv_0_0__1 * jac_affine_inv_0_0__1 +
                                                                       jac_affine_inv_0_1__1 * jac_affine_inv_0_1__1 +
                                                                       jac_affine_inv_0_2__1 * jac_affine_inv_0_2__1 );
                  const double elMatDiag_2__1 = tmp_kernel_op_0__1 * ( jac_affine_inv_1_0__1 * jac_affine_inv_1_0__1 +
                                                                       jac_affine_inv_1_1__1 * jac_affine_inv_1_1__1 +
                                                                       jac_affine_inv_1_2__1 * jac_affine_inv_1_2__1 );
                  const double elMatDiag_3__1 = tmp_kernel_op_0__1 * ( jac_affine_inv_2_0__1 * jac_affine_inv_2_0__1 +
                                                                       jac_affine_inv_2_1__1 * jac_affine_inv_2_1__1 +
                                                                       jac_affine_inv_2_2__1 * jac_affine_inv_2_2__1 );
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                 ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__1 ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_0__1 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                     ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__1 ) * ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL ) -
                                 ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_1__1 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL ) -
                                     ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL -
                                 ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_2__1 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__1 ) * ctr_1__1 / 2LL -
                                     ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) * ctr_1__1 + ctr_0__1 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                 ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__1 ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_3__1 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__1 ) * ( 2LL + ctr_1__1 ) / 2LL ) -
                                     ( -1LL * ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__1 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__1 ) * ( 1LL - ctr_2__1 + micro_edges_per_macro_edge ) + ctr_0__1 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.GREEN_UP */
         const double tmp_coords_jac_0__0   = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__0   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__0   = tmp_coords_jac_1__0 * 0.0;
         const double tmp_coords_jac_3__0   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_0__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_4__0   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__0   = tmp_coords_jac_1__0 * 1.0;
         const double tmp_coords_jac_6__0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_7__0   = tmp_coords_jac_2__0 * tmp_coords_jac_6__0;
         const double tmp_coords_jac_8__0   = tmp_coords_jac_7__0 + tmp_coords_jac_4__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_9__0   = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10__0  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_2__0 * tmp_coords_jac_9__0;
         const double tmp_coords_jac_11__0  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12__0  = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_13__0  = tmp_coords_jac_12__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_14__0  = tmp_coords_jac_13__0 + tmp_coords_jac_11__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_15__0  = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16__0  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_17__0  = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18__0  = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_19__0  = tmp_coords_jac_18__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_20__0  = tmp_coords_jac_19__0 + tmp_coords_jac_17__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_21__0  = tmp_coords_jac_2__0 * tmp_coords_jac_4__0;
         const double tmp_coords_jac_22__0  = tmp_coords_jac_11__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_23__0  = tmp_coords_jac_17__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_24__0  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_0__0 * tmp_coords_jac_5__0;
         const double tmp_coords_jac_25__0  = macro_vertex_coord_id_0comp1 + tmp_coords_jac_5__0 * tmp_coords_jac_9__0;
         const double tmp_coords_jac_26__0  = macro_vertex_coord_id_0comp2 + tmp_coords_jac_15__0 * tmp_coords_jac_5__0;
         const double p_affine_const_0_0__0 = tmp_coords_jac_3__0 + tmp_coords_jac_8__0;
         const double p_affine_const_0_1__0 = tmp_coords_jac_10__0 + tmp_coords_jac_14__0;
         const double p_affine_const_0_2__0 = tmp_coords_jac_16__0 + tmp_coords_jac_20__0;
         const double p_affine_const_1_0__0 =
             tmp_coords_jac_21__0 + tmp_coords_jac_3__0 + tmp_coords_jac_5__0 * tmp_coords_jac_6__0;
         const double p_affine_const_1_1__0 =
             tmp_coords_jac_10__0 + tmp_coords_jac_22__0 + tmp_coords_jac_12__0 * tmp_coords_jac_5__0;
         const double p_affine_const_1_2__0 =
             tmp_coords_jac_16__0 + tmp_coords_jac_23__0 + tmp_coords_jac_18__0 * tmp_coords_jac_5__0;
         const double p_affine_const_2_0__0 = tmp_coords_jac_21__0 + tmp_coords_jac_24__0 + tmp_coords_jac_7__0;
         const double p_affine_const_2_1__0 = tmp_coords_jac_13__0 + tmp_coords_jac_22__0 + tmp_coords_jac_25__0;
         const double p_affine_const_2_2__0 = tmp_coords_jac_19__0 + tmp_coords_jac_23__0 + tmp_coords_jac_26__0;
         const double p_affine_const_3_0__0 = tmp_coords_jac_24__0 + tmp_coords_jac_8__0;
         const double p_affine_const_3_1__0 = tmp_coords_jac_14__0 + tmp_coords_jac_25__0;
         const double p_affine_const_3_2__0 = tmp_coords_jac_20__0 + tmp_coords_jac_26__0;
         const double jac_affine_0_0__0     = p_affine_const_1_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_1__0     = p_affine_const_2_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_2__0     = p_affine_const_3_0__0 - p_affine_const_0_0__0;
         const double jac_affine_1_0__0     = p_affine_const_1_1__0 - p_affine_const_0_1__0;
         const double jac_affine_1_1__0     = p_affine_const_2_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_31__0  = jac_affine_0_2__0 * jac_affine_1_1__0;
         const double jac_affine_1_2__0     = p_affine_const_3_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_29__0  = jac_affine_0_1__0 * jac_affine_1_2__0;
         const double jac_affine_2_0__0     = p_affine_const_1_2__0 - p_affine_const_0_2__0;
         const double jac_affine_2_1__0     = p_affine_const_2_2__0 - p_affine_const_0_2__0;
         const double tmp_coords_jac_28__0  = jac_affine_1_2__0 * jac_affine_2_1__0;
         const double jac_affine_2_2__0     = p_affine_const_3_2__0 - p_affine_const_0_2__0;
         const double tmp_coords_jac_27__0  = jac_affine_1_1__0 * jac_affine_2_2__0;
         const double tmp_coords_jac_30__0  = jac_affine_0_1__0 * jac_affine_2_2__0;
         const double tmp_coords_jac_32__0 = jac_affine_0_0__0 * tmp_coords_jac_27__0 + jac_affine_2_0__0 * tmp_coords_jac_29__0 -
                                             jac_affine_0_0__0 * tmp_coords_jac_28__0 - jac_affine_1_0__0 * tmp_coords_jac_30__0 -
                                             jac_affine_2_0__0 * tmp_coords_jac_31__0 +
                                             jac_affine_0_2__0 * jac_affine_1_0__0 * jac_affine_2_1__0;
         const double tmp_coords_jac_33__0  = 1.0 / tmp_coords_jac_32__0;
         const double jac_affine_inv_0_0__0 = tmp_coords_jac_33__0 * ( tmp_coords_jac_27__0 - tmp_coords_jac_28__0 );
         const double jac_affine_inv_0_1__0 =
             tmp_coords_jac_33__0 * ( -1.0 * tmp_coords_jac_30__0 + jac_affine_0_2__0 * jac_affine_2_1__0 );
         const double jac_affine_inv_0_2__0 = tmp_coords_jac_33__0 * ( tmp_coords_jac_29__0 - tmp_coords_jac_31__0 );
         const double jac_affine_inv_1_0__0 =
             tmp_coords_jac_33__0 * ( jac_affine_1_2__0 * jac_affine_2_0__0 - jac_affine_1_0__0 * jac_affine_2_2__0 );
         const double jac_affine_inv_1_1__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_0__0 * jac_affine_2_2__0 - jac_affine_0_2__0 * jac_affine_2_0__0 );
         const double jac_affine_inv_1_2__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_2__0 * jac_affine_1_0__0 - jac_affine_0_0__0 * jac_affine_1_2__0 );
         const double jac_affine_inv_2_0__0 =
             tmp_coords_jac_33__0 * ( jac_affine_1_0__0 * jac_affine_2_1__0 - jac_affine_1_1__0 * jac_affine_2_0__0 );
         const double jac_affine_inv_2_1__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_1__0 * jac_affine_2_0__0 - jac_affine_0_0__0 * jac_affine_2_1__0 );
         const double jac_affine_inv_2_2__0 =
             tmp_coords_jac_33__0 * ( jac_affine_0_0__0 * jac_affine_1_1__0 - jac_affine_0_1__0 * jac_affine_1_0__0 );
         const double abs_det_jac_affine__0 = abs( tmp_coords_jac_32__0 );
         for ( int64_t ctr_2__0 = 0LL; ctr_2__0 < micro_edges_per_macro_edge; ctr_2__0 += 1LL )
         {
            for ( int64_t ctr_1__0 = 0LL; ctr_1__0 < -1LL * ctr_2__0 + micro_edges_per_macro_edge; ctr_1__0 += 1LL )
            {
               for ( int64_t ctr_0__0 = 0LL; ctr_0__0 < -1LL - ctr_1__0 - ctr_2__0 + micro_edges_per_macro_edge; ctr_0__0 += 1LL )
               {
                  const double p_affine_0_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__0;
                  const double p_affine_0_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__0;
                  const double p_affine_0_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__0;
                  const double p_affine_1_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2__0;
                  const double p_affine_1_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2__0;
                  const double p_affine_1_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2__0;
                  const double p_affine_2_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_2_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_2_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_0__0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_1__0 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2__0 );
                  const double p_affine_3_2__0 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0__0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1__0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2__0 );
                  const double tmp_kernel_op_0__0 = 0.16666666666666663 * abs_det_jac_affine__0;
                  const double elMatDiag_0__0 =
                      tmp_kernel_op_0__0 *
                      ( ( -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0 - jac_affine_inv_2_0__0 ) *
                            ( -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0 - jac_affine_inv_2_0__0 ) +
                        ( -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0 - jac_affine_inv_2_1__0 ) *
                            ( -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0 - jac_affine_inv_2_1__0 ) +
                        ( -1.0 * jac_affine_inv_0_2__0 - jac_affine_inv_1_2__0 - jac_affine_inv_2_2__0 ) *
                            ( -1.0 * jac_affine_inv_0_2__0 - jac_affine_inv_1_2__0 - jac_affine_inv_2_2__0 ) );
                  const double elMatDiag_1__0 = tmp_kernel_op_0__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_0_0__0 +
                                                                       jac_affine_inv_0_1__0 * jac_affine_inv_0_1__0 +
                                                                       jac_affine_inv_0_2__0 * jac_affine_inv_0_2__0 );
                  const double elMatDiag_2__0 = tmp_kernel_op_0__0 * ( jac_affine_inv_1_0__0 * jac_affine_inv_1_0__0 +
                                                                       jac_affine_inv_1_1__0 * jac_affine_inv_1_1__0 +
                                                                       jac_affine_inv_1_2__0 * jac_affine_inv_1_2__0 );
                  const double elMatDiag_3__0 = tmp_kernel_op_0__0 * ( jac_affine_inv_2_0__0 * jac_affine_inv_2_0__0 +
                                                                       jac_affine_inv_2_1__0 * jac_affine_inv_2_1__0 +
                                                                       jac_affine_inv_2_2__0 * jac_affine_inv_2_2__0 );
                  _data_invDiag_[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                 ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_0__0 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                     ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) -
                                 ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1__0 ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) + ctr_0__0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_1__0 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) -
                                     ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1__0 ) * ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) + ctr_0__0 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) -
                                 ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_2__0 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) -
                                     ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                 ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_3__0 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL -
                                     ( -1LL * ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) *
                                         ( 2LL - ctr_2__0 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL - ctr_2__0 + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
      {
         /* CellType.GREEN_DOWN */
         const double tmp_coords_jac_0  = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1  = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2  = tmp_coords_jac_1 * 0.0;
         const double tmp_coords_jac_3  = tmp_coords_jac_0 * tmp_coords_jac_2;
         const double tmp_coords_jac_4  = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5  = tmp_coords_jac_1 * 1.0;
         const double tmp_coords_jac_6  = tmp_coords_jac_4 * tmp_coords_jac_5;
         const double tmp_coords_jac_7  = macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_8  = macro_vertex_coord_id_0comp0 + tmp_coords_jac_6 + tmp_coords_jac_2 * tmp_coords_jac_7;
         const double tmp_coords_jac_9  = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_10 = tmp_coords_jac_2 * tmp_coords_jac_9;
         const double tmp_coords_jac_11 = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_12 = tmp_coords_jac_11 * tmp_coords_jac_5;
         const double tmp_coords_jac_13 = macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_14 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_12 + tmp_coords_jac_13 * tmp_coords_jac_2;
         const double tmp_coords_jac_15 = macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_16 = tmp_coords_jac_15 * tmp_coords_jac_2;
         const double tmp_coords_jac_17 = macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_18 = tmp_coords_jac_17 * tmp_coords_jac_5;
         const double tmp_coords_jac_19 = macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2;
         const double tmp_coords_jac_20 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_18 + tmp_coords_jac_19 * tmp_coords_jac_2;
         const double tmp_coords_jac_21 = tmp_coords_jac_0 * tmp_coords_jac_5;
         const double tmp_coords_jac_22 = tmp_coords_jac_5 * tmp_coords_jac_9;
         const double tmp_coords_jac_23 = tmp_coords_jac_15 * tmp_coords_jac_5;
         const double tmp_coords_jac_24 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_5 * tmp_coords_jac_7;
         const double tmp_coords_jac_25 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_13 * tmp_coords_jac_5;
         const double tmp_coords_jac_26 = macro_vertex_coord_id_0comp2 + tmp_coords_jac_19 * tmp_coords_jac_5;
         const double p_affine_const_0_0 = tmp_coords_jac_3 + tmp_coords_jac_8;
         const double p_affine_const_0_1 = tmp_coords_jac_10 + tmp_coords_jac_14;
         const double p_affine_const_0_2 = tmp_coords_jac_16 + tmp_coords_jac_20;
         const double p_affine_const_1_0 = tmp_coords_jac_21 + tmp_coords_jac_8;
         const double p_affine_const_1_1 = tmp_coords_jac_14 + tmp_coords_jac_22;
         const double p_affine_const_1_2 = tmp_coords_jac_20 + tmp_coords_jac_23;
         const double p_affine_const_2_0 = tmp_coords_jac_21 + tmp_coords_jac_24 + tmp_coords_jac_2 * tmp_coords_jac_4;
         const double p_affine_const_2_1 = tmp_coords_jac_22 + tmp_coords_jac_25 + tmp_coords_jac_11 * tmp_coords_jac_2;
         const double p_affine_const_2_2 = tmp_coords_jac_23 + tmp_coords_jac_26 + tmp_coords_jac_17 * tmp_coords_jac_2;
         const double p_affine_const_3_0 = tmp_coords_jac_24 + tmp_coords_jac_3 + tmp_coords_jac_6;
         const double p_affine_const_3_1 = tmp_coords_jac_10 + tmp_coords_jac_12 + tmp_coords_jac_25;
         const double p_affine_const_3_2 = tmp_coords_jac_16 + tmp_coords_jac_18 + tmp_coords_jac_26;
         const double jac_affine_0_0     = p_affine_const_1_0 - p_affine_const_0_0;
         const double jac_affine_0_1     = p_affine_const_2_0 - p_affine_const_0_0;
         const double jac_affine_0_2     = p_affine_const_3_0 - p_affine_const_0_0;
         const double jac_affine_1_0     = p_affine_const_1_1 - p_affine_const_0_1;
         const double jac_affine_1_1     = p_affine_const_2_1 - p_affine_const_0_1;
         const double tmp_coords_jac_31  = jac_affine_0_2 * jac_affine_1_1;
         const double jac_affine_1_2     = p_affine_const_3_1 - p_affine_const_0_1;
         const double tmp_coords_jac_29  = jac_affine_0_1 * jac_affine_1_2;
         const double jac_affine_2_0     = p_affine_const_1_2 - p_affine_const_0_2;
         const double jac_affine_2_1     = p_affine_const_2_2 - p_affine_const_0_2;
         const double tmp_coords_jac_28  = jac_affine_1_2 * jac_affine_2_1;
         const double jac_affine_2_2     = p_affine_const_3_2 - p_affine_const_0_2;
         const double tmp_coords_jac_27  = jac_affine_1_1 * jac_affine_2_2;
         const double tmp_coords_jac_30  = jac_affine_0_1 * jac_affine_2_2;
         const double tmp_coords_jac_32  = jac_affine_0_0 * tmp_coords_jac_27 + jac_affine_2_0 * tmp_coords_jac_29 -
                                          jac_affine_0_0 * tmp_coords_jac_28 - jac_affine_1_0 * tmp_coords_jac_30 -
                                          jac_affine_2_0 * tmp_coords_jac_31 + jac_affine_0_2 * jac_affine_1_0 * jac_affine_2_1;
         const double tmp_coords_jac_33  = 1.0 / tmp_coords_jac_32;
         const double jac_affine_inv_0_0 = tmp_coords_jac_33 * ( tmp_coords_jac_27 - tmp_coords_jac_28 );
         const double jac_affine_inv_0_1 = tmp_coords_jac_33 * ( -1.0 * tmp_coords_jac_30 + jac_affine_0_2 * jac_affine_2_1 );
         const double jac_affine_inv_0_2 = tmp_coords_jac_33 * ( tmp_coords_jac_29 - tmp_coords_jac_31 );
         const double jac_affine_inv_1_0 =
             tmp_coords_jac_33 * ( jac_affine_1_2 * jac_affine_2_0 - jac_affine_1_0 * jac_affine_2_2 );
         const double jac_affine_inv_1_1 =
             tmp_coords_jac_33 * ( jac_affine_0_0 * jac_affine_2_2 - jac_affine_0_2 * jac_affine_2_0 );
         const double jac_affine_inv_1_2 =
             tmp_coords_jac_33 * ( jac_affine_0_2 * jac_affine_1_0 - jac_affine_0_0 * jac_affine_1_2 );
         const double jac_affine_inv_2_0 =
             tmp_coords_jac_33 * ( jac_affine_1_0 * jac_affine_2_1 - jac_affine_1_1 * jac_affine_2_0 );
         const double jac_affine_inv_2_1 =
             tmp_coords_jac_33 * ( jac_affine_0_1 * jac_affine_2_0 - jac_affine_0_0 * jac_affine_2_1 );
         const double jac_affine_inv_2_2 =
             tmp_coords_jac_33 * ( jac_affine_0_0 * jac_affine_1_1 - jac_affine_0_1 * jac_affine_1_0 );
         const double abs_det_jac_affine = abs( tmp_coords_jac_32 );
         for ( int64_t ctr_2 = 0LL; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1LL )
         {
            for ( int64_t ctr_1 = 0LL; ctr_1 < -1LL * ctr_2 + micro_edges_per_macro_edge; ctr_1 += 1LL )
            {
               for ( int64_t ctr_0 = 0LL; ctr_0 < -1LL - ctr_1 - ctr_2 + micro_edges_per_macro_edge; ctr_0 += 1LL )
               {
                  const double p_affine_0_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2;
                  const double p_affine_0_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2;
                  const double p_affine_0_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2;
                  const double p_affine_1_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_2;
                  const double p_affine_1_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_2;
                  const double p_affine_1_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_2;
                  const double p_affine_2_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_2_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_2_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_0 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_0 =
                      macro_vertex_coord_id_0comp0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_1 =
                      macro_vertex_coord_id_0comp1 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_2 );
                  const double p_affine_3_2 =
                      macro_vertex_coord_id_0comp2 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_1comp2 - macro_vertex_coord_id_0comp2 ) * (double) ctr_0 +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_2comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_1 ) +
                      1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                          ( macro_vertex_coord_id_3comp2 - macro_vertex_coord_id_0comp2 ) * (double) ( 1LL + ctr_2 );
                  const double tmp_kernel_op_0 = 0.16666666666666663 * abs_det_jac_affine;
                  const double elMatDiag_0 =
                      tmp_kernel_op_0 * ( ( -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0 - jac_affine_inv_2_0 ) *
                                              ( -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0 - jac_affine_inv_2_0 ) +
                                          ( -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1 - jac_affine_inv_2_1 ) *
                                              ( -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1 - jac_affine_inv_2_1 ) +
                                          ( -1.0 * jac_affine_inv_0_2 - jac_affine_inv_1_2 - jac_affine_inv_2_2 ) *
                                              ( -1.0 * jac_affine_inv_0_2 - jac_affine_inv_1_2 - jac_affine_inv_2_2 ) );
                  const double elMatDiag_1 =
                      tmp_kernel_op_0 * ( jac_affine_inv_0_0 * jac_affine_inv_0_0 + jac_affine_inv_0_1 * jac_affine_inv_0_1 +
                                          jac_affine_inv_0_2 * jac_affine_inv_0_2 );
                  const double elMatDiag_2 =
                      tmp_kernel_op_0 * ( jac_affine_inv_1_0 * jac_affine_inv_1_0 + jac_affine_inv_1_1 * jac_affine_inv_1_1 +
                                          jac_affine_inv_1_2 * jac_affine_inv_1_2 );
                  const double elMatDiag_3 =
                      tmp_kernel_op_0 * ( jac_affine_inv_2_0 * jac_affine_inv_2_0 + jac_affine_inv_2_1 * jac_affine_inv_2_1 +
                                          jac_affine_inv_2_2 * jac_affine_inv_2_2 );
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                 ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_0 +
                      _data_invDiag_[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                     ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL -
                                 ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                     ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_1 +
                      _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL -
                                     ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) *
                                         ( 3LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                     ( 1LL + ctr_1 ) * ( 2LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                     ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                         ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL -
                                 ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ctr_1 + ctr_0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_2 + _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL -
                                                   ( -1LL * ctr_2 + micro_edges_per_macro_edge ) *
                                                       ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                                       ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                                   ( 1LL - ctr_2 + micro_edges_per_macro_edge ) * ctr_1 + ctr_0 +
                                                   ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                       ( 3LL + micro_edges_per_macro_edge ) / 6LL];
                  _data_invDiag_[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                 ( -1LL * ctr_2 + micro_edges_per_macro_edge ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                     ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                 ( 1LL + ctr_1 ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                 ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                     ( 3LL + micro_edges_per_macro_edge ) / 6LL] =
                      elMatDiag_3 + _data_invDiag_[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) -
                                                   ( -1LL * ctr_2 + micro_edges_per_macro_edge ) *
                                                       ( 1LL - ctr_2 + micro_edges_per_macro_edge ) *
                                                       ( 2LL - ctr_2 + micro_edges_per_macro_edge ) / 6LL +
                                                   ( 1LL + ctr_1 ) * ( 1LL - ctr_2 + micro_edges_per_macro_edge ) + ctr_0 +
                                                   ( 1LL + micro_edges_per_macro_edge ) * ( 2LL + micro_edges_per_macro_edge ) *
                                                       ( 3LL + micro_edges_per_macro_edge ) / 6LL];
               }
            }
         }
      }
   }
}
void P1ElementwiseDiffusion::computeInverseDiagonalOperatorValues_P1ElementwiseDiffusion_macro_2D(
    double* RESTRICT const _data_invDiag_,
    const double           macro_vertex_coord_id_0comp0,
    const double           macro_vertex_coord_id_0comp1,
    const double           macro_vertex_coord_id_1comp0,
    const double           macro_vertex_coord_id_1comp1,
    const double           macro_vertex_coord_id_2comp0,
    const double           macro_vertex_coord_id_2comp1,
    const int64_t          micro_edges_per_macro_edge,
    const double           micro_edges_per_macro_edge_float ) const
{
   {
      {
         /* FaceType.GRAY */
         const double tmp_coords_jac_0__0   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1__0   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2__0   = tmp_coords_jac_1__0 * 0.0;
         const double tmp_coords_jac_3__0   = tmp_coords_jac_0__0 * tmp_coords_jac_2__0;
         const double tmp_coords_jac_4__0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_5__0   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_2__0 * tmp_coords_jac_4__0;
         const double tmp_coords_jac_6__0   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_7__0   = tmp_coords_jac_2__0 * tmp_coords_jac_6__0;
         const double tmp_coords_jac_8__0   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_9__0   = macro_vertex_coord_id_0comp1 + tmp_coords_jac_2__0 * tmp_coords_jac_8__0;
         const double tmp_coords_jac_10__0  = tmp_coords_jac_1__0 * 1.0;
         const double p_affine_const_0_0__0 = tmp_coords_jac_3__0 + tmp_coords_jac_5__0;
         const double p_affine_const_0_1__0 = tmp_coords_jac_7__0 + tmp_coords_jac_9__0;
         const double p_affine_const_1_0__0 = tmp_coords_jac_5__0 + tmp_coords_jac_0__0 * tmp_coords_jac_10__0;
         const double p_affine_const_1_1__0 = tmp_coords_jac_9__0 + tmp_coords_jac_10__0 * tmp_coords_jac_6__0;
         const double p_affine_const_2_0__0 =
             macro_vertex_coord_id_0comp0 + tmp_coords_jac_3__0 + tmp_coords_jac_10__0 * tmp_coords_jac_4__0;
         const double p_affine_const_2_1__0 =
             macro_vertex_coord_id_0comp1 + tmp_coords_jac_7__0 + tmp_coords_jac_10__0 * tmp_coords_jac_8__0;
         const double jac_affine_0_0__0     = p_affine_const_1_0__0 - p_affine_const_0_0__0;
         const double jac_affine_0_1__0     = p_affine_const_2_0__0 - p_affine_const_0_0__0;
         const double jac_affine_1_0__0     = p_affine_const_1_1__0 - p_affine_const_0_1__0;
         const double jac_affine_1_1__0     = p_affine_const_2_1__0 - p_affine_const_0_1__0;
         const double tmp_coords_jac_11__0  = jac_affine_0_0__0 * jac_affine_1_1__0 - jac_affine_0_1__0 * jac_affine_1_0__0;
         const double tmp_coords_jac_12__0  = 1.0 / tmp_coords_jac_11__0;
         const double jac_affine_inv_0_0__0 = jac_affine_1_1__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_0_1__0 = -1.0 * jac_affine_0_1__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_1_0__0 = -1.0 * jac_affine_1_0__0 * tmp_coords_jac_12__0;
         const double jac_affine_inv_1_1__0 = jac_affine_0_0__0 * tmp_coords_jac_12__0;
         const double abs_det_jac_affine__0 = abs( tmp_coords_jac_11__0 );
         for ( int64_t ctr_1__0 = 0LL; ctr_1__0 < micro_edges_per_macro_edge; ctr_1__0 += 1LL )
         {
            for ( int64_t ctr_0__0 = 0LL; ctr_0__0 < -1LL * ctr_1__0 + micro_edges_per_macro_edge; ctr_0__0 += 1LL )
            {
               const double p_affine_0_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0;
               const double p_affine_0_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0;
               const double p_affine_1_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0__0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1__0;
               const double p_affine_1_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0__0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1__0;
               const double p_affine_2_0__0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1__0 );
               const double p_affine_2_1__0 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0__0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1__0 );
               const double tmp_kernel_op_0__0 = 0.5 * abs_det_jac_affine__0;
               const double elMatDiag_0__0 =
                   tmp_kernel_op_0__0 * ( ( -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0 ) *
                                              ( -1.0 * jac_affine_inv_0_0__0 - jac_affine_inv_1_0__0 ) +
                                          ( -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0 ) *
                                              ( -1.0 * jac_affine_inv_0_1__0 - jac_affine_inv_1_1__0 ) );
               const double elMatDiag_1__0 = tmp_kernel_op_0__0 * ( jac_affine_inv_0_0__0 * jac_affine_inv_0_0__0 +
                                                                    jac_affine_inv_0_1__0 * jac_affine_inv_0_1__0 );
               const double elMatDiag_2__0 = tmp_kernel_op_0__0 * ( jac_affine_inv_1_0__0 * jac_affine_inv_1_0__0 +
                                                                    jac_affine_inv_1_1__0 * jac_affine_inv_1_1__0 );
               _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) + ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 +
                              ctr_0__0] =
                   elMatDiag_0__0 + _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL ) +
                                                   ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               _data_invDiag_[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 +
                              ctr_0__0] =
                   elMatDiag_1__0 + _data_invDiag_[1LL - ( 1LL + ctr_1__0 ) * ctr_1__0 / 2LL +
                                                   ( 2LL + micro_edges_per_macro_edge ) * ctr_1__0 + ctr_0__0];
               _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) +
                              ( 1LL + ctr_1__0 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0__0] =
                   elMatDiag_2__0 + _data_invDiag_[-1LL * ( ( 1LL + ctr_1__0 ) * ( 2LL + ctr_1__0 ) / 2LL ) +
                                                   ( 1LL + ctr_1__0 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0__0];
            }
         }
      }
      {
         /* FaceType.BLUE */
         const double tmp_coords_jac_0   = macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_1   = 1.0 * ( 1.0 / micro_edges_per_macro_edge_float );
         const double tmp_coords_jac_2   = tmp_coords_jac_1 * 0.0;
         const double tmp_coords_jac_3   = macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0;
         const double tmp_coords_jac_4   = tmp_coords_jac_1 * 1.0;
         const double tmp_coords_jac_5   = macro_vertex_coord_id_0comp0 + tmp_coords_jac_3 * tmp_coords_jac_4;
         const double tmp_coords_jac_6   = macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_7   = macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1;
         const double tmp_coords_jac_8   = macro_vertex_coord_id_0comp1 + tmp_coords_jac_4 * tmp_coords_jac_7;
         const double tmp_coords_jac_9   = tmp_coords_jac_0 * tmp_coords_jac_4;
         const double tmp_coords_jac_10  = tmp_coords_jac_4 * tmp_coords_jac_6;
         const double p_affine_const_0_0 = tmp_coords_jac_5 + tmp_coords_jac_0 * tmp_coords_jac_2;
         const double p_affine_const_0_1 = tmp_coords_jac_8 + tmp_coords_jac_2 * tmp_coords_jac_6;
         const double p_affine_const_1_0 = macro_vertex_coord_id_0comp0 + tmp_coords_jac_9 + tmp_coords_jac_2 * tmp_coords_jac_3;
         const double p_affine_const_1_1 = macro_vertex_coord_id_0comp1 + tmp_coords_jac_10 + tmp_coords_jac_2 * tmp_coords_jac_7;
         const double p_affine_const_2_0 = tmp_coords_jac_5 + tmp_coords_jac_9;
         const double p_affine_const_2_1 = tmp_coords_jac_10 + tmp_coords_jac_8;
         const double jac_affine_0_0     = p_affine_const_1_0 - p_affine_const_0_0;
         const double jac_affine_0_1     = p_affine_const_2_0 - p_affine_const_0_0;
         const double jac_affine_1_0     = p_affine_const_1_1 - p_affine_const_0_1;
         const double jac_affine_1_1     = p_affine_const_2_1 - p_affine_const_0_1;
         const double tmp_coords_jac_11  = jac_affine_0_0 * jac_affine_1_1 - jac_affine_0_1 * jac_affine_1_0;
         const double tmp_coords_jac_12  = 1.0 / tmp_coords_jac_11;
         const double jac_affine_inv_0_0 = jac_affine_1_1 * tmp_coords_jac_12;
         const double jac_affine_inv_0_1 = -1.0 * jac_affine_0_1 * tmp_coords_jac_12;
         const double jac_affine_inv_1_0 = -1.0 * jac_affine_1_0 * tmp_coords_jac_12;
         const double jac_affine_inv_1_1 = jac_affine_0_0 * tmp_coords_jac_12;
         const double abs_det_jac_affine = abs( tmp_coords_jac_11 );
         for ( int64_t ctr_1 = 0LL; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1LL )
         {
            for ( int64_t ctr_0 = 0LL; ctr_0 < -1LL - ctr_1 + micro_edges_per_macro_edge; ctr_0 += 1LL )
            {
               const double p_affine_0_0 = macro_vertex_coord_id_0comp0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) *
                                               (double) ( 1LL + ctr_0 ) +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_1;
               const double p_affine_0_1 = macro_vertex_coord_id_0comp1 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) *
                                               (double) ( 1LL + ctr_0 ) +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_1;
               const double p_affine_1_0 = macro_vertex_coord_id_0comp0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ctr_0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) *
                                               (double) ( 1LL + ctr_1 );
               const double p_affine_1_1 = macro_vertex_coord_id_0comp1 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ctr_0 +
                                           1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                                               ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) *
                                               (double) ( 1LL + ctr_1 );
               const double p_affine_2_0 =
                   macro_vertex_coord_id_0comp0 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp0 - macro_vertex_coord_id_0comp0 ) * (double) ( 1LL + ctr_1 );
               const double p_affine_2_1 =
                   macro_vertex_coord_id_0comp1 +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_1comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_0 ) +
                   1.0 * ( 1.0 / micro_edges_per_macro_edge_float ) *
                       ( macro_vertex_coord_id_2comp1 - macro_vertex_coord_id_0comp1 ) * (double) ( 1LL + ctr_1 );
               const double tmp_kernel_op_0 = 0.5 * abs_det_jac_affine;
               const double elMatDiag_0 =
                   tmp_kernel_op_0 *
                   ( ( -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0 ) * ( -1.0 * jac_affine_inv_0_0 - jac_affine_inv_1_0 ) +
                     ( -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1 ) * ( -1.0 * jac_affine_inv_0_1 - jac_affine_inv_1_1 ) );
               const double elMatDiag_1 =
                   tmp_kernel_op_0 * ( jac_affine_inv_0_0 * jac_affine_inv_0_0 + jac_affine_inv_0_1 * jac_affine_inv_0_1 );
               const double elMatDiag_2 =
                   tmp_kernel_op_0 * ( jac_affine_inv_1_0 * jac_affine_inv_1_0 + jac_affine_inv_1_1 * jac_affine_inv_1_1 );
               _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1 + ctr_0] =
                   elMatDiag_0 +
                   _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ctr_1 / 2LL + ( 2LL + micro_edges_per_macro_edge ) * ctr_1 + ctr_0];
               _data_invDiag_[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) +
                              ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0] =
                   elMatDiag_1 + _data_invDiag_[-1LL * ( ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL ) +
                                                ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
               _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL +
                              ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0] =
                   elMatDiag_2 + _data_invDiag_[1LL - ( 1LL + ctr_1 ) * ( 2LL + ctr_1 ) / 2LL +
                                                ( 1LL + ctr_1 ) * ( 2LL + micro_edges_per_macro_edge ) + ctr_0];
            }
         }
      }
   }
}

} // namespace operatorgeneration

} // namespace hyteg
