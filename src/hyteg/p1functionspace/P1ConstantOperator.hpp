/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
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

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/P1RowSumForm.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_affine_qe.hpp"
#include "hyteg/p1functionspace/P1Operator.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p1functionspace/generatedKernels/gaussseidel_3D_macrocell_P1.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_2D_macroface_vertexdof_to_vertexdof.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_2D_macroface_vertexdof_to_vertexdof_backwards.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macrocell_P1.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macrocell_P1_backwards.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1_backwards.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1_one_sided.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1_one_sided_backwards.hpp"

namespace hyteg {

using walberla::real_t;

// todo add lumped, diagonal and inverse template flags and corresponding implementations
template < class P1Form, bool Diagonal = false, bool Lumped = false, bool InvertDiagonal = false >
class P1ConstantOperator : public P1Operator< P1Form >
{
   using P1Operator< P1Form >::P1Operator;
   using P1Operator< P1Form >::storage_;
   using P1Operator< P1Form >::form_;
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
   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : P1ConstantOperator( storage, minLevel, maxLevel, P1Form() )
   {}

   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form )
   : P1Operator< P1Form >( storage, minLevel, maxLevel, form )
   {
      // pre-assemble edge, face and cell stencils
      assembleStencils();
   }

   void scale( real_t scalar ){}; // todo implement

 protected:
   /// stencil assembly: stencils are pre-assembled -> nothing to do here! ///////////

   /* Initialize assembly of variable edge stencil.
      Will be called before iterating over edge whenever the stencil is applied.
   */
   inline void assemble_stencil_edge_init( Edge& edge, const uint_t level ) const {}

   /* Assembly of edge stencil.
      Will be called before stencil is applied to a particuar edge-DoF.
   */
   inline void assemble_stencil_edge( real_t* edge_stencil, const uint_t i ) const {}

   /* Initialize assembly of face stencil.
      Will be called before iterating over face whenever the stencil is applied.
   */
   inline void assemble_stencil_face_init( Face& face, const uint_t level ) const {}

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 2d domain.
   */
   inline void assemble_stencil_face( real_t* face_stencil, const uint_t i, const uint_t j ) const {}

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 3D domain.
   */
   inline void assemble_stencil_face3D( vertexdof::macroface::StencilMap_T& face_stencil, const uint_t i, const uint_t j ) const
   {}

   /* Initialize assembly of cell stencil.
      Will be called before iterating over cell whenever the stencil is applied.
   */
   inline void assemble_stencil_cell_init( Cell& cell, const uint_t level ) const {}

   /* Assembly of cell stencil.
      Will be called before stencil is applied to a particuar cell-DoF.
   */
   inline void assemble_stencil_cell( vertexdof::macrocell::StencilMap_T& cell_stencil,
                                      const uint_t                        i,
                                      const uint_t                        j,
                                      const uint_t                        k ) const
   {}

   /////////////////////////////////////////////////////////////////////////////////////////////////

   inline void apply_face3D_generated( Face&                                                    face,
                                       const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                                       const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                       const uint_t&                                            level,
                                       UpdateType                                               update ) const
   {
      if ( face.getNumNeighborCells() == 2 )
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Two-sided" ); }
      }
      else
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "One-sided" ); }
      }

      auto&        opr_data    = face.getData( faceStencil3DID_ )->getData( level );
      auto         src_data    = face.getData( srcId )->getPointer( level );
      auto         dst_data    = face.getData( dstId )->getPointer( level );
      const uint_t offset_gl_0 = levelinfo::num_microvertices_per_face( level );

      auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );

      auto neighbor_cell_0_local_vertex_id_0 = static_cast< int32_t >(
          neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 0 ) );
      auto neighbor_cell_0_local_vertex_id_1 = static_cast< int32_t >(
          neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 1 ) );
      auto neighbor_cell_0_local_vertex_id_2 = static_cast< int32_t >(
          neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 2 ) );

      if ( update == Replace )
      {
         vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace(
             dst_data,
             src_data,
             &src_data[offset_gl_0],
             static_cast< int32_t >( level ),
             neighbor_cell_0_local_vertex_id_0,
             neighbor_cell_0_local_vertex_id_1,
             neighbor_cell_0_local_vertex_id_2,
             opr_data[0] );
      }
      else
      {
         vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add(
             dst_data,
             src_data,
             &src_data[offset_gl_0],
             static_cast< int32_t >( level ),
             neighbor_cell_0_local_vertex_id_0,
             neighbor_cell_0_local_vertex_id_1,
             neighbor_cell_0_local_vertex_id_2,
             opr_data[0] );
      }

      if ( face.getNumNeighborCells() == 2 )
      {
         const uint_t offset_gl_1 =
             offset_gl_0 + levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microvertices_per_edge( level ) - 1 );

         auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

         auto neighbor_cell_1_local_vertex_id_0 = static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                                              .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                                              .at( 0 ) );
         auto neighbor_cell_1_local_vertex_id_1 = static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                                              .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                                              .at( 1 ) );
         auto neighbor_cell_1_local_vertex_id_2 = static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                                              .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                                              .at( 2 ) );

         vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add(
             dst_data,
             src_data,
             &src_data[offset_gl_1],
             static_cast< int32_t >( level ),
             neighbor_cell_1_local_vertex_id_0,
             neighbor_cell_1_local_vertex_id_1,
             neighbor_cell_1_local_vertex_id_2,
             opr_data[1] );
      }

      if ( face.getNumNeighborCells() == 2 )
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Two-sided" ); }
      }
      else
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "One-sided" ); }
      }
   }

   inline void apply_face_generated( Face&                                                    face,
                                     const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                                     const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                     const uint_t&                                            level,
                                     UpdateType                                               update ) const
   {
      real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
      real_t* src_data = face.getData( srcId )->getPointer( level );
      real_t* dst_data = face.getData( dstId )->getPointer( level );

      if ( update == hyteg::Replace )
      {
         vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace(
             dst_data, src_data, opr_data, static_cast< int32_t >( level ) );
      }
      else if ( update == hyteg::Add )
      {
         vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
             dst_data, src_data, opr_data, static_cast< int32_t >( level ) );
      }
   }

   inline void apply_cell_generated( Cell&                                                    cell,
                                     const PrimitiveDataID< FunctionMemory< real_t >, Cell >& srcId,
                                     const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
                                     const uint_t&                                            level,
                                     UpdateType                                               update ) const
   {
      auto&   opr_data = cell.getData( cellStencilID_ )->getData( level );
      real_t* src_data = cell.getData( srcId )->getPointer( level );
      real_t* dst_data = cell.getData( dstId )->getPointer( level );

      if ( update == Replace )
      {
         vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
             dst_data, src_data, static_cast< int32_t >( level ), opr_data );
      }
      else if ( update == Add )
      {
         vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
             dst_data, src_data, static_cast< int32_t >( level ), opr_data );
      }
   }

   inline void smooth_sor_face3D_generated( Face&                                                    face,
                                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                            const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
                                            const uint_t&                                            level,
                                            real_t                                                   relax,
                                            const bool&                                              backwards = false ) const
   {
      auto  rhs_data = face.getData( rhsId )->getPointer( level );
      auto  dst_data = face.getData( dstId )->getPointer( level );
      auto& stencil  = face.getData( faceStencil3DID_ )->getData( level );

      auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );

      auto neighbor_cell_0_local_vertex_id_0 = static_cast< int32_t >(
          neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 0 ) );
      auto neighbor_cell_0_local_vertex_id_1 = static_cast< int32_t >(
          neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 1 ) );
      auto neighbor_cell_0_local_vertex_id_2 = static_cast< int32_t >(
          neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 2 ) );

      const uint_t vertex_offset_gl_0 = levelinfo::num_microvertices_per_face( level );

      if ( face.getNumNeighborCells() == 1 )
      {
         this->timingTree_->start( "One-sided" );

         if ( backwards )
         {
            vertexdof::macroface::generated::sor_3D_macroface_P1_one_sided_backwards( dst_data,
                                                                                      &dst_data[vertex_offset_gl_0],
                                                                                      rhs_data,
                                                                                      static_cast< int32_t >( level ),
                                                                                      neighbor_cell_0_local_vertex_id_0,
                                                                                      neighbor_cell_0_local_vertex_id_1,
                                                                                      neighbor_cell_0_local_vertex_id_2,
                                                                                      relax,
                                                                                      stencil[0] );
         }
         else
         {
            vertexdof::macroface::generated::sor_3D_macroface_P1_one_sided( dst_data,
                                                                            &dst_data[vertex_offset_gl_0],
                                                                            rhs_data,
                                                                            static_cast< int32_t >( level ),
                                                                            neighbor_cell_0_local_vertex_id_0,
                                                                            neighbor_cell_0_local_vertex_id_1,
                                                                            neighbor_cell_0_local_vertex_id_2,
                                                                            relax,
                                                                            stencil[0] );
         }

         this->timingTree_->stop( "One-sided" );
      }

      if ( face.getNumNeighborCells() == 2 )
      {
         this->timingTree_->start( "Two-sided" );

         auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

         auto neighbor_cell_1_local_vertex_id_0 = static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                                              .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                                              .at( 0 ) );
         auto neighbor_cell_1_local_vertex_id_1 = static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                                              .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                                              .at( 1 ) );
         auto neighbor_cell_1_local_vertex_id_2 = static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                                              .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                                              .at( 2 ) );
         const uint_t vertex_offset_gl_1        = vertex_offset_gl_0 + levelinfo::num_microvertices_per_face_from_width(
                                                                    levelinfo::num_microvertices_per_edge( level ) - 1 );

         if ( neighbor_cell_0_local_vertex_id_0 > neighbor_cell_1_local_vertex_id_0 ||
              ( neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                neighbor_cell_0_local_vertex_id_1 > neighbor_cell_1_local_vertex_id_1 ) ||
              ( neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                neighbor_cell_0_local_vertex_id_1 == neighbor_cell_1_local_vertex_id_1 &&
                neighbor_cell_0_local_vertex_id_2 > neighbor_cell_1_local_vertex_id_2 ) )
         {
            if ( backwards )
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1_backwards( dst_data,
                                                                               &dst_data[vertex_offset_gl_1],
                                                                               &dst_data[vertex_offset_gl_0],
                                                                               rhs_data,
                                                                               static_cast< int32_t >( level ),
                                                                               neighbor_cell_1_local_vertex_id_0,
                                                                               neighbor_cell_1_local_vertex_id_1,
                                                                               neighbor_cell_1_local_vertex_id_2,
                                                                               neighbor_cell_0_local_vertex_id_0,
                                                                               neighbor_cell_0_local_vertex_id_1,
                                                                               neighbor_cell_0_local_vertex_id_2,
                                                                               relax,
                                                                               stencil[1],
                                                                               stencil[0] );
            }
            else
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1( dst_data,
                                                                     &dst_data[vertex_offset_gl_1],
                                                                     &dst_data[vertex_offset_gl_0],
                                                                     rhs_data,
                                                                     static_cast< int32_t >( level ),
                                                                     neighbor_cell_1_local_vertex_id_0,
                                                                     neighbor_cell_1_local_vertex_id_1,
                                                                     neighbor_cell_1_local_vertex_id_2,
                                                                     neighbor_cell_0_local_vertex_id_0,
                                                                     neighbor_cell_0_local_vertex_id_1,
                                                                     neighbor_cell_0_local_vertex_id_2,
                                                                     relax,
                                                                     stencil[1],
                                                                     stencil[0] );
            }
         }
         else
         {
            if ( backwards )
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1_backwards( dst_data,
                                                                               &dst_data[vertex_offset_gl_0],
                                                                               &dst_data[vertex_offset_gl_1],
                                                                               rhs_data,
                                                                               static_cast< int32_t >( level ),
                                                                               neighbor_cell_0_local_vertex_id_0,
                                                                               neighbor_cell_0_local_vertex_id_1,
                                                                               neighbor_cell_0_local_vertex_id_2,
                                                                               neighbor_cell_1_local_vertex_id_0,
                                                                               neighbor_cell_1_local_vertex_id_1,
                                                                               neighbor_cell_1_local_vertex_id_2,
                                                                               relax,
                                                                               stencil[0],
                                                                               stencil[1] );
            }
            else
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1( dst_data,
                                                                     &dst_data[vertex_offset_gl_0],
                                                                     &dst_data[vertex_offset_gl_1],
                                                                     rhs_data,
                                                                     static_cast< int32_t >( level ),
                                                                     neighbor_cell_0_local_vertex_id_0,
                                                                     neighbor_cell_0_local_vertex_id_1,
                                                                     neighbor_cell_0_local_vertex_id_2,
                                                                     neighbor_cell_1_local_vertex_id_0,
                                                                     neighbor_cell_1_local_vertex_id_1,
                                                                     neighbor_cell_1_local_vertex_id_2,
                                                                     relax,
                                                                     stencil[0],
                                                                     stencil[1] );
            }
         }

         this->timingTree_->stop( "Two-sided" );
      }
   }

   inline void smooth_sor_face_generated( Face&                                                    face,
                                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                          const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
                                          const uint_t&                                            level,
                                          real_t                                                   relax,
                                          const bool&                                              backwards = false ) const
   {
      auto rhs_data = face.getData( rhsId )->getPointer( level );
      auto dst_data = face.getData( dstId )->getPointer( level );
      auto stencil  = face.getData( faceStencilID_ )->getPointer( level );

      if ( backwards )
      {
         vertexdof::macroface::generated::sor_2D_macroface_vertexdof_to_vertexdof_backwards(
             dst_data, rhs_data, stencil, static_cast< int32_t >( level ), relax );
      }
      else
      {
         vertexdof::macroface::generated::sor_2D_macroface_vertexdof_to_vertexdof(
             dst_data, rhs_data, stencil, static_cast< int32_t >( level ), relax );
      }
   }

   inline void smooth_sor_cell_generated( Cell&                                                    cell,
                                          const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
                                          const PrimitiveDataID< FunctionMemory< real_t >, Cell >& rhsId,
                                          const uint_t&                                            level,
                                          real_t                                                   relax,
                                          const bool&                                              backwards = false ) const
   {
      auto  rhs_data = cell.getData( rhsId )->getPointer( level );
      auto  dst_data = cell.getData( dstId )->getPointer( level );
      auto& stencil  = cell.getData( cellStencilID_ )->getData( level );

      if ( backwards )
      {
         vertexdof::macrocell::generated::sor_3D_macrocell_P1_backwards(
             dst_data, rhs_data, static_cast< int32_t >( level ), stencil, relax );
      }
      else
      {
         vertexdof::macrocell::generated::sor_3D_macrocell_P1(
             dst_data, rhs_data, static_cast< int32_t >( level ), stencil, relax );
      }
   }

   inline bool backwards_sor_available() const { return true; }
   inline bool variableStencil() const { return false; }

   // assemble stencils for macro-edges, -faces and -cells
   void assembleStencils()
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
         for ( const auto& it : storage_->getCells() )
         {
            auto  cell          = it.second;
            auto& stencilMemory = cell->getData( cellStencilID_ )->getData( level );
            assemble_variableStencil_cell_init( *cell, level );
            assemble_variableStencil_cell( stencilMemory, 1, 1, 1 );
         }

         for ( auto& it : storage_->getFaces() )
         {
            Face& face          = *it.second;
            auto  stencilMemory = face.getData( faceStencilID_ )->getPointer( level );
            auto& stencilMap    = face.getData( faceStencil3DID_ )->getData( level );

            assemble_variableStencil_face_init( face, level );

            if ( storage_->hasGlobalCells() )
            {
               assemble_variableStencil_face3D( stencilMap, 1, 1 );
               WALBERLA_ASSERT_GREATER( face.getNumNeighborCells(), 0 );
            }
            else
            {
               assemble_variableStencil_face( stencilMemory, 0, 0 );
            }
         }

         for ( auto& it : storage_->getEdges() )
         {
            Edge& edge          = *it.second;
            auto  stencilMemory = edge.getData( edgeStencilID_ )->getPointer( level );
            // auto& stencilMap    = edge.getData(edgeStencil3DID_)->getData(level);

            assemble_variableStencil_edge_init( edge, level );

            // if (storage_->hasGlobalCells()) //! not implemented yet
            // {
            //    assemble_variableStencil_edge3D(stencilMap, 1);
            // }

            assemble_variableStencil_edge( stencilMemory, 1 ); // also assemble old version of 3D stencil
         }
      }
   }
};

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, fenics::NoAssemble > > P1ZeroOperator;

typedef P1ConstantOperator< forms::p1_diffusion_affine_q2 > P1ConstantLaplaceOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::UndefinedAssembly >, true >
    P1DiagonalLaplaceOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_0_otherwise > > P1ConstantEpsilonOperator_11;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_1_otherwise > > P1ConstantEpsilonOperator_12;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_2_otherwise > > P1ConstantEpsilonOperator_21;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_3_otherwise > > P1ConstantEpsilonOperator_22;

typedef P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >
    P1DivxOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >
                                                                                                           P1DivyOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > > P1DivzOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >
    P1DivTxOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >
                                                                                                            P1DivTyOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > > P1DivTzOperator;

typedef P1ConstantOperator< forms::p1_mass_affine_qe >                     P1ConstantMassOperator;
typedef P1ConstantOperator< forms::p1_mass_affine_qe, false, true, false > P1LumpedMassOperator;

typedef P1ConstantOperator< forms::p1_mass_affine_qe, false, true, true > P1LumpedInvMassOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >
    P1PSPGOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise >,
                            true,
                            false,
                            true >
    P1PSPGInvDiagOperator;

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >
    P2ToP1DivxVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >
    P2ToP1DivyVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >
    P2ToP1DivzVertexToVertexConstantOperator;

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >
    P1ToP1DivTxVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >
    P1ToP1DivTyVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >
    P1ToP1DivTzVertexToVertexConstantOperator;

typedef P1ConstantOperator< P1LinearCombinationForm > P1ConstantLinearCombinationOperator;

} // namespace hyteg
