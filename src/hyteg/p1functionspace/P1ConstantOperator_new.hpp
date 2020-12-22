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

#include "hyteg/p1functionspace/P1Operator.hpp"

namespace hyteg {

using walberla::real_t;

template < class P1Form >
class P1ConstantOperator_new : public P1Operator<P1Form>
{
   // todo: remove unneccessary stuff
   using P1Operator<P1Form>::P1Operator;
   using P1Operator<P1Form>::diagonalValues_;
   using P1Operator<P1Form>::inverseDiagonalValues_;
   using P1Operator<P1Form>::x0_;
   using P1Operator<P1Form>::dx_;
   using P1Operator<P1Form>::dy_;
   using P1Operator<P1Form>::dz_;
   using P1Operator<P1Form>::stencil_directions_2D_;
   using P1Operator<P1Form>::formS_;
   using P1Operator<P1Form>::formN_;
   using P1Operator<P1Form>::form_;
   using P1Operator<P1Form>::vertexStencilID_;
   using P1Operator<P1Form>::edgeStencilID_;
   using P1Operator<P1Form>::faceStencilID_;
   using P1Operator<P1Form>::edgeStencil3DID_;
   using P1Operator<P1Form>::faceStencil3DID_;
   using P1Operator<P1Form>::cellStencilID_;

 public:
   // todo Ctor

 protected:

   // todo: implement virtual functions
   /// stencil assembly ///////////

   /* Initialize assembly of variable edge stencil.
      Will be called before iterating over edge whenever the stencil is applied.
   */
   inline void assemble_stencil_edge_init(Edge& edge, const uint_t level)
   {
   }

   /* Assembly of edge stencil.
      Will be called before stencil is applied to a particuar edge-DoF.
      @return true if edge_stencil has been modified, false otherwise
   */
   inline bool assemble_stencil_edge(real_t* edge_stencil, const uint_t i)
   {
   }

   /* Initialize assembly of face stencil.
      Will be called before iterating over face whenever the stencil is applied.
   */
   inline void assemble_stencil_face_init(Face& face, const uint_t level)
   {
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 2d domain.
      @return true if face_stencil has been modified, false otherwise
   */
   inline bool assemble_stencil_face(real_t* face_stencil, const uint_t i, const uint_t j)
   {
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 3D domain.
      @return true if face_stencil has been modified, false otherwise
   */
   inline bool assemble_stencil_face3D(vertexdof::macroface::StencilMap_T& face_stencil, const uint_t i, const uint_t j)
   {
   }

   /* Initialize assembly of cell stencil.
      Will be called before iterating over cell whenever the stencil is applied.
   */
   inline void assemble_stencil_cell_init(Cell& cell, const uint_t level)
   {
   }

   /* Assembly of cell stencil.
      Will be called before stencil is applied to a particuar cell-DoF.
      @return true if cell_stencil has been modified, false otherwise
   */
   inline bool assemble_stencil_cell(vertexdof::macrocell::StencilMap_T& cell_stencil, const uint_t i, const uint_t j, const uint_t k)
   {
   }

   /////////////////////////////////////////////////////////////////////////////////////////////////


   inline void apply_face3D_generated(Face& face, const PrimitiveDataID<FunctionMemory< real_t >, Face>& srcId,
                                      const PrimitiveDataID<FunctionMemory< real_t >, Face>& dstId, const uint_t& level, UpdateType update)
   {
      if (face.getNumNeighborCells() == 2)
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start("Two-sided"); }
      }
      else
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start("One-sided"); }
      }

      auto&         opr_data    = face.getData(faceStencil3DID_)->getData(level);
      auto         src_data    = face.getData(srcId)->getPointer(level);
      auto         dst_data    = face.getData(dstId)->getPointer(level);
      const uint_t offset_gl_0 = levelinfo::num_microvertices_per_face(level);

      auto neighborCell0 = storage_->getCell(face.neighborCells()[0]);

      auto neighbor_cell_0_local_vertex_id_0 =
         static_cast< int32_t >(neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                .at(neighborCell0->getLocalFaceID(face.getID()))
                                .at(0));
      auto neighbor_cell_0_local_vertex_id_1 =
         static_cast< int32_t >(neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                .at(neighborCell0->getLocalFaceID(face.getID()))
                                .at(1));
      auto neighbor_cell_0_local_vertex_id_2 =
         static_cast< int32_t >(neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                .at(neighborCell0->getLocalFaceID(face.getID()))
                                .at(2));

      if (updateType == Replace)
      {
         vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace(
            dst_data,
            src_data,
            &src_data[offset_gl_0],
            static_cast< int32_t >(level),
            neighbor_cell_0_local_vertex_id_0,
            neighbor_cell_0_local_vertex_id_1,
            neighbor_cell_0_local_vertex_id_2,
            opr_data[0]);
      }
      else
      {
         vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add(
            dst_data,
            src_data,
            &src_data[offset_gl_0],
            static_cast< int32_t >(level),
            neighbor_cell_0_local_vertex_id_0,
            neighbor_cell_0_local_vertex_id_1,
            neighbor_cell_0_local_vertex_id_2,
            opr_data[0]);
      }

      if (face.getNumNeighborCells() == 2)
      {
         const uint_t offset_gl_1 = offset_gl_0 + levelinfo::num_microvertices_per_face_from_width(
                                       levelinfo::num_microvertices_per_edge(level) - 1);

         auto neighborCell1 = storage_->getCell(face.neighborCells()[1]);

         auto neighbor_cell_1_local_vertex_id_0 =
            static_cast< int32_t >(neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                   .at(neighborCell1->getLocalFaceID(face.getID()))
                                   .at(0));
         auto neighbor_cell_1_local_vertex_id_1 =
            static_cast< int32_t >(neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                   .at(neighborCell1->getLocalFaceID(face.getID()))
                                   .at(1));
         auto neighbor_cell_1_local_vertex_id_2 =
            static_cast< int32_t >(neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                   .at(neighborCell1->getLocalFaceID(face.getID()))
                                   .at(2));

         vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add(
            dst_data,
            src_data,
            &src_data[offset_gl_1],
            static_cast< int32_t >(level),
            neighbor_cell_1_local_vertex_id_0,
            neighbor_cell_1_local_vertex_id_1,
            neighbor_cell_1_local_vertex_id_2,
            opr_data[1]);
      }

      if (face.getNumNeighborCells() == 2)
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop("Two-sided"); }
      }
      else
      {
         WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop("One-sided"); }
      }
   }

   inline void apply_face_generated(Face& face, const PrimitiveDataID<FunctionMemory< real_t >, Face>& srcId,
                                    const PrimitiveDataID<FunctionMemory< real_t >, Face>& dstId, const uint_t& level, UpdateType update)
   {
      real_t* opr_data = face.getData(faceStencilID_)->getPointer(level);
      real_t* src_data = face.getData(srcId)->getPointer(level);
      real_t* dst_data = face.getData(dstId)->getPointer(level);

      if (updateType == hyteg::Replace)
      {
         vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace(
            dst_data, src_data, opr_data, static_cast< int32_t >(level));
      }
      else if (updateType == hyteg::Add)
      {
         vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
            dst_data, src_data, opr_data, static_cast< int32_t >(level));
      }
   }

   inline void apply_cell_generated(Cell& cell, const PrimitiveDataID<FunctionMemory< real_t >, Cell>& srcId,
                                    const PrimitiveDataID<FunctionMemory< real_t >, Cell>& dstId, const uint_t& level, UpdateType update)
   {
      auto&    opr_data = cell.getData(cellStencilID_)->getData(level);
      real_t* src_data = cell.getData(srcId)->getPointer(level);
      real_t* dst_data = cell.getData(dstId)->getPointer(level);

      if (updateType == Replace)
      {
         vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
            dst_data, src_data, static_cast< int32_t >(level), opr_data);
      }
      else if (updateType == Add)
      {
         vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
            dst_data, src_data, static_cast< int32_t >(level), opr_data);
      }
   }

   inline void smooth_sor_face3D_generated(Face& face, const PrimitiveDataID<FunctionMemory< real_t >, Face>& dstId,
                                           const PrimitiveDataID<FunctionMemory< real_t >, Face>& rhsId, const uint_t& level, real_t relax,
                                           const bool& backwards = false)
   {
      auto rhs_data = face.getData(rhsId)->getPointer(level);
      auto dst_data = face.getData(dstId)->getPointer(level);
      auto& stencil  = face.getData(faceStencil3DID_)->getData(level);

      auto neighborCell0 = storage_->getCell(face.neighborCells()[0]);

      auto neighbor_cell_0_local_vertex_id_0 =
         static_cast< int32_t >(neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                .at(neighborCell0->getLocalFaceID(face.getID()))
                                .at(0));
      auto neighbor_cell_0_local_vertex_id_1 =
         static_cast< int32_t >(neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                .at(neighborCell0->getLocalFaceID(face.getID()))
                                .at(1));
      auto neighbor_cell_0_local_vertex_id_2 =
         static_cast< int32_t >(neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                .at(neighborCell0->getLocalFaceID(face.getID()))
                                .at(2));

      const uint_t vertex_offset_gl_0 = levelinfo::num_microvertices_per_face(level);

      if (face.getNumNeighborCells() == 1)
      {
         this->timingTree_->start("One-sided");

         if (backwards)
         {
            vertexdof::macroface::generated::sor_3D_macroface_P1_one_sided_backwards(dst_data,
                                                                                     &dst_data[vertex_offset_gl_0],
                                                                                     rhs_data,
                                                                                     static_cast< int32_t >(level),
                                                                                     neighbor_cell_0_local_vertex_id_0,
                                                                                     neighbor_cell_0_local_vertex_id_1,
                                                                                     neighbor_cell_0_local_vertex_id_2,
                                                                                     relax,
                                                                                     stencil[0]);
         }
         else
         {
            vertexdof::macroface::generated::sor_3D_macroface_P1_one_sided(dst_data,
                                                                           &dst_data[vertex_offset_gl_0],
                                                                           rhs_data,
                                                                           static_cast< int32_t >(level),
                                                                           neighbor_cell_0_local_vertex_id_0,
                                                                           neighbor_cell_0_local_vertex_id_1,
                                                                           neighbor_cell_0_local_vertex_id_2,
                                                                           relax,
                                                                           stencil[0]);
         }

         this->timingTree_->stop("One-sided");
      }

      if (face.getNumNeighborCells() == 2)
      {
         this->timingTree_->start("Two-sided");

         auto neighborCell1 = storage_->getCell(face.neighborCells()[1]);

         auto neighbor_cell_1_local_vertex_id_0 =
            static_cast< int32_t >(neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                   .at(neighborCell1->getLocalFaceID(face.getID()))
                                   .at(0));
         auto neighbor_cell_1_local_vertex_id_1 =
            static_cast< int32_t >(neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                   .at(neighborCell1->getLocalFaceID(face.getID()))
                                   .at(1));
         auto neighbor_cell_1_local_vertex_id_2 =
            static_cast< int32_t >(neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                   .at(neighborCell1->getLocalFaceID(face.getID()))
                                   .at(2));
         const uint_t vertex_offset_gl_1 = vertex_offset_gl_0 + levelinfo::num_microvertices_per_face_from_width(
                                              levelinfo::num_microvertices_per_edge(level) - 1);

         if (neighbor_cell_0_local_vertex_id_0 > neighbor_cell_1_local_vertex_id_0 ||
               (neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                neighbor_cell_0_local_vertex_id_1 > neighbor_cell_1_local_vertex_id_1) ||
               (neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                neighbor_cell_0_local_vertex_id_1 == neighbor_cell_1_local_vertex_id_1 &&
                neighbor_cell_0_local_vertex_id_2 > neighbor_cell_1_local_vertex_id_2))
         {
            if (backwards)
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1_backwards(dst_data,
                                                                              &dst_data[vertex_offset_gl_1],
                                                                              &dst_data[vertex_offset_gl_0],
                                                                              rhs_data,
                                                                              static_cast< int32_t >(level),
                                                                              neighbor_cell_1_local_vertex_id_0,
                                                                              neighbor_cell_1_local_vertex_id_1,
                                                                              neighbor_cell_1_local_vertex_id_2,
                                                                              neighbor_cell_0_local_vertex_id_0,
                                                                              neighbor_cell_0_local_vertex_id_1,
                                                                              neighbor_cell_0_local_vertex_id_2,
                                                                              relax,
                                                                              stencil[1],
                                                                              stencil[0]);
            }
            else
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1(dst_data,
                                                                    &dst_data[vertex_offset_gl_1],
                                                                    &dst_data[vertex_offset_gl_0],
                                                                    rhs_data,
                                                                    static_cast< int32_t >(level),
                                                                    neighbor_cell_1_local_vertex_id_0,
                                                                    neighbor_cell_1_local_vertex_id_1,
                                                                    neighbor_cell_1_local_vertex_id_2,
                                                                    neighbor_cell_0_local_vertex_id_0,
                                                                    neighbor_cell_0_local_vertex_id_1,
                                                                    neighbor_cell_0_local_vertex_id_2,
                                                                    relax,
                                                                    stencil[1],
                                                                    stencil[0]);
            }
         }
         else
         {
            if (backwards)
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1_backwards(dst_data,
                                                                              &dst_data[vertex_offset_gl_0],
                                                                              &dst_data[vertex_offset_gl_1],
                                                                              rhs_data,
                                                                              static_cast< int32_t >(level),
                                                                              neighbor_cell_0_local_vertex_id_0,
                                                                              neighbor_cell_0_local_vertex_id_1,
                                                                              neighbor_cell_0_local_vertex_id_2,
                                                                              neighbor_cell_1_local_vertex_id_0,
                                                                              neighbor_cell_1_local_vertex_id_1,
                                                                              neighbor_cell_1_local_vertex_id_2,
                                                                              relax,
                                                                              stencil[0],
                                                                              stencil[1]);
            }
            else
            {
               vertexdof::macroface::generated::sor_3D_macroface_P1(dst_data,
                                                                    &dst_data[vertex_offset_gl_0],
                                                                    &dst_data[vertex_offset_gl_1],
                                                                    rhs_data,
                                                                    static_cast< int32_t >(level),
                                                                    neighbor_cell_0_local_vertex_id_0,
                                                                    neighbor_cell_0_local_vertex_id_1,
                                                                    neighbor_cell_0_local_vertex_id_2,
                                                                    neighbor_cell_1_local_vertex_id_0,
                                                                    neighbor_cell_1_local_vertex_id_1,
                                                                    neighbor_cell_1_local_vertex_id_2,
                                                                    relax,
                                                                    stencil[0],
                                                                    stencil[1]);
            }
         }

         this->timingTree_->stop("Two-sided");
      }

   }

   inline void smooth_sor_face_generated(Face& face, const PrimitiveDataID<FunctionMemory< real_t >, Face>& dstId,
                                         const PrimitiveDataID<FunctionMemory< real_t >, Face>& rhsId, const uint_t& level, real_t relax,
                                         const bool& backwards = false)
   {
      auto rhs_data = face.getData(rhsId)->getPointer(level);
      auto dst_data = face.getData(dstId)->getPointer(level);
      auto stencil  = face.getData(faceStencilID_)->getPointer(level);

      if (backwards)
      {
         vertexdof::macroface::generated::sor_2D_macroface_vertexdof_to_vertexdof_backwards(
            dst_data, rhs_data, stencil, static_cast< int32_t >(level), relax);
      }
      else
      {
         vertexdof::macroface::generated::sor_2D_macroface_vertexdof_to_vertexdof(
            dst_data, rhs_data, stencil, static_cast< int32_t >(level), relax);
      }
   }

   inline void smooth_sor_cell_generated(Cell& cell, const PrimitiveDataID<FunctionMemory< real_t >, Cell>& dstId,
                                         const PrimitiveDataID<FunctionMemory< real_t >, Cell>& rhsId, const uint_t& level, real_t relax,
                                         const bool& backwards = false)
   {

      auto rhs_data = cell.getData(rhsId)->getPointer(level);
      auto dst_data = cell.getData(dstId)->getPointer(level);
      auto& stencil  = cell.getData(cellStencilID_)->getData(level);

      if (backwards)
      {
         vertexdof::macrocell::generated::sor_3D_macrocell_P1_backwards(
            dst_data, rhs_data, static_cast< int32_t >(level), stencil, relax);
      }
      else
      {
         vertexdof::macrocell::generated::sor_3D_macrocell_P1(
            dst_data, rhs_data, static_cast< int32_t >(level), stencil, relax);
      }
   }

   // todo move to P1ConstantOperator and implement
   void assembleStencils();
   void assembleStencils3D();
};

} // namespace hyteg
