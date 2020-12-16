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

#include <array>

#include "hyteg/Operator.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"

#include "hyteg/p1functionspace/variablestencil/VertexDoFVariableStencil.hpp"
#include "P1Elements.hpp"
#include "hyteg/Stencil.hpp"

#include "core/OpenMP.h"

namespace hyteg {

using walberla::real_t;

template < class P1Form >
class P1Operator : public Operator< P1Function< real_t >, P1Function< real_t >>
{
 public:
   P1Operator(const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel)
      : P1Operator<P1Form>(storage, minLvel, maxLevel, P1Form())
   {}

   P1Operator(const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form)
      : Operator(storage, minLevel, maxLevel), form_(form)
   {}

   ~P1Operator() override = default;

   void apply(const P1Function< real_t >& src, const P1Function< real_t >& dst,
              size_t level, DoFType flag, UpdateType updateType = Replace) const
   {
      WALBERLA_ASSERT_NOT_IDENTICAL(std::addressof(src), std::addressof(dst));

      this->startTiming("Apply");
      src.communicate< Vertex, Edge >(level);
      src.communicate< Edge, Face >(level);
      src.communicate< Face, Cell >(level);

      src.communicate< Cell, Face >(level);
      src.communicate< Face, Edge >(level);
      src.communicate< Edge, Vertex >(level);

      this->timingTree_->start("Macro-Vertex");

      std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();

      for (int i = 0; i < int_c(vertexIDs.size()); i++)
      {
         Vertex& vertex = *this->getStorage()->getVertex(vertexIDs[uint_c(i)]);

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType(vertex.getMeshBoundaryFlag());

         if (testFlag(vertexBC, flag))
         {
            vertexdof::macrovertex::apply< real_t >(
               vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), level, updateType);
         }
      }

      this->timingTree_->stop("Macro-Vertex");

      this->timingTree_->start("Macro-Edge");

      if (level >= 1)
      {
         std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();

         for (int i = 0; i < int_c(edgeIDs.size()); i++)
         {
            Edge& edge = *this->getStorage()->getEdge(edgeIDs[uint_c(i)]);

            const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType(edge.getMeshBoundaryFlag());

            if (testFlag(edgeBC, flag))
            {
               apply_edge(edge, src.getEdgeDataID(), dst.getEdgeDataID(), level, updateType);
            }
         }
      }

      this->timingTree_->stop("Macro-Edge");

      this->timingTree_->start("Macro-Face");

      if (level >= 2)
      {
         std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();

         for (int i = 0; i < int_c(faceIDs.size()); i++)
         {
            Face& face = *this->getStorage()->getFace(faceIDs[uint_c(i)]);

            const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType(face.getMeshBoundaryFlag());

            if (testFlag(faceBC, flag))
            {
               if (storage_->hasGlobalCells())
               {
                  // todo generated kernels for variable operators
                  if (support_generated_kernels_ && hyteg::globalDefines::useGeneratedKernels)
                  {
                     if (face.getNumNeighborCells() == 2)
                     {
                        WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start("Two-sided"); }
                     }
                     else
                     {
                        WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start("One-sided"); }
                     }

                     auto         opr_data    = face.getData(faceStencil3DID_)->getData(level);
                     auto         src_data    = face.getData(src.getFaceDataID())->getPointer(level);
                     auto         dst_data    = face.getData(dst.getFaceDataID())->getPointer(level);
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
                  else
                  {
                     apply_face3D(face, src.getFaceDataID(), dst.getFaceDataID(), level, updateType);
                  }
               }
               else
               {
                  // todo generated kernels for variable operators
                  if (support_generated_kernels_ && hyteg::globalDefines::useGeneratedKernels)
                  {
                     real_t* opr_data = face.getData(faceStencilID_)->getPointer(level);
                     real_t* src_data = face.getData(src.getFaceDataID())->getPointer(level);
                     real_t* dst_data = face.getData(dst.getFaceDataID())->getPointer(level);

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
                  else
                  {
                     apply_face(face, src.getFaceDataID(), dst.getFaceDataID(), level, updateType);
                  }
               }
            }
         }
      }

      this->timingTree_->stop("Macro-Face");

      this->timingTree_->start("Macro-Cell");

      if (level >= 2)
      {
         std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();

         for (int i = 0; i < int_c(cellIDs.size()); i++)
         {
            Cell& cell = *this->getStorage()->getCell(cellIDs[uint_c(i)]);

            const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType(cell.getMeshBoundaryFlag());

            if (testFlag(cellBC, flag))
            {
               // todo generated kernels for variable operators
               if (support_generated_kernels_ && hyteg::globalDefines::useGeneratedKernels)
               {
                  auto    opr_data = cell.getData(cellStencilID_)->getData(level);
                  real_t* src_data = cell.getData(src.getCellDataID())->getPointer(level);
                  real_t* dst_data = cell.getData(dst.getCellDataID())->getPointer(level);

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
               else
               {
                  apply_cell(cell, src.getCellDataID(), dst.getCellDataID(), level, updateType);
               }
            }
         }
      }

      this->timingTree_->stop("Macro-Cell");

      this->stopTiming("Apply");
   }


   // todo implement
   void smooth_gs(const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag) const;

   // todo implement
   void smooth_gs_backwards(const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag) const
   {
      smooth_sor_backwards(dst, rhs, 1.0, level, flag);
   }

   // todo implement
   void smooth_sor(const P1Function< real_t >& dst,
                   const P1Function< real_t >& rhs,
                   real_t                      relax,
                   size_t                      level,
                   DoFType                     flag,
                   const bool&                 backwards = false) const;

   // todo implement
   void smooth_sor_backwards(const P1Function< real_t >& dst,
                             const P1Function< real_t >& rhs,
                             real_t                      relax,
                             size_t                      level,
                             DoFType                     flag) const
   {
      smooth_sor(dst, rhs, relax, level, flag, true);
   }


   // todo implement
   void smooth_jac(const P1Function< real_t >& dst,
                   const P1Function< real_t >& rhs,
                   const P1Function< real_t >& tmp,
                   const real_t&               relax,
                   size_t                      level,
                   DoFType                     flag) const;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues(false); }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() { computeDiagonalOperatorValues(true); }

   std::shared_ptr< P1Function< real_t >> getDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
         diagonalValues_,
         "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function.")
      return diagonalValues_;
   };

   std::shared_ptr< P1Function< real_t >> getInverseDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
         inverseDiagonalValues_,
         "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function.")
      return inverseDiagonalValues_;
   };

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >& getVertexStencilID() const { return vertexStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Edge >& getEdgeStencilID() const { return edgeStencilID_; }

   const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Edge >& getEdgeStencil3DID() const { return edgeStencil3DID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Face >& getFaceStencilID() const { return faceStencilID_; }

   const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face >& getFaceStencil3DID() const { return faceStencil3DID_; }

   const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >&  () const { return cellStencilID_; }

 private:
   void assembleStencils();

   void assembleStencils3D();

   inline void apply_edge(Edge& edge, const PrimitiveDataID<FunctionMemory< real_t >, Edge>& srcId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Edge>& dstId, const uint_t& level, UpdateType update)
   {

      using sD = stencilDirection;
      size_t rowsize = levelinfo::num_microvertices_per_edge(level);

      auto opr_data = edge.getData(edgeStencilID_)->getPointer(level);
      auto src = edge.getData(srcId)->getPointer(level);
      auto dst = edge.getData(dstId)->getPointer(level);

      assemble_stencil_edge_init(edge, level);

      real_t tmp;

      for (size_t i = 1; i < rowsize - 1; ++i)
      {
         assemble_stencil_edge(i);

         const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge(sD::VERTEX_W);
         const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge(sD::VERTEX_C);
         const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge(sD::VERTEX_E);

         const auto dofIdxW = vertexdof::macroedge::indexFromVertex(level, i, sD::VERTEX_W);
         const auto dofIdxC = vertexdof::macroedge::indexFromVertex(level, i, sD::VERTEX_C);
         const auto dofIdxE = vertexdof::macroedge::indexFromVertex(level, i, sD::VERTEX_E);

         tmp = opr_data[ stencilIdxW ] * src[ dofIdxW ]
               + opr_data[ stencilIdxC ] * src[ dofIdxC ]
               + opr_data[ stencilIdxE ] * src[ dofIdxE ];

         for (uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++)
         {
            const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace(sD::VERTEX_W, neighborFace);
            const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace(sD::VERTEX_E, neighborFace);
            const auto stencilWeightW = opr_data[ stencilIdxWNeighborFace ];
            const auto stencilWeightE = opr_data[ stencilIdxENeighborFace ];
            const auto dofIdxWNeighborFace = vertexdof::macroedge::indexFromVertexOnNeighborFace(level, i, neighborFace, sD::VERTEX_W);
            const auto dofIdxENeighborFace = vertexdof::macroedge::indexFromVertexOnNeighborFace(level, i, neighborFace, sD::VERTEX_E);
            tmp += stencilWeightW * src[dofIdxWNeighborFace] + stencilWeightE * src[dofIdxENeighborFace];
         }

         for (uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++)
         {
            const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell(neighborCell, edge.getNumNeighborFaces());
            const auto dofIdx = vertexdof::macroedge::indexFromVertexOnNeighborCell(level, i, neighborCell, edge.getNumNeighborFaces());
            tmp += opr_data[ stencilIdx ] * src[ dofIdx ];
         }

         if (update == Replace)
         {
            dst[vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C)] = tmp;
         }
         else if (update == Add)
         {
            dst[vertexdof::macroedge::indexFromVertex(level, i, stencilDirection::VERTEX_C)] += tmp;
         }
      }
   }

   inline void apply_face3D(Face& face, const PrimitiveDataID<FunctionMemory< real_t >, Face>& srcId,
                            const PrimitiveDataID<FunctionMemory< real_t >, Face>& dstId, const uint_t& level, UpdateType update)
   {
      auto&       opr_data = face.getData(faceStencil3DID_)->getData(level);
      real_t* src      = face.getData(srcId)->getPointer(level);
      real_t* dst      = face.getData(dstId)->getPointer(level);

      assemble_stencil_face_init(face, level);

      for (const auto& idxIt : vertexdof::macroface::Iterator(level, 1))
      {
         assemble_stencil_face(idxIt.x(), idxIt.y());

         real_t tmp = real_c(0);

         for (uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++)
         {
            auto neighborCell = storage.getCell(face.neighborCells().at(neighborCellIdx));
            auto centerIndexInCell =
               vertexdof::macroface::getIndexInNeighboringMacroCell(idxIt, face, neighborCellIdx, storage, level);

            for (auto stencilIt : opr_data[neighborCellIdx])
            {
               auto weight               = stencilIt.second;
               auto leafIndexInMacroCell = centerIndexInCell + stencilIt.first;
               auto leafIndexInMacroFace = vertexdof::macrocell::getIndexInNeighboringMacroFace(
                                              leafIndexInMacroCell, *neighborCell, neighborCell->getLocalFaceID(face.getID()), storage, level);

               uint_t leafArrayIndexInMacroFace;

               if (leafIndexInMacroFace.z() == 0)
               {
                  leafArrayIndexInMacroFace = vertexdof::macroface::index(level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y());
               }
               else
               {
                  WALBERLA_ASSERT_EQUAL(leafIndexInMacroFace.z(), 1);
                  leafArrayIndexInMacroFace = vertexdof::macroface::index(level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y(), neighborCellIdx);
               }

               tmp += weight * src[ leafArrayIndexInMacroFace ];
            }
         }

         if (update == Replace)
         {
            dst[ vertexdof::macroface::index(level, idxIt.x(), idxIt.y()) ] = tmp;
         }
         else if (update == Add)
         {
            dst[ vertexdof::macroface::index(level, idxIt.x(), idxIt.y()) ] += tmp;
         }
      }
   }

   inline void apply_face(Face& face, const PrimitiveDataID<FunctionMemory< real_t >, Face>& srcId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Face>& dstId, const uint_t& level, UpdateType update)
   {
      uint_t rowsize       = levelinfo::num_microvertices_per_edge(level);
      uint_t inner_rowsize = rowsize;

      real_t* opr_data = face.getData(faceStencilID_)->getPointer(level);
      real_t* src      = face.getData(srcId)->getPointer(level);
      real_t* dst      = face.getData(dstId)->getPointer(level);

      assemble_stencil_face_init(face, level);

      real_t tmp = real_c(0);

      for (uint_t j = 1; j < rowsize - 2; ++j)
      {
         for (uint_t i = 1; i < inner_rowsize - 2; ++i)
         {
            assemble_stencil_face(i, j);

            if (face.getNumNeighborCells() == 0)
            {
               static_assert(vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size");
               tmp = real_c(0);

               for (const auto direction : vertexdof::macroface::neighborsWithCenter)
               {
                  tmp += opr_data[vertexdof::stencilIndexFromVertex(direction)] *
                         src[vertexdof::macroface::indexFromVertex(level, i, j, direction)];
               }
            }
            else if (face.getNumNeighborCells() == 1)
            {
               tmp = real_c(0);

               for (const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter)
               {
                  tmp += opr_data[vertexdof::stencilIndexFromVertex(direction)] *
                         src[vertexdof::macroface::indexFromVertex(level, i, j, direction)];
               }
            }
            else if (face.getNumNeighborCells() == 2)
            {
               tmp = real_c(0);

               for (const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter)
               {
                  tmp += opr_data[vertexdof::stencilIndexFromVertex(direction)] *
                         src[vertexdof::macroface::indexFromVertex(level, i, j, direction)];
               }
            }

            WALBERLA_ASSERT_LESS(face.getNumNeighborCells(), 3);

            if (update == Replace)
            {
               dst[vertexdof::macroface::indexFromVertex(level, i, j, stencilDirection::VERTEX_C)] = tmp;
            }
            else
            {
               dst[vertexdof::macroface::indexFromVertex(level, i, j, stencilDirection::VERTEX_C)] += tmp;
            }
         }

         --inner_rowsize;
      }
   }


   inline void apply_cell(Cell& cell, const PrimitiveDataID<FunctionMemory< real_t >, Cell>& srcId,
                          const PrimitiveDataID<FunctionMemory< real_t >, Cell>& dstId, const uint_t& level, UpdateType update)
   {
      typedef stencilDirection sd;

      auto& operatorData     = cell.getData(cellStencilID_)->getData(level);
      const real_t* src = cell.getData(srcId)->getPointer(level);
      real_t* dst = cell.getData(dstId)->getPointer(level);

      assemble_stencil_cell_init(cell, level);

      real_t tmp;

      for (const auto& it : vertexdof::macrocell::Iterator(level, 1))
      {
         const uint_t x = it.x();
         const uint_t y = it.y();
         const uint_t z = it.z();

         assemble_stencil_cell(x, y, z);

         const uint_t centerIdx = vertexdof::macrocell::indexFromVertex(level, x, y, z, sd::VERTEX_C);

         tmp = operatorData.at({0, 0, 0}) * src[ centerIdx ];

         for (const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter)
         {
            const uint_t idx        = vertexdof::macrocell::indexFromVertex(level, x, y, z, neighbor);
            WALBERLA_ASSERT_GREATER(operatorData.count(vertexdof::logicalIndexOffsetFromVertex(neighbor)), 0);
            tmp += operatorData.at(vertexdof::logicalIndexOffsetFromVertex(neighbor)) * src[ idx ];
         }

         if (update == Replace)
         {
            dst[ centerIdx ] = tmp;
         }
         else
         {
            dst[ centerIdx ] += tmp;

         }
      }

   }

 protected:
   /// functions for variable stencil assembly. To be used in ,e.g., Ctor of constant operator, callback functions of variable operator, etc.

   /* Initialize assembly of variable edge stencil.
   */
   inline void assemble_variableStencil_edge_init(Edge& edge, const uint_t level)
   {
      real_t h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));


      // 3D version
      if (storage_->hasGlobalCells())
      {
         //todo
      }
      // 2D version
      else
      {
         // todo factor out common stuff
         x0_    = edge.getCoordinates()[0];
         dx_   = h * edge.getDirection();

         Face*  faceS   = storage_->getFace(edge.neighborFaces()[0]);
         Face*  faceN   = nullptr;

         if (edge.getNumNeighborFaces() == 2)
         {
            faceN   = storage_->getFace(edge.neighborFaces()[1]);
         }

         formS_.setGeometryMap(faceS->getGeometryMap());

         if (faceN) formN_.setGeometryMap(faceN->getGeometryMap());

         stencil_directions_2D_ = stencil::Directions2D(h, edge, faceS, faceN);
      }
   }

   /* assembly of variable edge stencil (requires assemble_variableStencil_edge_init() for appropriate edge and level).
   */
   inline void assemble_variableStencil_edge(real_t* edge_stencil, const uint_t i)
   {
      Point3D x = x0_ + i * dx_;

      // 3D version
      if (storage_->hasGlobalCells())
      {
         //todo
      }
      // 2D version
      else
      {
         // south face
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            formS_, {x, x + stencil_directions_2D_.W, x + stencil_directions_2D_.S}, P1Elements::P1Elements2D::elementSW, edge_stencil);
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            formS_, {x, x + stencil_directions_2D_.S, x + stencil_directions_2D_.SE}, P1Elements::P1Elements2D::elementS, edge_stencil);
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            formS_, {x, x + stencil_directions_2D_.SE, x + stencil_directions_2D_.E}, P1Elements::P1Elements2D::elementSE, edge_stencil);

         // north face
         if (north)
         {
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
               formN_, {x, x + stencil_directions_2D_.E, x + stencil_directions_2D_.N}, P1Elements::P1Elements2D::elementNE, edge_stencil);
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
               formN_, {x, x + stencil_directions_2D_.N, x + stencil_directions_2D_.NW}, P1Elements::P1Elements2D::elementN, edge_stencil);
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
               formN_, {x, x + stencil_directions_2D_.NW, x + stencil_directions_2D_.W}, P1Elements::P1Elements2D::elementNW, edge_stencil);
         }
      }
   }

   /* Initialize assembly of variable face stencil.
   */
   inline void assemble_variableStencil_face_init(Face& face, const uint_t level)
   {
      real_t h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

      // 3D version
      if (storage_->hasGlobalCells())
      {
         //todo
      }
      // 2D version
      else
      {
         // todo factor out common stuff
         x0_ = face.coords[0];
         dx_ = h * (face.coords[1] - face.coords[0]);
         dy_ = h * (face.coords[2] - face.coords[0]);

         form_.setGeometryMap(face.getGeometryMap());

         stencil_directions_2D_ = stencil::Directions2D(h, face);
      }
   }

   /* assembly of variable face stencil (requires assemble_variableStencil_face_init() for appropriate face and level).
   */
   inline void assemble_variableStencil_face(real_t* face_stencil, const unit_t i, const unit_t j)
   {
      Point3D x = x0_ + i * dx_ + j * dy_;

      // 3D version
      if (storage_->hasGlobalCells())
      {
         //todo
      }
      // 2D version
      else
      {
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            form, {x, x + stencil_directions_2D_.W, x + stencil_directions_2D_.S}, P1Elements::P1Elements2D::elementSW, face_stencil);
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            form, {x, x + stencil_directions_2D_.S, x + stencil_directions_2D_.SE}, P1Elements::P1Elements2D::elementS, face_stencil);
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            form, {x, x + stencil_directions_2D_.SE, x + stencil_directions_2D_.E}, P1Elements::P1Elements2D::elementSE, face_stencil);
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            form, {x, x + stencil_directions_2D_.E, x + stencil_directions_2D_.N}, P1Elements::P1Elements2D::elementNE, face_stencil);
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            form, {x, x + stencil_directions_2D_.N, x + stencil_directions_2D_.NW}, P1Elements::P1Elements2D::elementN, face_stencil);
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
            form, {x, x + stencil_directions_2D_.NW, x + stencil_directions_2D_.W}, P1Elements::P1Elements2D::elementNW, face_stencil);
      }
   }

   /* Initialize assembly of variable cell stencil.
   */
   inline void assemble_variableStencil_cell_init(Cell& cell, const uint_t level)
   {
      //todo
   }

   /* assembly of variable cell stencil (requires assemble_variableStencil_cell_init() for appropriate cell and level).
   */
   inline void assemble_variableStencil_cell(const unit_t i, const unit_t j, const unit_t k)
   {
      //todo
   }

   ////////////////////////////////////////////////////////////////////////////////////////////////////////

   /// callback functions for different stencil variants (constant, variable, surrogate, ... ) ///////////

   /* Initialize assembly of variable edge stencil.
      Will be called before iterating over edge whenever the stencil is applied.
   */
   virtual void assemble_stencil_edge_init(Edge& edge, const uint_t level) = 0;

   /* Assembly of edge stencil.
      Will be called before stencil is applied to a particuar DoF.
   */
   virtual void assemble_stencil_edge(const unit_t i) = 0;

   /* Initialize assembly of face stencil.
      Will be called before iterating over face whenever the stencil is applied.
   */
   virtual void assemble_stencil_face_init(Face& face, const uint_t level) = 0;

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar DoF.
   */
   virtual void assemble_stencil_face(const unit_t i, const unit_t j) = 0;

   /* Initialize assembly of cell stencil.
      Will be called before iterating over cell whenever the stencil is applied.
   */
   virtual void assemble_stencil_cell_init(Cell& cell, const uint_t level) = 0;

   /* Assembly of cell stencil.
      Will be called before stencil is applied to a particuar DoF.
   */
   virtual void assemble_stencil_cell(const unit_t i, const unit_t j, const unit_t k) = 0;

   /////////////////////////////////////////////////////////////////////////////////////////////////

 private:

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   // todo implement
   void computeDiagonalOperatorValues(bool invert);

   std::shared_ptr< P1Function< real_t >> diagonalValues_;
   std::shared_ptr< P1Function< real_t >> inverseDiagonalValues_;

   // todo implement
   void smooth_sor_macro_vertices(const P1Function< real_t >& dst,
                                  const P1Function< real_t >& rhs,
                                  real_t                      relax,
                                  size_t                      level,
                                  DoFType                     flag,
                                  const bool&                 backwards = false) const;

   // todo implement
   void smooth_sor_macro_edges(const P1Function< real_t >& dst,
                               const P1Function< real_t >& rhs,
                               real_t                      relax,
                               size_t                      level,
                               DoFType                     flag,
                               const bool&                 backwards = false) const;

   // todo implement
   void smooth_sor_macro_faces(const P1Function< real_t >& dst,
                               const P1Function< real_t >& rhs,
                               real_t                      relax,
                               size_t                      level,
                               DoFType                     flag,
                               const bool&                 backwards = false) const;

   // todo implement
   void smooth_sor_macro_cells(const P1Function< real_t >& dst,
                               const P1Function< real_t >& rhs,
                               real_t                      relax,
                               size_t                      level,
                               DoFType                     flag,
                               const bool&                 backwards = false) const;

 protected:
   // general data for stencil assembly
   Point3D x0_, dx_, dy_, dz_;

   // data for edge stencil assembly
   // todo different namespace?
   stencil::macroedge::Directions2D stencil_directions_2D_;
   P1Form formS_, formN_;

   PrimitiveDataID< StencilMemory< real_t >, Vertex > vertexStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >   edgeStencilID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macroedge::StencilMap_T >, Edge > edgeStencil3DID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >   faceStencilID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face > faceStencil3DID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell > cellStencilID_;

   bool support_generated_kernels_;
   P1Form form_;
};

} // namespace hyteg
