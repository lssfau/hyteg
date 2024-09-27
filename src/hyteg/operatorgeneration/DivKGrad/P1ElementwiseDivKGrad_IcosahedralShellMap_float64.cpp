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

#include "P1ElementwiseDivKGrad_IcosahedralShellMap_float64.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P1ElementwiseDivKGrad_IcosahedralShellMap_float64::
    P1ElementwiseDivKGrad_IcosahedralShellMap_float64(
        const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
        size_t maxLevel, const P0Function<walberla::float64> &_k)
    : Operator(storage, minLevel, maxLevel), k(_k) {}

void P1ElementwiseDivKGrad_IcosahedralShellMap_float64::apply(
    const P1Function<walberla::float64> &src,
    const P1Function<walberla::float64> &dst, uint_t level, DoFType flag,
    UpdateType updateType) const {
  this->startTiming("apply");

  // Make sure that halos are up-to-date
  this->timingTree_->start("pre-communication");
  if (this->storage_->hasGlobalCells()) {
    // Note that the order of communication is important, since the face -> cell
    // communication may overwrite parts of the halos that carry the
    // macro-vertex and macro-edge unknowns.
    src.communicate<Face, Cell>(level);
    src.communicate<Edge, Cell>(level);
    src.communicate<Vertex, Cell>(level);

  } else {
    WALBERLA_ABORT("Not implemented.");
  }
  this->timingTree_->stop("pre-communication");

  if (updateType == Replace) {
    // We need to zero the destination array (including halos).
    // However, we must not zero out anything that is not flagged with the
    // specified BCs. Therefore, we first zero out everything that flagged, and
    // then, later, the halos of the highest dim primitives.
    dst.interpolate(walberla::numeric_cast<walberla::float64>(0), level, flag);
  }

  if (storage_->hasGlobalCells()) {
    for (auto &it : storage_->getCells()) {
      Cell &cell = *it.second;

      // get hold of the actual numerical data in the functions
      walberla::float64 *_data_src =
          cell.getData(src.getCellDataID())->getPointer(level);
      walberla::float64 *_data_dst =
          cell.getData(dst.getCellDataID())->getPointer(level);
      walberla::float64 *_data_k =
          k.getDGFunction()->volumeDoFFunction()->dofMemory(it.first, level);

      // Zero out dst halos only
      //
      // This is also necessary when using update type == Add.
      // During additive comm we then skip zeroing the data on the lower-dim
      // primitives.
      for (const auto &idx : vertexdof::macrocell::Iterator(level)) {
        if (!vertexdof::macrocell::isOnCellFace(idx, level).empty()) {
          auto arrayIdx =
              vertexdof::macrocell::index(level, idx.x(), idx.y(), idx.z());
          _data_dst[arrayIdx] = walberla::float64(0);
        }
      }

      const auto micro_edges_per_macro_edge =
          (int64_t)levelinfo::num_microedges_per_edge(level);
      const auto num_microfaces_per_face =
          (int64_t)levelinfo::num_microfaces_per_face(level);
      const auto micro_edges_per_macro_edge_float =
          (walberla::float64)levelinfo::num_microedges_per_edge(level);
      const walberla::float64 macro_vertex_coord_id_0comp0 =
          (walberla::float64)cell.getCoordinates()[0][0];
      const walberla::float64 macro_vertex_coord_id_0comp1 =
          (walberla::float64)cell.getCoordinates()[0][1];
      const walberla::float64 macro_vertex_coord_id_0comp2 =
          (walberla::float64)cell.getCoordinates()[0][2];
      const walberla::float64 macro_vertex_coord_id_1comp0 =
          (walberla::float64)cell.getCoordinates()[1][0];
      const walberla::float64 macro_vertex_coord_id_1comp1 =
          (walberla::float64)cell.getCoordinates()[1][1];
      const walberla::float64 macro_vertex_coord_id_1comp2 =
          (walberla::float64)cell.getCoordinates()[1][2];
      const walberla::float64 macro_vertex_coord_id_2comp0 =
          (walberla::float64)cell.getCoordinates()[2][0];
      const walberla::float64 macro_vertex_coord_id_2comp1 =
          (walberla::float64)cell.getCoordinates()[2][1];
      const walberla::float64 macro_vertex_coord_id_2comp2 =
          (walberla::float64)cell.getCoordinates()[2][2];
      const walberla::float64 macro_vertex_coord_id_3comp0 =
          (walberla::float64)cell.getCoordinates()[3][0];
      const walberla::float64 macro_vertex_coord_id_3comp1 =
          (walberla::float64)cell.getCoordinates()[3][1];
      const walberla::float64 macro_vertex_coord_id_3comp2 =
          (walberla::float64)cell.getCoordinates()[3][2];
      WALBERLA_CHECK_NOT_NULLPTR(
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap()),
          "This operator requires the IcosahedralShellMap to be registered as "
          "GeometryMap on every macro-cell.")
      real_t radRefVertex =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->radRefVertex();
      real_t radRayVertex =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->radRayVertex();
      real_t refVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->refVertex()[0];
      real_t rayVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->rayVertex()[0];
      real_t thrVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->thrVertex()[0];
      real_t forVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->forVertex()[0];
      real_t refVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->refVertex()[1];
      real_t rayVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->rayVertex()[1];
      real_t thrVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->thrVertex()[1];
      real_t forVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->forVertex()[1];
      real_t refVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->refVertex()[2];
      real_t rayVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->rayVertex()[2];
      real_t thrVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->thrVertex()[2];
      real_t forVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->forVertex()[2];

      this->timingTree_->start("kernel");

      apply_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(

          _data_dst, _data_k, _data_src, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_0comp2,
          macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
          macro_vertex_coord_id_1comp2, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, macro_vertex_coord_id_2comp2,
          macro_vertex_coord_id_3comp0, macro_vertex_coord_id_3comp1,
          macro_vertex_coord_id_3comp2, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float);

      this->timingTree_->stop("kernel");
    }

    // Push result to lower-dimensional primitives
    //
    this->timingTree_->start("post-communication");
    // Note: We could avoid communication here by implementing the apply() also
    // for the respective
    //       lower dimensional primitives!
    dst.communicateAdditively<Cell, Face>(level, DoFType::All ^ flag, *storage_,
                                          updateType == Replace);
    dst.communicateAdditively<Cell, Edge>(level, DoFType::All ^ flag, *storage_,
                                          updateType == Replace);
    dst.communicateAdditively<Cell, Vertex>(level, DoFType::All ^ flag,
                                            *storage_, updateType == Replace);
    this->timingTree_->stop("post-communication");
  } else {
    WALBERLA_ABORT("Not implemented.");
  }

  this->stopTiming("apply");
}
void P1ElementwiseDivKGrad_IcosahedralShellMap_float64::toMatrix(
    const std::shared_ptr<SparseMatrixProxy> &mat, const P1Function<idx_t> &src,
    const P1Function<idx_t> &dst, uint_t level, DoFType flag) const {
  this->startTiming("toMatrix");

  // We currently ignore the flag provided!
  if (flag != All) {
    WALBERLA_LOG_WARNING_ON_ROOT(
        "Input flag ignored in toMatrix; using flag = All");
  }

  if (storage_->hasGlobalCells()) {
    this->timingTree_->start("pre-communication");

    this->timingTree_->stop("pre-communication");

    for (auto &it : storage_->getCells()) {
      Cell &cell = *it.second;

      // get hold of the actual numerical data
      idx_t *_data_src = cell.getData(src.getCellDataID())->getPointer(level);
      idx_t *_data_dst = cell.getData(dst.getCellDataID())->getPointer(level);
      walberla::float64 *_data_k =
          k.getDGFunction()->volumeDoFFunction()->dofMemory(it.first, level);

      const auto micro_edges_per_macro_edge =
          (int64_t)levelinfo::num_microedges_per_edge(level);
      const auto num_microfaces_per_face =
          (int64_t)levelinfo::num_microfaces_per_face(level);
      const auto micro_edges_per_macro_edge_float =
          (walberla::float64)levelinfo::num_microedges_per_edge(level);
      const walberla::float64 macro_vertex_coord_id_0comp0 =
          (walberla::float64)cell.getCoordinates()[0][0];
      const walberla::float64 macro_vertex_coord_id_0comp1 =
          (walberla::float64)cell.getCoordinates()[0][1];
      const walberla::float64 macro_vertex_coord_id_0comp2 =
          (walberla::float64)cell.getCoordinates()[0][2];
      const walberla::float64 macro_vertex_coord_id_1comp0 =
          (walberla::float64)cell.getCoordinates()[1][0];
      const walberla::float64 macro_vertex_coord_id_1comp1 =
          (walberla::float64)cell.getCoordinates()[1][1];
      const walberla::float64 macro_vertex_coord_id_1comp2 =
          (walberla::float64)cell.getCoordinates()[1][2];
      const walberla::float64 macro_vertex_coord_id_2comp0 =
          (walberla::float64)cell.getCoordinates()[2][0];
      const walberla::float64 macro_vertex_coord_id_2comp1 =
          (walberla::float64)cell.getCoordinates()[2][1];
      const walberla::float64 macro_vertex_coord_id_2comp2 =
          (walberla::float64)cell.getCoordinates()[2][2];
      const walberla::float64 macro_vertex_coord_id_3comp0 =
          (walberla::float64)cell.getCoordinates()[3][0];
      const walberla::float64 macro_vertex_coord_id_3comp1 =
          (walberla::float64)cell.getCoordinates()[3][1];
      const walberla::float64 macro_vertex_coord_id_3comp2 =
          (walberla::float64)cell.getCoordinates()[3][2];
      WALBERLA_CHECK_NOT_NULLPTR(
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap()),
          "This operator requires the IcosahedralShellMap to be registered as "
          "GeometryMap on every macro-cell.")
      real_t radRefVertex =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->radRefVertex();
      real_t radRayVertex =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->radRayVertex();
      real_t refVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->refVertex()[0];
      real_t rayVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->rayVertex()[0];
      real_t thrVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->thrVertex()[0];
      real_t forVertex_0 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->forVertex()[0];
      real_t refVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->refVertex()[1];
      real_t rayVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->rayVertex()[1];
      real_t thrVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->thrVertex()[1];
      real_t forVertex_1 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->forVertex()[1];
      real_t refVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->refVertex()[2];
      real_t rayVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->rayVertex()[2];
      real_t thrVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->thrVertex()[2];
      real_t forVertex_2 =
          std::dynamic_pointer_cast<IcosahedralShellMap>(cell.getGeometryMap())
              ->forVertex()[2];

      this->timingTree_->start("kernel");

      toMatrix_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(

          _data_dst, _data_k, _data_src, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_0comp2,
          macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
          macro_vertex_coord_id_1comp2, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, macro_vertex_coord_id_2comp2,
          macro_vertex_coord_id_3comp0, macro_vertex_coord_id_3comp1,
          macro_vertex_coord_id_3comp2, mat, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float);

      this->timingTree_->stop("kernel");
    }
  } else {
    this->timingTree_->start("pre-communication");

    this->timingTree_->stop("pre-communication");

    WALBERLA_ABORT("Not implemented.");
  }
  this->stopTiming("toMatrix");
}
void P1ElementwiseDivKGrad_IcosahedralShellMap_float64::
    computeInverseDiagonalOperatorValues() {
  this->startTiming("computeInverseDiagonalOperatorValues");

  if (invDiag_ == nullptr) {
    invDiag_ = std::make_shared<P1Function<walberla::float64>>(
        "inverse diagonal entries", storage_, minLevel_, maxLevel_);
  }

  for (uint_t level = minLevel_; level <= maxLevel_; level++) {
    invDiag_->setToZero(level);

    if (storage_->hasGlobalCells()) {
      this->timingTree_->start("pre-communication");

      this->timingTree_->stop("pre-communication");

      for (auto &it : storage_->getCells()) {
        Cell &cell = *it.second;

        // get hold of the actual numerical data
        walberla::float64 *_data_invDiag_ =
            cell.getData((*invDiag_).getCellDataID())->getPointer(level);
        walberla::float64 *_data_k =
            k.getDGFunction()->volumeDoFFunction()->dofMemory(it.first, level);

        const auto micro_edges_per_macro_edge =
            (int64_t)levelinfo::num_microedges_per_edge(level);
        const auto num_microfaces_per_face =
            (int64_t)levelinfo::num_microfaces_per_face(level);
        const auto micro_edges_per_macro_edge_float =
            (walberla::float64)levelinfo::num_microedges_per_edge(level);
        const walberla::float64 macro_vertex_coord_id_0comp0 =
            (walberla::float64)cell.getCoordinates()[0][0];
        const walberla::float64 macro_vertex_coord_id_0comp1 =
            (walberla::float64)cell.getCoordinates()[0][1];
        const walberla::float64 macro_vertex_coord_id_0comp2 =
            (walberla::float64)cell.getCoordinates()[0][2];
        const walberla::float64 macro_vertex_coord_id_1comp0 =
            (walberla::float64)cell.getCoordinates()[1][0];
        const walberla::float64 macro_vertex_coord_id_1comp1 =
            (walberla::float64)cell.getCoordinates()[1][1];
        const walberla::float64 macro_vertex_coord_id_1comp2 =
            (walberla::float64)cell.getCoordinates()[1][2];
        const walberla::float64 macro_vertex_coord_id_2comp0 =
            (walberla::float64)cell.getCoordinates()[2][0];
        const walberla::float64 macro_vertex_coord_id_2comp1 =
            (walberla::float64)cell.getCoordinates()[2][1];
        const walberla::float64 macro_vertex_coord_id_2comp2 =
            (walberla::float64)cell.getCoordinates()[2][2];
        const walberla::float64 macro_vertex_coord_id_3comp0 =
            (walberla::float64)cell.getCoordinates()[3][0];
        const walberla::float64 macro_vertex_coord_id_3comp1 =
            (walberla::float64)cell.getCoordinates()[3][1];
        const walberla::float64 macro_vertex_coord_id_3comp2 =
            (walberla::float64)cell.getCoordinates()[3][2];
        WALBERLA_CHECK_NOT_NULLPTR(
            std::dynamic_pointer_cast<IcosahedralShellMap>(
                cell.getGeometryMap()),
            "This operator requires the IcosahedralShellMap to be registered "
            "as GeometryMap on every macro-cell.")
        real_t radRefVertex = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                  cell.getGeometryMap())
                                  ->radRefVertex();
        real_t radRayVertex = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                  cell.getGeometryMap())
                                  ->radRayVertex();
        real_t refVertex_0 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->refVertex()[0];
        real_t rayVertex_0 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->rayVertex()[0];
        real_t thrVertex_0 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->thrVertex()[0];
        real_t forVertex_0 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->forVertex()[0];
        real_t refVertex_1 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->refVertex()[1];
        real_t rayVertex_1 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->rayVertex()[1];
        real_t thrVertex_1 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->thrVertex()[1];
        real_t forVertex_1 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->forVertex()[1];
        real_t refVertex_2 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->refVertex()[2];
        real_t rayVertex_2 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->rayVertex()[2];
        real_t thrVertex_2 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->thrVertex()[2];
        real_t forVertex_2 = std::dynamic_pointer_cast<IcosahedralShellMap>(
                                 cell.getGeometryMap())
                                 ->forVertex()[2];

        this->timingTree_->start("kernel");

        computeInverseDiagonalOperatorValues_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(

            _data_invDiag_, _data_k, macro_vertex_coord_id_0comp0,
            macro_vertex_coord_id_0comp1, macro_vertex_coord_id_0comp2,
            macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
            macro_vertex_coord_id_1comp2, macro_vertex_coord_id_2comp0,
            macro_vertex_coord_id_2comp1, macro_vertex_coord_id_2comp2,
            macro_vertex_coord_id_3comp0, macro_vertex_coord_id_3comp1,
            macro_vertex_coord_id_3comp2, micro_edges_per_macro_edge,
            micro_edges_per_macro_edge_float);

        this->timingTree_->stop("kernel");
      }

      // Push result to lower-dimensional primitives
      //
      this->timingTree_->start("post-communication");
      // Note: We could avoid communication here by implementing the apply()
      // also for the respective
      //       lower dimensional primitives!
      (*invDiag_).communicateAdditively<Cell, Face>(level);
      (*invDiag_).communicateAdditively<Cell, Edge>(level);
      (*invDiag_).communicateAdditively<Cell, Vertex>(level);
      this->timingTree_->stop("post-communication");
      (*invDiag_).invertElementwise(level);
    } else {
      this->timingTree_->start("pre-communication");

      this->timingTree_->stop("pre-communication");

      WALBERLA_ABORT("Not implemented.");
      (*invDiag_).invertElementwise(level);
    }
  }

  this->stopTiming("computeInverseDiagonalOperatorValues");
}
std::shared_ptr<P1Function<walberla::float64>>
P1ElementwiseDivKGrad_IcosahedralShellMap_float64::getInverseDiagonalValues()
    const {
  return invDiag_;
}
void P1ElementwiseDivKGrad_IcosahedralShellMap_float64::
    apply_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(
        walberla::float64 *RESTRICT _data_dst,
        walberla::float64 *RESTRICT _data_k,
        walberla::float64 *RESTRICT _data_src,
        walberla::float64 macro_vertex_coord_id_0comp0,
        walberla::float64 macro_vertex_coord_id_0comp1,
        walberla::float64 macro_vertex_coord_id_0comp2,
        walberla::float64 macro_vertex_coord_id_1comp0,
        walberla::float64 macro_vertex_coord_id_1comp1,
        walberla::float64 macro_vertex_coord_id_1comp2,
        walberla::float64 macro_vertex_coord_id_2comp0,
        walberla::float64 macro_vertex_coord_id_2comp1,
        walberla::float64 macro_vertex_coord_id_2comp2,
        walberla::float64 macro_vertex_coord_id_3comp0,
        walberla::float64 macro_vertex_coord_id_3comp1,
        walberla::float64 macro_vertex_coord_id_3comp2,
        int64_t micro_edges_per_macro_edge,
        walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    {
      /* CellType.WHITE_UP */
      const walberla::float64 tmp_coords_jac_0_WHITE_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 p_affine_const_0_0_WHITE_UP =
          macro_vertex_coord_id_0comp0;
      const walberla::float64 p_affine_const_0_1_WHITE_UP =
          macro_vertex_coord_id_0comp1;
      const walberla::float64 p_affine_const_0_2_WHITE_UP =
          macro_vertex_coord_id_0comp2;
      const walberla::float64 p_affine_const_1_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 p_affine_const_1_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 p_affine_const_1_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 p_affine_const_2_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 p_affine_const_2_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 p_affine_const_2_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_3_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 p_affine_const_3_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 p_affine_const_3_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 src_dof_0 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_1 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_2 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_3 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 +
                        ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_3 = src_dof_3 * tmp_kernel_op_0;
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_0 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_1 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_2 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_3 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
          }
    }
    {
      /* CellType.WHITE_DOWN */
      const walberla::float64 tmp_coords_jac_0_WHITE_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_WHITE_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_3_WHITE_DOWN =
          tmp_coords_jac_1_WHITE_DOWN + tmp_coords_jac_2_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_4_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_5_WHITE_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_6_WHITE_DOWN =
          tmp_coords_jac_4_WHITE_DOWN + tmp_coords_jac_5_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_7_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_8_WHITE_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_9_WHITE_DOWN =
          tmp_coords_jac_7_WHITE_DOWN + tmp_coords_jac_8_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_10_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_11_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_12_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_WHITE_DOWN =
          tmp_coords_jac_3_WHITE_DOWN;
      const walberla::float64 p_affine_const_0_1_WHITE_DOWN =
          tmp_coords_jac_6_WHITE_DOWN;
      const walberla::float64 p_affine_const_0_2_WHITE_DOWN =
          tmp_coords_jac_9_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_0_WHITE_DOWN =
          tmp_coords_jac_10_WHITE_DOWN + tmp_coords_jac_2_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_1_WHITE_DOWN =
          tmp_coords_jac_11_WHITE_DOWN + tmp_coords_jac_5_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_2_WHITE_DOWN =
          tmp_coords_jac_12_WHITE_DOWN + tmp_coords_jac_8_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_0_WHITE_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_10_WHITE_DOWN +
          tmp_coords_jac_1_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_1_WHITE_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_11_WHITE_DOWN +
          tmp_coords_jac_4_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_2_WHITE_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_12_WHITE_DOWN +
          tmp_coords_jac_7_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_0_WHITE_DOWN =
          tmp_coords_jac_10_WHITE_DOWN + tmp_coords_jac_3_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_1_WHITE_DOWN =
          tmp_coords_jac_11_WHITE_DOWN + tmp_coords_jac_6_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_2_WHITE_DOWN =
          tmp_coords_jac_12_WHITE_DOWN + tmp_coords_jac_9_WHITE_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 src_dof_0 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_1 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_2 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_3 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 +
                        ctr_1 * (-ctr_2 + micro_edges_per_macro_edge - 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 2) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_3 = src_dof_3 * tmp_kernel_op_0;
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_0 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_1 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_2 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_3 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
          }
    }
    {
      /* CellType.BLUE_UP */
      const walberla::float64 tmp_coords_jac_0_BLUE_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_BLUE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_2_BLUE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_3_BLUE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_4_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_5_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_6_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_0_0_BLUE_UP =
          tmp_coords_jac_1_BLUE_UP;
      const walberla::float64 p_affine_const_0_1_BLUE_UP =
          tmp_coords_jac_2_BLUE_UP;
      const walberla::float64 p_affine_const_0_2_BLUE_UP =
          tmp_coords_jac_3_BLUE_UP;
      const walberla::float64 p_affine_const_1_0_BLUE_UP =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_BLUE_UP;
      const walberla::float64 p_affine_const_1_1_BLUE_UP =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_BLUE_UP;
      const walberla::float64 p_affine_const_1_2_BLUE_UP =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_BLUE_UP;
      const walberla::float64 p_affine_const_2_0_BLUE_UP =
          tmp_coords_jac_1_BLUE_UP + tmp_coords_jac_4_BLUE_UP;
      const walberla::float64 p_affine_const_2_1_BLUE_UP =
          tmp_coords_jac_2_BLUE_UP + tmp_coords_jac_5_BLUE_UP;
      const walberla::float64 p_affine_const_2_2_BLUE_UP =
          tmp_coords_jac_3_BLUE_UP + tmp_coords_jac_6_BLUE_UP;
      const walberla::float64 p_affine_const_3_0_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0) +
          tmp_coords_jac_1_BLUE_UP;
      const walberla::float64 p_affine_const_3_1_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1) +
          tmp_coords_jac_2_BLUE_UP;
      const walberla::float64 p_affine_const_3_2_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2) +
          tmp_coords_jac_3_BLUE_UP;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 src_dof_0 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_1 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_2 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_3 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 1) *
                          (micro_edges_per_macro_edge + 1)) /
                         (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_3 = src_dof_3 * tmp_kernel_op_0;
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_0 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_1 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_2 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_3 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
          }
    }
    {
      /* CellType.BLUE_DOWN */
      const walberla::float64 tmp_coords_jac_0_BLUE_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_BLUE_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_BLUE_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_3_BLUE_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_4_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_5_BLUE_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_BLUE_DOWN;
      const walberla::float64 tmp_coords_jac_6_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_7_BLUE_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_6_BLUE_DOWN;
      const walberla::float64 tmp_coords_jac_8_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 tmp_coords_jac_9_BLUE_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_8_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_0_BLUE_DOWN =
          tmp_coords_jac_1_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_1_BLUE_DOWN =
          tmp_coords_jac_2_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_2_BLUE_DOWN =
          tmp_coords_jac_3_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_0_BLUE_DOWN =
          tmp_coords_jac_5_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_1_BLUE_DOWN =
          tmp_coords_jac_7_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_2_BLUE_DOWN =
          tmp_coords_jac_9_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_0_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0) +
          tmp_coords_jac_5_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_1_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1) +
          tmp_coords_jac_7_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_2_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2) +
          tmp_coords_jac_9_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_0_BLUE_DOWN =
          tmp_coords_jac_1_BLUE_DOWN + tmp_coords_jac_4_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_1_BLUE_DOWN =
          tmp_coords_jac_2_BLUE_DOWN + tmp_coords_jac_6_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_2_BLUE_DOWN =
          tmp_coords_jac_3_BLUE_DOWN + tmp_coords_jac_8_BLUE_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 src_dof_0 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_1 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_2 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_3 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        3 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_3 = src_dof_3 * tmp_kernel_op_0;
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_0 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_1 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_2 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_3 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
          }
    }
    {
      /* CellType.GREEN_UP */
      const walberla::float64 tmp_coords_jac_0_GREEN_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_GREEN_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_2_GREEN_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_3_GREEN_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_4_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_5_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_6_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_GREEN_UP =
          tmp_coords_jac_1_GREEN_UP;
      const walberla::float64 p_affine_const_0_1_GREEN_UP =
          tmp_coords_jac_2_GREEN_UP;
      const walberla::float64 p_affine_const_0_2_GREEN_UP =
          tmp_coords_jac_3_GREEN_UP;
      const walberla::float64 p_affine_const_1_0_GREEN_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 p_affine_const_1_1_GREEN_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 p_affine_const_1_2_GREEN_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_2_0_GREEN_UP =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_GREEN_UP;
      const walberla::float64 p_affine_const_2_1_GREEN_UP =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_GREEN_UP;
      const walberla::float64 p_affine_const_2_2_GREEN_UP =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_GREEN_UP;
      const walberla::float64 p_affine_const_3_0_GREEN_UP =
          tmp_coords_jac_1_GREEN_UP + tmp_coords_jac_4_GREEN_UP;
      const walberla::float64 p_affine_const_3_1_GREEN_UP =
          tmp_coords_jac_2_GREEN_UP + tmp_coords_jac_5_GREEN_UP;
      const walberla::float64 p_affine_const_3_2_GREEN_UP =
          tmp_coords_jac_3_GREEN_UP + tmp_coords_jac_6_GREEN_UP;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 src_dof_0 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_1 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_2 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_3 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_3 = src_dof_3 * tmp_kernel_op_0;
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_0 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_1 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_2 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_3 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
          }
    }
    {
      /* CellType.GREEN_DOWN */
      const walberla::float64 tmp_coords_jac_0_GREEN_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_GREEN_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_GREEN_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_3_GREEN_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_4_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_5_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_6_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_7_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_8_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_9_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN;
      const walberla::float64 p_affine_const_0_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN;
      const walberla::float64 p_affine_const_0_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN + tmp_coords_jac_4_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN + tmp_coords_jac_5_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN + tmp_coords_jac_6_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_0_GREEN_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_GREEN_DOWN +
          tmp_coords_jac_7_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_1_GREEN_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_GREEN_DOWN +
          tmp_coords_jac_8_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_2_GREEN_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_GREEN_DOWN +
          tmp_coords_jac_9_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN + tmp_coords_jac_7_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN + tmp_coords_jac_8_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN + tmp_coords_jac_9_GREEN_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 src_dof_0 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 src_dof_1 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_2 =
                _data_src[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            const walberla::float64 src_dof_3 =
                _data_src[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        4 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
            const walberla::float64 elMatVec_3 = src_dof_3 * tmp_kernel_op_0;
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_0 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6))];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_1 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) -
                          (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2) *
                            (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) +
                      1] =
                elMatVec_2 +
                _data_dst[ctr_0 +
                          ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6)) +
                          1];
            _data_dst[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6))] =
                elMatVec_3 +
                _data_dst[ctr_0 +
                          (ctr_1 + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                          (((-ctr_2 + micro_edges_per_macro_edge) *
                            (-ctr_2 + micro_edges_per_macro_edge + 1) *
                            (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                           (6)) +
                          (((micro_edges_per_macro_edge + 1) *
                            (micro_edges_per_macro_edge + 2) *
                            (micro_edges_per_macro_edge + 3)) /
                           (6))];
          }
    }
  }
}
void P1ElementwiseDivKGrad_IcosahedralShellMap_float64::
    toMatrix_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(
        idx_t *RESTRICT _data_dst, walberla::float64 *RESTRICT _data_k,
        idx_t *RESTRICT _data_src,
        walberla::float64 macro_vertex_coord_id_0comp0,
        walberla::float64 macro_vertex_coord_id_0comp1,
        walberla::float64 macro_vertex_coord_id_0comp2,
        walberla::float64 macro_vertex_coord_id_1comp0,
        walberla::float64 macro_vertex_coord_id_1comp1,
        walberla::float64 macro_vertex_coord_id_1comp2,
        walberla::float64 macro_vertex_coord_id_2comp0,
        walberla::float64 macro_vertex_coord_id_2comp1,
        walberla::float64 macro_vertex_coord_id_2comp2,
        walberla::float64 macro_vertex_coord_id_3comp0,
        walberla::float64 macro_vertex_coord_id_3comp1,
        walberla::float64 macro_vertex_coord_id_3comp2,
        std::shared_ptr<SparseMatrixProxy> mat,
        int64_t micro_edges_per_macro_edge,
        walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    {
      /* CellType.WHITE_UP */
      const walberla::float64 tmp_coords_jac_0_WHITE_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 p_affine_const_0_0_WHITE_UP =
          macro_vertex_coord_id_0comp0;
      const walberla::float64 p_affine_const_0_1_WHITE_UP =
          macro_vertex_coord_id_0comp1;
      const walberla::float64 p_affine_const_0_2_WHITE_UP =
          macro_vertex_coord_id_0comp2;
      const walberla::float64 p_affine_const_1_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 p_affine_const_1_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 p_affine_const_1_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 p_affine_const_2_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 p_affine_const_2_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 p_affine_const_2_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_3_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 p_affine_const_3_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 p_affine_const_3_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 +
                        ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
            const int64_t elMat_0_1 = 0;
            const int64_t elMat_0_2 = 0;
            const int64_t elMat_0_3 = 0;
            const int64_t elMat_1_0 = 0;
            const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
            const int64_t elMat_1_2 = 0;
            const int64_t elMat_1_3 = 0;
            const int64_t elMat_2_0 = 0;
            const int64_t elMat_2_1 = 0;
            const walberla::float64 elMat_2_2 = tmp_kernel_op_0;
            const int64_t elMat_2_3 = 0;
            const int64_t elMat_3_0 = 0;
            const int64_t elMat_3_1 = 0;
            const int64_t elMat_3_2 = 0;
            const walberla::float64 elMat_3_3 = tmp_kernel_op_0;

            std::vector<uint_t> _data_rowIdx(4);
            std::vector<uint_t> _data_colIdx(4);
            std::vector<real_t> _data_mat(16);

            _data_rowIdx[0] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[1] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[2] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[3] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[0] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[1] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[2] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[3] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));

            /* Apply basis transformation */

            _data_mat[0] = ((real_t)(elMat_0_0));
            _data_mat[1] = ((real_t)(elMat_0_1));
            _data_mat[2] = ((real_t)(elMat_0_2));
            _data_mat[3] = ((real_t)(elMat_0_3));
            _data_mat[4] = ((real_t)(elMat_1_0));
            _data_mat[5] = ((real_t)(elMat_1_1));
            _data_mat[6] = ((real_t)(elMat_1_2));
            _data_mat[7] = ((real_t)(elMat_1_3));
            _data_mat[8] = ((real_t)(elMat_2_0));
            _data_mat[9] = ((real_t)(elMat_2_1));
            _data_mat[10] = ((real_t)(elMat_2_2));
            _data_mat[11] = ((real_t)(elMat_2_3));
            _data_mat[12] = ((real_t)(elMat_3_0));
            _data_mat[13] = ((real_t)(elMat_3_1));
            _data_mat[14] = ((real_t)(elMat_3_2));
            _data_mat[15] = ((real_t)(elMat_3_3));

            mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
          }
    }
    {
      /* CellType.WHITE_DOWN */
      const walberla::float64 tmp_coords_jac_0_WHITE_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_WHITE_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_3_WHITE_DOWN =
          tmp_coords_jac_1_WHITE_DOWN + tmp_coords_jac_2_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_4_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_5_WHITE_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_6_WHITE_DOWN =
          tmp_coords_jac_4_WHITE_DOWN + tmp_coords_jac_5_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_7_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_8_WHITE_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_9_WHITE_DOWN =
          tmp_coords_jac_7_WHITE_DOWN + tmp_coords_jac_8_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_10_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_11_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_12_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_WHITE_DOWN =
          tmp_coords_jac_3_WHITE_DOWN;
      const walberla::float64 p_affine_const_0_1_WHITE_DOWN =
          tmp_coords_jac_6_WHITE_DOWN;
      const walberla::float64 p_affine_const_0_2_WHITE_DOWN =
          tmp_coords_jac_9_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_0_WHITE_DOWN =
          tmp_coords_jac_10_WHITE_DOWN + tmp_coords_jac_2_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_1_WHITE_DOWN =
          tmp_coords_jac_11_WHITE_DOWN + tmp_coords_jac_5_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_2_WHITE_DOWN =
          tmp_coords_jac_12_WHITE_DOWN + tmp_coords_jac_8_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_0_WHITE_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_10_WHITE_DOWN +
          tmp_coords_jac_1_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_1_WHITE_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_11_WHITE_DOWN +
          tmp_coords_jac_4_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_2_WHITE_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_12_WHITE_DOWN +
          tmp_coords_jac_7_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_0_WHITE_DOWN =
          tmp_coords_jac_10_WHITE_DOWN + tmp_coords_jac_3_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_1_WHITE_DOWN =
          tmp_coords_jac_11_WHITE_DOWN + tmp_coords_jac_6_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_2_WHITE_DOWN =
          tmp_coords_jac_12_WHITE_DOWN + tmp_coords_jac_9_WHITE_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 +
                        ctr_1 * (-ctr_2 + micro_edges_per_macro_edge - 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 2) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
            const int64_t elMat_0_1 = 0;
            const int64_t elMat_0_2 = 0;
            const int64_t elMat_0_3 = 0;
            const int64_t elMat_1_0 = 0;
            const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
            const int64_t elMat_1_2 = 0;
            const int64_t elMat_1_3 = 0;
            const int64_t elMat_2_0 = 0;
            const int64_t elMat_2_1 = 0;
            const walberla::float64 elMat_2_2 = tmp_kernel_op_0;
            const int64_t elMat_2_3 = 0;
            const int64_t elMat_3_0 = 0;
            const int64_t elMat_3_1 = 0;
            const int64_t elMat_3_2 = 0;
            const walberla::float64 elMat_3_3 = tmp_kernel_op_0;

            std::vector<uint_t> _data_rowIdx(4);
            std::vector<uint_t> _data_colIdx(4);
            std::vector<real_t> _data_mat(16);

            _data_rowIdx[0] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[1] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[2] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[3] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[0] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[1] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[2] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[3] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));

            /* Apply basis transformation */

            _data_mat[0] = ((real_t)(elMat_0_0));
            _data_mat[1] = ((real_t)(elMat_0_1));
            _data_mat[2] = ((real_t)(elMat_0_2));
            _data_mat[3] = ((real_t)(elMat_0_3));
            _data_mat[4] = ((real_t)(elMat_1_0));
            _data_mat[5] = ((real_t)(elMat_1_1));
            _data_mat[6] = ((real_t)(elMat_1_2));
            _data_mat[7] = ((real_t)(elMat_1_3));
            _data_mat[8] = ((real_t)(elMat_2_0));
            _data_mat[9] = ((real_t)(elMat_2_1));
            _data_mat[10] = ((real_t)(elMat_2_2));
            _data_mat[11] = ((real_t)(elMat_2_3));
            _data_mat[12] = ((real_t)(elMat_3_0));
            _data_mat[13] = ((real_t)(elMat_3_1));
            _data_mat[14] = ((real_t)(elMat_3_2));
            _data_mat[15] = ((real_t)(elMat_3_3));

            mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
          }
    }
    {
      /* CellType.BLUE_UP */
      const walberla::float64 tmp_coords_jac_0_BLUE_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_BLUE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_2_BLUE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_3_BLUE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_4_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_5_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_6_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_0_0_BLUE_UP =
          tmp_coords_jac_1_BLUE_UP;
      const walberla::float64 p_affine_const_0_1_BLUE_UP =
          tmp_coords_jac_2_BLUE_UP;
      const walberla::float64 p_affine_const_0_2_BLUE_UP =
          tmp_coords_jac_3_BLUE_UP;
      const walberla::float64 p_affine_const_1_0_BLUE_UP =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_BLUE_UP;
      const walberla::float64 p_affine_const_1_1_BLUE_UP =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_BLUE_UP;
      const walberla::float64 p_affine_const_1_2_BLUE_UP =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_BLUE_UP;
      const walberla::float64 p_affine_const_2_0_BLUE_UP =
          tmp_coords_jac_1_BLUE_UP + tmp_coords_jac_4_BLUE_UP;
      const walberla::float64 p_affine_const_2_1_BLUE_UP =
          tmp_coords_jac_2_BLUE_UP + tmp_coords_jac_5_BLUE_UP;
      const walberla::float64 p_affine_const_2_2_BLUE_UP =
          tmp_coords_jac_3_BLUE_UP + tmp_coords_jac_6_BLUE_UP;
      const walberla::float64 p_affine_const_3_0_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0) +
          tmp_coords_jac_1_BLUE_UP;
      const walberla::float64 p_affine_const_3_1_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1) +
          tmp_coords_jac_2_BLUE_UP;
      const walberla::float64 p_affine_const_3_2_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2) +
          tmp_coords_jac_3_BLUE_UP;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 1) *
                          (micro_edges_per_macro_edge + 1)) /
                         (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
            const int64_t elMat_0_1 = 0;
            const int64_t elMat_0_2 = 0;
            const int64_t elMat_0_3 = 0;
            const int64_t elMat_1_0 = 0;
            const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
            const int64_t elMat_1_2 = 0;
            const int64_t elMat_1_3 = 0;
            const int64_t elMat_2_0 = 0;
            const int64_t elMat_2_1 = 0;
            const walberla::float64 elMat_2_2 = tmp_kernel_op_0;
            const int64_t elMat_2_3 = 0;
            const int64_t elMat_3_0 = 0;
            const int64_t elMat_3_1 = 0;
            const int64_t elMat_3_2 = 0;
            const walberla::float64 elMat_3_3 = tmp_kernel_op_0;

            std::vector<uint_t> _data_rowIdx(4);
            std::vector<uint_t> _data_colIdx(4);
            std::vector<real_t> _data_mat(16);

            _data_rowIdx[0] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[1] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[2] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[3] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[0] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[1] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[2] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[3] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));

            /* Apply basis transformation */

            _data_mat[0] = ((real_t)(elMat_0_0));
            _data_mat[1] = ((real_t)(elMat_0_1));
            _data_mat[2] = ((real_t)(elMat_0_2));
            _data_mat[3] = ((real_t)(elMat_0_3));
            _data_mat[4] = ((real_t)(elMat_1_0));
            _data_mat[5] = ((real_t)(elMat_1_1));
            _data_mat[6] = ((real_t)(elMat_1_2));
            _data_mat[7] = ((real_t)(elMat_1_3));
            _data_mat[8] = ((real_t)(elMat_2_0));
            _data_mat[9] = ((real_t)(elMat_2_1));
            _data_mat[10] = ((real_t)(elMat_2_2));
            _data_mat[11] = ((real_t)(elMat_2_3));
            _data_mat[12] = ((real_t)(elMat_3_0));
            _data_mat[13] = ((real_t)(elMat_3_1));
            _data_mat[14] = ((real_t)(elMat_3_2));
            _data_mat[15] = ((real_t)(elMat_3_3));

            mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
          }
    }
    {
      /* CellType.BLUE_DOWN */
      const walberla::float64 tmp_coords_jac_0_BLUE_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_BLUE_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_BLUE_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_3_BLUE_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_4_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_5_BLUE_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_BLUE_DOWN;
      const walberla::float64 tmp_coords_jac_6_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_7_BLUE_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_6_BLUE_DOWN;
      const walberla::float64 tmp_coords_jac_8_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 tmp_coords_jac_9_BLUE_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_8_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_0_BLUE_DOWN =
          tmp_coords_jac_1_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_1_BLUE_DOWN =
          tmp_coords_jac_2_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_2_BLUE_DOWN =
          tmp_coords_jac_3_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_0_BLUE_DOWN =
          tmp_coords_jac_5_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_1_BLUE_DOWN =
          tmp_coords_jac_7_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_2_BLUE_DOWN =
          tmp_coords_jac_9_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_0_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0) +
          tmp_coords_jac_5_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_1_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1) +
          tmp_coords_jac_7_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_2_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2) +
          tmp_coords_jac_9_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_0_BLUE_DOWN =
          tmp_coords_jac_1_BLUE_DOWN + tmp_coords_jac_4_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_1_BLUE_DOWN =
          tmp_coords_jac_2_BLUE_DOWN + tmp_coords_jac_6_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_2_BLUE_DOWN =
          tmp_coords_jac_3_BLUE_DOWN + tmp_coords_jac_8_BLUE_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        3 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
            const int64_t elMat_0_1 = 0;
            const int64_t elMat_0_2 = 0;
            const int64_t elMat_0_3 = 0;
            const int64_t elMat_1_0 = 0;
            const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
            const int64_t elMat_1_2 = 0;
            const int64_t elMat_1_3 = 0;
            const int64_t elMat_2_0 = 0;
            const int64_t elMat_2_1 = 0;
            const walberla::float64 elMat_2_2 = tmp_kernel_op_0;
            const int64_t elMat_2_3 = 0;
            const int64_t elMat_3_0 = 0;
            const int64_t elMat_3_1 = 0;
            const int64_t elMat_3_2 = 0;
            const walberla::float64 elMat_3_3 = tmp_kernel_op_0;

            std::vector<uint_t> _data_rowIdx(4);
            std::vector<uint_t> _data_colIdx(4);
            std::vector<real_t> _data_mat(16);

            _data_rowIdx[0] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[1] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[2] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[3] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[0] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[1] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[2] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[3] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));

            /* Apply basis transformation */

            _data_mat[0] = ((real_t)(elMat_0_0));
            _data_mat[1] = ((real_t)(elMat_0_1));
            _data_mat[2] = ((real_t)(elMat_0_2));
            _data_mat[3] = ((real_t)(elMat_0_3));
            _data_mat[4] = ((real_t)(elMat_1_0));
            _data_mat[5] = ((real_t)(elMat_1_1));
            _data_mat[6] = ((real_t)(elMat_1_2));
            _data_mat[7] = ((real_t)(elMat_1_3));
            _data_mat[8] = ((real_t)(elMat_2_0));
            _data_mat[9] = ((real_t)(elMat_2_1));
            _data_mat[10] = ((real_t)(elMat_2_2));
            _data_mat[11] = ((real_t)(elMat_2_3));
            _data_mat[12] = ((real_t)(elMat_3_0));
            _data_mat[13] = ((real_t)(elMat_3_1));
            _data_mat[14] = ((real_t)(elMat_3_2));
            _data_mat[15] = ((real_t)(elMat_3_3));

            mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
          }
    }
    {
      /* CellType.GREEN_UP */
      const walberla::float64 tmp_coords_jac_0_GREEN_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_GREEN_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_2_GREEN_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_3_GREEN_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_4_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_5_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_6_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_GREEN_UP =
          tmp_coords_jac_1_GREEN_UP;
      const walberla::float64 p_affine_const_0_1_GREEN_UP =
          tmp_coords_jac_2_GREEN_UP;
      const walberla::float64 p_affine_const_0_2_GREEN_UP =
          tmp_coords_jac_3_GREEN_UP;
      const walberla::float64 p_affine_const_1_0_GREEN_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 p_affine_const_1_1_GREEN_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 p_affine_const_1_2_GREEN_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_2_0_GREEN_UP =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_GREEN_UP;
      const walberla::float64 p_affine_const_2_1_GREEN_UP =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_GREEN_UP;
      const walberla::float64 p_affine_const_2_2_GREEN_UP =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_GREEN_UP;
      const walberla::float64 p_affine_const_3_0_GREEN_UP =
          tmp_coords_jac_1_GREEN_UP + tmp_coords_jac_4_GREEN_UP;
      const walberla::float64 p_affine_const_3_1_GREEN_UP =
          tmp_coords_jac_2_GREEN_UP + tmp_coords_jac_5_GREEN_UP;
      const walberla::float64 p_affine_const_3_2_GREEN_UP =
          tmp_coords_jac_3_GREEN_UP + tmp_coords_jac_6_GREEN_UP;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
            const int64_t elMat_0_1 = 0;
            const int64_t elMat_0_2 = 0;
            const int64_t elMat_0_3 = 0;
            const int64_t elMat_1_0 = 0;
            const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
            const int64_t elMat_1_2 = 0;
            const int64_t elMat_1_3 = 0;
            const int64_t elMat_2_0 = 0;
            const int64_t elMat_2_1 = 0;
            const walberla::float64 elMat_2_2 = tmp_kernel_op_0;
            const int64_t elMat_2_3 = 0;
            const int64_t elMat_3_0 = 0;
            const int64_t elMat_3_1 = 0;
            const int64_t elMat_3_2 = 0;
            const walberla::float64 elMat_3_3 = tmp_kernel_op_0;

            std::vector<uint_t> _data_rowIdx(4);
            std::vector<uint_t> _data_colIdx(4);
            std::vector<real_t> _data_mat(16);

            _data_rowIdx[0] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[1] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[2] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[3] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[0] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[1] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[2] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[3] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));

            /* Apply basis transformation */

            _data_mat[0] = ((real_t)(elMat_0_0));
            _data_mat[1] = ((real_t)(elMat_0_1));
            _data_mat[2] = ((real_t)(elMat_0_2));
            _data_mat[3] = ((real_t)(elMat_0_3));
            _data_mat[4] = ((real_t)(elMat_1_0));
            _data_mat[5] = ((real_t)(elMat_1_1));
            _data_mat[6] = ((real_t)(elMat_1_2));
            _data_mat[7] = ((real_t)(elMat_1_3));
            _data_mat[8] = ((real_t)(elMat_2_0));
            _data_mat[9] = ((real_t)(elMat_2_1));
            _data_mat[10] = ((real_t)(elMat_2_2));
            _data_mat[11] = ((real_t)(elMat_2_3));
            _data_mat[12] = ((real_t)(elMat_3_0));
            _data_mat[13] = ((real_t)(elMat_3_1));
            _data_mat[14] = ((real_t)(elMat_3_2));
            _data_mat[15] = ((real_t)(elMat_3_3));

            mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
          }
    }
    {
      /* CellType.GREEN_DOWN */
      const walberla::float64 tmp_coords_jac_0_GREEN_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_GREEN_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_GREEN_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_3_GREEN_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_4_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_5_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_6_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_7_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_8_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_9_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN;
      const walberla::float64 p_affine_const_0_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN;
      const walberla::float64 p_affine_const_0_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN + tmp_coords_jac_4_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN + tmp_coords_jac_5_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN + tmp_coords_jac_6_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_0_GREEN_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_GREEN_DOWN +
          tmp_coords_jac_7_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_1_GREEN_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_GREEN_DOWN +
          tmp_coords_jac_8_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_2_GREEN_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_GREEN_DOWN +
          tmp_coords_jac_9_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN + tmp_coords_jac_7_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN + tmp_coords_jac_8_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN + tmp_coords_jac_9_GREEN_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        4 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
            const int64_t elMat_0_1 = 0;
            const int64_t elMat_0_2 = 0;
            const int64_t elMat_0_3 = 0;
            const int64_t elMat_1_0 = 0;
            const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
            const int64_t elMat_1_2 = 0;
            const int64_t elMat_1_3 = 0;
            const int64_t elMat_2_0 = 0;
            const int64_t elMat_2_1 = 0;
            const walberla::float64 elMat_2_2 = tmp_kernel_op_0;
            const int64_t elMat_2_3 = 0;
            const int64_t elMat_3_0 = 0;
            const int64_t elMat_3_1 = 0;
            const int64_t elMat_3_2 = 0;
            const walberla::float64 elMat_3_3 = tmp_kernel_op_0;

            std::vector<uint_t> _data_rowIdx(4);
            std::vector<uint_t> _data_colIdx(4);
            std::vector<real_t> _data_mat(16);

            _data_rowIdx[0] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_rowIdx[1] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[2] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_rowIdx[3] =
                ((uint64_t)(_data_dst
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[0] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6))]));
            _data_colIdx[1] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 2) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[2] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 ctr_1 *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6)) +
                                 1]));
            _data_colIdx[3] =
                ((uint64_t)(_data_src
                                [ctr_0 +
                                 (ctr_1 + 1) *
                                     (-ctr_2 + micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                                 (((-ctr_2 + micro_edges_per_macro_edge) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                  (6)) +
                                 (((micro_edges_per_macro_edge + 1) *
                                   (micro_edges_per_macro_edge + 2) *
                                   (micro_edges_per_macro_edge + 3)) /
                                  (6))]));

            /* Apply basis transformation */

            _data_mat[0] = ((real_t)(elMat_0_0));
            _data_mat[1] = ((real_t)(elMat_0_1));
            _data_mat[2] = ((real_t)(elMat_0_2));
            _data_mat[3] = ((real_t)(elMat_0_3));
            _data_mat[4] = ((real_t)(elMat_1_0));
            _data_mat[5] = ((real_t)(elMat_1_1));
            _data_mat[6] = ((real_t)(elMat_1_2));
            _data_mat[7] = ((real_t)(elMat_1_3));
            _data_mat[8] = ((real_t)(elMat_2_0));
            _data_mat[9] = ((real_t)(elMat_2_1));
            _data_mat[10] = ((real_t)(elMat_2_2));
            _data_mat[11] = ((real_t)(elMat_2_3));
            _data_mat[12] = ((real_t)(elMat_3_0));
            _data_mat[13] = ((real_t)(elMat_3_1));
            _data_mat[14] = ((real_t)(elMat_3_2));
            _data_mat[15] = ((real_t)(elMat_3_3));

            mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
          }
    }
  }
}
void P1ElementwiseDivKGrad_IcosahedralShellMap_float64::
    computeInverseDiagonalOperatorValues_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(
        walberla::float64 *RESTRICT _data_invDiag_,
        walberla::float64 *RESTRICT _data_k,
        walberla::float64 macro_vertex_coord_id_0comp0,
        walberla::float64 macro_vertex_coord_id_0comp1,
        walberla::float64 macro_vertex_coord_id_0comp2,
        walberla::float64 macro_vertex_coord_id_1comp0,
        walberla::float64 macro_vertex_coord_id_1comp1,
        walberla::float64 macro_vertex_coord_id_1comp2,
        walberla::float64 macro_vertex_coord_id_2comp0,
        walberla::float64 macro_vertex_coord_id_2comp1,
        walberla::float64 macro_vertex_coord_id_2comp2,
        walberla::float64 macro_vertex_coord_id_3comp0,
        walberla::float64 macro_vertex_coord_id_3comp1,
        walberla::float64 macro_vertex_coord_id_3comp2,
        int64_t micro_edges_per_macro_edge,
        walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    {
      /* CellType.WHITE_UP */
      const walberla::float64 tmp_coords_jac_0_WHITE_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 p_affine_const_0_0_WHITE_UP =
          macro_vertex_coord_id_0comp0;
      const walberla::float64 p_affine_const_0_1_WHITE_UP =
          macro_vertex_coord_id_0comp1;
      const walberla::float64 p_affine_const_0_2_WHITE_UP =
          macro_vertex_coord_id_0comp2;
      const walberla::float64 p_affine_const_1_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 p_affine_const_1_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 p_affine_const_1_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 p_affine_const_2_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 p_affine_const_2_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 p_affine_const_2_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_3_0_WHITE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 p_affine_const_3_1_WHITE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 p_affine_const_3_2_WHITE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 +
                        ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_3 = tmp_kernel_op_0;
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_0 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_1 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_2 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_3 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6))];
          }
    }
    {
      /* CellType.WHITE_DOWN */
      const walberla::float64 tmp_coords_jac_0_WHITE_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_WHITE_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_3_WHITE_DOWN =
          tmp_coords_jac_1_WHITE_DOWN + tmp_coords_jac_2_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_4_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_5_WHITE_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_6_WHITE_DOWN =
          tmp_coords_jac_4_WHITE_DOWN + tmp_coords_jac_5_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_7_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_8_WHITE_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_WHITE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_9_WHITE_DOWN =
          tmp_coords_jac_7_WHITE_DOWN + tmp_coords_jac_8_WHITE_DOWN;
      const walberla::float64 tmp_coords_jac_10_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_11_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_12_WHITE_DOWN =
          tmp_coords_jac_0_WHITE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_WHITE_DOWN =
          tmp_coords_jac_3_WHITE_DOWN;
      const walberla::float64 p_affine_const_0_1_WHITE_DOWN =
          tmp_coords_jac_6_WHITE_DOWN;
      const walberla::float64 p_affine_const_0_2_WHITE_DOWN =
          tmp_coords_jac_9_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_0_WHITE_DOWN =
          tmp_coords_jac_10_WHITE_DOWN + tmp_coords_jac_2_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_1_WHITE_DOWN =
          tmp_coords_jac_11_WHITE_DOWN + tmp_coords_jac_5_WHITE_DOWN;
      const walberla::float64 p_affine_const_1_2_WHITE_DOWN =
          tmp_coords_jac_12_WHITE_DOWN + tmp_coords_jac_8_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_0_WHITE_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_10_WHITE_DOWN +
          tmp_coords_jac_1_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_1_WHITE_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_11_WHITE_DOWN +
          tmp_coords_jac_4_WHITE_DOWN;
      const walberla::float64 p_affine_const_2_2_WHITE_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_12_WHITE_DOWN +
          tmp_coords_jac_7_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_0_WHITE_DOWN =
          tmp_coords_jac_10_WHITE_DOWN + tmp_coords_jac_3_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_1_WHITE_DOWN =
          tmp_coords_jac_11_WHITE_DOWN + tmp_coords_jac_6_WHITE_DOWN;
      const walberla::float64 p_affine_const_3_2_WHITE_DOWN =
          tmp_coords_jac_12_WHITE_DOWN + tmp_coords_jac_9_WHITE_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 +
                        ctr_1 * (-ctr_2 + micro_edges_per_macro_edge - 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 2) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_3 = tmp_kernel_op_0;
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_0 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_1 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_2 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_3 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
          }
    }
    {
      /* CellType.BLUE_UP */
      const walberla::float64 tmp_coords_jac_0_BLUE_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_BLUE_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_2_BLUE_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_3_BLUE_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_4_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_5_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_6_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_0_0_BLUE_UP =
          tmp_coords_jac_1_BLUE_UP;
      const walberla::float64 p_affine_const_0_1_BLUE_UP =
          tmp_coords_jac_2_BLUE_UP;
      const walberla::float64 p_affine_const_0_2_BLUE_UP =
          tmp_coords_jac_3_BLUE_UP;
      const walberla::float64 p_affine_const_1_0_BLUE_UP =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_BLUE_UP;
      const walberla::float64 p_affine_const_1_1_BLUE_UP =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_BLUE_UP;
      const walberla::float64 p_affine_const_1_2_BLUE_UP =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_BLUE_UP;
      const walberla::float64 p_affine_const_2_0_BLUE_UP =
          tmp_coords_jac_1_BLUE_UP + tmp_coords_jac_4_BLUE_UP;
      const walberla::float64 p_affine_const_2_1_BLUE_UP =
          tmp_coords_jac_2_BLUE_UP + tmp_coords_jac_5_BLUE_UP;
      const walberla::float64 p_affine_const_2_2_BLUE_UP =
          tmp_coords_jac_3_BLUE_UP + tmp_coords_jac_6_BLUE_UP;
      const walberla::float64 p_affine_const_3_0_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0) +
          tmp_coords_jac_1_BLUE_UP;
      const walberla::float64 p_affine_const_3_1_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1) +
          tmp_coords_jac_2_BLUE_UP;
      const walberla::float64 p_affine_const_3_2_BLUE_UP =
          tmp_coords_jac_0_BLUE_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2) +
          tmp_coords_jac_3_BLUE_UP;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 1) *
                          (micro_edges_per_macro_edge + 1)) /
                         (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_3 = tmp_kernel_op_0;
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_0 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_1 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_2 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_3 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
          }
    }
    {
      /* CellType.BLUE_DOWN */
      const walberla::float64 tmp_coords_jac_0_BLUE_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_BLUE_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_BLUE_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_3_BLUE_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_4_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_5_BLUE_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_BLUE_DOWN;
      const walberla::float64 tmp_coords_jac_6_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_7_BLUE_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_6_BLUE_DOWN;
      const walberla::float64 tmp_coords_jac_8_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 tmp_coords_jac_9_BLUE_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_8_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_0_BLUE_DOWN =
          tmp_coords_jac_1_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_1_BLUE_DOWN =
          tmp_coords_jac_2_BLUE_DOWN;
      const walberla::float64 p_affine_const_0_2_BLUE_DOWN =
          tmp_coords_jac_3_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_0_BLUE_DOWN =
          tmp_coords_jac_5_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_1_BLUE_DOWN =
          tmp_coords_jac_7_BLUE_DOWN;
      const walberla::float64 p_affine_const_1_2_BLUE_DOWN =
          tmp_coords_jac_9_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_0_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0) +
          tmp_coords_jac_5_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_1_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1) +
          tmp_coords_jac_7_BLUE_DOWN;
      const walberla::float64 p_affine_const_2_2_BLUE_DOWN =
          tmp_coords_jac_0_BLUE_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2) +
          tmp_coords_jac_9_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_0_BLUE_DOWN =
          tmp_coords_jac_1_BLUE_DOWN + tmp_coords_jac_4_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_1_BLUE_DOWN =
          tmp_coords_jac_2_BLUE_DOWN + tmp_coords_jac_6_BLUE_DOWN;
      const walberla::float64 p_affine_const_3_2_BLUE_DOWN =
          tmp_coords_jac_3_BLUE_DOWN + tmp_coords_jac_8_BLUE_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        3 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_3 = tmp_kernel_op_0;
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_0 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_1 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_2 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_3 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6))];
          }
    }
    {
      /* CellType.GREEN_UP */
      const walberla::float64 tmp_coords_jac_0_GREEN_UP =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_GREEN_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_2_GREEN_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_3_GREEN_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_4_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_5_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_6_GREEN_UP =
          tmp_coords_jac_0_GREEN_UP *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_GREEN_UP =
          tmp_coords_jac_1_GREEN_UP;
      const walberla::float64 p_affine_const_0_1_GREEN_UP =
          tmp_coords_jac_2_GREEN_UP;
      const walberla::float64 p_affine_const_0_2_GREEN_UP =
          tmp_coords_jac_3_GREEN_UP;
      const walberla::float64 p_affine_const_1_0_GREEN_UP =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 p_affine_const_1_1_GREEN_UP =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 p_affine_const_1_2_GREEN_UP =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_UP *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 p_affine_const_2_0_GREEN_UP =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_GREEN_UP;
      const walberla::float64 p_affine_const_2_1_GREEN_UP =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_GREEN_UP;
      const walberla::float64 p_affine_const_2_2_GREEN_UP =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_GREEN_UP;
      const walberla::float64 p_affine_const_3_0_GREEN_UP =
          tmp_coords_jac_1_GREEN_UP + tmp_coords_jac_4_GREEN_UP;
      const walberla::float64 p_affine_const_3_1_GREEN_UP =
          tmp_coords_jac_2_GREEN_UP + tmp_coords_jac_5_GREEN_UP;
      const walberla::float64 p_affine_const_3_2_GREEN_UP =
          tmp_coords_jac_3_GREEN_UP + tmp_coords_jac_6_GREEN_UP;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_3 = tmp_kernel_op_0;
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_0 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_1 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_2 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_3 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
          }
    }
    {
      /* CellType.GREEN_DOWN */
      const walberla::float64 tmp_coords_jac_0_GREEN_DOWN =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
      const walberla::float64 tmp_coords_jac_1_GREEN_DOWN =
          macro_vertex_coord_id_0comp0 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
      const walberla::float64 tmp_coords_jac_2_GREEN_DOWN =
          macro_vertex_coord_id_0comp1 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
      const walberla::float64 tmp_coords_jac_3_GREEN_DOWN =
          macro_vertex_coord_id_0comp2 +
          tmp_coords_jac_0_GREEN_DOWN *
              (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_2comp2);
      const walberla::float64 tmp_coords_jac_4_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
      const walberla::float64 tmp_coords_jac_5_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
      const walberla::float64 tmp_coords_jac_6_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_1comp2);
      const walberla::float64 tmp_coords_jac_7_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_3comp0);
      const walberla::float64 tmp_coords_jac_8_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_3comp1);
      const walberla::float64 tmp_coords_jac_9_GREEN_DOWN =
          tmp_coords_jac_0_GREEN_DOWN *
          (-macro_vertex_coord_id_0comp2 + macro_vertex_coord_id_3comp2);
      const walberla::float64 p_affine_const_0_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN;
      const walberla::float64 p_affine_const_0_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN;
      const walberla::float64 p_affine_const_0_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN + tmp_coords_jac_4_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN + tmp_coords_jac_5_GREEN_DOWN;
      const walberla::float64 p_affine_const_1_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN + tmp_coords_jac_6_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_0_GREEN_DOWN =
          macro_vertex_coord_id_0comp0 + tmp_coords_jac_4_GREEN_DOWN +
          tmp_coords_jac_7_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_1_GREEN_DOWN =
          macro_vertex_coord_id_0comp1 + tmp_coords_jac_5_GREEN_DOWN +
          tmp_coords_jac_8_GREEN_DOWN;
      const walberla::float64 p_affine_const_2_2_GREEN_DOWN =
          macro_vertex_coord_id_0comp2 + tmp_coords_jac_6_GREEN_DOWN +
          tmp_coords_jac_9_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_0_GREEN_DOWN =
          tmp_coords_jac_1_GREEN_DOWN + tmp_coords_jac_7_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_1_GREEN_DOWN =
          tmp_coords_jac_2_GREEN_DOWN + tmp_coords_jac_8_GREEN_DOWN;
      const walberla::float64 p_affine_const_3_2_GREEN_DOWN =
          tmp_coords_jac_3_GREEN_DOWN + tmp_coords_jac_9_GREEN_DOWN;
      for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
        for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
             ctr_1 += 1)
          for (int64_t ctr_0 = 0;
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 1;
               ctr_0 += 1) {
            const walberla::float64 p_affine_0_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_0_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_1_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2));
            const walberla::float64 p_affine_2_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_2_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_0 =
                macro_vertex_coord_id_0comp0 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_1comp0) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_2comp0) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp0 +
                     macro_vertex_coord_id_3comp0) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_1 =
                macro_vertex_coord_id_0comp1 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_1comp1) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_2comp1) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp1 +
                     macro_vertex_coord_id_3comp1) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 p_affine_3_2 =
                macro_vertex_coord_id_0comp2 +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_1comp2) *
                    1.0 * ((walberla::float64)(ctr_0)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_2comp2) *
                    1.0 * ((walberla::float64)(ctr_1 + 1)) +
                1.0 / (micro_edges_per_macro_edge_float) *
                    (-macro_vertex_coord_id_0comp2 +
                     macro_vertex_coord_id_3comp2) *
                    1.0 * ((walberla::float64)(ctr_2 + 1));
            const walberla::float64 k_dof_0 =
                _data_k[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge - 2) *
                          (micro_edges_per_macro_edge - 1)) /
                         (6)) +
                        4 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge - 1) *
                              (micro_edges_per_macro_edge + 1)) /
                             (6)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge - 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1)) /
                         (6))];
            const walberla::float64 tmp_kernel_op_0 =
                k_dof_0 * 0.16666666666666663;
            const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
            const walberla::float64 elMatDiag_3 = tmp_kernel_op_0;
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_0 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6))];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_1 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) -
                               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) +
                           1] =
                elMatDiag_2 +
                _data_invDiag_[ctr_0 +
                               ctr_1 *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6)) +
                               1];
            _data_invDiag_[ctr_0 +
                           (ctr_1 + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))] =
                elMatDiag_3 +
                _data_invDiag_[ctr_0 +
                               (ctr_1 + 1) *
                                   (-ctr_2 + micro_edges_per_macro_edge + 1) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                               (((-ctr_2 + micro_edges_per_macro_edge) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                                (6)) +
                               (((micro_edges_per_macro_edge + 1) *
                                 (micro_edges_per_macro_edge + 2) *
                                 (micro_edges_per_macro_edge + 3)) /
                                (6))];
          }
    }
  }
}

} // namespace operatorgeneration

} // namespace hyteg
