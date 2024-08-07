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
 * The entire file was generated with the HyTeG form generator.
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

#include "P1ElementwiseDiffusion_cubes_const_vect_float64.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P1ElementwiseDiffusion_cubes_const_vect_float64::
    P1ElementwiseDiffusion_cubes_const_vect_float64(
        const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
        size_t maxLevel)
    : Operator(storage, minLevel, maxLevel) {}

void P1ElementwiseDiffusion_cubes_const_vect_float64::apply(
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
    communication::syncFunctionBetweenPrimitives(
        src, level, communication::syncDirection_t::LOW2HIGH);
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

      this->timingTree_->start("kernel" + std::to_string(level));
      LIKWID_MARKER_START( "kernel"  + std::to_string(level) );
      apply_macro_3D(

          _data_dst, _data_src, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_0comp2,
          macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
          macro_vertex_coord_id_1comp2, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, macro_vertex_coord_id_2comp2,
          macro_vertex_coord_id_3comp0, macro_vertex_coord_id_3comp1,
          macro_vertex_coord_id_3comp2, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float);
      LIKWID_MARKER_STOP( "kernel" + std::to_string(level) );
      this->timingTree_->stop("kernel" + std::to_string(level));
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
    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      // get hold of the actual numerical data in the functions
      walberla::float64 *_data_src =
          face.getData(src.getFaceDataID())->getPointer(level);
      walberla::float64 *_data_dst =
          face.getData(dst.getFaceDataID())->getPointer(level);

      // Zero out dst halos only
      //
      // This is also necessary when using update type == Add.
      // During additive comm we then skip zeroing the data on the lower-dim
      // primitives.
      for (const auto &idx : vertexdof::macroface::Iterator(level)) {
        if (vertexdof::macroface::isVertexOnBoundary(level, idx)) {
          auto arrayIdx = vertexdof::macroface::index(level, idx.x(), idx.y());
          _data_dst[arrayIdx] = walberla::float64(0);
        }
      }

      const auto micro_edges_per_macro_edge =
          (int64_t)levelinfo::num_microedges_per_edge(level);
      const auto micro_edges_per_macro_edge_float =
          (walberla::float64)levelinfo::num_microedges_per_edge(level);
      const walberla::float64 macro_vertex_coord_id_0comp0 =
          (walberla::float64)face.getCoordinates()[0][0];
      const walberla::float64 macro_vertex_coord_id_0comp1 =
          (walberla::float64)face.getCoordinates()[0][1];
      const walberla::float64 macro_vertex_coord_id_1comp0 =
          (walberla::float64)face.getCoordinates()[1][0];
      const walberla::float64 macro_vertex_coord_id_1comp1 =
          (walberla::float64)face.getCoordinates()[1][1];
      const walberla::float64 macro_vertex_coord_id_2comp0 =
          (walberla::float64)face.getCoordinates()[2][0];
      const walberla::float64 macro_vertex_coord_id_2comp1 =
          (walberla::float64)face.getCoordinates()[2][1];

      this->timingTree_->start("kernel");

      apply_macro_2D(

          _data_dst, _data_src, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_1comp0,
          macro_vertex_coord_id_1comp1, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float);
      this->timingTree_->stop("kernel");
    }

    // Push result to lower-dimensional primitives
    //
    this->timingTree_->start("post-communication");
    // Note: We could avoid communication here by implementing the apply() also
    // for the respective
    //       lower dimensional primitives!
    dst.communicateAdditively<Face, Edge>(level, DoFType::All ^ flag, *storage_,
                                          updateType == Replace);
    dst.communicateAdditively<Face, Vertex>(level, DoFType::All ^ flag,
                                            *storage_, updateType == Replace);
    this->timingTree_->stop("post-communication");
  }

  this->stopTiming("apply");
}
void P1ElementwiseDiffusion_cubes_const_vect_float64::
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

        const auto micro_edges_per_macro_edge =
            (int64_t)levelinfo::num_microedges_per_edge(level);
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

        this->timingTree_->start("kernel");

        computeInverseDiagonalOperatorValues_macro_3D(

            _data_invDiag_, macro_vertex_coord_id_0comp0,
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
    } else {
      this->timingTree_->start("pre-communication");

      this->timingTree_->stop("pre-communication");

      for (auto &it : storage_->getFaces()) {
        Face &face = *it.second;

        // get hold of the actual numerical data
        walberla::float64 *_data_invDiag_ =
            face.getData((*invDiag_).getFaceDataID())->getPointer(level);

        const auto micro_edges_per_macro_edge =
            (int64_t)levelinfo::num_microedges_per_edge(level);
        const auto micro_edges_per_macro_edge_float =
            (walberla::float64)levelinfo::num_microedges_per_edge(level);
        const walberla::float64 macro_vertex_coord_id_0comp0 =
            (walberla::float64)face.getCoordinates()[0][0];
        const walberla::float64 macro_vertex_coord_id_0comp1 =
            (walberla::float64)face.getCoordinates()[0][1];
        const walberla::float64 macro_vertex_coord_id_1comp0 =
            (walberla::float64)face.getCoordinates()[1][0];
        const walberla::float64 macro_vertex_coord_id_1comp1 =
            (walberla::float64)face.getCoordinates()[1][1];
        const walberla::float64 macro_vertex_coord_id_2comp0 =
            (walberla::float64)face.getCoordinates()[2][0];
        const walberla::float64 macro_vertex_coord_id_2comp1 =
            (walberla::float64)face.getCoordinates()[2][1];

        this->timingTree_->start("kernel");

        computeInverseDiagonalOperatorValues_macro_2D(

            _data_invDiag_, macro_vertex_coord_id_0comp0,
            macro_vertex_coord_id_0comp1, macro_vertex_coord_id_1comp0,
            macro_vertex_coord_id_1comp1, macro_vertex_coord_id_2comp0,
            macro_vertex_coord_id_2comp1, micro_edges_per_macro_edge,
            micro_edges_per_macro_edge_float);
        this->timingTree_->stop("kernel");
      }

      // Push result to lower-dimensional primitives
      //
      this->timingTree_->start("post-communication");
      // Note: We could avoid communication here by implementing the apply()
      // also for the respective
      //       lower dimensional primitives!
      (*invDiag_).communicateAdditively<Face, Edge>(level);
      (*invDiag_).communicateAdditively<Face, Vertex>(level);
      this->timingTree_->stop("post-communication");
    }

    (*invDiag_).invertElementwise(level);
  }

  this->stopTiming("computeInverseDiagonalOperatorValues");
}
std::shared_ptr<P1Function<walberla::float64>>
P1ElementwiseDiffusion_cubes_const_vect_float64::getInverseDiagonalValues()
    const {
  return invDiag_;
}
void P1ElementwiseDiffusion_cubes_const_vect_float64::apply_macro_2D(
    walberla::float64 *RESTRICT _data_dst,
    walberla::float64 *RESTRICT _data_src,
    walberla::float64 macro_vertex_coord_id_0comp0,
    walberla::float64 macro_vertex_coord_id_0comp1,
    walberla::float64 macro_vertex_coord_id_1comp0,
    walberla::float64 macro_vertex_coord_id_1comp1,
    walberla::float64 macro_vertex_coord_id_2comp0,
    walberla::float64 macro_vertex_coord_id_2comp1,
    int64_t micro_edges_per_macro_edge,
    walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    const walberla::float64 tmp_coords_jac_0_BLUE =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
    const walberla::float64 tmp_coords_jac_1_BLUE =
        macro_vertex_coord_id_0comp0 +
        tmp_coords_jac_0_BLUE *
            (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
    const walberla::float64 tmp_coords_jac_2_BLUE =
        macro_vertex_coord_id_0comp1 +
        tmp_coords_jac_0_BLUE *
            (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
    const walberla::float64 tmp_coords_jac_3_BLUE =
        tmp_coords_jac_0_BLUE *
        (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
    const walberla::float64 tmp_coords_jac_4_BLUE =
        tmp_coords_jac_0_BLUE *
        (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
    const walberla::float64 p_affine_const_0_0_BLUE = tmp_coords_jac_1_BLUE;
    const walberla::float64 p_affine_const_0_1_BLUE = tmp_coords_jac_2_BLUE;
    const walberla::float64 p_affine_const_1_0_BLUE =
        macro_vertex_coord_id_0comp0 + tmp_coords_jac_3_BLUE;
    const walberla::float64 p_affine_const_1_1_BLUE =
        macro_vertex_coord_id_0comp1 + tmp_coords_jac_4_BLUE;
    const walberla::float64 p_affine_const_2_0_BLUE =
        tmp_coords_jac_1_BLUE + tmp_coords_jac_3_BLUE;
    const walberla::float64 p_affine_const_2_1_BLUE =
        tmp_coords_jac_2_BLUE + tmp_coords_jac_4_BLUE;
    const walberla::float64 jac_affine_0_0_BLUE =
        -p_affine_const_0_0_BLUE + p_affine_const_1_0_BLUE;
    const walberla::float64 jac_affine_0_1_BLUE =
        -p_affine_const_0_0_BLUE + p_affine_const_2_0_BLUE;
    const walberla::float64 jac_affine_1_0_BLUE =
        -p_affine_const_0_1_BLUE + p_affine_const_1_1_BLUE;
    const walberla::float64 jac_affine_1_1_BLUE =
        -p_affine_const_0_1_BLUE + p_affine_const_2_1_BLUE;
    const walberla::float64 tmp_coords_jac_5_BLUE =
        jac_affine_0_0_BLUE * jac_affine_1_1_BLUE -
        jac_affine_0_1_BLUE * jac_affine_1_0_BLUE;
    const walberla::float64 tmp_coords_jac_6_BLUE =
        1.0 / (tmp_coords_jac_5_BLUE);
    const walberla::float64 jac_affine_inv_0_0_BLUE =
        jac_affine_1_1_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 jac_affine_inv_0_1_BLUE =
        -jac_affine_0_1_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 jac_affine_inv_1_0_BLUE =
        -jac_affine_1_0_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 jac_affine_inv_1_1_BLUE =
        jac_affine_0_0_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 abs_det_jac_affine_BLUE =
        abs(tmp_coords_jac_5_BLUE);
    const walberla::float64 tmp_coords_jac_0_GRAY =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
    const walberla::float64 p_affine_const_0_0_GRAY =
        macro_vertex_coord_id_0comp0;
    const walberla::float64 p_affine_const_0_1_GRAY =
        macro_vertex_coord_id_0comp1;
    const walberla::float64 p_affine_const_1_0_GRAY =
        macro_vertex_coord_id_0comp0 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
    const walberla::float64 p_affine_const_1_1_GRAY =
        macro_vertex_coord_id_0comp1 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
    const walberla::float64 p_affine_const_2_0_GRAY =
        macro_vertex_coord_id_0comp0 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
    const walberla::float64 p_affine_const_2_1_GRAY =
        macro_vertex_coord_id_0comp1 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
    const walberla::float64 jac_affine_0_0_GRAY =
        -p_affine_const_0_0_GRAY + p_affine_const_1_0_GRAY;
    const walberla::float64 jac_affine_0_1_GRAY =
        -p_affine_const_0_0_GRAY + p_affine_const_2_0_GRAY;
    const walberla::float64 jac_affine_1_0_GRAY =
        -p_affine_const_0_1_GRAY + p_affine_const_1_1_GRAY;
    const walberla::float64 jac_affine_1_1_GRAY =
        -p_affine_const_0_1_GRAY + p_affine_const_2_1_GRAY;
    const walberla::float64 tmp_coords_jac_1_GRAY =
        jac_affine_0_0_GRAY * jac_affine_1_1_GRAY -
        jac_affine_0_1_GRAY * jac_affine_1_0_GRAY;
    const walberla::float64 tmp_coords_jac_2_GRAY =
        1.0 / (tmp_coords_jac_1_GRAY);
    const walberla::float64 jac_affine_inv_0_0_GRAY =
        jac_affine_1_1_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 jac_affine_inv_0_1_GRAY =
        -jac_affine_0_1_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 jac_affine_inv_1_0_GRAY =
        -jac_affine_1_0_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 jac_affine_inv_1_1_GRAY =
        jac_affine_0_0_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 abs_det_jac_affine_GRAY =
        abs(tmp_coords_jac_1_GRAY);
    const walberla::float64 tmp_kernel_op_0 =
        -jac_affine_inv_0_0_GRAY - jac_affine_inv_1_0_GRAY;
    const walberla::float64 tmp_kernel_op_1 =
        -jac_affine_inv_0_1_GRAY - jac_affine_inv_1_1_GRAY;
    const walberla::float64 tmp_kernel_op_2 = abs_det_jac_affine_GRAY * 0.5;
    const walberla::float64 tmp_kernel_op_4 =
        jac_affine_inv_0_0_GRAY * tmp_kernel_op_0 +
        jac_affine_inv_0_1_GRAY * tmp_kernel_op_1;
    const walberla::float64 tmp_kernel_op_6 =
        jac_affine_inv_1_0_GRAY * tmp_kernel_op_0 +
        jac_affine_inv_1_1_GRAY * tmp_kernel_op_1;
    const walberla::float64 tmp_kernel_op_8 =
        jac_affine_inv_0_0_GRAY * jac_affine_inv_1_0_GRAY +
        jac_affine_inv_0_1_GRAY * jac_affine_inv_1_1_GRAY;
    const walberla::float64 tmp_moved_constant_3 =
        -jac_affine_inv_0_0_BLUE - jac_affine_inv_1_0_BLUE;
    const walberla::float64 tmp_moved_constant_4 =
        -jac_affine_inv_0_1_BLUE - jac_affine_inv_1_1_BLUE;
    const walberla::float64 tmp_moved_constant_5 =
        abs_det_jac_affine_BLUE * 0.5;
    const walberla::float64 tmp_moved_constant_7 =
        jac_affine_inv_0_0_BLUE * tmp_moved_constant_3 +
        jac_affine_inv_0_1_BLUE * tmp_moved_constant_4;
    const walberla::float64 tmp_moved_constant_9 =
        jac_affine_inv_1_0_BLUE * tmp_moved_constant_3 +
        jac_affine_inv_1_1_BLUE * tmp_moved_constant_4;
    const walberla::float64 tmp_moved_constant_11 =
        jac_affine_inv_0_0_BLUE * jac_affine_inv_1_0_BLUE +
        jac_affine_inv_0_1_BLUE * jac_affine_inv_1_1_BLUE;
    for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1) {
      {
        for (int64_t ctr_0 = 0;
             ctr_0 <
             (int64_t)((-ctr_1 + micro_edges_per_macro_edge - 1) / (4)) * (4);
             ctr_0 += 4) {
          const __m256d src_dof_0 = _mm256_loadu_pd(
              &_data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2))]);
          const __m256d src_dof_1 = _mm256_loadu_pd(
              &_data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]);
          const __m256d src_dof_2 = _mm256_loadu_pd(
              &_data_src[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]);
          const __m256d tmp_kernel_op_3 = _mm256_mul_pd(
              src_dof_0, _mm256_set_pd(tmp_kernel_op_2, tmp_kernel_op_2,
                                       tmp_kernel_op_2, tmp_kernel_op_2));
          const __m256d tmp_kernel_op_5 = _mm256_mul_pd(
              src_dof_1, _mm256_set_pd(tmp_kernel_op_2, tmp_kernel_op_2,
                                       tmp_kernel_op_2, tmp_kernel_op_2));
          const __m256d tmp_kernel_op_7 = _mm256_mul_pd(
              src_dof_2, _mm256_set_pd(tmp_kernel_op_2, tmp_kernel_op_2,
                                       tmp_kernel_op_2, tmp_kernel_op_2));
          const __m256d elMatVec_0 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(
                      tmp_kernel_op_3,
                      _mm256_add_pd(
                          _mm256_mul_pd(
                              _mm256_set_pd(tmp_kernel_op_0, tmp_kernel_op_0,
                                            tmp_kernel_op_0, tmp_kernel_op_0),
                              _mm256_set_pd(tmp_kernel_op_0, tmp_kernel_op_0,
                                            tmp_kernel_op_0, tmp_kernel_op_0)),
                          _mm256_mul_pd(
                              _mm256_set_pd(tmp_kernel_op_1, tmp_kernel_op_1,
                                            tmp_kernel_op_1, tmp_kernel_op_1),
                              _mm256_set_pd(tmp_kernel_op_1, tmp_kernel_op_1,
                                            tmp_kernel_op_1,
                                            tmp_kernel_op_1)))),
                  _mm256_mul_pd(tmp_kernel_op_5,
                                _mm256_set_pd(tmp_kernel_op_4, tmp_kernel_op_4,
                                              tmp_kernel_op_4,
                                              tmp_kernel_op_4))),
              _mm256_mul_pd(tmp_kernel_op_7,
                            _mm256_set_pd(tmp_kernel_op_6, tmp_kernel_op_6,
                                          tmp_kernel_op_6, tmp_kernel_op_6)));
          const __m256d elMatVec_1 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(
                      tmp_kernel_op_5,
                      _mm256_add_pd(
                          _mm256_mul_pd(_mm256_set_pd(jac_affine_inv_0_0_GRAY,
                                                      jac_affine_inv_0_0_GRAY,
                                                      jac_affine_inv_0_0_GRAY,
                                                      jac_affine_inv_0_0_GRAY),
                                        _mm256_set_pd(jac_affine_inv_0_0_GRAY,
                                                      jac_affine_inv_0_0_GRAY,
                                                      jac_affine_inv_0_0_GRAY,
                                                      jac_affine_inv_0_0_GRAY)),
                          _mm256_mul_pd(
                              _mm256_set_pd(jac_affine_inv_0_1_GRAY,
                                            jac_affine_inv_0_1_GRAY,
                                            jac_affine_inv_0_1_GRAY,
                                            jac_affine_inv_0_1_GRAY),
                              _mm256_set_pd(jac_affine_inv_0_1_GRAY,
                                            jac_affine_inv_0_1_GRAY,
                                            jac_affine_inv_0_1_GRAY,
                                            jac_affine_inv_0_1_GRAY)))),
                  _mm256_mul_pd(tmp_kernel_op_3,
                                _mm256_set_pd(tmp_kernel_op_4, tmp_kernel_op_4,
                                              tmp_kernel_op_4,
                                              tmp_kernel_op_4))),
              _mm256_mul_pd(tmp_kernel_op_7,
                            _mm256_set_pd(tmp_kernel_op_8, tmp_kernel_op_8,
                                          tmp_kernel_op_8, tmp_kernel_op_8)));
          const __m256d elMatVec_2 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(
                      tmp_kernel_op_7,
                      _mm256_add_pd(
                          _mm256_mul_pd(_mm256_set_pd(jac_affine_inv_1_0_GRAY,
                                                      jac_affine_inv_1_0_GRAY,
                                                      jac_affine_inv_1_0_GRAY,
                                                      jac_affine_inv_1_0_GRAY),
                                        _mm256_set_pd(jac_affine_inv_1_0_GRAY,
                                                      jac_affine_inv_1_0_GRAY,
                                                      jac_affine_inv_1_0_GRAY,
                                                      jac_affine_inv_1_0_GRAY)),
                          _mm256_mul_pd(
                              _mm256_set_pd(jac_affine_inv_1_1_GRAY,
                                            jac_affine_inv_1_1_GRAY,
                                            jac_affine_inv_1_1_GRAY,
                                            jac_affine_inv_1_1_GRAY),
                              _mm256_set_pd(jac_affine_inv_1_1_GRAY,
                                            jac_affine_inv_1_1_GRAY,
                                            jac_affine_inv_1_1_GRAY,
                                            jac_affine_inv_1_1_GRAY)))),
                  _mm256_mul_pd(tmp_kernel_op_3,
                                _mm256_set_pd(tmp_kernel_op_6, tmp_kernel_op_6,
                                              tmp_kernel_op_6,
                                              tmp_kernel_op_6))),
              _mm256_mul_pd(tmp_kernel_op_5,
                            _mm256_set_pd(tmp_kernel_op_8, tmp_kernel_op_8,
                                          tmp_kernel_op_8, tmp_kernel_op_8)));
          {
            {
              _mm256_storeu_pd(
                  &_data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                             ((ctr_1 * (ctr_1 + 1)) / (2))],
                  _mm256_add_pd(
                      elMatVec_0,
                      _mm256_loadu_pd(
                          &_data_dst[ctr_0 +
                                     ctr_1 * (micro_edges_per_macro_edge + 2) -
                                     ((ctr_1 * (ctr_1 + 1)) / (2))])));
              _mm256_storeu_pd(
                  &_data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) + 1],
                  _mm256_add_pd(
                      elMatVec_1,
                      _mm256_loadu_pd(
                          &_data_dst[ctr_0 +
                                     ctr_1 * (micro_edges_per_macro_edge + 2) -
                                     ((ctr_1 * (ctr_1 + 1)) / (2)) + 1])));
              _mm256_storeu_pd(
                  &_data_dst[ctr_0 +
                             (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                             (((ctr_1 + 1) * (ctr_1 + 2)) / (2))],
                  _mm256_add_pd(
                      elMatVec_2,
                      _mm256_loadu_pd(
                          &_data_dst[ctr_0 +
                                     (ctr_1 + 1) *
                                         (micro_edges_per_macro_edge + 2) -
                                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2))])));
            }
          }
          const __m256d tmp_moved_constant_0 = _mm256_loadu_pd(
              &_data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]);
          const __m256d tmp_moved_constant_1 = _mm256_loadu_pd(
              &_data_src[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]);
          const __m256d tmp_moved_constant_2 = _mm256_loadu_pd(
              &_data_src[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1]);
          const __m256d tmp_moved_constant_6 = _mm256_mul_pd(
              tmp_moved_constant_0,
              _mm256_set_pd(tmp_moved_constant_5, tmp_moved_constant_5,
                            tmp_moved_constant_5, tmp_moved_constant_5));
          const __m256d tmp_moved_constant_8 = _mm256_mul_pd(
              tmp_moved_constant_1,
              _mm256_set_pd(tmp_moved_constant_5, tmp_moved_constant_5,
                            tmp_moved_constant_5, tmp_moved_constant_5));
          const __m256d tmp_moved_constant_10 = _mm256_mul_pd(
              tmp_moved_constant_2,
              _mm256_set_pd(tmp_moved_constant_5, tmp_moved_constant_5,
                            tmp_moved_constant_5, tmp_moved_constant_5));
          const __m256d tmp_moved_constant_12 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(
                      tmp_moved_constant_6,
                      _mm256_add_pd(
                          _mm256_mul_pd(
                              _mm256_set_pd(
                                  tmp_moved_constant_3, tmp_moved_constant_3,
                                  tmp_moved_constant_3, tmp_moved_constant_3),
                              _mm256_set_pd(
                                  tmp_moved_constant_3, tmp_moved_constant_3,
                                  tmp_moved_constant_3, tmp_moved_constant_3)),
                          _mm256_mul_pd(_mm256_set_pd(tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4),
                                        _mm256_set_pd(tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4)))),
                  _mm256_mul_pd(tmp_moved_constant_8,
                                _mm256_set_pd(tmp_moved_constant_7,
                                              tmp_moved_constant_7,
                                              tmp_moved_constant_7,
                                              tmp_moved_constant_7))),
              _mm256_mul_pd(
                  tmp_moved_constant_10,
                  _mm256_set_pd(tmp_moved_constant_9, tmp_moved_constant_9,
                                tmp_moved_constant_9, tmp_moved_constant_9)));
          const __m256d tmp_moved_constant_13 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(
                      tmp_moved_constant_8,
                      _mm256_add_pd(
                          _mm256_mul_pd(_mm256_set_pd(jac_affine_inv_0_0_BLUE,
                                                      jac_affine_inv_0_0_BLUE,
                                                      jac_affine_inv_0_0_BLUE,
                                                      jac_affine_inv_0_0_BLUE),
                                        _mm256_set_pd(jac_affine_inv_0_0_BLUE,
                                                      jac_affine_inv_0_0_BLUE,
                                                      jac_affine_inv_0_0_BLUE,
                                                      jac_affine_inv_0_0_BLUE)),
                          _mm256_mul_pd(
                              _mm256_set_pd(jac_affine_inv_0_1_BLUE,
                                            jac_affine_inv_0_1_BLUE,
                                            jac_affine_inv_0_1_BLUE,
                                            jac_affine_inv_0_1_BLUE),
                              _mm256_set_pd(jac_affine_inv_0_1_BLUE,
                                            jac_affine_inv_0_1_BLUE,
                                            jac_affine_inv_0_1_BLUE,
                                            jac_affine_inv_0_1_BLUE)))),
                  _mm256_mul_pd(tmp_moved_constant_10,
                                _mm256_set_pd(tmp_moved_constant_11,
                                              tmp_moved_constant_11,
                                              tmp_moved_constant_11,
                                              tmp_moved_constant_11))),
              _mm256_mul_pd(
                  tmp_moved_constant_6,
                  _mm256_set_pd(tmp_moved_constant_7, tmp_moved_constant_7,
                                tmp_moved_constant_7, tmp_moved_constant_7)));
          const __m256d tmp_moved_constant_14 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(
                      tmp_moved_constant_10,
                      _mm256_add_pd(
                          _mm256_mul_pd(_mm256_set_pd(jac_affine_inv_1_0_BLUE,
                                                      jac_affine_inv_1_0_BLUE,
                                                      jac_affine_inv_1_0_BLUE,
                                                      jac_affine_inv_1_0_BLUE),
                                        _mm256_set_pd(jac_affine_inv_1_0_BLUE,
                                                      jac_affine_inv_1_0_BLUE,
                                                      jac_affine_inv_1_0_BLUE,
                                                      jac_affine_inv_1_0_BLUE)),
                          _mm256_mul_pd(
                              _mm256_set_pd(jac_affine_inv_1_1_BLUE,
                                            jac_affine_inv_1_1_BLUE,
                                            jac_affine_inv_1_1_BLUE,
                                            jac_affine_inv_1_1_BLUE),
                              _mm256_set_pd(jac_affine_inv_1_1_BLUE,
                                            jac_affine_inv_1_1_BLUE,
                                            jac_affine_inv_1_1_BLUE,
                                            jac_affine_inv_1_1_BLUE)))),
                  _mm256_mul_pd(tmp_moved_constant_8,
                                _mm256_set_pd(tmp_moved_constant_11,
                                              tmp_moved_constant_11,
                                              tmp_moved_constant_11,
                                              tmp_moved_constant_11))),
              _mm256_mul_pd(
                  tmp_moved_constant_6,
                  _mm256_set_pd(tmp_moved_constant_9, tmp_moved_constant_9,
                                tmp_moved_constant_9, tmp_moved_constant_9)));
          {
            {
              _mm256_storeu_pd(
                  &_data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) + 1],
                  _mm256_add_pd(
                      tmp_moved_constant_12,
                      _mm256_loadu_pd(
                          &_data_dst[ctr_0 +
                                     ctr_1 * (micro_edges_per_macro_edge + 2) -
                                     ((ctr_1 * (ctr_1 + 1)) / (2)) + 1])));
              _mm256_storeu_pd(
                  &_data_dst[ctr_0 +
                             (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                             (((ctr_1 + 1) * (ctr_1 + 2)) / (2))],
                  _mm256_add_pd(
                      tmp_moved_constant_13,
                      _mm256_loadu_pd(
                          &_data_dst[ctr_0 +
                                     (ctr_1 + 1) *
                                         (micro_edges_per_macro_edge + 2) -
                                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2))])));
              _mm256_storeu_pd(
                  &_data_dst[ctr_0 +
                             (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                             (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1],
                  _mm256_add_pd(
                      tmp_moved_constant_14,
                      _mm256_loadu_pd(
                          &_data_dst
                              [ctr_0 +
                               (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1])));
            }
          }
        }
        for (int64_t ctr_0 =
                 (int64_t)((-ctr_1 + micro_edges_per_macro_edge - 1) / (4)) *
                 (4);
             ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1; ctr_0 += 1) {
          const walberla::float64 src_dof_0 =
              _data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 src_dof_1 =
              _data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          const walberla::float64 src_dof_2 =
              _data_src[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          const walberla::float64 tmp_kernel_op_3 = src_dof_0 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_5 = src_dof_1 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_7 = src_dof_2 * tmp_kernel_op_2;
          const walberla::float64 elMatVec_0 =
              tmp_kernel_op_3 * ((tmp_kernel_op_0 * tmp_kernel_op_0) +
                                 (tmp_kernel_op_1 * tmp_kernel_op_1)) +
              tmp_kernel_op_4 * tmp_kernel_op_5 +
              tmp_kernel_op_6 * tmp_kernel_op_7;
          const walberla::float64 elMatVec_1 =
              tmp_kernel_op_3 * tmp_kernel_op_4 +
              tmp_kernel_op_5 *
                  ((jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY) +
                   (jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY)) +
              tmp_kernel_op_7 * tmp_kernel_op_8;
          const walberla::float64 elMatVec_2 =
              tmp_kernel_op_3 * tmp_kernel_op_6 +
              tmp_kernel_op_5 * tmp_kernel_op_8 +
              tmp_kernel_op_7 *
                  ((jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY) +
                   (jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY));
          {
            {
              _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2))] =
                  elMatVec_0 +
                  _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                            ((ctr_1 * (ctr_1 + 1)) / (2))];
              _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
                  elMatVec_1 +
                  _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
              _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
                  elMatVec_2 +
                  _data_dst[ctr_0 +
                            (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                            (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
            }
          }
          const walberla::float64 tmp_moved_constant_0 =
              _data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          const walberla::float64 tmp_moved_constant_1 =
              _data_src[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          const walberla::float64 tmp_moved_constant_2 =
              _data_src[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
          const walberla::float64 tmp_moved_constant_6 =
              tmp_moved_constant_0 * tmp_moved_constant_5;
          const walberla::float64 tmp_moved_constant_8 =
              tmp_moved_constant_1 * tmp_moved_constant_5;
          const walberla::float64 tmp_moved_constant_10 =
              tmp_moved_constant_2 * tmp_moved_constant_5;
          const walberla::float64 tmp_moved_constant_12 =
              tmp_moved_constant_10 * tmp_moved_constant_9 +
              tmp_moved_constant_6 *
                  ((tmp_moved_constant_3 * tmp_moved_constant_3) +
                   (tmp_moved_constant_4 * tmp_moved_constant_4)) +
              tmp_moved_constant_7 * tmp_moved_constant_8;
          const walberla::float64 tmp_moved_constant_13 =
              tmp_moved_constant_10 * tmp_moved_constant_11 +
              tmp_moved_constant_6 * tmp_moved_constant_7 +
              tmp_moved_constant_8 *
                  ((jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE) +
                   (jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE));
          const walberla::float64 tmp_moved_constant_14 =
              tmp_moved_constant_10 *
                  ((jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE) +
                   (jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE)) +
              tmp_moved_constant_11 * tmp_moved_constant_8 +
              tmp_moved_constant_6 * tmp_moved_constant_9;
          {
            {
              _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
                  tmp_moved_constant_12 +
                  _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
              _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
                  tmp_moved_constant_13 +
                  _data_dst[ctr_0 +
                            (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                            (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
              _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1] =
                  tmp_moved_constant_14 +
                  _data_dst[ctr_0 +
                            (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                            (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
            }
          }
        }
      }
      const walberla::float64 src_dof_0 =
          _data_src[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                    micro_edges_per_macro_edge - ((ctr_1 * (ctr_1 + 1)) / (2)) -
                    1];
      const walberla::float64 src_dof_1 =
          _data_src[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                    micro_edges_per_macro_edge - ((ctr_1 * (ctr_1 + 1)) / (2))];
      const walberla::float64 src_dof_2 =
          _data_src[-ctr_1 + micro_edges_per_macro_edge +
                    (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) - 1];
      const walberla::float64 tmp_kernel_op_3 = src_dof_0 * tmp_kernel_op_2;
      const walberla::float64 tmp_kernel_op_5 = src_dof_1 * tmp_kernel_op_2;
      const walberla::float64 tmp_kernel_op_7 = src_dof_2 * tmp_kernel_op_2;
      const walberla::float64 elMatVec_0 =
          tmp_kernel_op_3 * ((tmp_kernel_op_0 * tmp_kernel_op_0) +
                             (tmp_kernel_op_1 * tmp_kernel_op_1)) +
          tmp_kernel_op_4 * tmp_kernel_op_5 + tmp_kernel_op_6 * tmp_kernel_op_7;
      const walberla::float64 elMatVec_1 =
          tmp_kernel_op_3 * tmp_kernel_op_4 +
          tmp_kernel_op_5 *
              ((jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY) +
               (jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY)) +
          tmp_kernel_op_7 * tmp_kernel_op_8;
      const walberla::float64 elMatVec_2 =
          tmp_kernel_op_3 * tmp_kernel_op_6 +
          tmp_kernel_op_5 * tmp_kernel_op_8 +
          tmp_kernel_op_7 *
              ((jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY) +
               (jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY));
      {
        {
          {
            _data_dst[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                      micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) - 1] =
                elMatVec_0 +
                _data_dst[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                          micro_edges_per_macro_edge -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) - 1];
            _data_dst[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                      micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2))] =
                elMatVec_1 +
                _data_dst[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                          micro_edges_per_macro_edge -
                          ((ctr_1 * (ctr_1 + 1)) / (2))];
            _data_dst[-ctr_1 + micro_edges_per_macro_edge +
                      (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) - 1] =
                elMatVec_2 +
                _data_dst[-ctr_1 + micro_edges_per_macro_edge +
                          (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) - 1];
          }
        }
      }
    }
  }
}
void P1ElementwiseDiffusion_cubes_const_vect_float64::apply_macro_3D(
    walberla::float64 *RESTRICT _data_dst,
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
    const walberla::float64 tmp_coords_jac_0_GREEN_DOWN =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_GREEN_DOWN =
        -p_affine_const_0_0_GREEN_DOWN + p_affine_const_1_0_GREEN_DOWN;
    const walberla::float64 jac_affine_0_1_GREEN_DOWN =
        -p_affine_const_0_0_GREEN_DOWN + p_affine_const_2_0_GREEN_DOWN;
    const walberla::float64 jac_affine_0_2_GREEN_DOWN =
        -p_affine_const_0_0_GREEN_DOWN + p_affine_const_3_0_GREEN_DOWN;
    const walberla::float64 jac_affine_1_0_GREEN_DOWN =
        -p_affine_const_0_1_GREEN_DOWN + p_affine_const_1_1_GREEN_DOWN;
    const walberla::float64 jac_affine_1_1_GREEN_DOWN =
        -p_affine_const_0_1_GREEN_DOWN + p_affine_const_2_1_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_14_GREEN_DOWN =
        jac_affine_0_2_GREEN_DOWN * jac_affine_1_1_GREEN_DOWN;
    const walberla::float64 jac_affine_1_2_GREEN_DOWN =
        -p_affine_const_0_1_GREEN_DOWN + p_affine_const_3_1_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_12_GREEN_DOWN =
        jac_affine_0_1_GREEN_DOWN * jac_affine_1_2_GREEN_DOWN;
    const walberla::float64 jac_affine_2_0_GREEN_DOWN =
        -p_affine_const_0_2_GREEN_DOWN + p_affine_const_1_2_GREEN_DOWN;
    const walberla::float64 jac_affine_2_1_GREEN_DOWN =
        -p_affine_const_0_2_GREEN_DOWN + p_affine_const_2_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_11_GREEN_DOWN =
        jac_affine_1_2_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN;
    const walberla::float64 jac_affine_2_2_GREEN_DOWN =
        -p_affine_const_0_2_GREEN_DOWN + p_affine_const_3_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_10_GREEN_DOWN =
        jac_affine_1_1_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_13_GREEN_DOWN =
        jac_affine_0_1_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_15_GREEN_DOWN =
        jac_affine_0_0_GREEN_DOWN * tmp_coords_jac_10_GREEN_DOWN -
        jac_affine_0_0_GREEN_DOWN * tmp_coords_jac_11_GREEN_DOWN +
        jac_affine_0_2_GREEN_DOWN * jac_affine_1_0_GREEN_DOWN *
            jac_affine_2_1_GREEN_DOWN -
        jac_affine_1_0_GREEN_DOWN * tmp_coords_jac_13_GREEN_DOWN +
        jac_affine_2_0_GREEN_DOWN * tmp_coords_jac_12_GREEN_DOWN -
        jac_affine_2_0_GREEN_DOWN * tmp_coords_jac_14_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_16_GREEN_DOWN =
        1.0 / (tmp_coords_jac_15_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_0_0_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (tmp_coords_jac_10_GREEN_DOWN - tmp_coords_jac_11_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_0_1_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_0_2_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN -
         tmp_coords_jac_13_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_0_2_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (tmp_coords_jac_12_GREEN_DOWN - tmp_coords_jac_14_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_1_0_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (-jac_affine_1_0_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN +
         jac_affine_1_2_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_1_1_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_0_0_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN -
         jac_affine_0_2_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_1_2_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (-jac_affine_0_0_GREEN_DOWN * jac_affine_1_2_GREEN_DOWN +
         jac_affine_0_2_GREEN_DOWN * jac_affine_1_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_2_0_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_1_0_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN -
         jac_affine_1_1_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_2_1_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (-jac_affine_0_0_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN +
         jac_affine_0_1_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_2_2_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_0_0_GREEN_DOWN * jac_affine_1_1_GREEN_DOWN -
         jac_affine_0_1_GREEN_DOWN * jac_affine_1_0_GREEN_DOWN);
    const walberla::float64 abs_det_jac_affine_GREEN_DOWN =
        abs(tmp_coords_jac_15_GREEN_DOWN);
    const walberla::float64 tmp_coords_jac_0_GREEN_UP =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_GREEN_UP =
        -p_affine_const_0_0_GREEN_UP + p_affine_const_1_0_GREEN_UP;
    const walberla::float64 jac_affine_0_1_GREEN_UP =
        -p_affine_const_0_0_GREEN_UP + p_affine_const_2_0_GREEN_UP;
    const walberla::float64 jac_affine_0_2_GREEN_UP =
        -p_affine_const_0_0_GREEN_UP + p_affine_const_3_0_GREEN_UP;
    const walberla::float64 jac_affine_1_0_GREEN_UP =
        -p_affine_const_0_1_GREEN_UP + p_affine_const_1_1_GREEN_UP;
    const walberla::float64 jac_affine_1_1_GREEN_UP =
        -p_affine_const_0_1_GREEN_UP + p_affine_const_2_1_GREEN_UP;
    const walberla::float64 tmp_coords_jac_11_GREEN_UP =
        jac_affine_0_2_GREEN_UP * jac_affine_1_1_GREEN_UP;
    const walberla::float64 jac_affine_1_2_GREEN_UP =
        -p_affine_const_0_1_GREEN_UP + p_affine_const_3_1_GREEN_UP;
    const walberla::float64 tmp_coords_jac_9_GREEN_UP =
        jac_affine_0_1_GREEN_UP * jac_affine_1_2_GREEN_UP;
    const walberla::float64 jac_affine_2_0_GREEN_UP =
        -p_affine_const_0_2_GREEN_UP + p_affine_const_1_2_GREEN_UP;
    const walberla::float64 jac_affine_2_1_GREEN_UP =
        -p_affine_const_0_2_GREEN_UP + p_affine_const_2_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_8_GREEN_UP =
        jac_affine_1_2_GREEN_UP * jac_affine_2_1_GREEN_UP;
    const walberla::float64 jac_affine_2_2_GREEN_UP =
        -p_affine_const_0_2_GREEN_UP + p_affine_const_3_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_7_GREEN_UP =
        jac_affine_1_1_GREEN_UP * jac_affine_2_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_10_GREEN_UP =
        jac_affine_0_1_GREEN_UP * jac_affine_2_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_12_GREEN_UP =
        jac_affine_0_0_GREEN_UP * tmp_coords_jac_7_GREEN_UP -
        jac_affine_0_0_GREEN_UP * tmp_coords_jac_8_GREEN_UP +
        jac_affine_0_2_GREEN_UP * jac_affine_1_0_GREEN_UP *
            jac_affine_2_1_GREEN_UP -
        jac_affine_1_0_GREEN_UP * tmp_coords_jac_10_GREEN_UP -
        jac_affine_2_0_GREEN_UP * tmp_coords_jac_11_GREEN_UP +
        jac_affine_2_0_GREEN_UP * tmp_coords_jac_9_GREEN_UP;
    const walberla::float64 tmp_coords_jac_13_GREEN_UP =
        1.0 / (tmp_coords_jac_12_GREEN_UP);
    const walberla::float64 jac_affine_inv_0_0_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (tmp_coords_jac_7_GREEN_UP - tmp_coords_jac_8_GREEN_UP);
    const walberla::float64 jac_affine_inv_0_1_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_0_2_GREEN_UP * jac_affine_2_1_GREEN_UP -
         tmp_coords_jac_10_GREEN_UP);
    const walberla::float64 jac_affine_inv_0_2_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-tmp_coords_jac_11_GREEN_UP + tmp_coords_jac_9_GREEN_UP);
    const walberla::float64 jac_affine_inv_1_0_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-jac_affine_1_0_GREEN_UP * jac_affine_2_2_GREEN_UP +
         jac_affine_1_2_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_1_1_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_0_0_GREEN_UP * jac_affine_2_2_GREEN_UP -
         jac_affine_0_2_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_1_2_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-jac_affine_0_0_GREEN_UP * jac_affine_1_2_GREEN_UP +
         jac_affine_0_2_GREEN_UP * jac_affine_1_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_2_0_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_1_0_GREEN_UP * jac_affine_2_1_GREEN_UP -
         jac_affine_1_1_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_2_1_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-jac_affine_0_0_GREEN_UP * jac_affine_2_1_GREEN_UP +
         jac_affine_0_1_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_2_2_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_0_0_GREEN_UP * jac_affine_1_1_GREEN_UP -
         jac_affine_0_1_GREEN_UP * jac_affine_1_0_GREEN_UP);
    const walberla::float64 abs_det_jac_affine_GREEN_UP =
        abs(tmp_coords_jac_12_GREEN_UP);
    const walberla::float64 tmp_coords_jac_0_BLUE_DOWN =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_BLUE_DOWN =
        -p_affine_const_0_0_BLUE_DOWN + p_affine_const_1_0_BLUE_DOWN;
    const walberla::float64 jac_affine_0_1_BLUE_DOWN =
        -p_affine_const_0_0_BLUE_DOWN + p_affine_const_2_0_BLUE_DOWN;
    const walberla::float64 jac_affine_0_2_BLUE_DOWN =
        -p_affine_const_0_0_BLUE_DOWN + p_affine_const_3_0_BLUE_DOWN;
    const walberla::float64 jac_affine_1_0_BLUE_DOWN =
        -p_affine_const_0_1_BLUE_DOWN + p_affine_const_1_1_BLUE_DOWN;
    const walberla::float64 jac_affine_1_1_BLUE_DOWN =
        -p_affine_const_0_1_BLUE_DOWN + p_affine_const_2_1_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_14_BLUE_DOWN =
        jac_affine_0_2_BLUE_DOWN * jac_affine_1_1_BLUE_DOWN;
    const walberla::float64 jac_affine_1_2_BLUE_DOWN =
        -p_affine_const_0_1_BLUE_DOWN + p_affine_const_3_1_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_12_BLUE_DOWN =
        jac_affine_0_1_BLUE_DOWN * jac_affine_1_2_BLUE_DOWN;
    const walberla::float64 jac_affine_2_0_BLUE_DOWN =
        -p_affine_const_0_2_BLUE_DOWN + p_affine_const_1_2_BLUE_DOWN;
    const walberla::float64 jac_affine_2_1_BLUE_DOWN =
        -p_affine_const_0_2_BLUE_DOWN + p_affine_const_2_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_11_BLUE_DOWN =
        jac_affine_1_2_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN;
    const walberla::float64 jac_affine_2_2_BLUE_DOWN =
        -p_affine_const_0_2_BLUE_DOWN + p_affine_const_3_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_10_BLUE_DOWN =
        jac_affine_1_1_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_13_BLUE_DOWN =
        jac_affine_0_1_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_15_BLUE_DOWN =
        jac_affine_0_0_BLUE_DOWN * tmp_coords_jac_10_BLUE_DOWN -
        jac_affine_0_0_BLUE_DOWN * tmp_coords_jac_11_BLUE_DOWN +
        jac_affine_0_2_BLUE_DOWN * jac_affine_1_0_BLUE_DOWN *
            jac_affine_2_1_BLUE_DOWN -
        jac_affine_1_0_BLUE_DOWN * tmp_coords_jac_13_BLUE_DOWN +
        jac_affine_2_0_BLUE_DOWN * tmp_coords_jac_12_BLUE_DOWN -
        jac_affine_2_0_BLUE_DOWN * tmp_coords_jac_14_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_16_BLUE_DOWN =
        1.0 / (tmp_coords_jac_15_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_0_0_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (tmp_coords_jac_10_BLUE_DOWN - tmp_coords_jac_11_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_0_1_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_0_2_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN -
         tmp_coords_jac_13_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_0_2_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (tmp_coords_jac_12_BLUE_DOWN - tmp_coords_jac_14_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_1_0_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (-jac_affine_1_0_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN +
         jac_affine_1_2_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_1_1_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_0_0_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN -
         jac_affine_0_2_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_1_2_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (-jac_affine_0_0_BLUE_DOWN * jac_affine_1_2_BLUE_DOWN +
         jac_affine_0_2_BLUE_DOWN * jac_affine_1_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_2_0_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_1_0_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN -
         jac_affine_1_1_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_2_1_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (-jac_affine_0_0_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN +
         jac_affine_0_1_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_2_2_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_0_0_BLUE_DOWN * jac_affine_1_1_BLUE_DOWN -
         jac_affine_0_1_BLUE_DOWN * jac_affine_1_0_BLUE_DOWN);
    const walberla::float64 abs_det_jac_affine_BLUE_DOWN =
        abs(tmp_coords_jac_15_BLUE_DOWN);
    const walberla::float64 tmp_coords_jac_0_BLUE_UP =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_BLUE_UP =
        -p_affine_const_0_0_BLUE_UP + p_affine_const_1_0_BLUE_UP;
    const walberla::float64 jac_affine_0_1_BLUE_UP =
        -p_affine_const_0_0_BLUE_UP + p_affine_const_2_0_BLUE_UP;
    const walberla::float64 jac_affine_0_2_BLUE_UP =
        -p_affine_const_0_0_BLUE_UP + p_affine_const_3_0_BLUE_UP;
    const walberla::float64 jac_affine_1_0_BLUE_UP =
        -p_affine_const_0_1_BLUE_UP + p_affine_const_1_1_BLUE_UP;
    const walberla::float64 jac_affine_1_1_BLUE_UP =
        -p_affine_const_0_1_BLUE_UP + p_affine_const_2_1_BLUE_UP;
    const walberla::float64 tmp_coords_jac_11_BLUE_UP =
        jac_affine_0_2_BLUE_UP * jac_affine_1_1_BLUE_UP;
    const walberla::float64 jac_affine_1_2_BLUE_UP =
        -p_affine_const_0_1_BLUE_UP + p_affine_const_3_1_BLUE_UP;
    const walberla::float64 tmp_coords_jac_9_BLUE_UP =
        jac_affine_0_1_BLUE_UP * jac_affine_1_2_BLUE_UP;
    const walberla::float64 jac_affine_2_0_BLUE_UP =
        -p_affine_const_0_2_BLUE_UP + p_affine_const_1_2_BLUE_UP;
    const walberla::float64 jac_affine_2_1_BLUE_UP =
        -p_affine_const_0_2_BLUE_UP + p_affine_const_2_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_8_BLUE_UP =
        jac_affine_1_2_BLUE_UP * jac_affine_2_1_BLUE_UP;
    const walberla::float64 jac_affine_2_2_BLUE_UP =
        -p_affine_const_0_2_BLUE_UP + p_affine_const_3_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_7_BLUE_UP =
        jac_affine_1_1_BLUE_UP * jac_affine_2_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_10_BLUE_UP =
        jac_affine_0_1_BLUE_UP * jac_affine_2_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_12_BLUE_UP =
        jac_affine_0_0_BLUE_UP * tmp_coords_jac_7_BLUE_UP -
        jac_affine_0_0_BLUE_UP * tmp_coords_jac_8_BLUE_UP +
        jac_affine_0_2_BLUE_UP * jac_affine_1_0_BLUE_UP *
            jac_affine_2_1_BLUE_UP -
        jac_affine_1_0_BLUE_UP * tmp_coords_jac_10_BLUE_UP -
        jac_affine_2_0_BLUE_UP * tmp_coords_jac_11_BLUE_UP +
        jac_affine_2_0_BLUE_UP * tmp_coords_jac_9_BLUE_UP;
    const walberla::float64 tmp_coords_jac_13_BLUE_UP =
        1.0 / (tmp_coords_jac_12_BLUE_UP);
    const walberla::float64 jac_affine_inv_0_0_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (tmp_coords_jac_7_BLUE_UP - tmp_coords_jac_8_BLUE_UP);
    const walberla::float64 jac_affine_inv_0_1_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_0_2_BLUE_UP * jac_affine_2_1_BLUE_UP -
         tmp_coords_jac_10_BLUE_UP);
    const walberla::float64 jac_affine_inv_0_2_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-tmp_coords_jac_11_BLUE_UP + tmp_coords_jac_9_BLUE_UP);
    const walberla::float64 jac_affine_inv_1_0_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-jac_affine_1_0_BLUE_UP * jac_affine_2_2_BLUE_UP +
         jac_affine_1_2_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_1_1_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_0_0_BLUE_UP * jac_affine_2_2_BLUE_UP -
         jac_affine_0_2_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_1_2_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-jac_affine_0_0_BLUE_UP * jac_affine_1_2_BLUE_UP +
         jac_affine_0_2_BLUE_UP * jac_affine_1_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_2_0_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_1_0_BLUE_UP * jac_affine_2_1_BLUE_UP -
         jac_affine_1_1_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_2_1_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-jac_affine_0_0_BLUE_UP * jac_affine_2_1_BLUE_UP +
         jac_affine_0_1_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_2_2_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_0_0_BLUE_UP * jac_affine_1_1_BLUE_UP -
         jac_affine_0_1_BLUE_UP * jac_affine_1_0_BLUE_UP);
    const walberla::float64 abs_det_jac_affine_BLUE_UP =
        abs(tmp_coords_jac_12_BLUE_UP);
    const walberla::float64 tmp_coords_jac_0_WHITE_DOWN =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_WHITE_DOWN =
        -p_affine_const_0_0_WHITE_DOWN + p_affine_const_1_0_WHITE_DOWN;
    const walberla::float64 jac_affine_0_1_WHITE_DOWN =
        -p_affine_const_0_0_WHITE_DOWN + p_affine_const_2_0_WHITE_DOWN;
    const walberla::float64 jac_affine_0_2_WHITE_DOWN =
        -p_affine_const_0_0_WHITE_DOWN + p_affine_const_3_0_WHITE_DOWN;
    const walberla::float64 jac_affine_1_0_WHITE_DOWN =
        -p_affine_const_0_1_WHITE_DOWN + p_affine_const_1_1_WHITE_DOWN;
    const walberla::float64 jac_affine_1_1_WHITE_DOWN =
        -p_affine_const_0_1_WHITE_DOWN + p_affine_const_2_1_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_17_WHITE_DOWN =
        jac_affine_0_2_WHITE_DOWN * jac_affine_1_1_WHITE_DOWN;
    const walberla::float64 jac_affine_1_2_WHITE_DOWN =
        -p_affine_const_0_1_WHITE_DOWN + p_affine_const_3_1_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_15_WHITE_DOWN =
        jac_affine_0_1_WHITE_DOWN * jac_affine_1_2_WHITE_DOWN;
    const walberla::float64 jac_affine_2_0_WHITE_DOWN =
        -p_affine_const_0_2_WHITE_DOWN + p_affine_const_1_2_WHITE_DOWN;
    const walberla::float64 jac_affine_2_1_WHITE_DOWN =
        -p_affine_const_0_2_WHITE_DOWN + p_affine_const_2_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_14_WHITE_DOWN =
        jac_affine_1_2_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN;
    const walberla::float64 jac_affine_2_2_WHITE_DOWN =
        -p_affine_const_0_2_WHITE_DOWN + p_affine_const_3_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_13_WHITE_DOWN =
        jac_affine_1_1_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_16_WHITE_DOWN =
        jac_affine_0_1_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_18_WHITE_DOWN =
        jac_affine_0_0_WHITE_DOWN * tmp_coords_jac_13_WHITE_DOWN -
        jac_affine_0_0_WHITE_DOWN * tmp_coords_jac_14_WHITE_DOWN +
        jac_affine_0_2_WHITE_DOWN * jac_affine_1_0_WHITE_DOWN *
            jac_affine_2_1_WHITE_DOWN -
        jac_affine_1_0_WHITE_DOWN * tmp_coords_jac_16_WHITE_DOWN +
        jac_affine_2_0_WHITE_DOWN * tmp_coords_jac_15_WHITE_DOWN -
        jac_affine_2_0_WHITE_DOWN * tmp_coords_jac_17_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_19_WHITE_DOWN =
        1.0 / (tmp_coords_jac_18_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_0_0_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (tmp_coords_jac_13_WHITE_DOWN - tmp_coords_jac_14_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_0_1_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_0_2_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN -
         tmp_coords_jac_16_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_0_2_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (tmp_coords_jac_15_WHITE_DOWN - tmp_coords_jac_17_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_1_0_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (-jac_affine_1_0_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN +
         jac_affine_1_2_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_1_1_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_0_0_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN -
         jac_affine_0_2_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_1_2_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (-jac_affine_0_0_WHITE_DOWN * jac_affine_1_2_WHITE_DOWN +
         jac_affine_0_2_WHITE_DOWN * jac_affine_1_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_2_0_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_1_0_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN -
         jac_affine_1_1_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_2_1_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (-jac_affine_0_0_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN +
         jac_affine_0_1_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_2_2_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_0_0_WHITE_DOWN * jac_affine_1_1_WHITE_DOWN -
         jac_affine_0_1_WHITE_DOWN * jac_affine_1_0_WHITE_DOWN);
    const walberla::float64 abs_det_jac_affine_WHITE_DOWN =
        abs(tmp_coords_jac_18_WHITE_DOWN);
    const walberla::float64 tmp_coords_jac_0_WHITE_UP =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_WHITE_UP =
        -p_affine_const_0_0_WHITE_UP + p_affine_const_1_0_WHITE_UP;
    const walberla::float64 jac_affine_0_1_WHITE_UP =
        -p_affine_const_0_0_WHITE_UP + p_affine_const_2_0_WHITE_UP;
    const walberla::float64 jac_affine_0_2_WHITE_UP =
        -p_affine_const_0_0_WHITE_UP + p_affine_const_3_0_WHITE_UP;
    const walberla::float64 jac_affine_1_0_WHITE_UP =
        -p_affine_const_0_1_WHITE_UP + p_affine_const_1_1_WHITE_UP;
    const walberla::float64 jac_affine_1_1_WHITE_UP =
        -p_affine_const_0_1_WHITE_UP + p_affine_const_2_1_WHITE_UP;
    const walberla::float64 tmp_coords_jac_5_WHITE_UP =
        jac_affine_0_2_WHITE_UP * jac_affine_1_1_WHITE_UP;
    const walberla::float64 jac_affine_1_2_WHITE_UP =
        -p_affine_const_0_1_WHITE_UP + p_affine_const_3_1_WHITE_UP;
    const walberla::float64 tmp_coords_jac_3_WHITE_UP =
        jac_affine_0_1_WHITE_UP * jac_affine_1_2_WHITE_UP;
    const walberla::float64 jac_affine_2_0_WHITE_UP =
        -p_affine_const_0_2_WHITE_UP + p_affine_const_1_2_WHITE_UP;
    const walberla::float64 jac_affine_2_1_WHITE_UP =
        -p_affine_const_0_2_WHITE_UP + p_affine_const_2_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_2_WHITE_UP =
        jac_affine_1_2_WHITE_UP * jac_affine_2_1_WHITE_UP;
    const walberla::float64 jac_affine_2_2_WHITE_UP =
        -p_affine_const_0_2_WHITE_UP + p_affine_const_3_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_1_WHITE_UP =
        jac_affine_1_1_WHITE_UP * jac_affine_2_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_4_WHITE_UP =
        jac_affine_0_1_WHITE_UP * jac_affine_2_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_6_WHITE_UP =
        jac_affine_0_0_WHITE_UP * tmp_coords_jac_1_WHITE_UP -
        jac_affine_0_0_WHITE_UP * tmp_coords_jac_2_WHITE_UP +
        jac_affine_0_2_WHITE_UP * jac_affine_1_0_WHITE_UP *
            jac_affine_2_1_WHITE_UP -
        jac_affine_1_0_WHITE_UP * tmp_coords_jac_4_WHITE_UP +
        jac_affine_2_0_WHITE_UP * tmp_coords_jac_3_WHITE_UP -
        jac_affine_2_0_WHITE_UP * tmp_coords_jac_5_WHITE_UP;
    const walberla::float64 tmp_coords_jac_7_WHITE_UP =
        1.0 / (tmp_coords_jac_6_WHITE_UP);
    const walberla::float64 jac_affine_inv_0_0_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (tmp_coords_jac_1_WHITE_UP - tmp_coords_jac_2_WHITE_UP);
    const walberla::float64 jac_affine_inv_0_1_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_0_2_WHITE_UP * jac_affine_2_1_WHITE_UP -
         tmp_coords_jac_4_WHITE_UP);
    const walberla::float64 jac_affine_inv_0_2_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (tmp_coords_jac_3_WHITE_UP - tmp_coords_jac_5_WHITE_UP);
    const walberla::float64 jac_affine_inv_1_0_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (-jac_affine_1_0_WHITE_UP * jac_affine_2_2_WHITE_UP +
         jac_affine_1_2_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_1_1_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_0_0_WHITE_UP * jac_affine_2_2_WHITE_UP -
         jac_affine_0_2_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_1_2_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (-jac_affine_0_0_WHITE_UP * jac_affine_1_2_WHITE_UP +
         jac_affine_0_2_WHITE_UP * jac_affine_1_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_2_0_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_1_0_WHITE_UP * jac_affine_2_1_WHITE_UP -
         jac_affine_1_1_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_2_1_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (-jac_affine_0_0_WHITE_UP * jac_affine_2_1_WHITE_UP +
         jac_affine_0_1_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_2_2_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_0_0_WHITE_UP * jac_affine_1_1_WHITE_UP -
         jac_affine_0_1_WHITE_UP * jac_affine_1_0_WHITE_UP);
    const walberla::float64 abs_det_jac_affine_WHITE_UP =
        abs(tmp_coords_jac_6_WHITE_UP);
    const walberla::float64 tmp_kernel_op_0 = -jac_affine_inv_0_0_WHITE_UP -
                                              jac_affine_inv_1_0_WHITE_UP -
                                              jac_affine_inv_2_0_WHITE_UP;
    const walberla::float64 tmp_kernel_op_1 = -jac_affine_inv_0_1_WHITE_UP -
                                              jac_affine_inv_1_1_WHITE_UP -
                                              jac_affine_inv_2_1_WHITE_UP;
    const walberla::float64 tmp_kernel_op_2 = -jac_affine_inv_0_2_WHITE_UP -
                                              jac_affine_inv_1_2_WHITE_UP -
                                              jac_affine_inv_2_2_WHITE_UP;
    const walberla::float64 tmp_kernel_op_3 =
        abs_det_jac_affine_WHITE_UP * 0.16666666666666663;
    const walberla::float64 tmp_kernel_op_5 =
        jac_affine_inv_0_0_WHITE_UP * tmp_kernel_op_0 +
        jac_affine_inv_0_1_WHITE_UP * tmp_kernel_op_1 +
        jac_affine_inv_0_2_WHITE_UP * tmp_kernel_op_2;
    const walberla::float64 tmp_kernel_op_7 =
        jac_affine_inv_1_0_WHITE_UP * tmp_kernel_op_0 +
        jac_affine_inv_1_1_WHITE_UP * tmp_kernel_op_1 +
        jac_affine_inv_1_2_WHITE_UP * tmp_kernel_op_2;
    const walberla::float64 tmp_kernel_op_9 =
        jac_affine_inv_2_0_WHITE_UP * tmp_kernel_op_0 +
        jac_affine_inv_2_1_WHITE_UP * tmp_kernel_op_1 +
        jac_affine_inv_2_2_WHITE_UP * tmp_kernel_op_2;
    const walberla::float64 tmp_kernel_op_11 =
        jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_1_0_WHITE_UP +
        jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_1_1_WHITE_UP +
        jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_1_2_WHITE_UP;
    const walberla::float64 tmp_kernel_op_12 =
        jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP +
        jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP +
        jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_2_2_WHITE_UP;
    const walberla::float64 tmp_kernel_op_13 =
        jac_affine_inv_1_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP +
        jac_affine_inv_1_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP +
        jac_affine_inv_1_2_WHITE_UP * jac_affine_inv_2_2_WHITE_UP;
    const walberla::float64 tmp_moved_constant_4 =
        -jac_affine_inv_0_0_WHITE_DOWN - jac_affine_inv_1_0_WHITE_DOWN -
        jac_affine_inv_2_0_WHITE_DOWN;
    const walberla::float64 tmp_moved_constant_5 =
        -jac_affine_inv_0_1_WHITE_DOWN - jac_affine_inv_1_1_WHITE_DOWN -
        jac_affine_inv_2_1_WHITE_DOWN;
    const walberla::float64 tmp_moved_constant_6 =
        -jac_affine_inv_0_2_WHITE_DOWN - jac_affine_inv_1_2_WHITE_DOWN -
        jac_affine_inv_2_2_WHITE_DOWN;
    const walberla::float64 tmp_moved_constant_7 =
        abs_det_jac_affine_WHITE_DOWN * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_9 =
        jac_affine_inv_0_0_WHITE_DOWN * tmp_moved_constant_4 +
        jac_affine_inv_0_1_WHITE_DOWN * tmp_moved_constant_5 +
        jac_affine_inv_0_2_WHITE_DOWN * tmp_moved_constant_6;
    const walberla::float64 tmp_moved_constant_11 =
        jac_affine_inv_1_0_WHITE_DOWN * tmp_moved_constant_4 +
        jac_affine_inv_1_1_WHITE_DOWN * tmp_moved_constant_5 +
        jac_affine_inv_1_2_WHITE_DOWN * tmp_moved_constant_6;
    const walberla::float64 tmp_moved_constant_13 =
        jac_affine_inv_2_0_WHITE_DOWN * tmp_moved_constant_4 +
        jac_affine_inv_2_1_WHITE_DOWN * tmp_moved_constant_5 +
        jac_affine_inv_2_2_WHITE_DOWN * tmp_moved_constant_6;
    const walberla::float64 tmp_moved_constant_15 =
        jac_affine_inv_0_0_WHITE_DOWN * jac_affine_inv_1_0_WHITE_DOWN +
        jac_affine_inv_0_1_WHITE_DOWN * jac_affine_inv_1_1_WHITE_DOWN +
        jac_affine_inv_0_2_WHITE_DOWN * jac_affine_inv_1_2_WHITE_DOWN;
    const walberla::float64 tmp_moved_constant_16 =
        jac_affine_inv_0_0_WHITE_DOWN * jac_affine_inv_2_0_WHITE_DOWN +
        jac_affine_inv_0_1_WHITE_DOWN * jac_affine_inv_2_1_WHITE_DOWN +
        jac_affine_inv_0_2_WHITE_DOWN * jac_affine_inv_2_2_WHITE_DOWN;
    const walberla::float64 tmp_moved_constant_17 =
        jac_affine_inv_1_0_WHITE_DOWN * jac_affine_inv_2_0_WHITE_DOWN +
        jac_affine_inv_1_1_WHITE_DOWN * jac_affine_inv_2_1_WHITE_DOWN +
        jac_affine_inv_1_2_WHITE_DOWN * jac_affine_inv_2_2_WHITE_DOWN;
    const walberla::float64 tmp_moved_constant_26 =
        -jac_affine_inv_0_0_BLUE_UP - jac_affine_inv_1_0_BLUE_UP -
        jac_affine_inv_2_0_BLUE_UP;
    const walberla::float64 tmp_moved_constant_27 =
        -jac_affine_inv_0_1_BLUE_UP - jac_affine_inv_1_1_BLUE_UP -
        jac_affine_inv_2_1_BLUE_UP;
    const walberla::float64 tmp_moved_constant_28 =
        -jac_affine_inv_0_2_BLUE_UP - jac_affine_inv_1_2_BLUE_UP -
        jac_affine_inv_2_2_BLUE_UP;
    const walberla::float64 tmp_moved_constant_29 =
        abs_det_jac_affine_BLUE_UP * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_31 =
        jac_affine_inv_0_0_BLUE_UP * tmp_moved_constant_26 +
        jac_affine_inv_0_1_BLUE_UP * tmp_moved_constant_27 +
        jac_affine_inv_0_2_BLUE_UP * tmp_moved_constant_28;
    const walberla::float64 tmp_moved_constant_33 =
        jac_affine_inv_1_0_BLUE_UP * tmp_moved_constant_26 +
        jac_affine_inv_1_1_BLUE_UP * tmp_moved_constant_27 +
        jac_affine_inv_1_2_BLUE_UP * tmp_moved_constant_28;
    const walberla::float64 tmp_moved_constant_35 =
        jac_affine_inv_2_0_BLUE_UP * tmp_moved_constant_26 +
        jac_affine_inv_2_1_BLUE_UP * tmp_moved_constant_27 +
        jac_affine_inv_2_2_BLUE_UP * tmp_moved_constant_28;
    const walberla::float64 tmp_moved_constant_37 =
        jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_1_0_BLUE_UP +
        jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_1_1_BLUE_UP +
        jac_affine_inv_0_2_BLUE_UP * jac_affine_inv_1_2_BLUE_UP;
    const walberla::float64 tmp_moved_constant_38 =
        jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP +
        jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP +
        jac_affine_inv_0_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP;
    const walberla::float64 tmp_moved_constant_39 =
        jac_affine_inv_1_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP +
        jac_affine_inv_1_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP +
        jac_affine_inv_1_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP;
    const walberla::float64 tmp_moved_constant_48 =
        -jac_affine_inv_0_0_BLUE_DOWN - jac_affine_inv_1_0_BLUE_DOWN -
        jac_affine_inv_2_0_BLUE_DOWN;
    const walberla::float64 tmp_moved_constant_49 =
        -jac_affine_inv_0_1_BLUE_DOWN - jac_affine_inv_1_1_BLUE_DOWN -
        jac_affine_inv_2_1_BLUE_DOWN;
    const walberla::float64 tmp_moved_constant_50 =
        -jac_affine_inv_0_2_BLUE_DOWN - jac_affine_inv_1_2_BLUE_DOWN -
        jac_affine_inv_2_2_BLUE_DOWN;
    const walberla::float64 tmp_moved_constant_51 =
        abs_det_jac_affine_BLUE_DOWN * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_53 =
        jac_affine_inv_0_0_BLUE_DOWN * tmp_moved_constant_48 +
        jac_affine_inv_0_1_BLUE_DOWN * tmp_moved_constant_49 +
        jac_affine_inv_0_2_BLUE_DOWN * tmp_moved_constant_50;
    const walberla::float64 tmp_moved_constant_55 =
        jac_affine_inv_1_0_BLUE_DOWN * tmp_moved_constant_48 +
        jac_affine_inv_1_1_BLUE_DOWN * tmp_moved_constant_49 +
        jac_affine_inv_1_2_BLUE_DOWN * tmp_moved_constant_50;
    const walberla::float64 tmp_moved_constant_57 =
        jac_affine_inv_2_0_BLUE_DOWN * tmp_moved_constant_48 +
        jac_affine_inv_2_1_BLUE_DOWN * tmp_moved_constant_49 +
        jac_affine_inv_2_2_BLUE_DOWN * tmp_moved_constant_50;
    const walberla::float64 tmp_moved_constant_59 =
        jac_affine_inv_0_0_BLUE_DOWN * jac_affine_inv_1_0_BLUE_DOWN +
        jac_affine_inv_0_1_BLUE_DOWN * jac_affine_inv_1_1_BLUE_DOWN +
        jac_affine_inv_0_2_BLUE_DOWN * jac_affine_inv_1_2_BLUE_DOWN;
    const walberla::float64 tmp_moved_constant_60 =
        jac_affine_inv_0_0_BLUE_DOWN * jac_affine_inv_2_0_BLUE_DOWN +
        jac_affine_inv_0_1_BLUE_DOWN * jac_affine_inv_2_1_BLUE_DOWN +
        jac_affine_inv_0_2_BLUE_DOWN * jac_affine_inv_2_2_BLUE_DOWN;
    const walberla::float64 tmp_moved_constant_61 =
        jac_affine_inv_1_0_BLUE_DOWN * jac_affine_inv_2_0_BLUE_DOWN +
        jac_affine_inv_1_1_BLUE_DOWN * jac_affine_inv_2_1_BLUE_DOWN +
        jac_affine_inv_1_2_BLUE_DOWN * jac_affine_inv_2_2_BLUE_DOWN;
    const walberla::float64 tmp_moved_constant_70 =
        -jac_affine_inv_0_0_GREEN_UP - jac_affine_inv_1_0_GREEN_UP -
        jac_affine_inv_2_0_GREEN_UP;
    const walberla::float64 tmp_moved_constant_71 =
        -jac_affine_inv_0_1_GREEN_UP - jac_affine_inv_1_1_GREEN_UP -
        jac_affine_inv_2_1_GREEN_UP;
    const walberla::float64 tmp_moved_constant_72 =
        -jac_affine_inv_0_2_GREEN_UP - jac_affine_inv_1_2_GREEN_UP -
        jac_affine_inv_2_2_GREEN_UP;
    const walberla::float64 tmp_moved_constant_73 =
        abs_det_jac_affine_GREEN_UP * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_75 =
        jac_affine_inv_0_0_GREEN_UP * tmp_moved_constant_70 +
        jac_affine_inv_0_1_GREEN_UP * tmp_moved_constant_71 +
        jac_affine_inv_0_2_GREEN_UP * tmp_moved_constant_72;
    const walberla::float64 tmp_moved_constant_77 =
        jac_affine_inv_1_0_GREEN_UP * tmp_moved_constant_70 +
        jac_affine_inv_1_1_GREEN_UP * tmp_moved_constant_71 +
        jac_affine_inv_1_2_GREEN_UP * tmp_moved_constant_72;
    const walberla::float64 tmp_moved_constant_79 =
        jac_affine_inv_2_0_GREEN_UP * tmp_moved_constant_70 +
        jac_affine_inv_2_1_GREEN_UP * tmp_moved_constant_71 +
        jac_affine_inv_2_2_GREEN_UP * tmp_moved_constant_72;
    const walberla::float64 tmp_moved_constant_81 =
        jac_affine_inv_0_0_GREEN_UP * jac_affine_inv_1_0_GREEN_UP +
        jac_affine_inv_0_1_GREEN_UP * jac_affine_inv_1_1_GREEN_UP +
        jac_affine_inv_0_2_GREEN_UP * jac_affine_inv_1_2_GREEN_UP;
    const walberla::float64 tmp_moved_constant_82 =
        jac_affine_inv_0_0_GREEN_UP * jac_affine_inv_2_0_GREEN_UP +
        jac_affine_inv_0_1_GREEN_UP * jac_affine_inv_2_1_GREEN_UP +
        jac_affine_inv_0_2_GREEN_UP * jac_affine_inv_2_2_GREEN_UP;
    const walberla::float64 tmp_moved_constant_83 =
        jac_affine_inv_1_0_GREEN_UP * jac_affine_inv_2_0_GREEN_UP +
        jac_affine_inv_1_1_GREEN_UP * jac_affine_inv_2_1_GREEN_UP +
        jac_affine_inv_1_2_GREEN_UP * jac_affine_inv_2_2_GREEN_UP;
    const walberla::float64 tmp_moved_constant_92 =
        -jac_affine_inv_0_0_GREEN_DOWN - jac_affine_inv_1_0_GREEN_DOWN -
        jac_affine_inv_2_0_GREEN_DOWN;
    const walberla::float64 tmp_moved_constant_93 =
        -jac_affine_inv_0_1_GREEN_DOWN - jac_affine_inv_1_1_GREEN_DOWN -
        jac_affine_inv_2_1_GREEN_DOWN;
    const walberla::float64 tmp_moved_constant_94 =
        -jac_affine_inv_0_2_GREEN_DOWN - jac_affine_inv_1_2_GREEN_DOWN -
        jac_affine_inv_2_2_GREEN_DOWN;
    const walberla::float64 tmp_moved_constant_95 =
        abs_det_jac_affine_GREEN_DOWN * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_97 =
        jac_affine_inv_0_0_GREEN_DOWN * tmp_moved_constant_92 +
        jac_affine_inv_0_1_GREEN_DOWN * tmp_moved_constant_93 +
        jac_affine_inv_0_2_GREEN_DOWN * tmp_moved_constant_94;
    const walberla::float64 tmp_moved_constant_99 =
        jac_affine_inv_1_0_GREEN_DOWN * tmp_moved_constant_92 +
        jac_affine_inv_1_1_GREEN_DOWN * tmp_moved_constant_93 +
        jac_affine_inv_1_2_GREEN_DOWN * tmp_moved_constant_94;
    const walberla::float64 tmp_moved_constant_101 =
        jac_affine_inv_2_0_GREEN_DOWN * tmp_moved_constant_92 +
        jac_affine_inv_2_1_GREEN_DOWN * tmp_moved_constant_93 +
        jac_affine_inv_2_2_GREEN_DOWN * tmp_moved_constant_94;
    const walberla::float64 tmp_moved_constant_103 =
        jac_affine_inv_0_0_GREEN_DOWN * jac_affine_inv_1_0_GREEN_DOWN +
        jac_affine_inv_0_1_GREEN_DOWN * jac_affine_inv_1_1_GREEN_DOWN +
        jac_affine_inv_0_2_GREEN_DOWN * jac_affine_inv_1_2_GREEN_DOWN;
    const walberla::float64 tmp_moved_constant_104 =
        jac_affine_inv_0_0_GREEN_DOWN * jac_affine_inv_2_0_GREEN_DOWN +
        jac_affine_inv_0_1_GREEN_DOWN * jac_affine_inv_2_1_GREEN_DOWN +
        jac_affine_inv_0_2_GREEN_DOWN * jac_affine_inv_2_2_GREEN_DOWN;
    const walberla::float64 tmp_moved_constant_105 =
        jac_affine_inv_1_0_GREEN_DOWN * jac_affine_inv_2_0_GREEN_DOWN +
        jac_affine_inv_1_1_GREEN_DOWN * jac_affine_inv_2_1_GREEN_DOWN +
        jac_affine_inv_1_2_GREEN_DOWN * jac_affine_inv_2_2_GREEN_DOWN;
    for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
      for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
           ctr_1 += 1) {
        {
          for (int64_t ctr_0 = 0;
               ctr_0 <
               (int64_t)((-ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2) /
                         (4)) *
                   (4);
               ctr_0 += 4) {
            const __m256d src_dof_0 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6)) -
                           (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2) *
                             (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                            (6))]);
            const __m256d src_dof_1 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d src_dof_2 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d src_dof_3 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))]);
            const __m256d tmp_kernel_op_4 = _mm256_mul_pd(
                src_dof_0, _mm256_set_pd(tmp_kernel_op_3, tmp_kernel_op_3,
                                         tmp_kernel_op_3, tmp_kernel_op_3));
            const __m256d tmp_kernel_op_6 = _mm256_mul_pd(
                src_dof_1, _mm256_set_pd(tmp_kernel_op_3, tmp_kernel_op_3,
                                         tmp_kernel_op_3, tmp_kernel_op_3));
            const __m256d tmp_kernel_op_8 = _mm256_mul_pd(
                src_dof_2, _mm256_set_pd(tmp_kernel_op_3, tmp_kernel_op_3,
                                         tmp_kernel_op_3, tmp_kernel_op_3));
            const __m256d tmp_kernel_op_10 = _mm256_mul_pd(
                src_dof_3, _mm256_set_pd(tmp_kernel_op_3, tmp_kernel_op_3,
                                         tmp_kernel_op_3, tmp_kernel_op_3));
            const __m256d elMatVec_0 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_kernel_op_4,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            tmp_kernel_op_0, tmp_kernel_op_0,
                                            tmp_kernel_op_0, tmp_kernel_op_0),
                                        _mm256_set_pd(
                                            tmp_kernel_op_0, tmp_kernel_op_0,
                                            tmp_kernel_op_0, tmp_kernel_op_0)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            tmp_kernel_op_1, tmp_kernel_op_1,
                                            tmp_kernel_op_1, tmp_kernel_op_1),
                                        _mm256_set_pd(
                                            tmp_kernel_op_1, tmp_kernel_op_1,
                                            tmp_kernel_op_1, tmp_kernel_op_1))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(
                                        tmp_kernel_op_2, tmp_kernel_op_2,
                                        tmp_kernel_op_2, tmp_kernel_op_2),
                                    _mm256_set_pd(
                                        tmp_kernel_op_2, tmp_kernel_op_2,
                                        tmp_kernel_op_2, tmp_kernel_op_2)))),
                        _mm256_mul_pd(
                            tmp_kernel_op_6,
                            _mm256_set_pd(tmp_kernel_op_5, tmp_kernel_op_5,
                                          tmp_kernel_op_5, tmp_kernel_op_5))),
                    _mm256_mul_pd(
                        tmp_kernel_op_8,
                        _mm256_set_pd(tmp_kernel_op_7, tmp_kernel_op_7,
                                      tmp_kernel_op_7, tmp_kernel_op_7))),
                _mm256_mul_pd(tmp_kernel_op_10,
                              _mm256_set_pd(tmp_kernel_op_9, tmp_kernel_op_9,
                                            tmp_kernel_op_9, tmp_kernel_op_9)));
            const __m256d elMatVec_1 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_kernel_op_6,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_WHITE_UP,
                                            jac_affine_inv_0_0_WHITE_UP,
                                            jac_affine_inv_0_0_WHITE_UP,
                                            jac_affine_inv_0_0_WHITE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_WHITE_UP,
                                            jac_affine_inv_0_0_WHITE_UP,
                                            jac_affine_inv_0_0_WHITE_UP,
                                            jac_affine_inv_0_0_WHITE_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_WHITE_UP,
                                            jac_affine_inv_0_1_WHITE_UP,
                                            jac_affine_inv_0_1_WHITE_UP,
                                            jac_affine_inv_0_1_WHITE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_WHITE_UP,
                                            jac_affine_inv_0_1_WHITE_UP,
                                            jac_affine_inv_0_1_WHITE_UP,
                                            jac_affine_inv_0_1_WHITE_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_0_2_WHITE_UP,
                                                  jac_affine_inv_0_2_WHITE_UP,
                                                  jac_affine_inv_0_2_WHITE_UP,
                                                  jac_affine_inv_0_2_WHITE_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_WHITE_UP,
                                        jac_affine_inv_0_2_WHITE_UP,
                                        jac_affine_inv_0_2_WHITE_UP,
                                        jac_affine_inv_0_2_WHITE_UP)))),
                        _mm256_mul_pd(
                            tmp_kernel_op_8,
                            _mm256_set_pd(tmp_kernel_op_11, tmp_kernel_op_11,
                                          tmp_kernel_op_11, tmp_kernel_op_11))),
                    _mm256_mul_pd(
                        tmp_kernel_op_10,
                        _mm256_set_pd(tmp_kernel_op_12, tmp_kernel_op_12,
                                      tmp_kernel_op_12, tmp_kernel_op_12))),
                _mm256_mul_pd(tmp_kernel_op_4,
                              _mm256_set_pd(tmp_kernel_op_5, tmp_kernel_op_5,
                                            tmp_kernel_op_5, tmp_kernel_op_5)));
            const __m256d elMatVec_2 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_kernel_op_8,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_WHITE_UP,
                                            jac_affine_inv_1_0_WHITE_UP,
                                            jac_affine_inv_1_0_WHITE_UP,
                                            jac_affine_inv_1_0_WHITE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_WHITE_UP,
                                            jac_affine_inv_1_0_WHITE_UP,
                                            jac_affine_inv_1_0_WHITE_UP,
                                            jac_affine_inv_1_0_WHITE_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_WHITE_UP,
                                            jac_affine_inv_1_1_WHITE_UP,
                                            jac_affine_inv_1_1_WHITE_UP,
                                            jac_affine_inv_1_1_WHITE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_WHITE_UP,
                                            jac_affine_inv_1_1_WHITE_UP,
                                            jac_affine_inv_1_1_WHITE_UP,
                                            jac_affine_inv_1_1_WHITE_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_1_2_WHITE_UP,
                                                  jac_affine_inv_1_2_WHITE_UP,
                                                  jac_affine_inv_1_2_WHITE_UP,
                                                  jac_affine_inv_1_2_WHITE_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_WHITE_UP,
                                        jac_affine_inv_1_2_WHITE_UP,
                                        jac_affine_inv_1_2_WHITE_UP,
                                        jac_affine_inv_1_2_WHITE_UP)))),
                        _mm256_mul_pd(
                            tmp_kernel_op_6,
                            _mm256_set_pd(tmp_kernel_op_11, tmp_kernel_op_11,
                                          tmp_kernel_op_11, tmp_kernel_op_11))),
                    _mm256_mul_pd(
                        tmp_kernel_op_10,
                        _mm256_set_pd(tmp_kernel_op_13, tmp_kernel_op_13,
                                      tmp_kernel_op_13, tmp_kernel_op_13))),
                _mm256_mul_pd(tmp_kernel_op_4,
                              _mm256_set_pd(tmp_kernel_op_7, tmp_kernel_op_7,
                                            tmp_kernel_op_7, tmp_kernel_op_7)));
            const __m256d elMatVec_3 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_kernel_op_10,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_WHITE_UP,
                                            jac_affine_inv_2_0_WHITE_UP,
                                            jac_affine_inv_2_0_WHITE_UP,
                                            jac_affine_inv_2_0_WHITE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_WHITE_UP,
                                            jac_affine_inv_2_0_WHITE_UP,
                                            jac_affine_inv_2_0_WHITE_UP,
                                            jac_affine_inv_2_0_WHITE_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_WHITE_UP,
                                            jac_affine_inv_2_1_WHITE_UP,
                                            jac_affine_inv_2_1_WHITE_UP,
                                            jac_affine_inv_2_1_WHITE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_WHITE_UP,
                                            jac_affine_inv_2_1_WHITE_UP,
                                            jac_affine_inv_2_1_WHITE_UP,
                                            jac_affine_inv_2_1_WHITE_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_2_2_WHITE_UP,
                                                  jac_affine_inv_2_2_WHITE_UP,
                                                  jac_affine_inv_2_2_WHITE_UP,
                                                  jac_affine_inv_2_2_WHITE_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_WHITE_UP,
                                        jac_affine_inv_2_2_WHITE_UP,
                                        jac_affine_inv_2_2_WHITE_UP,
                                        jac_affine_inv_2_2_WHITE_UP)))),
                        _mm256_mul_pd(
                            tmp_kernel_op_6,
                            _mm256_set_pd(tmp_kernel_op_12, tmp_kernel_op_12,
                                          tmp_kernel_op_12, tmp_kernel_op_12))),
                    _mm256_mul_pd(
                        tmp_kernel_op_8,
                        _mm256_set_pd(tmp_kernel_op_13, tmp_kernel_op_13,
                                      tmp_kernel_op_13, tmp_kernel_op_13))),
                _mm256_mul_pd(tmp_kernel_op_4,
                              _mm256_set_pd(tmp_kernel_op_9, tmp_kernel_op_9,
                                            tmp_kernel_op_9, tmp_kernel_op_9)));
            {
              {
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        elMatVec_0,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        elMatVec_1,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        elMatVec_2,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        elMatVec_3,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
              }
            }
            const __m256d tmp_moved_constant_0 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_1 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_2 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d tmp_moved_constant_3 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_8 = _mm256_mul_pd(
                tmp_moved_constant_0,
                _mm256_set_pd(tmp_moved_constant_7, tmp_moved_constant_7,
                              tmp_moved_constant_7, tmp_moved_constant_7));
            const __m256d tmp_moved_constant_10 = _mm256_mul_pd(
                tmp_moved_constant_1,
                _mm256_set_pd(tmp_moved_constant_7, tmp_moved_constant_7,
                              tmp_moved_constant_7, tmp_moved_constant_7));
            const __m256d tmp_moved_constant_12 = _mm256_mul_pd(
                tmp_moved_constant_2,
                _mm256_set_pd(tmp_moved_constant_7, tmp_moved_constant_7,
                              tmp_moved_constant_7, tmp_moved_constant_7));
            const __m256d tmp_moved_constant_14 = _mm256_mul_pd(
                tmp_moved_constant_3,
                _mm256_set_pd(tmp_moved_constant_7, tmp_moved_constant_7,
                              tmp_moved_constant_7, tmp_moved_constant_7));
            const __m256d tmp_moved_constant_18 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_8,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4),
                                        _mm256_set_pd(tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4,
                                                      tmp_moved_constant_4)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_5,
                                                      tmp_moved_constant_5,
                                                      tmp_moved_constant_5,
                                                      tmp_moved_constant_5),
                                        _mm256_set_pd(tmp_moved_constant_5,
                                                      tmp_moved_constant_5,
                                                      tmp_moved_constant_5,
                                                      tmp_moved_constant_5))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(tmp_moved_constant_6,
                                                  tmp_moved_constant_6,
                                                  tmp_moved_constant_6,
                                                  tmp_moved_constant_6),
                                    _mm256_set_pd(tmp_moved_constant_6,
                                                  tmp_moved_constant_6,
                                                  tmp_moved_constant_6,
                                                  tmp_moved_constant_6)))),
                        _mm256_mul_pd(tmp_moved_constant_12,
                                      _mm256_set_pd(tmp_moved_constant_11,
                                                    tmp_moved_constant_11,
                                                    tmp_moved_constant_11,
                                                    tmp_moved_constant_11))),
                    _mm256_mul_pd(tmp_moved_constant_14,
                                  _mm256_set_pd(tmp_moved_constant_13,
                                                tmp_moved_constant_13,
                                                tmp_moved_constant_13,
                                                tmp_moved_constant_13))),
                _mm256_mul_pd(
                    tmp_moved_constant_10,
                    _mm256_set_pd(tmp_moved_constant_9, tmp_moved_constant_9,
                                  tmp_moved_constant_9, tmp_moved_constant_9)));
            const __m256d tmp_moved_constant_19 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_10,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_WHITE_DOWN,
                                            jac_affine_inv_0_0_WHITE_DOWN,
                                            jac_affine_inv_0_0_WHITE_DOWN,
                                            jac_affine_inv_0_0_WHITE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_WHITE_DOWN,
                                            jac_affine_inv_0_0_WHITE_DOWN,
                                            jac_affine_inv_0_0_WHITE_DOWN,
                                            jac_affine_inv_0_0_WHITE_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_WHITE_DOWN,
                                            jac_affine_inv_0_1_WHITE_DOWN,
                                            jac_affine_inv_0_1_WHITE_DOWN,
                                            jac_affine_inv_0_1_WHITE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_WHITE_DOWN,
                                            jac_affine_inv_0_1_WHITE_DOWN,
                                            jac_affine_inv_0_1_WHITE_DOWN,
                                            jac_affine_inv_0_1_WHITE_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_WHITE_DOWN,
                                        jac_affine_inv_0_2_WHITE_DOWN,
                                        jac_affine_inv_0_2_WHITE_DOWN,
                                        jac_affine_inv_0_2_WHITE_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_WHITE_DOWN,
                                        jac_affine_inv_0_2_WHITE_DOWN,
                                        jac_affine_inv_0_2_WHITE_DOWN,
                                        jac_affine_inv_0_2_WHITE_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_12,
                                      _mm256_set_pd(tmp_moved_constant_15,
                                                    tmp_moved_constant_15,
                                                    tmp_moved_constant_15,
                                                    tmp_moved_constant_15))),
                    _mm256_mul_pd(tmp_moved_constant_14,
                                  _mm256_set_pd(tmp_moved_constant_16,
                                                tmp_moved_constant_16,
                                                tmp_moved_constant_16,
                                                tmp_moved_constant_16))),
                _mm256_mul_pd(
                    tmp_moved_constant_8,
                    _mm256_set_pd(tmp_moved_constant_9, tmp_moved_constant_9,
                                  tmp_moved_constant_9, tmp_moved_constant_9)));
            const __m256d tmp_moved_constant_20 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_12,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_WHITE_DOWN,
                                            jac_affine_inv_1_0_WHITE_DOWN,
                                            jac_affine_inv_1_0_WHITE_DOWN,
                                            jac_affine_inv_1_0_WHITE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_WHITE_DOWN,
                                            jac_affine_inv_1_0_WHITE_DOWN,
                                            jac_affine_inv_1_0_WHITE_DOWN,
                                            jac_affine_inv_1_0_WHITE_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_WHITE_DOWN,
                                            jac_affine_inv_1_1_WHITE_DOWN,
                                            jac_affine_inv_1_1_WHITE_DOWN,
                                            jac_affine_inv_1_1_WHITE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_WHITE_DOWN,
                                            jac_affine_inv_1_1_WHITE_DOWN,
                                            jac_affine_inv_1_1_WHITE_DOWN,
                                            jac_affine_inv_1_1_WHITE_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_WHITE_DOWN,
                                        jac_affine_inv_1_2_WHITE_DOWN,
                                        jac_affine_inv_1_2_WHITE_DOWN,
                                        jac_affine_inv_1_2_WHITE_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_WHITE_DOWN,
                                        jac_affine_inv_1_2_WHITE_DOWN,
                                        jac_affine_inv_1_2_WHITE_DOWN,
                                        jac_affine_inv_1_2_WHITE_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_8,
                                      _mm256_set_pd(tmp_moved_constant_11,
                                                    tmp_moved_constant_11,
                                                    tmp_moved_constant_11,
                                                    tmp_moved_constant_11))),
                    _mm256_mul_pd(tmp_moved_constant_10,
                                  _mm256_set_pd(tmp_moved_constant_15,
                                                tmp_moved_constant_15,
                                                tmp_moved_constant_15,
                                                tmp_moved_constant_15))),
                _mm256_mul_pd(tmp_moved_constant_14,
                              _mm256_set_pd(tmp_moved_constant_17,
                                            tmp_moved_constant_17,
                                            tmp_moved_constant_17,
                                            tmp_moved_constant_17)));
            const __m256d tmp_moved_constant_21 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_14,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_WHITE_DOWN,
                                            jac_affine_inv_2_0_WHITE_DOWN,
                                            jac_affine_inv_2_0_WHITE_DOWN,
                                            jac_affine_inv_2_0_WHITE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_WHITE_DOWN,
                                            jac_affine_inv_2_0_WHITE_DOWN,
                                            jac_affine_inv_2_0_WHITE_DOWN,
                                            jac_affine_inv_2_0_WHITE_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_WHITE_DOWN,
                                            jac_affine_inv_2_1_WHITE_DOWN,
                                            jac_affine_inv_2_1_WHITE_DOWN,
                                            jac_affine_inv_2_1_WHITE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_WHITE_DOWN,
                                            jac_affine_inv_2_1_WHITE_DOWN,
                                            jac_affine_inv_2_1_WHITE_DOWN,
                                            jac_affine_inv_2_1_WHITE_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_WHITE_DOWN,
                                        jac_affine_inv_2_2_WHITE_DOWN,
                                        jac_affine_inv_2_2_WHITE_DOWN,
                                        jac_affine_inv_2_2_WHITE_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_WHITE_DOWN,
                                        jac_affine_inv_2_2_WHITE_DOWN,
                                        jac_affine_inv_2_2_WHITE_DOWN,
                                        jac_affine_inv_2_2_WHITE_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_8,
                                      _mm256_set_pd(tmp_moved_constant_13,
                                                    tmp_moved_constant_13,
                                                    tmp_moved_constant_13,
                                                    tmp_moved_constant_13))),
                    _mm256_mul_pd(tmp_moved_constant_10,
                                  _mm256_set_pd(tmp_moved_constant_16,
                                                tmp_moved_constant_16,
                                                tmp_moved_constant_16,
                                                tmp_moved_constant_16))),
                _mm256_mul_pd(tmp_moved_constant_12,
                              _mm256_set_pd(tmp_moved_constant_17,
                                            tmp_moved_constant_17,
                                            tmp_moved_constant_17,
                                            tmp_moved_constant_17)));
            {
              {
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_18,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_19,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_20,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_21,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
              }
            }
            const __m256d tmp_moved_constant_22 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_23 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d tmp_moved_constant_24 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_25 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_30 = _mm256_mul_pd(
                tmp_moved_constant_22,
                _mm256_set_pd(tmp_moved_constant_29, tmp_moved_constant_29,
                              tmp_moved_constant_29, tmp_moved_constant_29));
            const __m256d tmp_moved_constant_32 = _mm256_mul_pd(
                tmp_moved_constant_23,
                _mm256_set_pd(tmp_moved_constant_29, tmp_moved_constant_29,
                              tmp_moved_constant_29, tmp_moved_constant_29));
            const __m256d tmp_moved_constant_34 = _mm256_mul_pd(
                tmp_moved_constant_24,
                _mm256_set_pd(tmp_moved_constant_29, tmp_moved_constant_29,
                              tmp_moved_constant_29, tmp_moved_constant_29));
            const __m256d tmp_moved_constant_36 = _mm256_mul_pd(
                tmp_moved_constant_25,
                _mm256_set_pd(tmp_moved_constant_29, tmp_moved_constant_29,
                              tmp_moved_constant_29, tmp_moved_constant_29));
            const __m256d tmp_moved_constant_40 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_30,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_26,
                                                      tmp_moved_constant_26,
                                                      tmp_moved_constant_26,
                                                      tmp_moved_constant_26),
                                        _mm256_set_pd(tmp_moved_constant_26,
                                                      tmp_moved_constant_26,
                                                      tmp_moved_constant_26,
                                                      tmp_moved_constant_26)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_27,
                                                      tmp_moved_constant_27,
                                                      tmp_moved_constant_27,
                                                      tmp_moved_constant_27),
                                        _mm256_set_pd(tmp_moved_constant_27,
                                                      tmp_moved_constant_27,
                                                      tmp_moved_constant_27,
                                                      tmp_moved_constant_27))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(tmp_moved_constant_28,
                                                  tmp_moved_constant_28,
                                                  tmp_moved_constant_28,
                                                  tmp_moved_constant_28),
                                    _mm256_set_pd(tmp_moved_constant_28,
                                                  tmp_moved_constant_28,
                                                  tmp_moved_constant_28,
                                                  tmp_moved_constant_28)))),
                        _mm256_mul_pd(tmp_moved_constant_32,
                                      _mm256_set_pd(tmp_moved_constant_31,
                                                    tmp_moved_constant_31,
                                                    tmp_moved_constant_31,
                                                    tmp_moved_constant_31))),
                    _mm256_mul_pd(tmp_moved_constant_34,
                                  _mm256_set_pd(tmp_moved_constant_33,
                                                tmp_moved_constant_33,
                                                tmp_moved_constant_33,
                                                tmp_moved_constant_33))),
                _mm256_mul_pd(tmp_moved_constant_36,
                              _mm256_set_pd(tmp_moved_constant_35,
                                            tmp_moved_constant_35,
                                            tmp_moved_constant_35,
                                            tmp_moved_constant_35)));
            const __m256d tmp_moved_constant_41 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_32,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_BLUE_UP,
                                            jac_affine_inv_0_0_BLUE_UP,
                                            jac_affine_inv_0_0_BLUE_UP,
                                            jac_affine_inv_0_0_BLUE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_BLUE_UP,
                                            jac_affine_inv_0_0_BLUE_UP,
                                            jac_affine_inv_0_0_BLUE_UP,
                                            jac_affine_inv_0_0_BLUE_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_BLUE_UP,
                                            jac_affine_inv_0_1_BLUE_UP,
                                            jac_affine_inv_0_1_BLUE_UP,
                                            jac_affine_inv_0_1_BLUE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_BLUE_UP,
                                            jac_affine_inv_0_1_BLUE_UP,
                                            jac_affine_inv_0_1_BLUE_UP,
                                            jac_affine_inv_0_1_BLUE_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_0_2_BLUE_UP,
                                                  jac_affine_inv_0_2_BLUE_UP,
                                                  jac_affine_inv_0_2_BLUE_UP,
                                                  jac_affine_inv_0_2_BLUE_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_BLUE_UP,
                                        jac_affine_inv_0_2_BLUE_UP,
                                        jac_affine_inv_0_2_BLUE_UP,
                                        jac_affine_inv_0_2_BLUE_UP)))),
                        _mm256_mul_pd(tmp_moved_constant_30,
                                      _mm256_set_pd(tmp_moved_constant_31,
                                                    tmp_moved_constant_31,
                                                    tmp_moved_constant_31,
                                                    tmp_moved_constant_31))),
                    _mm256_mul_pd(tmp_moved_constant_34,
                                  _mm256_set_pd(tmp_moved_constant_37,
                                                tmp_moved_constant_37,
                                                tmp_moved_constant_37,
                                                tmp_moved_constant_37))),
                _mm256_mul_pd(tmp_moved_constant_36,
                              _mm256_set_pd(tmp_moved_constant_38,
                                            tmp_moved_constant_38,
                                            tmp_moved_constant_38,
                                            tmp_moved_constant_38)));
            const __m256d tmp_moved_constant_42 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_34,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_BLUE_UP,
                                            jac_affine_inv_1_0_BLUE_UP,
                                            jac_affine_inv_1_0_BLUE_UP,
                                            jac_affine_inv_1_0_BLUE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_BLUE_UP,
                                            jac_affine_inv_1_0_BLUE_UP,
                                            jac_affine_inv_1_0_BLUE_UP,
                                            jac_affine_inv_1_0_BLUE_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_BLUE_UP,
                                            jac_affine_inv_1_1_BLUE_UP,
                                            jac_affine_inv_1_1_BLUE_UP,
                                            jac_affine_inv_1_1_BLUE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_BLUE_UP,
                                            jac_affine_inv_1_1_BLUE_UP,
                                            jac_affine_inv_1_1_BLUE_UP,
                                            jac_affine_inv_1_1_BLUE_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_1_2_BLUE_UP,
                                                  jac_affine_inv_1_2_BLUE_UP,
                                                  jac_affine_inv_1_2_BLUE_UP,
                                                  jac_affine_inv_1_2_BLUE_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_BLUE_UP,
                                        jac_affine_inv_1_2_BLUE_UP,
                                        jac_affine_inv_1_2_BLUE_UP,
                                        jac_affine_inv_1_2_BLUE_UP)))),
                        _mm256_mul_pd(tmp_moved_constant_30,
                                      _mm256_set_pd(tmp_moved_constant_33,
                                                    tmp_moved_constant_33,
                                                    tmp_moved_constant_33,
                                                    tmp_moved_constant_33))),
                    _mm256_mul_pd(tmp_moved_constant_32,
                                  _mm256_set_pd(tmp_moved_constant_37,
                                                tmp_moved_constant_37,
                                                tmp_moved_constant_37,
                                                tmp_moved_constant_37))),
                _mm256_mul_pd(tmp_moved_constant_36,
                              _mm256_set_pd(tmp_moved_constant_39,
                                            tmp_moved_constant_39,
                                            tmp_moved_constant_39,
                                            tmp_moved_constant_39)));
            const __m256d tmp_moved_constant_43 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_36,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_BLUE_UP,
                                            jac_affine_inv_2_0_BLUE_UP,
                                            jac_affine_inv_2_0_BLUE_UP,
                                            jac_affine_inv_2_0_BLUE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_BLUE_UP,
                                            jac_affine_inv_2_0_BLUE_UP,
                                            jac_affine_inv_2_0_BLUE_UP,
                                            jac_affine_inv_2_0_BLUE_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_BLUE_UP,
                                            jac_affine_inv_2_1_BLUE_UP,
                                            jac_affine_inv_2_1_BLUE_UP,
                                            jac_affine_inv_2_1_BLUE_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_BLUE_UP,
                                            jac_affine_inv_2_1_BLUE_UP,
                                            jac_affine_inv_2_1_BLUE_UP,
                                            jac_affine_inv_2_1_BLUE_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_2_2_BLUE_UP,
                                                  jac_affine_inv_2_2_BLUE_UP,
                                                  jac_affine_inv_2_2_BLUE_UP,
                                                  jac_affine_inv_2_2_BLUE_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_BLUE_UP,
                                        jac_affine_inv_2_2_BLUE_UP,
                                        jac_affine_inv_2_2_BLUE_UP,
                                        jac_affine_inv_2_2_BLUE_UP)))),
                        _mm256_mul_pd(tmp_moved_constant_30,
                                      _mm256_set_pd(tmp_moved_constant_35,
                                                    tmp_moved_constant_35,
                                                    tmp_moved_constant_35,
                                                    tmp_moved_constant_35))),
                    _mm256_mul_pd(tmp_moved_constant_32,
                                  _mm256_set_pd(tmp_moved_constant_38,
                                                tmp_moved_constant_38,
                                                tmp_moved_constant_38,
                                                tmp_moved_constant_38))),
                _mm256_mul_pd(tmp_moved_constant_34,
                              _mm256_set_pd(tmp_moved_constant_39,
                                            tmp_moved_constant_39,
                                            tmp_moved_constant_39,
                                            tmp_moved_constant_39)));
            {
              {
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_40,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_41,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_42,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_43,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
              }
            }
            const __m256d tmp_moved_constant_44 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d tmp_moved_constant_45 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))]);
            const __m256d tmp_moved_constant_46 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_47 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d tmp_moved_constant_52 = _mm256_mul_pd(
                tmp_moved_constant_44,
                _mm256_set_pd(tmp_moved_constant_51, tmp_moved_constant_51,
                              tmp_moved_constant_51, tmp_moved_constant_51));
            const __m256d tmp_moved_constant_54 = _mm256_mul_pd(
                tmp_moved_constant_45,
                _mm256_set_pd(tmp_moved_constant_51, tmp_moved_constant_51,
                              tmp_moved_constant_51, tmp_moved_constant_51));
            const __m256d tmp_moved_constant_56 = _mm256_mul_pd(
                tmp_moved_constant_46,
                _mm256_set_pd(tmp_moved_constant_51, tmp_moved_constant_51,
                              tmp_moved_constant_51, tmp_moved_constant_51));
            const __m256d tmp_moved_constant_58 = _mm256_mul_pd(
                tmp_moved_constant_47,
                _mm256_set_pd(tmp_moved_constant_51, tmp_moved_constant_51,
                              tmp_moved_constant_51, tmp_moved_constant_51));
            const __m256d tmp_moved_constant_62 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_52,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_48,
                                                      tmp_moved_constant_48,
                                                      tmp_moved_constant_48,
                                                      tmp_moved_constant_48),
                                        _mm256_set_pd(tmp_moved_constant_48,
                                                      tmp_moved_constant_48,
                                                      tmp_moved_constant_48,
                                                      tmp_moved_constant_48)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_49,
                                                      tmp_moved_constant_49,
                                                      tmp_moved_constant_49,
                                                      tmp_moved_constant_49),
                                        _mm256_set_pd(tmp_moved_constant_49,
                                                      tmp_moved_constant_49,
                                                      tmp_moved_constant_49,
                                                      tmp_moved_constant_49))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(tmp_moved_constant_50,
                                                  tmp_moved_constant_50,
                                                  tmp_moved_constant_50,
                                                  tmp_moved_constant_50),
                                    _mm256_set_pd(tmp_moved_constant_50,
                                                  tmp_moved_constant_50,
                                                  tmp_moved_constant_50,
                                                  tmp_moved_constant_50)))),
                        _mm256_mul_pd(tmp_moved_constant_54,
                                      _mm256_set_pd(tmp_moved_constant_53,
                                                    tmp_moved_constant_53,
                                                    tmp_moved_constant_53,
                                                    tmp_moved_constant_53))),
                    _mm256_mul_pd(tmp_moved_constant_56,
                                  _mm256_set_pd(tmp_moved_constant_55,
                                                tmp_moved_constant_55,
                                                tmp_moved_constant_55,
                                                tmp_moved_constant_55))),
                _mm256_mul_pd(tmp_moved_constant_58,
                              _mm256_set_pd(tmp_moved_constant_57,
                                            tmp_moved_constant_57,
                                            tmp_moved_constant_57,
                                            tmp_moved_constant_57)));
            const __m256d tmp_moved_constant_63 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_54,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_BLUE_DOWN,
                                            jac_affine_inv_0_0_BLUE_DOWN,
                                            jac_affine_inv_0_0_BLUE_DOWN,
                                            jac_affine_inv_0_0_BLUE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_BLUE_DOWN,
                                            jac_affine_inv_0_0_BLUE_DOWN,
                                            jac_affine_inv_0_0_BLUE_DOWN,
                                            jac_affine_inv_0_0_BLUE_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_BLUE_DOWN,
                                            jac_affine_inv_0_1_BLUE_DOWN,
                                            jac_affine_inv_0_1_BLUE_DOWN,
                                            jac_affine_inv_0_1_BLUE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_BLUE_DOWN,
                                            jac_affine_inv_0_1_BLUE_DOWN,
                                            jac_affine_inv_0_1_BLUE_DOWN,
                                            jac_affine_inv_0_1_BLUE_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_0_2_BLUE_DOWN,
                                                  jac_affine_inv_0_2_BLUE_DOWN,
                                                  jac_affine_inv_0_2_BLUE_DOWN,
                                                  jac_affine_inv_0_2_BLUE_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_BLUE_DOWN,
                                        jac_affine_inv_0_2_BLUE_DOWN,
                                        jac_affine_inv_0_2_BLUE_DOWN,
                                        jac_affine_inv_0_2_BLUE_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_52,
                                      _mm256_set_pd(tmp_moved_constant_53,
                                                    tmp_moved_constant_53,
                                                    tmp_moved_constant_53,
                                                    tmp_moved_constant_53))),
                    _mm256_mul_pd(tmp_moved_constant_56,
                                  _mm256_set_pd(tmp_moved_constant_59,
                                                tmp_moved_constant_59,
                                                tmp_moved_constant_59,
                                                tmp_moved_constant_59))),
                _mm256_mul_pd(tmp_moved_constant_58,
                              _mm256_set_pd(tmp_moved_constant_60,
                                            tmp_moved_constant_60,
                                            tmp_moved_constant_60,
                                            tmp_moved_constant_60)));
            const __m256d tmp_moved_constant_64 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_56,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_BLUE_DOWN,
                                            jac_affine_inv_1_0_BLUE_DOWN,
                                            jac_affine_inv_1_0_BLUE_DOWN,
                                            jac_affine_inv_1_0_BLUE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_BLUE_DOWN,
                                            jac_affine_inv_1_0_BLUE_DOWN,
                                            jac_affine_inv_1_0_BLUE_DOWN,
                                            jac_affine_inv_1_0_BLUE_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_BLUE_DOWN,
                                            jac_affine_inv_1_1_BLUE_DOWN,
                                            jac_affine_inv_1_1_BLUE_DOWN,
                                            jac_affine_inv_1_1_BLUE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_BLUE_DOWN,
                                            jac_affine_inv_1_1_BLUE_DOWN,
                                            jac_affine_inv_1_1_BLUE_DOWN,
                                            jac_affine_inv_1_1_BLUE_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_1_2_BLUE_DOWN,
                                                  jac_affine_inv_1_2_BLUE_DOWN,
                                                  jac_affine_inv_1_2_BLUE_DOWN,
                                                  jac_affine_inv_1_2_BLUE_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_BLUE_DOWN,
                                        jac_affine_inv_1_2_BLUE_DOWN,
                                        jac_affine_inv_1_2_BLUE_DOWN,
                                        jac_affine_inv_1_2_BLUE_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_52,
                                      _mm256_set_pd(tmp_moved_constant_55,
                                                    tmp_moved_constant_55,
                                                    tmp_moved_constant_55,
                                                    tmp_moved_constant_55))),
                    _mm256_mul_pd(tmp_moved_constant_54,
                                  _mm256_set_pd(tmp_moved_constant_59,
                                                tmp_moved_constant_59,
                                                tmp_moved_constant_59,
                                                tmp_moved_constant_59))),
                _mm256_mul_pd(tmp_moved_constant_58,
                              _mm256_set_pd(tmp_moved_constant_61,
                                            tmp_moved_constant_61,
                                            tmp_moved_constant_61,
                                            tmp_moved_constant_61)));
            const __m256d tmp_moved_constant_65 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_58,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_BLUE_DOWN,
                                            jac_affine_inv_2_0_BLUE_DOWN,
                                            jac_affine_inv_2_0_BLUE_DOWN,
                                            jac_affine_inv_2_0_BLUE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_BLUE_DOWN,
                                            jac_affine_inv_2_0_BLUE_DOWN,
                                            jac_affine_inv_2_0_BLUE_DOWN,
                                            jac_affine_inv_2_0_BLUE_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_BLUE_DOWN,
                                            jac_affine_inv_2_1_BLUE_DOWN,
                                            jac_affine_inv_2_1_BLUE_DOWN,
                                            jac_affine_inv_2_1_BLUE_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_BLUE_DOWN,
                                            jac_affine_inv_2_1_BLUE_DOWN,
                                            jac_affine_inv_2_1_BLUE_DOWN,
                                            jac_affine_inv_2_1_BLUE_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_2_2_BLUE_DOWN,
                                                  jac_affine_inv_2_2_BLUE_DOWN,
                                                  jac_affine_inv_2_2_BLUE_DOWN,
                                                  jac_affine_inv_2_2_BLUE_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_BLUE_DOWN,
                                        jac_affine_inv_2_2_BLUE_DOWN,
                                        jac_affine_inv_2_2_BLUE_DOWN,
                                        jac_affine_inv_2_2_BLUE_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_52,
                                      _mm256_set_pd(tmp_moved_constant_57,
                                                    tmp_moved_constant_57,
                                                    tmp_moved_constant_57,
                                                    tmp_moved_constant_57))),
                    _mm256_mul_pd(tmp_moved_constant_54,
                                  _mm256_set_pd(tmp_moved_constant_60,
                                                tmp_moved_constant_60,
                                                tmp_moved_constant_60,
                                                tmp_moved_constant_60))),
                _mm256_mul_pd(tmp_moved_constant_56,
                              _mm256_set_pd(tmp_moved_constant_61,
                                            tmp_moved_constant_61,
                                            tmp_moved_constant_61,
                                            tmp_moved_constant_61)));
            {
              {
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_62,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_63,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_64,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_65,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
              }
            }
            const __m256d tmp_moved_constant_66 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_67 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d tmp_moved_constant_68 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
                           ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                           ((ctr_1 * (ctr_1 + 1)) / (2)) -
                           (((-ctr_2 + micro_edges_per_macro_edge) *
                             (-ctr_2 + micro_edges_per_macro_edge + 1) *
                             (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                            (6)) +
                           (((micro_edges_per_macro_edge + 1) *
                             (micro_edges_per_macro_edge + 2) *
                             (micro_edges_per_macro_edge + 3)) /
                            (6))]);
            const __m256d tmp_moved_constant_69 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_74 = _mm256_mul_pd(
                tmp_moved_constant_66,
                _mm256_set_pd(tmp_moved_constant_73, tmp_moved_constant_73,
                              tmp_moved_constant_73, tmp_moved_constant_73));
            const __m256d tmp_moved_constant_76 = _mm256_mul_pd(
                tmp_moved_constant_67,
                _mm256_set_pd(tmp_moved_constant_73, tmp_moved_constant_73,
                              tmp_moved_constant_73, tmp_moved_constant_73));
            const __m256d tmp_moved_constant_78 = _mm256_mul_pd(
                tmp_moved_constant_68,
                _mm256_set_pd(tmp_moved_constant_73, tmp_moved_constant_73,
                              tmp_moved_constant_73, tmp_moved_constant_73));
            const __m256d tmp_moved_constant_80 = _mm256_mul_pd(
                tmp_moved_constant_69,
                _mm256_set_pd(tmp_moved_constant_73, tmp_moved_constant_73,
                              tmp_moved_constant_73, tmp_moved_constant_73));
            const __m256d tmp_moved_constant_84 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_74,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_70,
                                                      tmp_moved_constant_70,
                                                      tmp_moved_constant_70,
                                                      tmp_moved_constant_70),
                                        _mm256_set_pd(tmp_moved_constant_70,
                                                      tmp_moved_constant_70,
                                                      tmp_moved_constant_70,
                                                      tmp_moved_constant_70)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_71,
                                                      tmp_moved_constant_71,
                                                      tmp_moved_constant_71,
                                                      tmp_moved_constant_71),
                                        _mm256_set_pd(tmp_moved_constant_71,
                                                      tmp_moved_constant_71,
                                                      tmp_moved_constant_71,
                                                      tmp_moved_constant_71))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(tmp_moved_constant_72,
                                                  tmp_moved_constant_72,
                                                  tmp_moved_constant_72,
                                                  tmp_moved_constant_72),
                                    _mm256_set_pd(tmp_moved_constant_72,
                                                  tmp_moved_constant_72,
                                                  tmp_moved_constant_72,
                                                  tmp_moved_constant_72)))),
                        _mm256_mul_pd(tmp_moved_constant_76,
                                      _mm256_set_pd(tmp_moved_constant_75,
                                                    tmp_moved_constant_75,
                                                    tmp_moved_constant_75,
                                                    tmp_moved_constant_75))),
                    _mm256_mul_pd(tmp_moved_constant_78,
                                  _mm256_set_pd(tmp_moved_constant_77,
                                                tmp_moved_constant_77,
                                                tmp_moved_constant_77,
                                                tmp_moved_constant_77))),
                _mm256_mul_pd(tmp_moved_constant_80,
                              _mm256_set_pd(tmp_moved_constant_79,
                                            tmp_moved_constant_79,
                                            tmp_moved_constant_79,
                                            tmp_moved_constant_79)));
            const __m256d tmp_moved_constant_85 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_76,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_GREEN_UP,
                                            jac_affine_inv_0_0_GREEN_UP,
                                            jac_affine_inv_0_0_GREEN_UP,
                                            jac_affine_inv_0_0_GREEN_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_GREEN_UP,
                                            jac_affine_inv_0_0_GREEN_UP,
                                            jac_affine_inv_0_0_GREEN_UP,
                                            jac_affine_inv_0_0_GREEN_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_GREEN_UP,
                                            jac_affine_inv_0_1_GREEN_UP,
                                            jac_affine_inv_0_1_GREEN_UP,
                                            jac_affine_inv_0_1_GREEN_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_GREEN_UP,
                                            jac_affine_inv_0_1_GREEN_UP,
                                            jac_affine_inv_0_1_GREEN_UP,
                                            jac_affine_inv_0_1_GREEN_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_0_2_GREEN_UP,
                                                  jac_affine_inv_0_2_GREEN_UP,
                                                  jac_affine_inv_0_2_GREEN_UP,
                                                  jac_affine_inv_0_2_GREEN_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_GREEN_UP,
                                        jac_affine_inv_0_2_GREEN_UP,
                                        jac_affine_inv_0_2_GREEN_UP,
                                        jac_affine_inv_0_2_GREEN_UP)))),
                        _mm256_mul_pd(tmp_moved_constant_74,
                                      _mm256_set_pd(tmp_moved_constant_75,
                                                    tmp_moved_constant_75,
                                                    tmp_moved_constant_75,
                                                    tmp_moved_constant_75))),
                    _mm256_mul_pd(tmp_moved_constant_78,
                                  _mm256_set_pd(tmp_moved_constant_81,
                                                tmp_moved_constant_81,
                                                tmp_moved_constant_81,
                                                tmp_moved_constant_81))),
                _mm256_mul_pd(tmp_moved_constant_80,
                              _mm256_set_pd(tmp_moved_constant_82,
                                            tmp_moved_constant_82,
                                            tmp_moved_constant_82,
                                            tmp_moved_constant_82)));
            const __m256d tmp_moved_constant_86 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_78,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_GREEN_UP,
                                            jac_affine_inv_1_0_GREEN_UP,
                                            jac_affine_inv_1_0_GREEN_UP,
                                            jac_affine_inv_1_0_GREEN_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_GREEN_UP,
                                            jac_affine_inv_1_0_GREEN_UP,
                                            jac_affine_inv_1_0_GREEN_UP,
                                            jac_affine_inv_1_0_GREEN_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_GREEN_UP,
                                            jac_affine_inv_1_1_GREEN_UP,
                                            jac_affine_inv_1_1_GREEN_UP,
                                            jac_affine_inv_1_1_GREEN_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_GREEN_UP,
                                            jac_affine_inv_1_1_GREEN_UP,
                                            jac_affine_inv_1_1_GREEN_UP,
                                            jac_affine_inv_1_1_GREEN_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_1_2_GREEN_UP,
                                                  jac_affine_inv_1_2_GREEN_UP,
                                                  jac_affine_inv_1_2_GREEN_UP,
                                                  jac_affine_inv_1_2_GREEN_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_GREEN_UP,
                                        jac_affine_inv_1_2_GREEN_UP,
                                        jac_affine_inv_1_2_GREEN_UP,
                                        jac_affine_inv_1_2_GREEN_UP)))),
                        _mm256_mul_pd(tmp_moved_constant_74,
                                      _mm256_set_pd(tmp_moved_constant_77,
                                                    tmp_moved_constant_77,
                                                    tmp_moved_constant_77,
                                                    tmp_moved_constant_77))),
                    _mm256_mul_pd(tmp_moved_constant_76,
                                  _mm256_set_pd(tmp_moved_constant_81,
                                                tmp_moved_constant_81,
                                                tmp_moved_constant_81,
                                                tmp_moved_constant_81))),
                _mm256_mul_pd(tmp_moved_constant_80,
                              _mm256_set_pd(tmp_moved_constant_83,
                                            tmp_moved_constant_83,
                                            tmp_moved_constant_83,
                                            tmp_moved_constant_83)));
            const __m256d tmp_moved_constant_87 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_80,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_GREEN_UP,
                                            jac_affine_inv_2_0_GREEN_UP,
                                            jac_affine_inv_2_0_GREEN_UP,
                                            jac_affine_inv_2_0_GREEN_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_GREEN_UP,
                                            jac_affine_inv_2_0_GREEN_UP,
                                            jac_affine_inv_2_0_GREEN_UP,
                                            jac_affine_inv_2_0_GREEN_UP)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_GREEN_UP,
                                            jac_affine_inv_2_1_GREEN_UP,
                                            jac_affine_inv_2_1_GREEN_UP,
                                            jac_affine_inv_2_1_GREEN_UP),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_GREEN_UP,
                                            jac_affine_inv_2_1_GREEN_UP,
                                            jac_affine_inv_2_1_GREEN_UP,
                                            jac_affine_inv_2_1_GREEN_UP))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(jac_affine_inv_2_2_GREEN_UP,
                                                  jac_affine_inv_2_2_GREEN_UP,
                                                  jac_affine_inv_2_2_GREEN_UP,
                                                  jac_affine_inv_2_2_GREEN_UP),
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_GREEN_UP,
                                        jac_affine_inv_2_2_GREEN_UP,
                                        jac_affine_inv_2_2_GREEN_UP,
                                        jac_affine_inv_2_2_GREEN_UP)))),
                        _mm256_mul_pd(tmp_moved_constant_74,
                                      _mm256_set_pd(tmp_moved_constant_79,
                                                    tmp_moved_constant_79,
                                                    tmp_moved_constant_79,
                                                    tmp_moved_constant_79))),
                    _mm256_mul_pd(tmp_moved_constant_76,
                                  _mm256_set_pd(tmp_moved_constant_82,
                                                tmp_moved_constant_82,
                                                tmp_moved_constant_82,
                                                tmp_moved_constant_82))),
                _mm256_mul_pd(tmp_moved_constant_78,
                              _mm256_set_pd(tmp_moved_constant_83,
                                            tmp_moved_constant_83,
                                            tmp_moved_constant_83,
                                            tmp_moved_constant_83)));
            {
              {
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_84,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_85,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_86,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_87,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
              }
            }
            const __m256d tmp_moved_constant_88 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d tmp_moved_constant_89 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_90 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                           1]);
            const __m256d tmp_moved_constant_91 = _mm256_loadu_pd(
                &_data_src[ctr_0 +
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
                            (6))]);
            const __m256d tmp_moved_constant_96 = _mm256_mul_pd(
                tmp_moved_constant_88,
                _mm256_set_pd(tmp_moved_constant_95, tmp_moved_constant_95,
                              tmp_moved_constant_95, tmp_moved_constant_95));
            const __m256d tmp_moved_constant_98 = _mm256_mul_pd(
                tmp_moved_constant_89,
                _mm256_set_pd(tmp_moved_constant_95, tmp_moved_constant_95,
                              tmp_moved_constant_95, tmp_moved_constant_95));
            const __m256d tmp_moved_constant_100 = _mm256_mul_pd(
                tmp_moved_constant_90,
                _mm256_set_pd(tmp_moved_constant_95, tmp_moved_constant_95,
                              tmp_moved_constant_95, tmp_moved_constant_95));
            const __m256d tmp_moved_constant_102 = _mm256_mul_pd(
                tmp_moved_constant_91,
                _mm256_set_pd(tmp_moved_constant_95, tmp_moved_constant_95,
                              tmp_moved_constant_95, tmp_moved_constant_95));
            const __m256d tmp_moved_constant_106 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_96,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_92,
                                                      tmp_moved_constant_92,
                                                      tmp_moved_constant_92,
                                                      tmp_moved_constant_92),
                                        _mm256_set_pd(tmp_moved_constant_92,
                                                      tmp_moved_constant_92,
                                                      tmp_moved_constant_92,
                                                      tmp_moved_constant_92)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(tmp_moved_constant_93,
                                                      tmp_moved_constant_93,
                                                      tmp_moved_constant_93,
                                                      tmp_moved_constant_93),
                                        _mm256_set_pd(tmp_moved_constant_93,
                                                      tmp_moved_constant_93,
                                                      tmp_moved_constant_93,
                                                      tmp_moved_constant_93))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(tmp_moved_constant_94,
                                                  tmp_moved_constant_94,
                                                  tmp_moved_constant_94,
                                                  tmp_moved_constant_94),
                                    _mm256_set_pd(tmp_moved_constant_94,
                                                  tmp_moved_constant_94,
                                                  tmp_moved_constant_94,
                                                  tmp_moved_constant_94)))),
                        _mm256_mul_pd(tmp_moved_constant_102,
                                      _mm256_set_pd(tmp_moved_constant_101,
                                                    tmp_moved_constant_101,
                                                    tmp_moved_constant_101,
                                                    tmp_moved_constant_101))),
                    _mm256_mul_pd(tmp_moved_constant_98,
                                  _mm256_set_pd(tmp_moved_constant_97,
                                                tmp_moved_constant_97,
                                                tmp_moved_constant_97,
                                                tmp_moved_constant_97))),
                _mm256_mul_pd(tmp_moved_constant_100,
                              _mm256_set_pd(tmp_moved_constant_99,
                                            tmp_moved_constant_99,
                                            tmp_moved_constant_99,
                                            tmp_moved_constant_99)));
            const __m256d tmp_moved_constant_107 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_98,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_GREEN_DOWN,
                                            jac_affine_inv_0_0_GREEN_DOWN,
                                            jac_affine_inv_0_0_GREEN_DOWN,
                                            jac_affine_inv_0_0_GREEN_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_0_GREEN_DOWN,
                                            jac_affine_inv_0_0_GREEN_DOWN,
                                            jac_affine_inv_0_0_GREEN_DOWN,
                                            jac_affine_inv_0_0_GREEN_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_GREEN_DOWN,
                                            jac_affine_inv_0_1_GREEN_DOWN,
                                            jac_affine_inv_0_1_GREEN_DOWN,
                                            jac_affine_inv_0_1_GREEN_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_0_1_GREEN_DOWN,
                                            jac_affine_inv_0_1_GREEN_DOWN,
                                            jac_affine_inv_0_1_GREEN_DOWN,
                                            jac_affine_inv_0_1_GREEN_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_GREEN_DOWN,
                                        jac_affine_inv_0_2_GREEN_DOWN,
                                        jac_affine_inv_0_2_GREEN_DOWN,
                                        jac_affine_inv_0_2_GREEN_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_0_2_GREEN_DOWN,
                                        jac_affine_inv_0_2_GREEN_DOWN,
                                        jac_affine_inv_0_2_GREEN_DOWN,
                                        jac_affine_inv_0_2_GREEN_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_100,
                                      _mm256_set_pd(tmp_moved_constant_103,
                                                    tmp_moved_constant_103,
                                                    tmp_moved_constant_103,
                                                    tmp_moved_constant_103))),
                    _mm256_mul_pd(tmp_moved_constant_102,
                                  _mm256_set_pd(tmp_moved_constant_104,
                                                tmp_moved_constant_104,
                                                tmp_moved_constant_104,
                                                tmp_moved_constant_104))),
                _mm256_mul_pd(tmp_moved_constant_96,
                              _mm256_set_pd(tmp_moved_constant_97,
                                            tmp_moved_constant_97,
                                            tmp_moved_constant_97,
                                            tmp_moved_constant_97)));
            const __m256d tmp_moved_constant_108 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_100,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_GREEN_DOWN,
                                            jac_affine_inv_1_0_GREEN_DOWN,
                                            jac_affine_inv_1_0_GREEN_DOWN,
                                            jac_affine_inv_1_0_GREEN_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_0_GREEN_DOWN,
                                            jac_affine_inv_1_0_GREEN_DOWN,
                                            jac_affine_inv_1_0_GREEN_DOWN,
                                            jac_affine_inv_1_0_GREEN_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_GREEN_DOWN,
                                            jac_affine_inv_1_1_GREEN_DOWN,
                                            jac_affine_inv_1_1_GREEN_DOWN,
                                            jac_affine_inv_1_1_GREEN_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_1_1_GREEN_DOWN,
                                            jac_affine_inv_1_1_GREEN_DOWN,
                                            jac_affine_inv_1_1_GREEN_DOWN,
                                            jac_affine_inv_1_1_GREEN_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_GREEN_DOWN,
                                        jac_affine_inv_1_2_GREEN_DOWN,
                                        jac_affine_inv_1_2_GREEN_DOWN,
                                        jac_affine_inv_1_2_GREEN_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_1_2_GREEN_DOWN,
                                        jac_affine_inv_1_2_GREEN_DOWN,
                                        jac_affine_inv_1_2_GREEN_DOWN,
                                        jac_affine_inv_1_2_GREEN_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_98,
                                      _mm256_set_pd(tmp_moved_constant_103,
                                                    tmp_moved_constant_103,
                                                    tmp_moved_constant_103,
                                                    tmp_moved_constant_103))),
                    _mm256_mul_pd(tmp_moved_constant_102,
                                  _mm256_set_pd(tmp_moved_constant_105,
                                                tmp_moved_constant_105,
                                                tmp_moved_constant_105,
                                                tmp_moved_constant_105))),
                _mm256_mul_pd(tmp_moved_constant_96,
                              _mm256_set_pd(tmp_moved_constant_99,
                                            tmp_moved_constant_99,
                                            tmp_moved_constant_99,
                                            tmp_moved_constant_99)));
            const __m256d tmp_moved_constant_109 = _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(
                            tmp_moved_constant_102,
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_GREEN_DOWN,
                                            jac_affine_inv_2_0_GREEN_DOWN,
                                            jac_affine_inv_2_0_GREEN_DOWN,
                                            jac_affine_inv_2_0_GREEN_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_0_GREEN_DOWN,
                                            jac_affine_inv_2_0_GREEN_DOWN,
                                            jac_affine_inv_2_0_GREEN_DOWN,
                                            jac_affine_inv_2_0_GREEN_DOWN)),
                                    _mm256_mul_pd(
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_GREEN_DOWN,
                                            jac_affine_inv_2_1_GREEN_DOWN,
                                            jac_affine_inv_2_1_GREEN_DOWN,
                                            jac_affine_inv_2_1_GREEN_DOWN),
                                        _mm256_set_pd(
                                            jac_affine_inv_2_1_GREEN_DOWN,
                                            jac_affine_inv_2_1_GREEN_DOWN,
                                            jac_affine_inv_2_1_GREEN_DOWN,
                                            jac_affine_inv_2_1_GREEN_DOWN))),
                                _mm256_mul_pd(
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_GREEN_DOWN,
                                        jac_affine_inv_2_2_GREEN_DOWN,
                                        jac_affine_inv_2_2_GREEN_DOWN,
                                        jac_affine_inv_2_2_GREEN_DOWN),
                                    _mm256_set_pd(
                                        jac_affine_inv_2_2_GREEN_DOWN,
                                        jac_affine_inv_2_2_GREEN_DOWN,
                                        jac_affine_inv_2_2_GREEN_DOWN,
                                        jac_affine_inv_2_2_GREEN_DOWN)))),
                        _mm256_mul_pd(tmp_moved_constant_96,
                                      _mm256_set_pd(tmp_moved_constant_101,
                                                    tmp_moved_constant_101,
                                                    tmp_moved_constant_101,
                                                    tmp_moved_constant_101))),
                    _mm256_mul_pd(tmp_moved_constant_98,
                                  _mm256_set_pd(tmp_moved_constant_104,
                                                tmp_moved_constant_104,
                                                tmp_moved_constant_104,
                                                tmp_moved_constant_104))),
                _mm256_mul_pd(tmp_moved_constant_100,
                              _mm256_set_pd(tmp_moved_constant_105,
                                            tmp_moved_constant_105,
                                            tmp_moved_constant_105,
                                            tmp_moved_constant_105)));
            {
              {
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_106,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_107,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                               1],
                    _mm256_add_pd(
                        tmp_moved_constant_108,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                 1])));
                _mm256_storeu_pd(
                    &_data_dst[ctr_0 +
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
                                (6))],
                    _mm256_add_pd(
                        tmp_moved_constant_109,
                        _mm256_loadu_pd(
                            &_data_dst
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
                                  (6))])));
              }
            }
          }
          for (int64_t ctr_0 =
                   (int64_t)((-ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2) /
                             (4)) *
                   (4);
               ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2;
               ctr_0 += 1) {
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
            const walberla::float64 tmp_kernel_op_4 =
                src_dof_0 * tmp_kernel_op_3;
            const walberla::float64 tmp_kernel_op_6 =
                src_dof_1 * tmp_kernel_op_3;
            const walberla::float64 tmp_kernel_op_8 =
                src_dof_2 * tmp_kernel_op_3;
            const walberla::float64 tmp_kernel_op_10 =
                src_dof_3 * tmp_kernel_op_3;
            const walberla::float64 elMatVec_0 =
                tmp_kernel_op_10 * tmp_kernel_op_9 +
                tmp_kernel_op_4 * ((tmp_kernel_op_0 * tmp_kernel_op_0) +
                                   (tmp_kernel_op_1 * tmp_kernel_op_1) +
                                   (tmp_kernel_op_2 * tmp_kernel_op_2)) +
                tmp_kernel_op_5 * tmp_kernel_op_6 +
                tmp_kernel_op_7 * tmp_kernel_op_8;
            const walberla::float64 elMatVec_1 =
                tmp_kernel_op_10 * tmp_kernel_op_12 +
                tmp_kernel_op_11 * tmp_kernel_op_8 +
                tmp_kernel_op_4 * tmp_kernel_op_5 +
                tmp_kernel_op_6 * ((jac_affine_inv_0_0_WHITE_UP *
                                    jac_affine_inv_0_0_WHITE_UP) +
                                   (jac_affine_inv_0_1_WHITE_UP *
                                    jac_affine_inv_0_1_WHITE_UP) +
                                   (jac_affine_inv_0_2_WHITE_UP *
                                    jac_affine_inv_0_2_WHITE_UP));
            const walberla::float64 elMatVec_2 =
                tmp_kernel_op_10 * tmp_kernel_op_13 +
                tmp_kernel_op_11 * tmp_kernel_op_6 +
                tmp_kernel_op_4 * tmp_kernel_op_7 +
                tmp_kernel_op_8 * ((jac_affine_inv_1_0_WHITE_UP *
                                    jac_affine_inv_1_0_WHITE_UP) +
                                   (jac_affine_inv_1_1_WHITE_UP *
                                    jac_affine_inv_1_1_WHITE_UP) +
                                   (jac_affine_inv_1_2_WHITE_UP *
                                    jac_affine_inv_1_2_WHITE_UP));
            const walberla::float64 elMatVec_3 =
                tmp_kernel_op_10 * ((jac_affine_inv_2_0_WHITE_UP *
                                     jac_affine_inv_2_0_WHITE_UP) +
                                    (jac_affine_inv_2_1_WHITE_UP *
                                     jac_affine_inv_2_1_WHITE_UP) +
                                    (jac_affine_inv_2_2_WHITE_UP *
                                     jac_affine_inv_2_2_WHITE_UP)) +
                tmp_kernel_op_12 * tmp_kernel_op_6 +
                tmp_kernel_op_13 * tmp_kernel_op_8 +
                tmp_kernel_op_4 * tmp_kernel_op_9;
            {
              {
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
            const walberla::float64 tmp_moved_constant_0 =
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
            const walberla::float64 tmp_moved_constant_1 =
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
            const walberla::float64 tmp_moved_constant_2 =
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
            const walberla::float64 tmp_moved_constant_3 =
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
            const walberla::float64 tmp_moved_constant_8 =
                tmp_moved_constant_0 * tmp_moved_constant_7;
            const walberla::float64 tmp_moved_constant_10 =
                tmp_moved_constant_1 * tmp_moved_constant_7;
            const walberla::float64 tmp_moved_constant_12 =
                tmp_moved_constant_2 * tmp_moved_constant_7;
            const walberla::float64 tmp_moved_constant_14 =
                tmp_moved_constant_3 * tmp_moved_constant_7;
            const walberla::float64 tmp_moved_constant_18 =
                tmp_moved_constant_10 * tmp_moved_constant_9 +
                tmp_moved_constant_11 * tmp_moved_constant_12 +
                tmp_moved_constant_13 * tmp_moved_constant_14 +
                tmp_moved_constant_8 *
                    ((tmp_moved_constant_4 * tmp_moved_constant_4) +
                     (tmp_moved_constant_5 * tmp_moved_constant_5) +
                     (tmp_moved_constant_6 * tmp_moved_constant_6));
            const walberla::float64 tmp_moved_constant_19 =
                tmp_moved_constant_10 * ((jac_affine_inv_0_0_WHITE_DOWN *
                                          jac_affine_inv_0_0_WHITE_DOWN) +
                                         (jac_affine_inv_0_1_WHITE_DOWN *
                                          jac_affine_inv_0_1_WHITE_DOWN) +
                                         (jac_affine_inv_0_2_WHITE_DOWN *
                                          jac_affine_inv_0_2_WHITE_DOWN)) +
                tmp_moved_constant_12 * tmp_moved_constant_15 +
                tmp_moved_constant_14 * tmp_moved_constant_16 +
                tmp_moved_constant_8 * tmp_moved_constant_9;
            const walberla::float64 tmp_moved_constant_20 =
                tmp_moved_constant_10 * tmp_moved_constant_15 +
                tmp_moved_constant_11 * tmp_moved_constant_8 +
                tmp_moved_constant_12 * ((jac_affine_inv_1_0_WHITE_DOWN *
                                          jac_affine_inv_1_0_WHITE_DOWN) +
                                         (jac_affine_inv_1_1_WHITE_DOWN *
                                          jac_affine_inv_1_1_WHITE_DOWN) +
                                         (jac_affine_inv_1_2_WHITE_DOWN *
                                          jac_affine_inv_1_2_WHITE_DOWN)) +
                tmp_moved_constant_14 * tmp_moved_constant_17;
            const walberla::float64 tmp_moved_constant_21 =
                tmp_moved_constant_10 * tmp_moved_constant_16 +
                tmp_moved_constant_12 * tmp_moved_constant_17 +
                tmp_moved_constant_13 * tmp_moved_constant_8 +
                tmp_moved_constant_14 * ((jac_affine_inv_2_0_WHITE_DOWN *
                                          jac_affine_inv_2_0_WHITE_DOWN) +
                                         (jac_affine_inv_2_1_WHITE_DOWN *
                                          jac_affine_inv_2_1_WHITE_DOWN) +
                                         (jac_affine_inv_2_2_WHITE_DOWN *
                                          jac_affine_inv_2_2_WHITE_DOWN));
            {
              {
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
                          1] =
                    tmp_moved_constant_18 +
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
                    tmp_moved_constant_19 +
                    _data_dst[ctr_0 +
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
                           (6))] =
                    tmp_moved_constant_20 +
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
                    tmp_moved_constant_21 +
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
            const walberla::float64 tmp_moved_constant_22 =
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
            const walberla::float64 tmp_moved_constant_23 =
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
            const walberla::float64 tmp_moved_constant_24 =
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
            const walberla::float64 tmp_moved_constant_25 =
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
            const walberla::float64 tmp_moved_constant_30 =
                tmp_moved_constant_22 * tmp_moved_constant_29;
            const walberla::float64 tmp_moved_constant_32 =
                tmp_moved_constant_23 * tmp_moved_constant_29;
            const walberla::float64 tmp_moved_constant_34 =
                tmp_moved_constant_24 * tmp_moved_constant_29;
            const walberla::float64 tmp_moved_constant_36 =
                tmp_moved_constant_25 * tmp_moved_constant_29;
            const walberla::float64 tmp_moved_constant_40 =
                tmp_moved_constant_30 *
                    ((tmp_moved_constant_26 * tmp_moved_constant_26) +
                     (tmp_moved_constant_27 * tmp_moved_constant_27) +
                     (tmp_moved_constant_28 * tmp_moved_constant_28)) +
                tmp_moved_constant_31 * tmp_moved_constant_32 +
                tmp_moved_constant_33 * tmp_moved_constant_34 +
                tmp_moved_constant_35 * tmp_moved_constant_36;
            const walberla::float64 tmp_moved_constant_41 =
                tmp_moved_constant_30 * tmp_moved_constant_31 +
                tmp_moved_constant_32 *
                    ((jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_0_0_BLUE_UP) +
                     (jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_0_1_BLUE_UP) +
                     (jac_affine_inv_0_2_BLUE_UP *
                      jac_affine_inv_0_2_BLUE_UP)) +
                tmp_moved_constant_34 * tmp_moved_constant_37 +
                tmp_moved_constant_36 * tmp_moved_constant_38;
            const walberla::float64 tmp_moved_constant_42 =
                tmp_moved_constant_30 * tmp_moved_constant_33 +
                tmp_moved_constant_32 * tmp_moved_constant_37 +
                tmp_moved_constant_34 *
                    ((jac_affine_inv_1_0_BLUE_UP * jac_affine_inv_1_0_BLUE_UP) +
                     (jac_affine_inv_1_1_BLUE_UP * jac_affine_inv_1_1_BLUE_UP) +
                     (jac_affine_inv_1_2_BLUE_UP *
                      jac_affine_inv_1_2_BLUE_UP)) +
                tmp_moved_constant_36 * tmp_moved_constant_39;
            const walberla::float64 tmp_moved_constant_43 =
                tmp_moved_constant_30 * tmp_moved_constant_35 +
                tmp_moved_constant_32 * tmp_moved_constant_38 +
                tmp_moved_constant_34 * tmp_moved_constant_39 +
                tmp_moved_constant_36 *
                    ((jac_affine_inv_2_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP) +
                     (jac_affine_inv_2_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP) +
                     (jac_affine_inv_2_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP));
            {
              {
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
                    tmp_moved_constant_40 +
                    _data_dst[ctr_0 +
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
                           (6))] =
                    tmp_moved_constant_41 +
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
                    tmp_moved_constant_42 +
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
                    tmp_moved_constant_43 +
                    _data_dst[ctr_0 +
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
            const walberla::float64 tmp_moved_constant_44 =
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
            const walberla::float64 tmp_moved_constant_45 =
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
            const walberla::float64 tmp_moved_constant_46 =
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
            const walberla::float64 tmp_moved_constant_47 =
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
            const walberla::float64 tmp_moved_constant_52 =
                tmp_moved_constant_44 * tmp_moved_constant_51;
            const walberla::float64 tmp_moved_constant_54 =
                tmp_moved_constant_45 * tmp_moved_constant_51;
            const walberla::float64 tmp_moved_constant_56 =
                tmp_moved_constant_46 * tmp_moved_constant_51;
            const walberla::float64 tmp_moved_constant_58 =
                tmp_moved_constant_47 * tmp_moved_constant_51;
            const walberla::float64 tmp_moved_constant_62 =
                tmp_moved_constant_52 *
                    ((tmp_moved_constant_48 * tmp_moved_constant_48) +
                     (tmp_moved_constant_49 * tmp_moved_constant_49) +
                     (tmp_moved_constant_50 * tmp_moved_constant_50)) +
                tmp_moved_constant_53 * tmp_moved_constant_54 +
                tmp_moved_constant_55 * tmp_moved_constant_56 +
                tmp_moved_constant_57 * tmp_moved_constant_58;
            const walberla::float64 tmp_moved_constant_63 =
                tmp_moved_constant_52 * tmp_moved_constant_53 +
                tmp_moved_constant_54 * ((jac_affine_inv_0_0_BLUE_DOWN *
                                          jac_affine_inv_0_0_BLUE_DOWN) +
                                         (jac_affine_inv_0_1_BLUE_DOWN *
                                          jac_affine_inv_0_1_BLUE_DOWN) +
                                         (jac_affine_inv_0_2_BLUE_DOWN *
                                          jac_affine_inv_0_2_BLUE_DOWN)) +
                tmp_moved_constant_56 * tmp_moved_constant_59 +
                tmp_moved_constant_58 * tmp_moved_constant_60;
            const walberla::float64 tmp_moved_constant_64 =
                tmp_moved_constant_52 * tmp_moved_constant_55 +
                tmp_moved_constant_54 * tmp_moved_constant_59 +
                tmp_moved_constant_56 * ((jac_affine_inv_1_0_BLUE_DOWN *
                                          jac_affine_inv_1_0_BLUE_DOWN) +
                                         (jac_affine_inv_1_1_BLUE_DOWN *
                                          jac_affine_inv_1_1_BLUE_DOWN) +
                                         (jac_affine_inv_1_2_BLUE_DOWN *
                                          jac_affine_inv_1_2_BLUE_DOWN)) +
                tmp_moved_constant_58 * tmp_moved_constant_61;
            const walberla::float64 tmp_moved_constant_65 =
                tmp_moved_constant_52 * tmp_moved_constant_57 +
                tmp_moved_constant_54 * tmp_moved_constant_60 +
                tmp_moved_constant_56 * tmp_moved_constant_61 +
                tmp_moved_constant_58 * ((jac_affine_inv_2_0_BLUE_DOWN *
                                          jac_affine_inv_2_0_BLUE_DOWN) +
                                         (jac_affine_inv_2_1_BLUE_DOWN *
                                          jac_affine_inv_2_1_BLUE_DOWN) +
                                         (jac_affine_inv_2_2_BLUE_DOWN *
                                          jac_affine_inv_2_2_BLUE_DOWN));
            {
              {
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
                           (6))] =
                    tmp_moved_constant_62 +
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
                    tmp_moved_constant_63 +
                    _data_dst[ctr_0 +
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
                    tmp_moved_constant_64 +
                    _data_dst[ctr_0 +
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
                           (6))] =
                    tmp_moved_constant_65 +
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
            const walberla::float64 tmp_moved_constant_66 =
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
            const walberla::float64 tmp_moved_constant_67 =
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
            const walberla::float64 tmp_moved_constant_68 =
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
            const walberla::float64 tmp_moved_constant_69 =
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
            const walberla::float64 tmp_moved_constant_74 =
                tmp_moved_constant_66 * tmp_moved_constant_73;
            const walberla::float64 tmp_moved_constant_76 =
                tmp_moved_constant_67 * tmp_moved_constant_73;
            const walberla::float64 tmp_moved_constant_78 =
                tmp_moved_constant_68 * tmp_moved_constant_73;
            const walberla::float64 tmp_moved_constant_80 =
                tmp_moved_constant_69 * tmp_moved_constant_73;
            const walberla::float64 tmp_moved_constant_84 =
                tmp_moved_constant_74 *
                    ((tmp_moved_constant_70 * tmp_moved_constant_70) +
                     (tmp_moved_constant_71 * tmp_moved_constant_71) +
                     (tmp_moved_constant_72 * tmp_moved_constant_72)) +
                tmp_moved_constant_75 * tmp_moved_constant_76 +
                tmp_moved_constant_77 * tmp_moved_constant_78 +
                tmp_moved_constant_79 * tmp_moved_constant_80;
            const walberla::float64 tmp_moved_constant_85 =
                tmp_moved_constant_74 * tmp_moved_constant_75 +
                tmp_moved_constant_76 * ((jac_affine_inv_0_0_GREEN_UP *
                                          jac_affine_inv_0_0_GREEN_UP) +
                                         (jac_affine_inv_0_1_GREEN_UP *
                                          jac_affine_inv_0_1_GREEN_UP) +
                                         (jac_affine_inv_0_2_GREEN_UP *
                                          jac_affine_inv_0_2_GREEN_UP)) +
                tmp_moved_constant_78 * tmp_moved_constant_81 +
                tmp_moved_constant_80 * tmp_moved_constant_82;
            const walberla::float64 tmp_moved_constant_86 =
                tmp_moved_constant_74 * tmp_moved_constant_77 +
                tmp_moved_constant_76 * tmp_moved_constant_81 +
                tmp_moved_constant_78 * ((jac_affine_inv_1_0_GREEN_UP *
                                          jac_affine_inv_1_0_GREEN_UP) +
                                         (jac_affine_inv_1_1_GREEN_UP *
                                          jac_affine_inv_1_1_GREEN_UP) +
                                         (jac_affine_inv_1_2_GREEN_UP *
                                          jac_affine_inv_1_2_GREEN_UP)) +
                tmp_moved_constant_80 * tmp_moved_constant_83;
            const walberla::float64 tmp_moved_constant_87 =
                tmp_moved_constant_74 * tmp_moved_constant_79 +
                tmp_moved_constant_76 * tmp_moved_constant_82 +
                tmp_moved_constant_78 * tmp_moved_constant_83 +
                tmp_moved_constant_80 * ((jac_affine_inv_2_0_GREEN_UP *
                                          jac_affine_inv_2_0_GREEN_UP) +
                                         (jac_affine_inv_2_1_GREEN_UP *
                                          jac_affine_inv_2_1_GREEN_UP) +
                                         (jac_affine_inv_2_2_GREEN_UP *
                                          jac_affine_inv_2_2_GREEN_UP));
            {
              {
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
                    tmp_moved_constant_84 +
                    _data_dst[ctr_0 +
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
                           (6))] =
                    tmp_moved_constant_85 +
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
                    tmp_moved_constant_86 +
                    _data_dst[ctr_0 +
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
                    tmp_moved_constant_87 +
                    _data_dst[ctr_0 +
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
            const walberla::float64 tmp_moved_constant_88 =
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
            const walberla::float64 tmp_moved_constant_89 =
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
            const walberla::float64 tmp_moved_constant_90 =
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
            const walberla::float64 tmp_moved_constant_91 =
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
            const walberla::float64 tmp_moved_constant_96 =
                tmp_moved_constant_88 * tmp_moved_constant_95;
            const walberla::float64 tmp_moved_constant_98 =
                tmp_moved_constant_89 * tmp_moved_constant_95;
            const walberla::float64 tmp_moved_constant_100 =
                tmp_moved_constant_90 * tmp_moved_constant_95;
            const walberla::float64 tmp_moved_constant_102 =
                tmp_moved_constant_91 * tmp_moved_constant_95;
            const walberla::float64 tmp_moved_constant_106 =
                tmp_moved_constant_100 * tmp_moved_constant_99 +
                tmp_moved_constant_101 * tmp_moved_constant_102 +
                tmp_moved_constant_96 *
                    ((tmp_moved_constant_92 * tmp_moved_constant_92) +
                     (tmp_moved_constant_93 * tmp_moved_constant_93) +
                     (tmp_moved_constant_94 * tmp_moved_constant_94)) +
                tmp_moved_constant_97 * tmp_moved_constant_98;
            const walberla::float64 tmp_moved_constant_107 =
                tmp_moved_constant_100 * tmp_moved_constant_103 +
                tmp_moved_constant_102 * tmp_moved_constant_104 +
                tmp_moved_constant_96 * tmp_moved_constant_97 +
                tmp_moved_constant_98 * ((jac_affine_inv_0_0_GREEN_DOWN *
                                          jac_affine_inv_0_0_GREEN_DOWN) +
                                         (jac_affine_inv_0_1_GREEN_DOWN *
                                          jac_affine_inv_0_1_GREEN_DOWN) +
                                         (jac_affine_inv_0_2_GREEN_DOWN *
                                          jac_affine_inv_0_2_GREEN_DOWN));
            const walberla::float64 tmp_moved_constant_108 =
                tmp_moved_constant_100 * ((jac_affine_inv_1_0_GREEN_DOWN *
                                           jac_affine_inv_1_0_GREEN_DOWN) +
                                          (jac_affine_inv_1_1_GREEN_DOWN *
                                           jac_affine_inv_1_1_GREEN_DOWN) +
                                          (jac_affine_inv_1_2_GREEN_DOWN *
                                           jac_affine_inv_1_2_GREEN_DOWN)) +
                tmp_moved_constant_102 * tmp_moved_constant_105 +
                tmp_moved_constant_103 * tmp_moved_constant_98 +
                tmp_moved_constant_96 * tmp_moved_constant_99;
            const walberla::float64 tmp_moved_constant_109 =
                tmp_moved_constant_100 * tmp_moved_constant_105 +
                tmp_moved_constant_101 * tmp_moved_constant_96 +
                tmp_moved_constant_102 * ((jac_affine_inv_2_0_GREEN_DOWN *
                                           jac_affine_inv_2_0_GREEN_DOWN) +
                                          (jac_affine_inv_2_1_GREEN_DOWN *
                                           jac_affine_inv_2_1_GREEN_DOWN) +
                                          (jac_affine_inv_2_2_GREEN_DOWN *
                                           jac_affine_inv_2_2_GREEN_DOWN)) +
                tmp_moved_constant_104 * tmp_moved_constant_98;
            {
              {
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
                           (6))] =
                    tmp_moved_constant_106 +
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
                    tmp_moved_constant_107 +
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
                    tmp_moved_constant_108 +
                    _data_dst[ctr_0 +
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
                           (6))] =
                    tmp_moved_constant_109 +
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
        if (-ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2 >= 0) {
          const walberla::float64 src_dof_0 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2];
          const walberla::float64 src_dof_1 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1];
          const walberla::float64 src_dof_2 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2];
          const walberla::float64 src_dof_3 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2];
          const walberla::float64 tmp_kernel_op_0 =
              -jac_affine_inv_0_0_WHITE_UP - jac_affine_inv_1_0_WHITE_UP -
              jac_affine_inv_2_0_WHITE_UP;
          const walberla::float64 tmp_kernel_op_1 =
              -jac_affine_inv_0_1_WHITE_UP - jac_affine_inv_1_1_WHITE_UP -
              jac_affine_inv_2_1_WHITE_UP;
          const walberla::float64 tmp_kernel_op_2 =
              -jac_affine_inv_0_2_WHITE_UP - jac_affine_inv_1_2_WHITE_UP -
              jac_affine_inv_2_2_WHITE_UP;
          const walberla::float64 tmp_kernel_op_3 =
              abs_det_jac_affine_WHITE_UP * 0.16666666666666663;
          const walberla::float64 tmp_kernel_op_4 = src_dof_0 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_5 =
              jac_affine_inv_0_0_WHITE_UP * tmp_kernel_op_0 +
              jac_affine_inv_0_1_WHITE_UP * tmp_kernel_op_1 +
              jac_affine_inv_0_2_WHITE_UP * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_6 = src_dof_1 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_7 =
              jac_affine_inv_1_0_WHITE_UP * tmp_kernel_op_0 +
              jac_affine_inv_1_1_WHITE_UP * tmp_kernel_op_1 +
              jac_affine_inv_1_2_WHITE_UP * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_8 = src_dof_2 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_9 =
              jac_affine_inv_2_0_WHITE_UP * tmp_kernel_op_0 +
              jac_affine_inv_2_1_WHITE_UP * tmp_kernel_op_1 +
              jac_affine_inv_2_2_WHITE_UP * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_10 =
              src_dof_3 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_11 =
              jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_1_0_WHITE_UP +
              jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_1_1_WHITE_UP +
              jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_1_2_WHITE_UP;
          const walberla::float64 tmp_kernel_op_12 =
              jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP +
              jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP +
              jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_2_2_WHITE_UP;
          const walberla::float64 tmp_kernel_op_13 =
              jac_affine_inv_1_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP +
              jac_affine_inv_1_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP +
              jac_affine_inv_1_2_WHITE_UP * jac_affine_inv_2_2_WHITE_UP;
          const walberla::float64 elMatVec_0 =
              tmp_kernel_op_10 * tmp_kernel_op_9 +
              tmp_kernel_op_4 * ((tmp_kernel_op_0 * tmp_kernel_op_0) +
                                 (tmp_kernel_op_1 * tmp_kernel_op_1) +
                                 (tmp_kernel_op_2 * tmp_kernel_op_2)) +
              tmp_kernel_op_5 * tmp_kernel_op_6 +
              tmp_kernel_op_7 * tmp_kernel_op_8;
          const walberla::float64 elMatVec_1 =
              tmp_kernel_op_10 * tmp_kernel_op_12 +
              tmp_kernel_op_11 * tmp_kernel_op_8 +
              tmp_kernel_op_4 * tmp_kernel_op_5 +
              tmp_kernel_op_6 *
                  ((jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_0_0_WHITE_UP) +
                   (jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_0_1_WHITE_UP) +
                   (jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_0_2_WHITE_UP));
          const walberla::float64 elMatVec_2 =
              tmp_kernel_op_10 * tmp_kernel_op_13 +
              tmp_kernel_op_11 * tmp_kernel_op_6 +
              tmp_kernel_op_4 * tmp_kernel_op_7 +
              tmp_kernel_op_8 *
                  ((jac_affine_inv_1_0_WHITE_UP * jac_affine_inv_1_0_WHITE_UP) +
                   (jac_affine_inv_1_1_WHITE_UP * jac_affine_inv_1_1_WHITE_UP) +
                   (jac_affine_inv_1_2_WHITE_UP * jac_affine_inv_1_2_WHITE_UP));
          const walberla::float64 elMatVec_3 =
              tmp_kernel_op_10 *
                  ((jac_affine_inv_2_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP) +
                   (jac_affine_inv_2_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP) +
                   (jac_affine_inv_2_2_WHITE_UP *
                    jac_affine_inv_2_2_WHITE_UP)) +
              tmp_kernel_op_12 * tmp_kernel_op_6 +
              tmp_kernel_op_13 * tmp_kernel_op_8 +
              tmp_kernel_op_4 * tmp_kernel_op_9;
          {
            {
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2] =
                  elMatVec_0 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) *
                              (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            2];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  elMatVec_1 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) *
                              (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2] =
                  elMatVec_2 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            2];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2] =
                  elMatVec_3 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            2];
            }
          }
          const walberla::float64 tmp_moved_constant_110 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_111 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_112 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_113 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_114 =
              -jac_affine_inv_0_0_BLUE_UP - jac_affine_inv_1_0_BLUE_UP -
              jac_affine_inv_2_0_BLUE_UP;
          const walberla::float64 tmp_moved_constant_115 =
              -jac_affine_inv_0_1_BLUE_UP - jac_affine_inv_1_1_BLUE_UP -
              jac_affine_inv_2_1_BLUE_UP;
          const walberla::float64 tmp_moved_constant_116 =
              -jac_affine_inv_0_2_BLUE_UP - jac_affine_inv_1_2_BLUE_UP -
              jac_affine_inv_2_2_BLUE_UP;
          const walberla::float64 tmp_moved_constant_117 =
              abs_det_jac_affine_BLUE_UP * 0.16666666666666663;
          const walberla::float64 tmp_moved_constant_118 =
              tmp_moved_constant_110 * tmp_moved_constant_117;
          const walberla::float64 tmp_moved_constant_119 =
              jac_affine_inv_0_0_BLUE_UP * tmp_moved_constant_114 +
              jac_affine_inv_0_1_BLUE_UP * tmp_moved_constant_115 +
              jac_affine_inv_0_2_BLUE_UP * tmp_moved_constant_116;
          const walberla::float64 tmp_moved_constant_120 =
              tmp_moved_constant_111 * tmp_moved_constant_117;
          const walberla::float64 tmp_moved_constant_121 =
              jac_affine_inv_1_0_BLUE_UP * tmp_moved_constant_114 +
              jac_affine_inv_1_1_BLUE_UP * tmp_moved_constant_115 +
              jac_affine_inv_1_2_BLUE_UP * tmp_moved_constant_116;
          const walberla::float64 tmp_moved_constant_122 =
              tmp_moved_constant_112 * tmp_moved_constant_117;
          const walberla::float64 tmp_moved_constant_123 =
              jac_affine_inv_2_0_BLUE_UP * tmp_moved_constant_114 +
              jac_affine_inv_2_1_BLUE_UP * tmp_moved_constant_115 +
              jac_affine_inv_2_2_BLUE_UP * tmp_moved_constant_116;
          const walberla::float64 tmp_moved_constant_124 =
              tmp_moved_constant_113 * tmp_moved_constant_117;
          const walberla::float64 tmp_moved_constant_125 =
              jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_1_0_BLUE_UP +
              jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_1_1_BLUE_UP +
              jac_affine_inv_0_2_BLUE_UP * jac_affine_inv_1_2_BLUE_UP;
          const walberla::float64 tmp_moved_constant_126 =
              jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP +
              jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP +
              jac_affine_inv_0_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP;
          const walberla::float64 tmp_moved_constant_127 =
              jac_affine_inv_1_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP +
              jac_affine_inv_1_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP +
              jac_affine_inv_1_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP;
          const walberla::float64 tmp_moved_constant_128 =
              tmp_moved_constant_118 *
                  ((tmp_moved_constant_114 * tmp_moved_constant_114) +
                   (tmp_moved_constant_115 * tmp_moved_constant_115) +
                   (tmp_moved_constant_116 * tmp_moved_constant_116)) +
              tmp_moved_constant_119 * tmp_moved_constant_120 +
              tmp_moved_constant_121 * tmp_moved_constant_122 +
              tmp_moved_constant_123 * tmp_moved_constant_124;
          const walberla::float64 tmp_moved_constant_129 =
              tmp_moved_constant_118 * tmp_moved_constant_119 +
              tmp_moved_constant_120 *
                  ((jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_0_0_BLUE_UP) +
                   (jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_0_1_BLUE_UP) +
                   (jac_affine_inv_0_2_BLUE_UP * jac_affine_inv_0_2_BLUE_UP)) +
              tmp_moved_constant_122 * tmp_moved_constant_125 +
              tmp_moved_constant_124 * tmp_moved_constant_126;
          const walberla::float64 tmp_moved_constant_130 =
              tmp_moved_constant_118 * tmp_moved_constant_121 +
              tmp_moved_constant_120 * tmp_moved_constant_125 +
              tmp_moved_constant_122 *
                  ((jac_affine_inv_1_0_BLUE_UP * jac_affine_inv_1_0_BLUE_UP) +
                   (jac_affine_inv_1_1_BLUE_UP * jac_affine_inv_1_1_BLUE_UP) +
                   (jac_affine_inv_1_2_BLUE_UP * jac_affine_inv_1_2_BLUE_UP)) +
              tmp_moved_constant_124 * tmp_moved_constant_127;
          const walberla::float64 tmp_moved_constant_131 =
              tmp_moved_constant_118 * tmp_moved_constant_123 +
              tmp_moved_constant_120 * tmp_moved_constant_126 +
              tmp_moved_constant_122 * tmp_moved_constant_127 +
              tmp_moved_constant_124 *
                  ((jac_affine_inv_2_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP) +
                   (jac_affine_inv_2_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP) +
                   (jac_affine_inv_2_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP));
          {
            {
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  tmp_moved_constant_128 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) *
                              (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2] =
                  tmp_moved_constant_129 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            2];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        1] =
                  tmp_moved_constant_130 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            1];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  tmp_moved_constant_131 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
            }
          }
          const walberla::float64 tmp_moved_constant_132 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_133 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_134 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_135 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_136 =
              -jac_affine_inv_0_0_BLUE_DOWN - jac_affine_inv_1_0_BLUE_DOWN -
              jac_affine_inv_2_0_BLUE_DOWN;
          const walberla::float64 tmp_moved_constant_137 =
              -jac_affine_inv_0_1_BLUE_DOWN - jac_affine_inv_1_1_BLUE_DOWN -
              jac_affine_inv_2_1_BLUE_DOWN;
          const walberla::float64 tmp_moved_constant_138 =
              -jac_affine_inv_0_2_BLUE_DOWN - jac_affine_inv_1_2_BLUE_DOWN -
              jac_affine_inv_2_2_BLUE_DOWN;
          const walberla::float64 tmp_moved_constant_139 =
              abs_det_jac_affine_BLUE_DOWN * 0.16666666666666663;
          const walberla::float64 tmp_moved_constant_140 =
              tmp_moved_constant_132 * tmp_moved_constant_139;
          const walberla::float64 tmp_moved_constant_141 =
              jac_affine_inv_0_0_BLUE_DOWN * tmp_moved_constant_136 +
              jac_affine_inv_0_1_BLUE_DOWN * tmp_moved_constant_137 +
              jac_affine_inv_0_2_BLUE_DOWN * tmp_moved_constant_138;
          const walberla::float64 tmp_moved_constant_142 =
              tmp_moved_constant_133 * tmp_moved_constant_139;
          const walberla::float64 tmp_moved_constant_143 =
              jac_affine_inv_1_0_BLUE_DOWN * tmp_moved_constant_136 +
              jac_affine_inv_1_1_BLUE_DOWN * tmp_moved_constant_137 +
              jac_affine_inv_1_2_BLUE_DOWN * tmp_moved_constant_138;
          const walberla::float64 tmp_moved_constant_144 =
              tmp_moved_constant_134 * tmp_moved_constant_139;
          const walberla::float64 tmp_moved_constant_145 =
              jac_affine_inv_2_0_BLUE_DOWN * tmp_moved_constant_136 +
              jac_affine_inv_2_1_BLUE_DOWN * tmp_moved_constant_137 +
              jac_affine_inv_2_2_BLUE_DOWN * tmp_moved_constant_138;
          const walberla::float64 tmp_moved_constant_146 =
              tmp_moved_constant_135 * tmp_moved_constant_139;
          const walberla::float64 tmp_moved_constant_147 =
              jac_affine_inv_0_0_BLUE_DOWN * jac_affine_inv_1_0_BLUE_DOWN +
              jac_affine_inv_0_1_BLUE_DOWN * jac_affine_inv_1_1_BLUE_DOWN +
              jac_affine_inv_0_2_BLUE_DOWN * jac_affine_inv_1_2_BLUE_DOWN;
          const walberla::float64 tmp_moved_constant_148 =
              jac_affine_inv_0_0_BLUE_DOWN * jac_affine_inv_2_0_BLUE_DOWN +
              jac_affine_inv_0_1_BLUE_DOWN * jac_affine_inv_2_1_BLUE_DOWN +
              jac_affine_inv_0_2_BLUE_DOWN * jac_affine_inv_2_2_BLUE_DOWN;
          const walberla::float64 tmp_moved_constant_149 =
              jac_affine_inv_1_0_BLUE_DOWN * jac_affine_inv_2_0_BLUE_DOWN +
              jac_affine_inv_1_1_BLUE_DOWN * jac_affine_inv_2_1_BLUE_DOWN +
              jac_affine_inv_1_2_BLUE_DOWN * jac_affine_inv_2_2_BLUE_DOWN;
          const walberla::float64 tmp_moved_constant_150 =
              tmp_moved_constant_140 *
                  ((tmp_moved_constant_136 * tmp_moved_constant_136) +
                   (tmp_moved_constant_137 * tmp_moved_constant_137) +
                   (tmp_moved_constant_138 * tmp_moved_constant_138)) +
              tmp_moved_constant_141 * tmp_moved_constant_142 +
              tmp_moved_constant_143 * tmp_moved_constant_144 +
              tmp_moved_constant_145 * tmp_moved_constant_146;
          const walberla::float64 tmp_moved_constant_151 =
              tmp_moved_constant_140 * tmp_moved_constant_141 +
              tmp_moved_constant_142 * ((jac_affine_inv_0_0_BLUE_DOWN *
                                         jac_affine_inv_0_0_BLUE_DOWN) +
                                        (jac_affine_inv_0_1_BLUE_DOWN *
                                         jac_affine_inv_0_1_BLUE_DOWN) +
                                        (jac_affine_inv_0_2_BLUE_DOWN *
                                         jac_affine_inv_0_2_BLUE_DOWN)) +
              tmp_moved_constant_144 * tmp_moved_constant_147 +
              tmp_moved_constant_146 * tmp_moved_constant_148;
          const walberla::float64 tmp_moved_constant_152 =
              tmp_moved_constant_140 * tmp_moved_constant_143 +
              tmp_moved_constant_142 * tmp_moved_constant_147 +
              tmp_moved_constant_144 * ((jac_affine_inv_1_0_BLUE_DOWN *
                                         jac_affine_inv_1_0_BLUE_DOWN) +
                                        (jac_affine_inv_1_1_BLUE_DOWN *
                                         jac_affine_inv_1_1_BLUE_DOWN) +
                                        (jac_affine_inv_1_2_BLUE_DOWN *
                                         jac_affine_inv_1_2_BLUE_DOWN)) +
              tmp_moved_constant_146 * tmp_moved_constant_149;
          const walberla::float64 tmp_moved_constant_153 =
              tmp_moved_constant_140 * tmp_moved_constant_145 +
              tmp_moved_constant_142 * tmp_moved_constant_148 +
              tmp_moved_constant_144 * tmp_moved_constant_149 +
              tmp_moved_constant_146 * ((jac_affine_inv_2_0_BLUE_DOWN *
                                         jac_affine_inv_2_0_BLUE_DOWN) +
                                        (jac_affine_inv_2_1_BLUE_DOWN *
                                         jac_affine_inv_2_1_BLUE_DOWN) +
                                        (jac_affine_inv_2_2_BLUE_DOWN *
                                         jac_affine_inv_2_2_BLUE_DOWN));
          {
            {
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2] =
                  tmp_moved_constant_150 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            2];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2] =
                  tmp_moved_constant_151 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            2];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  tmp_moved_constant_152 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2] =
                  tmp_moved_constant_153 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            2];
            }
          }
          const walberla::float64 tmp_moved_constant_154 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_155 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_156 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_157 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_158 =
              -jac_affine_inv_0_0_GREEN_UP - jac_affine_inv_1_0_GREEN_UP -
              jac_affine_inv_2_0_GREEN_UP;
          const walberla::float64 tmp_moved_constant_159 =
              -jac_affine_inv_0_1_GREEN_UP - jac_affine_inv_1_1_GREEN_UP -
              jac_affine_inv_2_1_GREEN_UP;
          const walberla::float64 tmp_moved_constant_160 =
              -jac_affine_inv_0_2_GREEN_UP - jac_affine_inv_1_2_GREEN_UP -
              jac_affine_inv_2_2_GREEN_UP;
          const walberla::float64 tmp_moved_constant_161 =
              abs_det_jac_affine_GREEN_UP * 0.16666666666666663;
          const walberla::float64 tmp_moved_constant_162 =
              tmp_moved_constant_154 * tmp_moved_constant_161;
          const walberla::float64 tmp_moved_constant_163 =
              jac_affine_inv_0_0_GREEN_UP * tmp_moved_constant_158 +
              jac_affine_inv_0_1_GREEN_UP * tmp_moved_constant_159 +
              jac_affine_inv_0_2_GREEN_UP * tmp_moved_constant_160;
          const walberla::float64 tmp_moved_constant_164 =
              tmp_moved_constant_155 * tmp_moved_constant_161;
          const walberla::float64 tmp_moved_constant_165 =
              jac_affine_inv_1_0_GREEN_UP * tmp_moved_constant_158 +
              jac_affine_inv_1_1_GREEN_UP * tmp_moved_constant_159 +
              jac_affine_inv_1_2_GREEN_UP * tmp_moved_constant_160;
          const walberla::float64 tmp_moved_constant_166 =
              tmp_moved_constant_156 * tmp_moved_constant_161;
          const walberla::float64 tmp_moved_constant_167 =
              jac_affine_inv_2_0_GREEN_UP * tmp_moved_constant_158 +
              jac_affine_inv_2_1_GREEN_UP * tmp_moved_constant_159 +
              jac_affine_inv_2_2_GREEN_UP * tmp_moved_constant_160;
          const walberla::float64 tmp_moved_constant_168 =
              tmp_moved_constant_157 * tmp_moved_constant_161;
          const walberla::float64 tmp_moved_constant_169 =
              jac_affine_inv_0_0_GREEN_UP * jac_affine_inv_1_0_GREEN_UP +
              jac_affine_inv_0_1_GREEN_UP * jac_affine_inv_1_1_GREEN_UP +
              jac_affine_inv_0_2_GREEN_UP * jac_affine_inv_1_2_GREEN_UP;
          const walberla::float64 tmp_moved_constant_170 =
              jac_affine_inv_0_0_GREEN_UP * jac_affine_inv_2_0_GREEN_UP +
              jac_affine_inv_0_1_GREEN_UP * jac_affine_inv_2_1_GREEN_UP +
              jac_affine_inv_0_2_GREEN_UP * jac_affine_inv_2_2_GREEN_UP;
          const walberla::float64 tmp_moved_constant_171 =
              jac_affine_inv_1_0_GREEN_UP * jac_affine_inv_2_0_GREEN_UP +
              jac_affine_inv_1_1_GREEN_UP * jac_affine_inv_2_1_GREEN_UP +
              jac_affine_inv_1_2_GREEN_UP * jac_affine_inv_2_2_GREEN_UP;
          const walberla::float64 tmp_moved_constant_172 =
              tmp_moved_constant_162 *
                  ((tmp_moved_constant_158 * tmp_moved_constant_158) +
                   (tmp_moved_constant_159 * tmp_moved_constant_159) +
                   (tmp_moved_constant_160 * tmp_moved_constant_160)) +
              tmp_moved_constant_163 * tmp_moved_constant_164 +
              tmp_moved_constant_165 * tmp_moved_constant_166 +
              tmp_moved_constant_167 * tmp_moved_constant_168;
          const walberla::float64 tmp_moved_constant_173 =
              tmp_moved_constant_162 * tmp_moved_constant_163 +
              tmp_moved_constant_164 *
                  ((jac_affine_inv_0_0_GREEN_UP * jac_affine_inv_0_0_GREEN_UP) +
                   (jac_affine_inv_0_1_GREEN_UP * jac_affine_inv_0_1_GREEN_UP) +
                   (jac_affine_inv_0_2_GREEN_UP *
                    jac_affine_inv_0_2_GREEN_UP)) +
              tmp_moved_constant_166 * tmp_moved_constant_169 +
              tmp_moved_constant_168 * tmp_moved_constant_170;
          const walberla::float64 tmp_moved_constant_174 =
              tmp_moved_constant_162 * tmp_moved_constant_165 +
              tmp_moved_constant_164 * tmp_moved_constant_169 +
              tmp_moved_constant_166 *
                  ((jac_affine_inv_1_0_GREEN_UP * jac_affine_inv_1_0_GREEN_UP) +
                   (jac_affine_inv_1_1_GREEN_UP * jac_affine_inv_1_1_GREEN_UP) +
                   (jac_affine_inv_1_2_GREEN_UP *
                    jac_affine_inv_1_2_GREEN_UP)) +
              tmp_moved_constant_168 * tmp_moved_constant_171;
          const walberla::float64 tmp_moved_constant_175 =
              tmp_moved_constant_162 * tmp_moved_constant_167 +
              tmp_moved_constant_164 * tmp_moved_constant_170 +
              tmp_moved_constant_166 * tmp_moved_constant_171 +
              tmp_moved_constant_168 *
                  ((jac_affine_inv_2_0_GREEN_UP * jac_affine_inv_2_0_GREEN_UP) +
                   (jac_affine_inv_2_1_GREEN_UP * jac_affine_inv_2_1_GREEN_UP) +
                   (jac_affine_inv_2_2_GREEN_UP * jac_affine_inv_2_2_GREEN_UP));
          {
            {
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  tmp_moved_constant_172 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) *
                              (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2] =
                  tmp_moved_constant_173 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            2];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        2] =
                  tmp_moved_constant_174 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            2];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  tmp_moved_constant_175 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
            }
          }
          const walberla::float64 tmp_moved_constant_176 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_177 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_178 =
              _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1];
          const walberla::float64 tmp_moved_constant_179 =
              _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2];
          const walberla::float64 tmp_moved_constant_180 =
              -jac_affine_inv_0_0_GREEN_DOWN - jac_affine_inv_1_0_GREEN_DOWN -
              jac_affine_inv_2_0_GREEN_DOWN;
          const walberla::float64 tmp_moved_constant_181 =
              -jac_affine_inv_0_1_GREEN_DOWN - jac_affine_inv_1_1_GREEN_DOWN -
              jac_affine_inv_2_1_GREEN_DOWN;
          const walberla::float64 tmp_moved_constant_182 =
              -jac_affine_inv_0_2_GREEN_DOWN - jac_affine_inv_1_2_GREEN_DOWN -
              jac_affine_inv_2_2_GREEN_DOWN;
          const walberla::float64 tmp_moved_constant_183 =
              abs_det_jac_affine_GREEN_DOWN * 0.16666666666666663;
          const walberla::float64 tmp_moved_constant_184 =
              tmp_moved_constant_176 * tmp_moved_constant_183;
          const walberla::float64 tmp_moved_constant_185 =
              jac_affine_inv_0_0_GREEN_DOWN * tmp_moved_constant_180 +
              jac_affine_inv_0_1_GREEN_DOWN * tmp_moved_constant_181 +
              jac_affine_inv_0_2_GREEN_DOWN * tmp_moved_constant_182;
          const walberla::float64 tmp_moved_constant_186 =
              tmp_moved_constant_177 * tmp_moved_constant_183;
          const walberla::float64 tmp_moved_constant_187 =
              jac_affine_inv_1_0_GREEN_DOWN * tmp_moved_constant_180 +
              jac_affine_inv_1_1_GREEN_DOWN * tmp_moved_constant_181 +
              jac_affine_inv_1_2_GREEN_DOWN * tmp_moved_constant_182;
          const walberla::float64 tmp_moved_constant_188 =
              tmp_moved_constant_178 * tmp_moved_constant_183;
          const walberla::float64 tmp_moved_constant_189 =
              jac_affine_inv_2_0_GREEN_DOWN * tmp_moved_constant_180 +
              jac_affine_inv_2_1_GREEN_DOWN * tmp_moved_constant_181 +
              jac_affine_inv_2_2_GREEN_DOWN * tmp_moved_constant_182;
          const walberla::float64 tmp_moved_constant_190 =
              tmp_moved_constant_179 * tmp_moved_constant_183;
          const walberla::float64 tmp_moved_constant_191 =
              jac_affine_inv_0_0_GREEN_DOWN * jac_affine_inv_1_0_GREEN_DOWN +
              jac_affine_inv_0_1_GREEN_DOWN * jac_affine_inv_1_1_GREEN_DOWN +
              jac_affine_inv_0_2_GREEN_DOWN * jac_affine_inv_1_2_GREEN_DOWN;
          const walberla::float64 tmp_moved_constant_192 =
              jac_affine_inv_0_0_GREEN_DOWN * jac_affine_inv_2_0_GREEN_DOWN +
              jac_affine_inv_0_1_GREEN_DOWN * jac_affine_inv_2_1_GREEN_DOWN +
              jac_affine_inv_0_2_GREEN_DOWN * jac_affine_inv_2_2_GREEN_DOWN;
          const walberla::float64 tmp_moved_constant_193 =
              jac_affine_inv_1_0_GREEN_DOWN * jac_affine_inv_2_0_GREEN_DOWN +
              jac_affine_inv_1_1_GREEN_DOWN * jac_affine_inv_2_1_GREEN_DOWN +
              jac_affine_inv_1_2_GREEN_DOWN * jac_affine_inv_2_2_GREEN_DOWN;
          const walberla::float64 tmp_moved_constant_194 =
              tmp_moved_constant_184 *
                  ((tmp_moved_constant_180 * tmp_moved_constant_180) +
                   (tmp_moved_constant_181 * tmp_moved_constant_181) +
                   (tmp_moved_constant_182 * tmp_moved_constant_182)) +
              tmp_moved_constant_185 * tmp_moved_constant_186 +
              tmp_moved_constant_187 * tmp_moved_constant_188 +
              tmp_moved_constant_189 * tmp_moved_constant_190;
          const walberla::float64 tmp_moved_constant_195 =
              tmp_moved_constant_184 * tmp_moved_constant_185 +
              tmp_moved_constant_186 * ((jac_affine_inv_0_0_GREEN_DOWN *
                                         jac_affine_inv_0_0_GREEN_DOWN) +
                                        (jac_affine_inv_0_1_GREEN_DOWN *
                                         jac_affine_inv_0_1_GREEN_DOWN) +
                                        (jac_affine_inv_0_2_GREEN_DOWN *
                                         jac_affine_inv_0_2_GREEN_DOWN)) +
              tmp_moved_constant_188 * tmp_moved_constant_191 +
              tmp_moved_constant_190 * tmp_moved_constant_192;
          const walberla::float64 tmp_moved_constant_196 =
              tmp_moved_constant_184 * tmp_moved_constant_187 +
              tmp_moved_constant_186 * tmp_moved_constant_191 +
              tmp_moved_constant_188 * ((jac_affine_inv_1_0_GREEN_DOWN *
                                         jac_affine_inv_1_0_GREEN_DOWN) +
                                        (jac_affine_inv_1_1_GREEN_DOWN *
                                         jac_affine_inv_1_1_GREEN_DOWN) +
                                        (jac_affine_inv_1_2_GREEN_DOWN *
                                         jac_affine_inv_1_2_GREEN_DOWN)) +
              tmp_moved_constant_190 * tmp_moved_constant_193;
          const walberla::float64 tmp_moved_constant_197 =
              tmp_moved_constant_184 * tmp_moved_constant_189 +
              tmp_moved_constant_186 * tmp_moved_constant_192 +
              tmp_moved_constant_188 * tmp_moved_constant_193 +
              tmp_moved_constant_190 * ((jac_affine_inv_2_0_GREEN_DOWN *
                                         jac_affine_inv_2_0_GREEN_DOWN) +
                                        (jac_affine_inv_2_1_GREEN_DOWN *
                                         jac_affine_inv_2_1_GREEN_DOWN) +
                                        (jac_affine_inv_2_2_GREEN_DOWN *
                                         jac_affine_inv_2_2_GREEN_DOWN));
          {
            {
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2] =
                  tmp_moved_constant_194 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            2];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        1] =
                  tmp_moved_constant_195 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            1];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  tmp_moved_constant_196 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        2] =
                  tmp_moved_constant_197 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            2];
            }
          }
        }
        const walberla::float64 src_dof_0 =
            _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      1];
        const walberla::float64 src_dof_1 =
            _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))];
        const walberla::float64 src_dof_2 =
            _data_src[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      1];
        const walberla::float64 src_dof_3 =
            _data_src[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                      ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) -
                      (((-ctr_2 + micro_edges_per_macro_edge) *
                        (-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                       (6)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      1];
        const walberla::float64 tmp_kernel_op_4 = src_dof_0 * tmp_kernel_op_3;
        const walberla::float64 tmp_kernel_op_6 = src_dof_1 * tmp_kernel_op_3;
        const walberla::float64 tmp_kernel_op_8 = src_dof_2 * tmp_kernel_op_3;
        const walberla::float64 tmp_kernel_op_10 = src_dof_3 * tmp_kernel_op_3;
        const walberla::float64 elMatVec_0 =
            tmp_kernel_op_10 * tmp_kernel_op_9 +
            tmp_kernel_op_4 * ((tmp_kernel_op_0 * tmp_kernel_op_0) +
                               (tmp_kernel_op_1 * tmp_kernel_op_1) +
                               (tmp_kernel_op_2 * tmp_kernel_op_2)) +
            tmp_kernel_op_5 * tmp_kernel_op_6 +
            tmp_kernel_op_7 * tmp_kernel_op_8;
        const walberla::float64 elMatVec_1 =
            tmp_kernel_op_10 * tmp_kernel_op_12 +
            tmp_kernel_op_11 * tmp_kernel_op_8 +
            tmp_kernel_op_4 * tmp_kernel_op_5 +
            tmp_kernel_op_6 *
                ((jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_0_0_WHITE_UP) +
                 (jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_0_1_WHITE_UP) +
                 (jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_0_2_WHITE_UP));
        const walberla::float64 elMatVec_2 =
            tmp_kernel_op_10 * tmp_kernel_op_13 +
            tmp_kernel_op_11 * tmp_kernel_op_6 +
            tmp_kernel_op_4 * tmp_kernel_op_7 +
            tmp_kernel_op_8 *
                ((jac_affine_inv_1_0_WHITE_UP * jac_affine_inv_1_0_WHITE_UP) +
                 (jac_affine_inv_1_1_WHITE_UP * jac_affine_inv_1_1_WHITE_UP) +
                 (jac_affine_inv_1_2_WHITE_UP * jac_affine_inv_1_2_WHITE_UP));
        const walberla::float64 elMatVec_3 =
            tmp_kernel_op_10 *
                ((jac_affine_inv_2_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP) +
                 (jac_affine_inv_2_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP) +
                 (jac_affine_inv_2_2_WHITE_UP * jac_affine_inv_2_2_WHITE_UP)) +
            tmp_kernel_op_12 * tmp_kernel_op_6 +
            tmp_kernel_op_13 * tmp_kernel_op_8 +
            tmp_kernel_op_4 * tmp_kernel_op_9;
        {
          {
            {
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  elMatVec_0 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) *
                              (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2) *
                          (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                         (6))] =
                  elMatVec_1 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2) *
                              (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                             (6))];
              _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                         (6)) -
                        1] =
                  elMatVec_2 +
                  _data_dst[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                             (6)) -
                            1];
              _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                        ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) -
                        (((-ctr_2 + micro_edges_per_macro_edge) *
                          (-ctr_2 + micro_edges_per_macro_edge + 1) *
                          (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                         (6)) +
                        (((micro_edges_per_macro_edge + 1) *
                          (micro_edges_per_macro_edge + 2) *
                          (micro_edges_per_macro_edge + 3)) /
                         (6)) -
                        1] =
                  elMatVec_3 +
                  _data_dst[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6)) -
                            1];
            }
          }
        }
      }
  }
}
void P1ElementwiseDiffusion_cubes_const_vect_float64::
    computeInverseDiagonalOperatorValues_macro_2D(
        walberla::float64 *RESTRICT _data_invDiag_,
        walberla::float64 macro_vertex_coord_id_0comp0,
        walberla::float64 macro_vertex_coord_id_0comp1,
        walberla::float64 macro_vertex_coord_id_1comp0,
        walberla::float64 macro_vertex_coord_id_1comp1,
        walberla::float64 macro_vertex_coord_id_2comp0,
        walberla::float64 macro_vertex_coord_id_2comp1,
        int64_t micro_edges_per_macro_edge,
        walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    const walberla::float64 tmp_coords_jac_0_BLUE =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
    const walberla::float64 tmp_coords_jac_1_BLUE =
        macro_vertex_coord_id_0comp0 +
        tmp_coords_jac_0_BLUE *
            (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
    const walberla::float64 tmp_coords_jac_2_BLUE =
        macro_vertex_coord_id_0comp1 +
        tmp_coords_jac_0_BLUE *
            (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
    const walberla::float64 tmp_coords_jac_3_BLUE =
        tmp_coords_jac_0_BLUE *
        (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
    const walberla::float64 tmp_coords_jac_4_BLUE =
        tmp_coords_jac_0_BLUE *
        (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
    const walberla::float64 p_affine_const_0_0_BLUE = tmp_coords_jac_1_BLUE;
    const walberla::float64 p_affine_const_0_1_BLUE = tmp_coords_jac_2_BLUE;
    const walberla::float64 p_affine_const_1_0_BLUE =
        macro_vertex_coord_id_0comp0 + tmp_coords_jac_3_BLUE;
    const walberla::float64 p_affine_const_1_1_BLUE =
        macro_vertex_coord_id_0comp1 + tmp_coords_jac_4_BLUE;
    const walberla::float64 p_affine_const_2_0_BLUE =
        tmp_coords_jac_1_BLUE + tmp_coords_jac_3_BLUE;
    const walberla::float64 p_affine_const_2_1_BLUE =
        tmp_coords_jac_2_BLUE + tmp_coords_jac_4_BLUE;
    const walberla::float64 jac_affine_0_0_BLUE =
        -p_affine_const_0_0_BLUE + p_affine_const_1_0_BLUE;
    const walberla::float64 jac_affine_0_1_BLUE =
        -p_affine_const_0_0_BLUE + p_affine_const_2_0_BLUE;
    const walberla::float64 jac_affine_1_0_BLUE =
        -p_affine_const_0_1_BLUE + p_affine_const_1_1_BLUE;
    const walberla::float64 jac_affine_1_1_BLUE =
        -p_affine_const_0_1_BLUE + p_affine_const_2_1_BLUE;
    const walberla::float64 tmp_coords_jac_5_BLUE =
        jac_affine_0_0_BLUE * jac_affine_1_1_BLUE -
        jac_affine_0_1_BLUE * jac_affine_1_0_BLUE;
    const walberla::float64 tmp_coords_jac_6_BLUE =
        1.0 / (tmp_coords_jac_5_BLUE);
    const walberla::float64 jac_affine_inv_0_0_BLUE =
        jac_affine_1_1_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 jac_affine_inv_0_1_BLUE =
        -jac_affine_0_1_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 jac_affine_inv_1_0_BLUE =
        -jac_affine_1_0_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 jac_affine_inv_1_1_BLUE =
        jac_affine_0_0_BLUE * tmp_coords_jac_6_BLUE;
    const walberla::float64 abs_det_jac_affine_BLUE =
        abs(tmp_coords_jac_5_BLUE);
    const walberla::float64 tmp_coords_jac_0_GRAY =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
    const walberla::float64 p_affine_const_0_0_GRAY =
        macro_vertex_coord_id_0comp0;
    const walberla::float64 p_affine_const_0_1_GRAY =
        macro_vertex_coord_id_0comp1;
    const walberla::float64 p_affine_const_1_0_GRAY =
        macro_vertex_coord_id_0comp0 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_1comp0);
    const walberla::float64 p_affine_const_1_1_GRAY =
        macro_vertex_coord_id_0comp1 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_1comp1);
    const walberla::float64 p_affine_const_2_0_GRAY =
        macro_vertex_coord_id_0comp0 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp0 + macro_vertex_coord_id_2comp0);
    const walberla::float64 p_affine_const_2_1_GRAY =
        macro_vertex_coord_id_0comp1 +
        tmp_coords_jac_0_GRAY *
            (-macro_vertex_coord_id_0comp1 + macro_vertex_coord_id_2comp1);
    const walberla::float64 jac_affine_0_0_GRAY =
        -p_affine_const_0_0_GRAY + p_affine_const_1_0_GRAY;
    const walberla::float64 jac_affine_0_1_GRAY =
        -p_affine_const_0_0_GRAY + p_affine_const_2_0_GRAY;
    const walberla::float64 jac_affine_1_0_GRAY =
        -p_affine_const_0_1_GRAY + p_affine_const_1_1_GRAY;
    const walberla::float64 jac_affine_1_1_GRAY =
        -p_affine_const_0_1_GRAY + p_affine_const_2_1_GRAY;
    const walberla::float64 tmp_coords_jac_1_GRAY =
        jac_affine_0_0_GRAY * jac_affine_1_1_GRAY -
        jac_affine_0_1_GRAY * jac_affine_1_0_GRAY;
    const walberla::float64 tmp_coords_jac_2_GRAY =
        1.0 / (tmp_coords_jac_1_GRAY);
    const walberla::float64 jac_affine_inv_0_0_GRAY =
        jac_affine_1_1_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 jac_affine_inv_0_1_GRAY =
        -jac_affine_0_1_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 jac_affine_inv_1_0_GRAY =
        -jac_affine_1_0_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 jac_affine_inv_1_1_GRAY =
        jac_affine_0_0_GRAY * tmp_coords_jac_2_GRAY;
    const walberla::float64 abs_det_jac_affine_GRAY =
        abs(tmp_coords_jac_1_GRAY);
    const walberla::float64 tmp_kernel_op_0 = abs_det_jac_affine_GRAY * 0.5;
    const walberla::float64 elMatDiag_0 =
        tmp_kernel_op_0 *
        (((-jac_affine_inv_0_0_GRAY - jac_affine_inv_1_0_GRAY) *
          (-jac_affine_inv_0_0_GRAY - jac_affine_inv_1_0_GRAY)) +
         ((-jac_affine_inv_0_1_GRAY - jac_affine_inv_1_1_GRAY) *
          (-jac_affine_inv_0_1_GRAY - jac_affine_inv_1_1_GRAY)));
    const walberla::float64 elMatDiag_1 =
        tmp_kernel_op_0 * ((jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY) +
                           (jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY));
    const walberla::float64 elMatDiag_2 =
        tmp_kernel_op_0 * ((jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY) +
                           (jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY));
    const walberla::float64 tmp_moved_constant_0 =
        abs_det_jac_affine_BLUE * 0.5;
    const walberla::float64 tmp_moved_constant_1 =
        tmp_moved_constant_0 *
        (((-jac_affine_inv_0_0_BLUE - jac_affine_inv_1_0_BLUE) *
          (-jac_affine_inv_0_0_BLUE - jac_affine_inv_1_0_BLUE)) +
         ((-jac_affine_inv_0_1_BLUE - jac_affine_inv_1_1_BLUE) *
          (-jac_affine_inv_0_1_BLUE - jac_affine_inv_1_1_BLUE)));
    const walberla::float64 tmp_moved_constant_2 =
        tmp_moved_constant_0 *
        ((jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE) +
         (jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE));
    const walberla::float64 tmp_moved_constant_3 =
        tmp_moved_constant_0 *
        ((jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE) +
         (jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE));
    for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1) {
      {
        for (int64_t ctr_0 = 0;
             ctr_0 <
             (int64_t)((-ctr_1 + micro_edges_per_macro_edge - 1) / (4)) * (4);
             ctr_0 += 4) {
          {{_mm256_storeu_pd(
              &_data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2))],
              _mm256_add_pd(
                  _mm256_set_pd(elMatDiag_0, elMatDiag_0, elMatDiag_0,
                                elMatDiag_0),
                  _mm256_loadu_pd(
                      &_data_invDiag_[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 2) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2))])));
          _mm256_storeu_pd(
              &_data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2)) + 1],
              _mm256_add_pd(
                  _mm256_set_pd(elMatDiag_1, elMatDiag_1, elMatDiag_1,
                                elMatDiag_1),
                  _mm256_loadu_pd(
                      &_data_invDiag_[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 2) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) + 1])));
          _mm256_storeu_pd(
              &_data_invDiag_[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2))],
              _mm256_add_pd(
                  _mm256_set_pd(elMatDiag_2, elMatDiag_2, elMatDiag_2,
                                elMatDiag_2),
                  _mm256_loadu_pd(
                      &_data_invDiag_[ctr_0 +
                                      (ctr_1 + 1) *
                                          (micro_edges_per_macro_edge + 2) -
                                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2))])));
        }
      }
      {
        {
          _mm256_storeu_pd(
              &_data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2)) + 1],
              _mm256_add_pd(
                  _mm256_set_pd(tmp_moved_constant_1, tmp_moved_constant_1,
                                tmp_moved_constant_1, tmp_moved_constant_1),
                  _mm256_loadu_pd(
                      &_data_invDiag_[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 2) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) + 1])));
          _mm256_storeu_pd(
              &_data_invDiag_[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2))],
              _mm256_add_pd(
                  _mm256_set_pd(tmp_moved_constant_2, tmp_moved_constant_2,
                                tmp_moved_constant_2, tmp_moved_constant_2),
                  _mm256_loadu_pd(
                      &_data_invDiag_[ctr_0 +
                                      (ctr_1 + 1) *
                                          (micro_edges_per_macro_edge + 2) -
                                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2))])));
          _mm256_storeu_pd(
              &_data_invDiag_[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1],
              _mm256_add_pd(
                  _mm256_set_pd(tmp_moved_constant_3, tmp_moved_constant_3,
                                tmp_moved_constant_3, tmp_moved_constant_3),
                  _mm256_loadu_pd(
                      &_data_invDiag_
                          [ctr_0 +
                           (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                           (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1])));
        }
      }
    }
    for (int64_t ctr_0 =
             (int64_t)((-ctr_1 + micro_edges_per_macro_edge - 1) / (4)) * (4);
         ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1; ctr_0 += 1) {
      {{_data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                       ((ctr_1 * (ctr_1 + 1)) / (2))] =
            elMatDiag_0 +
            _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                           ((ctr_1 * (ctr_1 + 1)) / (2))];
      _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
          elMatDiag_1 +
          _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
      _data_invDiag_[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
          elMatDiag_2 +
          _data_invDiag_[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
    }
  }
  {
    {
      _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
          tmp_moved_constant_1 +
          _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
      _data_invDiag_[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
          tmp_moved_constant_2 +
          _data_invDiag_[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
      _data_invDiag_[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1] =
          tmp_moved_constant_3 +
          _data_invDiag_[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
    }
  }
}
} // namespace operatorgeneration
{
  {
    {
      _data_invDiag_[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                     micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) - 1] =
          elMatDiag_0 +
          _data_invDiag_[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                         micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) - 1];
      _data_invDiag_[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                     micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2))] =
          elMatDiag_1 +
          _data_invDiag_[ctr_1 * (micro_edges_per_macro_edge + 2) - ctr_1 +
                         micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2))];
      _data_invDiag_[-ctr_1 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) - 1] =
          elMatDiag_2 +
          _data_invDiag_[-ctr_1 + micro_edges_per_macro_edge +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) - 1];
    }
  }
}
} // namespace hyteg
}
}
void P1ElementwiseDiffusion_cubes_const_vect_float64::
    computeInverseDiagonalOperatorValues_macro_3D(
        walberla::float64 *RESTRICT _data_invDiag_,
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
    const walberla::float64 tmp_coords_jac_0_GREEN_DOWN =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_GREEN_DOWN =
        -p_affine_const_0_0_GREEN_DOWN + p_affine_const_1_0_GREEN_DOWN;
    const walberla::float64 jac_affine_0_1_GREEN_DOWN =
        -p_affine_const_0_0_GREEN_DOWN + p_affine_const_2_0_GREEN_DOWN;
    const walberla::float64 jac_affine_0_2_GREEN_DOWN =
        -p_affine_const_0_0_GREEN_DOWN + p_affine_const_3_0_GREEN_DOWN;
    const walberla::float64 jac_affine_1_0_GREEN_DOWN =
        -p_affine_const_0_1_GREEN_DOWN + p_affine_const_1_1_GREEN_DOWN;
    const walberla::float64 jac_affine_1_1_GREEN_DOWN =
        -p_affine_const_0_1_GREEN_DOWN + p_affine_const_2_1_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_14_GREEN_DOWN =
        jac_affine_0_2_GREEN_DOWN * jac_affine_1_1_GREEN_DOWN;
    const walberla::float64 jac_affine_1_2_GREEN_DOWN =
        -p_affine_const_0_1_GREEN_DOWN + p_affine_const_3_1_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_12_GREEN_DOWN =
        jac_affine_0_1_GREEN_DOWN * jac_affine_1_2_GREEN_DOWN;
    const walberla::float64 jac_affine_2_0_GREEN_DOWN =
        -p_affine_const_0_2_GREEN_DOWN + p_affine_const_1_2_GREEN_DOWN;
    const walberla::float64 jac_affine_2_1_GREEN_DOWN =
        -p_affine_const_0_2_GREEN_DOWN + p_affine_const_2_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_11_GREEN_DOWN =
        jac_affine_1_2_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN;
    const walberla::float64 jac_affine_2_2_GREEN_DOWN =
        -p_affine_const_0_2_GREEN_DOWN + p_affine_const_3_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_10_GREEN_DOWN =
        jac_affine_1_1_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_13_GREEN_DOWN =
        jac_affine_0_1_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_15_GREEN_DOWN =
        jac_affine_0_0_GREEN_DOWN * tmp_coords_jac_10_GREEN_DOWN -
        jac_affine_0_0_GREEN_DOWN * tmp_coords_jac_11_GREEN_DOWN +
        jac_affine_0_2_GREEN_DOWN * jac_affine_1_0_GREEN_DOWN *
            jac_affine_2_1_GREEN_DOWN -
        jac_affine_1_0_GREEN_DOWN * tmp_coords_jac_13_GREEN_DOWN +
        jac_affine_2_0_GREEN_DOWN * tmp_coords_jac_12_GREEN_DOWN -
        jac_affine_2_0_GREEN_DOWN * tmp_coords_jac_14_GREEN_DOWN;
    const walberla::float64 tmp_coords_jac_16_GREEN_DOWN =
        1.0 / (tmp_coords_jac_15_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_0_0_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (tmp_coords_jac_10_GREEN_DOWN - tmp_coords_jac_11_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_0_1_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_0_2_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN -
         tmp_coords_jac_13_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_0_2_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (tmp_coords_jac_12_GREEN_DOWN - tmp_coords_jac_14_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_1_0_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (-jac_affine_1_0_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN +
         jac_affine_1_2_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_1_1_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_0_0_GREEN_DOWN * jac_affine_2_2_GREEN_DOWN -
         jac_affine_0_2_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_1_2_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (-jac_affine_0_0_GREEN_DOWN * jac_affine_1_2_GREEN_DOWN +
         jac_affine_0_2_GREEN_DOWN * jac_affine_1_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_2_0_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_1_0_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN -
         jac_affine_1_1_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_2_1_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (-jac_affine_0_0_GREEN_DOWN * jac_affine_2_1_GREEN_DOWN +
         jac_affine_0_1_GREEN_DOWN * jac_affine_2_0_GREEN_DOWN);
    const walberla::float64 jac_affine_inv_2_2_GREEN_DOWN =
        tmp_coords_jac_16_GREEN_DOWN *
        (jac_affine_0_0_GREEN_DOWN * jac_affine_1_1_GREEN_DOWN -
         jac_affine_0_1_GREEN_DOWN * jac_affine_1_0_GREEN_DOWN);
    const walberla::float64 abs_det_jac_affine_GREEN_DOWN =
        abs(tmp_coords_jac_15_GREEN_DOWN);
    const walberla::float64 tmp_coords_jac_0_GREEN_UP =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_GREEN_UP =
        -p_affine_const_0_0_GREEN_UP + p_affine_const_1_0_GREEN_UP;
    const walberla::float64 jac_affine_0_1_GREEN_UP =
        -p_affine_const_0_0_GREEN_UP + p_affine_const_2_0_GREEN_UP;
    const walberla::float64 jac_affine_0_2_GREEN_UP =
        -p_affine_const_0_0_GREEN_UP + p_affine_const_3_0_GREEN_UP;
    const walberla::float64 jac_affine_1_0_GREEN_UP =
        -p_affine_const_0_1_GREEN_UP + p_affine_const_1_1_GREEN_UP;
    const walberla::float64 jac_affine_1_1_GREEN_UP =
        -p_affine_const_0_1_GREEN_UP + p_affine_const_2_1_GREEN_UP;
    const walberla::float64 tmp_coords_jac_11_GREEN_UP =
        jac_affine_0_2_GREEN_UP * jac_affine_1_1_GREEN_UP;
    const walberla::float64 jac_affine_1_2_GREEN_UP =
        -p_affine_const_0_1_GREEN_UP + p_affine_const_3_1_GREEN_UP;
    const walberla::float64 tmp_coords_jac_9_GREEN_UP =
        jac_affine_0_1_GREEN_UP * jac_affine_1_2_GREEN_UP;
    const walberla::float64 jac_affine_2_0_GREEN_UP =
        -p_affine_const_0_2_GREEN_UP + p_affine_const_1_2_GREEN_UP;
    const walberla::float64 jac_affine_2_1_GREEN_UP =
        -p_affine_const_0_2_GREEN_UP + p_affine_const_2_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_8_GREEN_UP =
        jac_affine_1_2_GREEN_UP * jac_affine_2_1_GREEN_UP;
    const walberla::float64 jac_affine_2_2_GREEN_UP =
        -p_affine_const_0_2_GREEN_UP + p_affine_const_3_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_7_GREEN_UP =
        jac_affine_1_1_GREEN_UP * jac_affine_2_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_10_GREEN_UP =
        jac_affine_0_1_GREEN_UP * jac_affine_2_2_GREEN_UP;
    const walberla::float64 tmp_coords_jac_12_GREEN_UP =
        jac_affine_0_0_GREEN_UP * tmp_coords_jac_7_GREEN_UP -
        jac_affine_0_0_GREEN_UP * tmp_coords_jac_8_GREEN_UP +
        jac_affine_0_2_GREEN_UP * jac_affine_1_0_GREEN_UP *
            jac_affine_2_1_GREEN_UP -
        jac_affine_1_0_GREEN_UP * tmp_coords_jac_10_GREEN_UP -
        jac_affine_2_0_GREEN_UP * tmp_coords_jac_11_GREEN_UP +
        jac_affine_2_0_GREEN_UP * tmp_coords_jac_9_GREEN_UP;
    const walberla::float64 tmp_coords_jac_13_GREEN_UP =
        1.0 / (tmp_coords_jac_12_GREEN_UP);
    const walberla::float64 jac_affine_inv_0_0_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (tmp_coords_jac_7_GREEN_UP - tmp_coords_jac_8_GREEN_UP);
    const walberla::float64 jac_affine_inv_0_1_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_0_2_GREEN_UP * jac_affine_2_1_GREEN_UP -
         tmp_coords_jac_10_GREEN_UP);
    const walberla::float64 jac_affine_inv_0_2_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-tmp_coords_jac_11_GREEN_UP + tmp_coords_jac_9_GREEN_UP);
    const walberla::float64 jac_affine_inv_1_0_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-jac_affine_1_0_GREEN_UP * jac_affine_2_2_GREEN_UP +
         jac_affine_1_2_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_1_1_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_0_0_GREEN_UP * jac_affine_2_2_GREEN_UP -
         jac_affine_0_2_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_1_2_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-jac_affine_0_0_GREEN_UP * jac_affine_1_2_GREEN_UP +
         jac_affine_0_2_GREEN_UP * jac_affine_1_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_2_0_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_1_0_GREEN_UP * jac_affine_2_1_GREEN_UP -
         jac_affine_1_1_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_2_1_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (-jac_affine_0_0_GREEN_UP * jac_affine_2_1_GREEN_UP +
         jac_affine_0_1_GREEN_UP * jac_affine_2_0_GREEN_UP);
    const walberla::float64 jac_affine_inv_2_2_GREEN_UP =
        tmp_coords_jac_13_GREEN_UP *
        (jac_affine_0_0_GREEN_UP * jac_affine_1_1_GREEN_UP -
         jac_affine_0_1_GREEN_UP * jac_affine_1_0_GREEN_UP);
    const walberla::float64 abs_det_jac_affine_GREEN_UP =
        abs(tmp_coords_jac_12_GREEN_UP);
    const walberla::float64 tmp_coords_jac_0_BLUE_DOWN =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_BLUE_DOWN =
        -p_affine_const_0_0_BLUE_DOWN + p_affine_const_1_0_BLUE_DOWN;
    const walberla::float64 jac_affine_0_1_BLUE_DOWN =
        -p_affine_const_0_0_BLUE_DOWN + p_affine_const_2_0_BLUE_DOWN;
    const walberla::float64 jac_affine_0_2_BLUE_DOWN =
        -p_affine_const_0_0_BLUE_DOWN + p_affine_const_3_0_BLUE_DOWN;
    const walberla::float64 jac_affine_1_0_BLUE_DOWN =
        -p_affine_const_0_1_BLUE_DOWN + p_affine_const_1_1_BLUE_DOWN;
    const walberla::float64 jac_affine_1_1_BLUE_DOWN =
        -p_affine_const_0_1_BLUE_DOWN + p_affine_const_2_1_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_14_BLUE_DOWN =
        jac_affine_0_2_BLUE_DOWN * jac_affine_1_1_BLUE_DOWN;
    const walberla::float64 jac_affine_1_2_BLUE_DOWN =
        -p_affine_const_0_1_BLUE_DOWN + p_affine_const_3_1_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_12_BLUE_DOWN =
        jac_affine_0_1_BLUE_DOWN * jac_affine_1_2_BLUE_DOWN;
    const walberla::float64 jac_affine_2_0_BLUE_DOWN =
        -p_affine_const_0_2_BLUE_DOWN + p_affine_const_1_2_BLUE_DOWN;
    const walberla::float64 jac_affine_2_1_BLUE_DOWN =
        -p_affine_const_0_2_BLUE_DOWN + p_affine_const_2_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_11_BLUE_DOWN =
        jac_affine_1_2_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN;
    const walberla::float64 jac_affine_2_2_BLUE_DOWN =
        -p_affine_const_0_2_BLUE_DOWN + p_affine_const_3_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_10_BLUE_DOWN =
        jac_affine_1_1_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_13_BLUE_DOWN =
        jac_affine_0_1_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_15_BLUE_DOWN =
        jac_affine_0_0_BLUE_DOWN * tmp_coords_jac_10_BLUE_DOWN -
        jac_affine_0_0_BLUE_DOWN * tmp_coords_jac_11_BLUE_DOWN +
        jac_affine_0_2_BLUE_DOWN * jac_affine_1_0_BLUE_DOWN *
            jac_affine_2_1_BLUE_DOWN -
        jac_affine_1_0_BLUE_DOWN * tmp_coords_jac_13_BLUE_DOWN +
        jac_affine_2_0_BLUE_DOWN * tmp_coords_jac_12_BLUE_DOWN -
        jac_affine_2_0_BLUE_DOWN * tmp_coords_jac_14_BLUE_DOWN;
    const walberla::float64 tmp_coords_jac_16_BLUE_DOWN =
        1.0 / (tmp_coords_jac_15_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_0_0_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (tmp_coords_jac_10_BLUE_DOWN - tmp_coords_jac_11_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_0_1_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_0_2_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN -
         tmp_coords_jac_13_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_0_2_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (tmp_coords_jac_12_BLUE_DOWN - tmp_coords_jac_14_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_1_0_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (-jac_affine_1_0_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN +
         jac_affine_1_2_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_1_1_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_0_0_BLUE_DOWN * jac_affine_2_2_BLUE_DOWN -
         jac_affine_0_2_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_1_2_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (-jac_affine_0_0_BLUE_DOWN * jac_affine_1_2_BLUE_DOWN +
         jac_affine_0_2_BLUE_DOWN * jac_affine_1_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_2_0_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_1_0_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN -
         jac_affine_1_1_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_2_1_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (-jac_affine_0_0_BLUE_DOWN * jac_affine_2_1_BLUE_DOWN +
         jac_affine_0_1_BLUE_DOWN * jac_affine_2_0_BLUE_DOWN);
    const walberla::float64 jac_affine_inv_2_2_BLUE_DOWN =
        tmp_coords_jac_16_BLUE_DOWN *
        (jac_affine_0_0_BLUE_DOWN * jac_affine_1_1_BLUE_DOWN -
         jac_affine_0_1_BLUE_DOWN * jac_affine_1_0_BLUE_DOWN);
    const walberla::float64 abs_det_jac_affine_BLUE_DOWN =
        abs(tmp_coords_jac_15_BLUE_DOWN);
    const walberla::float64 tmp_coords_jac_0_BLUE_UP =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_BLUE_UP =
        -p_affine_const_0_0_BLUE_UP + p_affine_const_1_0_BLUE_UP;
    const walberla::float64 jac_affine_0_1_BLUE_UP =
        -p_affine_const_0_0_BLUE_UP + p_affine_const_2_0_BLUE_UP;
    const walberla::float64 jac_affine_0_2_BLUE_UP =
        -p_affine_const_0_0_BLUE_UP + p_affine_const_3_0_BLUE_UP;
    const walberla::float64 jac_affine_1_0_BLUE_UP =
        -p_affine_const_0_1_BLUE_UP + p_affine_const_1_1_BLUE_UP;
    const walberla::float64 jac_affine_1_1_BLUE_UP =
        -p_affine_const_0_1_BLUE_UP + p_affine_const_2_1_BLUE_UP;
    const walberla::float64 tmp_coords_jac_11_BLUE_UP =
        jac_affine_0_2_BLUE_UP * jac_affine_1_1_BLUE_UP;
    const walberla::float64 jac_affine_1_2_BLUE_UP =
        -p_affine_const_0_1_BLUE_UP + p_affine_const_3_1_BLUE_UP;
    const walberla::float64 tmp_coords_jac_9_BLUE_UP =
        jac_affine_0_1_BLUE_UP * jac_affine_1_2_BLUE_UP;
    const walberla::float64 jac_affine_2_0_BLUE_UP =
        -p_affine_const_0_2_BLUE_UP + p_affine_const_1_2_BLUE_UP;
    const walberla::float64 jac_affine_2_1_BLUE_UP =
        -p_affine_const_0_2_BLUE_UP + p_affine_const_2_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_8_BLUE_UP =
        jac_affine_1_2_BLUE_UP * jac_affine_2_1_BLUE_UP;
    const walberla::float64 jac_affine_2_2_BLUE_UP =
        -p_affine_const_0_2_BLUE_UP + p_affine_const_3_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_7_BLUE_UP =
        jac_affine_1_1_BLUE_UP * jac_affine_2_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_10_BLUE_UP =
        jac_affine_0_1_BLUE_UP * jac_affine_2_2_BLUE_UP;
    const walberla::float64 tmp_coords_jac_12_BLUE_UP =
        jac_affine_0_0_BLUE_UP * tmp_coords_jac_7_BLUE_UP -
        jac_affine_0_0_BLUE_UP * tmp_coords_jac_8_BLUE_UP +
        jac_affine_0_2_BLUE_UP * jac_affine_1_0_BLUE_UP *
            jac_affine_2_1_BLUE_UP -
        jac_affine_1_0_BLUE_UP * tmp_coords_jac_10_BLUE_UP -
        jac_affine_2_0_BLUE_UP * tmp_coords_jac_11_BLUE_UP +
        jac_affine_2_0_BLUE_UP * tmp_coords_jac_9_BLUE_UP;
    const walberla::float64 tmp_coords_jac_13_BLUE_UP =
        1.0 / (tmp_coords_jac_12_BLUE_UP);
    const walberla::float64 jac_affine_inv_0_0_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (tmp_coords_jac_7_BLUE_UP - tmp_coords_jac_8_BLUE_UP);
    const walberla::float64 jac_affine_inv_0_1_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_0_2_BLUE_UP * jac_affine_2_1_BLUE_UP -
         tmp_coords_jac_10_BLUE_UP);
    const walberla::float64 jac_affine_inv_0_2_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-tmp_coords_jac_11_BLUE_UP + tmp_coords_jac_9_BLUE_UP);
    const walberla::float64 jac_affine_inv_1_0_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-jac_affine_1_0_BLUE_UP * jac_affine_2_2_BLUE_UP +
         jac_affine_1_2_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_1_1_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_0_0_BLUE_UP * jac_affine_2_2_BLUE_UP -
         jac_affine_0_2_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_1_2_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-jac_affine_0_0_BLUE_UP * jac_affine_1_2_BLUE_UP +
         jac_affine_0_2_BLUE_UP * jac_affine_1_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_2_0_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_1_0_BLUE_UP * jac_affine_2_1_BLUE_UP -
         jac_affine_1_1_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_2_1_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (-jac_affine_0_0_BLUE_UP * jac_affine_2_1_BLUE_UP +
         jac_affine_0_1_BLUE_UP * jac_affine_2_0_BLUE_UP);
    const walberla::float64 jac_affine_inv_2_2_BLUE_UP =
        tmp_coords_jac_13_BLUE_UP *
        (jac_affine_0_0_BLUE_UP * jac_affine_1_1_BLUE_UP -
         jac_affine_0_1_BLUE_UP * jac_affine_1_0_BLUE_UP);
    const walberla::float64 abs_det_jac_affine_BLUE_UP =
        abs(tmp_coords_jac_12_BLUE_UP);
    const walberla::float64 tmp_coords_jac_0_WHITE_DOWN =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_WHITE_DOWN =
        -p_affine_const_0_0_WHITE_DOWN + p_affine_const_1_0_WHITE_DOWN;
    const walberla::float64 jac_affine_0_1_WHITE_DOWN =
        -p_affine_const_0_0_WHITE_DOWN + p_affine_const_2_0_WHITE_DOWN;
    const walberla::float64 jac_affine_0_2_WHITE_DOWN =
        -p_affine_const_0_0_WHITE_DOWN + p_affine_const_3_0_WHITE_DOWN;
    const walberla::float64 jac_affine_1_0_WHITE_DOWN =
        -p_affine_const_0_1_WHITE_DOWN + p_affine_const_1_1_WHITE_DOWN;
    const walberla::float64 jac_affine_1_1_WHITE_DOWN =
        -p_affine_const_0_1_WHITE_DOWN + p_affine_const_2_1_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_17_WHITE_DOWN =
        jac_affine_0_2_WHITE_DOWN * jac_affine_1_1_WHITE_DOWN;
    const walberla::float64 jac_affine_1_2_WHITE_DOWN =
        -p_affine_const_0_1_WHITE_DOWN + p_affine_const_3_1_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_15_WHITE_DOWN =
        jac_affine_0_1_WHITE_DOWN * jac_affine_1_2_WHITE_DOWN;
    const walberla::float64 jac_affine_2_0_WHITE_DOWN =
        -p_affine_const_0_2_WHITE_DOWN + p_affine_const_1_2_WHITE_DOWN;
    const walberla::float64 jac_affine_2_1_WHITE_DOWN =
        -p_affine_const_0_2_WHITE_DOWN + p_affine_const_2_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_14_WHITE_DOWN =
        jac_affine_1_2_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN;
    const walberla::float64 jac_affine_2_2_WHITE_DOWN =
        -p_affine_const_0_2_WHITE_DOWN + p_affine_const_3_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_13_WHITE_DOWN =
        jac_affine_1_1_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_16_WHITE_DOWN =
        jac_affine_0_1_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_18_WHITE_DOWN =
        jac_affine_0_0_WHITE_DOWN * tmp_coords_jac_13_WHITE_DOWN -
        jac_affine_0_0_WHITE_DOWN * tmp_coords_jac_14_WHITE_DOWN +
        jac_affine_0_2_WHITE_DOWN * jac_affine_1_0_WHITE_DOWN *
            jac_affine_2_1_WHITE_DOWN -
        jac_affine_1_0_WHITE_DOWN * tmp_coords_jac_16_WHITE_DOWN +
        jac_affine_2_0_WHITE_DOWN * tmp_coords_jac_15_WHITE_DOWN -
        jac_affine_2_0_WHITE_DOWN * tmp_coords_jac_17_WHITE_DOWN;
    const walberla::float64 tmp_coords_jac_19_WHITE_DOWN =
        1.0 / (tmp_coords_jac_18_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_0_0_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (tmp_coords_jac_13_WHITE_DOWN - tmp_coords_jac_14_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_0_1_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_0_2_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN -
         tmp_coords_jac_16_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_0_2_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (tmp_coords_jac_15_WHITE_DOWN - tmp_coords_jac_17_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_1_0_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (-jac_affine_1_0_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN +
         jac_affine_1_2_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_1_1_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_0_0_WHITE_DOWN * jac_affine_2_2_WHITE_DOWN -
         jac_affine_0_2_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_1_2_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (-jac_affine_0_0_WHITE_DOWN * jac_affine_1_2_WHITE_DOWN +
         jac_affine_0_2_WHITE_DOWN * jac_affine_1_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_2_0_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_1_0_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN -
         jac_affine_1_1_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_2_1_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (-jac_affine_0_0_WHITE_DOWN * jac_affine_2_1_WHITE_DOWN +
         jac_affine_0_1_WHITE_DOWN * jac_affine_2_0_WHITE_DOWN);
    const walberla::float64 jac_affine_inv_2_2_WHITE_DOWN =
        tmp_coords_jac_19_WHITE_DOWN *
        (jac_affine_0_0_WHITE_DOWN * jac_affine_1_1_WHITE_DOWN -
         jac_affine_0_1_WHITE_DOWN * jac_affine_1_0_WHITE_DOWN);
    const walberla::float64 abs_det_jac_affine_WHITE_DOWN =
        abs(tmp_coords_jac_18_WHITE_DOWN);
    const walberla::float64 tmp_coords_jac_0_WHITE_UP =
        1.0 / (micro_edges_per_macro_edge_float)*1.0;
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
    const walberla::float64 jac_affine_0_0_WHITE_UP =
        -p_affine_const_0_0_WHITE_UP + p_affine_const_1_0_WHITE_UP;
    const walberla::float64 jac_affine_0_1_WHITE_UP =
        -p_affine_const_0_0_WHITE_UP + p_affine_const_2_0_WHITE_UP;
    const walberla::float64 jac_affine_0_2_WHITE_UP =
        -p_affine_const_0_0_WHITE_UP + p_affine_const_3_0_WHITE_UP;
    const walberla::float64 jac_affine_1_0_WHITE_UP =
        -p_affine_const_0_1_WHITE_UP + p_affine_const_1_1_WHITE_UP;
    const walberla::float64 jac_affine_1_1_WHITE_UP =
        -p_affine_const_0_1_WHITE_UP + p_affine_const_2_1_WHITE_UP;
    const walberla::float64 tmp_coords_jac_5_WHITE_UP =
        jac_affine_0_2_WHITE_UP * jac_affine_1_1_WHITE_UP;
    const walberla::float64 jac_affine_1_2_WHITE_UP =
        -p_affine_const_0_1_WHITE_UP + p_affine_const_3_1_WHITE_UP;
    const walberla::float64 tmp_coords_jac_3_WHITE_UP =
        jac_affine_0_1_WHITE_UP * jac_affine_1_2_WHITE_UP;
    const walberla::float64 jac_affine_2_0_WHITE_UP =
        -p_affine_const_0_2_WHITE_UP + p_affine_const_1_2_WHITE_UP;
    const walberla::float64 jac_affine_2_1_WHITE_UP =
        -p_affine_const_0_2_WHITE_UP + p_affine_const_2_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_2_WHITE_UP =
        jac_affine_1_2_WHITE_UP * jac_affine_2_1_WHITE_UP;
    const walberla::float64 jac_affine_2_2_WHITE_UP =
        -p_affine_const_0_2_WHITE_UP + p_affine_const_3_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_1_WHITE_UP =
        jac_affine_1_1_WHITE_UP * jac_affine_2_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_4_WHITE_UP =
        jac_affine_0_1_WHITE_UP * jac_affine_2_2_WHITE_UP;
    const walberla::float64 tmp_coords_jac_6_WHITE_UP =
        jac_affine_0_0_WHITE_UP * tmp_coords_jac_1_WHITE_UP -
        jac_affine_0_0_WHITE_UP * tmp_coords_jac_2_WHITE_UP +
        jac_affine_0_2_WHITE_UP * jac_affine_1_0_WHITE_UP *
            jac_affine_2_1_WHITE_UP -
        jac_affine_1_0_WHITE_UP * tmp_coords_jac_4_WHITE_UP +
        jac_affine_2_0_WHITE_UP * tmp_coords_jac_3_WHITE_UP -
        jac_affine_2_0_WHITE_UP * tmp_coords_jac_5_WHITE_UP;
    const walberla::float64 tmp_coords_jac_7_WHITE_UP =
        1.0 / (tmp_coords_jac_6_WHITE_UP);
    const walberla::float64 jac_affine_inv_0_0_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (tmp_coords_jac_1_WHITE_UP - tmp_coords_jac_2_WHITE_UP);
    const walberla::float64 jac_affine_inv_0_1_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_0_2_WHITE_UP * jac_affine_2_1_WHITE_UP -
         tmp_coords_jac_4_WHITE_UP);
    const walberla::float64 jac_affine_inv_0_2_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (tmp_coords_jac_3_WHITE_UP - tmp_coords_jac_5_WHITE_UP);
    const walberla::float64 jac_affine_inv_1_0_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (-jac_affine_1_0_WHITE_UP * jac_affine_2_2_WHITE_UP +
         jac_affine_1_2_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_1_1_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_0_0_WHITE_UP * jac_affine_2_2_WHITE_UP -
         jac_affine_0_2_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_1_2_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (-jac_affine_0_0_WHITE_UP * jac_affine_1_2_WHITE_UP +
         jac_affine_0_2_WHITE_UP * jac_affine_1_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_2_0_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_1_0_WHITE_UP * jac_affine_2_1_WHITE_UP -
         jac_affine_1_1_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_2_1_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (-jac_affine_0_0_WHITE_UP * jac_affine_2_1_WHITE_UP +
         jac_affine_0_1_WHITE_UP * jac_affine_2_0_WHITE_UP);
    const walberla::float64 jac_affine_inv_2_2_WHITE_UP =
        tmp_coords_jac_7_WHITE_UP *
        (jac_affine_0_0_WHITE_UP * jac_affine_1_1_WHITE_UP -
         jac_affine_0_1_WHITE_UP * jac_affine_1_0_WHITE_UP);
    const walberla::float64 abs_det_jac_affine_WHITE_UP =
        abs(tmp_coords_jac_6_WHITE_UP);
    const walberla::float64 tmp_kernel_op_0 =
        abs_det_jac_affine_WHITE_UP * 0.16666666666666663;
    const walberla::float64 elMatDiag_0 =
        tmp_kernel_op_0 *
        (((-jac_affine_inv_0_0_WHITE_UP - jac_affine_inv_1_0_WHITE_UP -
           jac_affine_inv_2_0_WHITE_UP) *
          (-jac_affine_inv_0_0_WHITE_UP - jac_affine_inv_1_0_WHITE_UP -
           jac_affine_inv_2_0_WHITE_UP)) +
         ((-jac_affine_inv_0_1_WHITE_UP - jac_affine_inv_1_1_WHITE_UP -
           jac_affine_inv_2_1_WHITE_UP) *
          (-jac_affine_inv_0_1_WHITE_UP - jac_affine_inv_1_1_WHITE_UP -
           jac_affine_inv_2_1_WHITE_UP)) +
         ((-jac_affine_inv_0_2_WHITE_UP - jac_affine_inv_1_2_WHITE_UP -
           jac_affine_inv_2_2_WHITE_UP) *
          (-jac_affine_inv_0_2_WHITE_UP - jac_affine_inv_1_2_WHITE_UP -
           jac_affine_inv_2_2_WHITE_UP)));
    const walberla::float64 elMatDiag_1 =
        tmp_kernel_op_0 *
        ((jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_0_0_WHITE_UP) +
         (jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_0_1_WHITE_UP) +
         (jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_0_2_WHITE_UP));
    const walberla::float64 elMatDiag_2 =
        tmp_kernel_op_0 *
        ((jac_affine_inv_1_0_WHITE_UP * jac_affine_inv_1_0_WHITE_UP) +
         (jac_affine_inv_1_1_WHITE_UP * jac_affine_inv_1_1_WHITE_UP) +
         (jac_affine_inv_1_2_WHITE_UP * jac_affine_inv_1_2_WHITE_UP));
    const walberla::float64 elMatDiag_3 =
        tmp_kernel_op_0 *
        ((jac_affine_inv_2_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP) +
         (jac_affine_inv_2_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP) +
         (jac_affine_inv_2_2_WHITE_UP * jac_affine_inv_2_2_WHITE_UP));
    const walberla::float64 tmp_moved_constant_0 =
        abs_det_jac_affine_WHITE_DOWN * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_1 =
        tmp_moved_constant_0 *
        (((-jac_affine_inv_0_0_WHITE_DOWN - jac_affine_inv_1_0_WHITE_DOWN -
           jac_affine_inv_2_0_WHITE_DOWN) *
          (-jac_affine_inv_0_0_WHITE_DOWN - jac_affine_inv_1_0_WHITE_DOWN -
           jac_affine_inv_2_0_WHITE_DOWN)) +
         ((-jac_affine_inv_0_1_WHITE_DOWN - jac_affine_inv_1_1_WHITE_DOWN -
           jac_affine_inv_2_1_WHITE_DOWN) *
          (-jac_affine_inv_0_1_WHITE_DOWN - jac_affine_inv_1_1_WHITE_DOWN -
           jac_affine_inv_2_1_WHITE_DOWN)) +
         ((-jac_affine_inv_0_2_WHITE_DOWN - jac_affine_inv_1_2_WHITE_DOWN -
           jac_affine_inv_2_2_WHITE_DOWN) *
          (-jac_affine_inv_0_2_WHITE_DOWN - jac_affine_inv_1_2_WHITE_DOWN -
           jac_affine_inv_2_2_WHITE_DOWN)));
    const walberla::float64 tmp_moved_constant_2 =
        tmp_moved_constant_0 *
        ((jac_affine_inv_0_0_WHITE_DOWN * jac_affine_inv_0_0_WHITE_DOWN) +
         (jac_affine_inv_0_1_WHITE_DOWN * jac_affine_inv_0_1_WHITE_DOWN) +
         (jac_affine_inv_0_2_WHITE_DOWN * jac_affine_inv_0_2_WHITE_DOWN));
    const walberla::float64 tmp_moved_constant_3 =
        tmp_moved_constant_0 *
        ((jac_affine_inv_1_0_WHITE_DOWN * jac_affine_inv_1_0_WHITE_DOWN) +
         (jac_affine_inv_1_1_WHITE_DOWN * jac_affine_inv_1_1_WHITE_DOWN) +
         (jac_affine_inv_1_2_WHITE_DOWN * jac_affine_inv_1_2_WHITE_DOWN));
    const walberla::float64 tmp_moved_constant_4 =
        tmp_moved_constant_0 *
        ((jac_affine_inv_2_0_WHITE_DOWN * jac_affine_inv_2_0_WHITE_DOWN) +
         (jac_affine_inv_2_1_WHITE_DOWN * jac_affine_inv_2_1_WHITE_DOWN) +
         (jac_affine_inv_2_2_WHITE_DOWN * jac_affine_inv_2_2_WHITE_DOWN));
    const walberla::float64 tmp_moved_constant_5 =
        abs_det_jac_affine_BLUE_UP * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_6 =
        tmp_moved_constant_5 *
        (((-jac_affine_inv_0_0_BLUE_UP - jac_affine_inv_1_0_BLUE_UP -
           jac_affine_inv_2_0_BLUE_UP) *
          (-jac_affine_inv_0_0_BLUE_UP - jac_affine_inv_1_0_BLUE_UP -
           jac_affine_inv_2_0_BLUE_UP)) +
         ((-jac_affine_inv_0_1_BLUE_UP - jac_affine_inv_1_1_BLUE_UP -
           jac_affine_inv_2_1_BLUE_UP) *
          (-jac_affine_inv_0_1_BLUE_UP - jac_affine_inv_1_1_BLUE_UP -
           jac_affine_inv_2_1_BLUE_UP)) +
         ((-jac_affine_inv_0_2_BLUE_UP - jac_affine_inv_1_2_BLUE_UP -
           jac_affine_inv_2_2_BLUE_UP) *
          (-jac_affine_inv_0_2_BLUE_UP - jac_affine_inv_1_2_BLUE_UP -
           jac_affine_inv_2_2_BLUE_UP)));
    const walberla::float64 tmp_moved_constant_7 =
        tmp_moved_constant_5 *
        ((jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_0_0_BLUE_UP) +
         (jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_0_1_BLUE_UP) +
         (jac_affine_inv_0_2_BLUE_UP * jac_affine_inv_0_2_BLUE_UP));
    const walberla::float64 tmp_moved_constant_8 =
        tmp_moved_constant_5 *
        ((jac_affine_inv_1_0_BLUE_UP * jac_affine_inv_1_0_BLUE_UP) +
         (jac_affine_inv_1_1_BLUE_UP * jac_affine_inv_1_1_BLUE_UP) +
         (jac_affine_inv_1_2_BLUE_UP * jac_affine_inv_1_2_BLUE_UP));
    const walberla::float64 tmp_moved_constant_9 =
        tmp_moved_constant_5 *
        ((jac_affine_inv_2_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP) +
         (jac_affine_inv_2_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP) +
         (jac_affine_inv_2_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP));
    const walberla::float64 tmp_moved_constant_10 =
        abs_det_jac_affine_BLUE_DOWN * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_11 =
        tmp_moved_constant_10 *
        (((-jac_affine_inv_0_0_BLUE_DOWN - jac_affine_inv_1_0_BLUE_DOWN -
           jac_affine_inv_2_0_BLUE_DOWN) *
          (-jac_affine_inv_0_0_BLUE_DOWN - jac_affine_inv_1_0_BLUE_DOWN -
           jac_affine_inv_2_0_BLUE_DOWN)) +
         ((-jac_affine_inv_0_1_BLUE_DOWN - jac_affine_inv_1_1_BLUE_DOWN -
           jac_affine_inv_2_1_BLUE_DOWN) *
          (-jac_affine_inv_0_1_BLUE_DOWN - jac_affine_inv_1_1_BLUE_DOWN -
           jac_affine_inv_2_1_BLUE_DOWN)) +
         ((-jac_affine_inv_0_2_BLUE_DOWN - jac_affine_inv_1_2_BLUE_DOWN -
           jac_affine_inv_2_2_BLUE_DOWN) *
          (-jac_affine_inv_0_2_BLUE_DOWN - jac_affine_inv_1_2_BLUE_DOWN -
           jac_affine_inv_2_2_BLUE_DOWN)));
    const walberla::float64 tmp_moved_constant_12 =
        tmp_moved_constant_10 *
        ((jac_affine_inv_0_0_BLUE_DOWN * jac_affine_inv_0_0_BLUE_DOWN) +
         (jac_affine_inv_0_1_BLUE_DOWN * jac_affine_inv_0_1_BLUE_DOWN) +
         (jac_affine_inv_0_2_BLUE_DOWN * jac_affine_inv_0_2_BLUE_DOWN));
    const walberla::float64 tmp_moved_constant_13 =
        tmp_moved_constant_10 *
        ((jac_affine_inv_1_0_BLUE_DOWN * jac_affine_inv_1_0_BLUE_DOWN) +
         (jac_affine_inv_1_1_BLUE_DOWN * jac_affine_inv_1_1_BLUE_DOWN) +
         (jac_affine_inv_1_2_BLUE_DOWN * jac_affine_inv_1_2_BLUE_DOWN));
    const walberla::float64 tmp_moved_constant_14 =
        tmp_moved_constant_10 *
        ((jac_affine_inv_2_0_BLUE_DOWN * jac_affine_inv_2_0_BLUE_DOWN) +
         (jac_affine_inv_2_1_BLUE_DOWN * jac_affine_inv_2_1_BLUE_DOWN) +
         (jac_affine_inv_2_2_BLUE_DOWN * jac_affine_inv_2_2_BLUE_DOWN));
    const walberla::float64 tmp_moved_constant_15 =
        abs_det_jac_affine_GREEN_UP * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_16 =
        tmp_moved_constant_15 *
        (((-jac_affine_inv_0_0_GREEN_UP - jac_affine_inv_1_0_GREEN_UP -
           jac_affine_inv_2_0_GREEN_UP) *
          (-jac_affine_inv_0_0_GREEN_UP - jac_affine_inv_1_0_GREEN_UP -
           jac_affine_inv_2_0_GREEN_UP)) +
         ((-jac_affine_inv_0_1_GREEN_UP - jac_affine_inv_1_1_GREEN_UP -
           jac_affine_inv_2_1_GREEN_UP) *
          (-jac_affine_inv_0_1_GREEN_UP - jac_affine_inv_1_1_GREEN_UP -
           jac_affine_inv_2_1_GREEN_UP)) +
         ((-jac_affine_inv_0_2_GREEN_UP - jac_affine_inv_1_2_GREEN_UP -
           jac_affine_inv_2_2_GREEN_UP) *
          (-jac_affine_inv_0_2_GREEN_UP - jac_affine_inv_1_2_GREEN_UP -
           jac_affine_inv_2_2_GREEN_UP)));
    const walberla::float64 tmp_moved_constant_17 =
        tmp_moved_constant_15 *
        ((jac_affine_inv_0_0_GREEN_UP * jac_affine_inv_0_0_GREEN_UP) +
         (jac_affine_inv_0_1_GREEN_UP * jac_affine_inv_0_1_GREEN_UP) +
         (jac_affine_inv_0_2_GREEN_UP * jac_affine_inv_0_2_GREEN_UP));
    const walberla::float64 tmp_moved_constant_18 =
        tmp_moved_constant_15 *
        ((jac_affine_inv_1_0_GREEN_UP * jac_affine_inv_1_0_GREEN_UP) +
         (jac_affine_inv_1_1_GREEN_UP * jac_affine_inv_1_1_GREEN_UP) +
         (jac_affine_inv_1_2_GREEN_UP * jac_affine_inv_1_2_GREEN_UP));
    const walberla::float64 tmp_moved_constant_19 =
        tmp_moved_constant_15 *
        ((jac_affine_inv_2_0_GREEN_UP * jac_affine_inv_2_0_GREEN_UP) +
         (jac_affine_inv_2_1_GREEN_UP * jac_affine_inv_2_1_GREEN_UP) +
         (jac_affine_inv_2_2_GREEN_UP * jac_affine_inv_2_2_GREEN_UP));
    const walberla::float64 tmp_moved_constant_20 =
        abs_det_jac_affine_GREEN_DOWN * 0.16666666666666663;
    const walberla::float64 tmp_moved_constant_21 =
        tmp_moved_constant_20 *
        (((-jac_affine_inv_0_0_GREEN_DOWN - jac_affine_inv_1_0_GREEN_DOWN -
           jac_affine_inv_2_0_GREEN_DOWN) *
          (-jac_affine_inv_0_0_GREEN_DOWN - jac_affine_inv_1_0_GREEN_DOWN -
           jac_affine_inv_2_0_GREEN_DOWN)) +
         ((-jac_affine_inv_0_1_GREEN_DOWN - jac_affine_inv_1_1_GREEN_DOWN -
           jac_affine_inv_2_1_GREEN_DOWN) *
          (-jac_affine_inv_0_1_GREEN_DOWN - jac_affine_inv_1_1_GREEN_DOWN -
           jac_affine_inv_2_1_GREEN_DOWN)) +
         ((-jac_affine_inv_0_2_GREEN_DOWN - jac_affine_inv_1_2_GREEN_DOWN -
           jac_affine_inv_2_2_GREEN_DOWN) *
          (-jac_affine_inv_0_2_GREEN_DOWN - jac_affine_inv_1_2_GREEN_DOWN -
           jac_affine_inv_2_2_GREEN_DOWN)));
    const walberla::float64 tmp_moved_constant_22 =
        tmp_moved_constant_20 *
        ((jac_affine_inv_0_0_GREEN_DOWN * jac_affine_inv_0_0_GREEN_DOWN) +
         (jac_affine_inv_0_1_GREEN_DOWN * jac_affine_inv_0_1_GREEN_DOWN) +
         (jac_affine_inv_0_2_GREEN_DOWN * jac_affine_inv_0_2_GREEN_DOWN));
    const walberla::float64 tmp_moved_constant_23 =
        tmp_moved_constant_20 *
        ((jac_affine_inv_1_0_GREEN_DOWN * jac_affine_inv_1_0_GREEN_DOWN) +
         (jac_affine_inv_1_1_GREEN_DOWN * jac_affine_inv_1_1_GREEN_DOWN) +
         (jac_affine_inv_1_2_GREEN_DOWN * jac_affine_inv_1_2_GREEN_DOWN));
    const walberla::float64 tmp_moved_constant_24 =
        tmp_moved_constant_20 *
        ((jac_affine_inv_2_0_GREEN_DOWN * jac_affine_inv_2_0_GREEN_DOWN) +
         (jac_affine_inv_2_1_GREEN_DOWN * jac_affine_inv_2_1_GREEN_DOWN) +
         (jac_affine_inv_2_2_GREEN_DOWN * jac_affine_inv_2_2_GREEN_DOWN));
    for (int64_t ctr_2 = 0; ctr_2 < micro_edges_per_macro_edge; ctr_2 += 1)
      for (int64_t ctr_1 = 0; ctr_1 < -ctr_2 + micro_edges_per_macro_edge;
           ctr_1 += 1) {
        {
          for (int64_t ctr_0 = 0;
               ctr_0 <
               (int64_t)((-ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2) /
                         (4)) *
                   (4);
               ctr_0 += 4) {
            {{_mm256_storeu_pd(
                &_data_invDiag_
                    [ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6))],
                _mm256_add_pd(
                    _mm256_set_pd(elMatDiag_0, elMatDiag_0, elMatDiag_0,
                                  elMatDiag_0),
                    _mm256_loadu_pd(
                        &_data_invDiag_
                            [ctr_0 +
                             ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) +
                             (((micro_edges_per_macro_edge + 1) *
                               (micro_edges_per_macro_edge + 2) *
                               (micro_edges_per_macro_edge + 3)) /
                              (6)) -
                             (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2) *
                               (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                              (6))])));
            _mm256_storeu_pd(
                &_data_invDiag_[ctr_0 +
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
                                1],
                _mm256_add_pd(
                    _mm256_set_pd(elMatDiag_1, elMatDiag_1, elMatDiag_1,
                                  elMatDiag_1),
                    _mm256_loadu_pd(
                        &_data_invDiag_
                            [ctr_0 +
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
                             1])));
            _mm256_storeu_pd(
                &_data_invDiag_[ctr_0 +
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
                                 (6))],
                _mm256_add_pd(
                    _mm256_set_pd(elMatDiag_2, elMatDiag_2, elMatDiag_2,
                                  elMatDiag_2),
                    _mm256_loadu_pd(
                        &_data_invDiag_
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
                              (6))])));
            _mm256_storeu_pd(
                &_data_invDiag_[ctr_0 +
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
                                 (6))],
                _mm256_add_pd(
                    _mm256_set_pd(elMatDiag_3, elMatDiag_3, elMatDiag_3,
                                  elMatDiag_3),
                    _mm256_loadu_pd(
                        &_data_invDiag_
                            [ctr_0 +
                             ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) -
                             (((-ctr_2 + micro_edges_per_macro_edge) *
                               (-ctr_2 + micro_edges_per_macro_edge + 1) *
                               (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                              (6)) +
                             (((micro_edges_per_macro_edge + 1) *
                               (micro_edges_per_macro_edge + 2) *
                               (micro_edges_per_macro_edge + 3)) /
                              (6))])));
          }
        }
        {{_mm256_storeu_pd(
            &_data_invDiag_[ctr_0 +
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
                            1],
            _mm256_add_pd(
                _mm256_set_pd(tmp_moved_constant_1, tmp_moved_constant_1,
                              tmp_moved_constant_1, tmp_moved_constant_1),
                _mm256_loadu_pd(
                    &_data_invDiag_
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
                         1])));
        _mm256_storeu_pd(
            &_data_invDiag_[ctr_0 +
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
                            1],
            _mm256_add_pd(
                _mm256_set_pd(tmp_moved_constant_2, tmp_moved_constant_2,
                              tmp_moved_constant_2, tmp_moved_constant_2),
                _mm256_loadu_pd(
                    &_data_invDiag_
                        [ctr_0 +
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
                         1])));
        _mm256_storeu_pd(
            &_data_invDiag_[ctr_0 +
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
                             (6))],
            _mm256_add_pd(
                _mm256_set_pd(tmp_moved_constant_3, tmp_moved_constant_3,
                              tmp_moved_constant_3, tmp_moved_constant_3),
                _mm256_loadu_pd(
                    &_data_invDiag_
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
                          (6))])));
        _mm256_storeu_pd(
            &_data_invDiag_[ctr_0 +
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
                            1],
            _mm256_add_pd(
                _mm256_set_pd(tmp_moved_constant_4, tmp_moved_constant_4,
                              tmp_moved_constant_4, tmp_moved_constant_4),
                _mm256_loadu_pd(
                    &_data_invDiag_
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
                         1])));
      }
  }
  {{_mm256_storeu_pd(
      &_data_invDiag_[ctr_0 +
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
                      1],
      _mm256_add_pd(
          _mm256_set_pd(tmp_moved_constant_6, tmp_moved_constant_6,
                        tmp_moved_constant_6, tmp_moved_constant_6),
          _mm256_loadu_pd(
              &_data_invDiag_
                  [ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                   ((ctr_1 * (ctr_1 + 1)) / (2)) +
                   (((micro_edges_per_macro_edge + 1) *
                     (micro_edges_per_macro_edge + 2) *
                     (micro_edges_per_macro_edge + 3)) /
                    (6)) -
                   (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                     (-ctr_2 + micro_edges_per_macro_edge + 2) *
                     (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                    (6)) +
                   1])));
  _mm256_storeu_pd(
      &_data_invDiag_[ctr_0 +
                      (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                      (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                      (((micro_edges_per_macro_edge + 1) *
                        (micro_edges_per_macro_edge + 2) *
                        (micro_edges_per_macro_edge + 3)) /
                       (6)) -
                      (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                        (-ctr_2 + micro_edges_per_macro_edge + 2) *
                        (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                       (6))],
      _mm256_add_pd(
          _mm256_set_pd(tmp_moved_constant_7, tmp_moved_constant_7,
                        tmp_moved_constant_7, tmp_moved_constant_7),
          _mm256_loadu_pd(
              &_data_invDiag_[ctr_0 +
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
                               (6))])));
  _mm256_storeu_pd(
      &_data_invDiag_[ctr_0 +
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
                      1],
      _mm256_add_pd(
          _mm256_set_pd(tmp_moved_constant_8, tmp_moved_constant_8,
                        tmp_moved_constant_8, tmp_moved_constant_8),
          _mm256_loadu_pd(
              &_data_invDiag_[ctr_0 +
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
                              1])));
  _mm256_storeu_pd(
      &_data_invDiag_[ctr_0 +
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
                      1],
      _mm256_add_pd(
          _mm256_set_pd(tmp_moved_constant_9, tmp_moved_constant_9,
                        tmp_moved_constant_9, tmp_moved_constant_9),
          _mm256_loadu_pd(
              &_data_invDiag_
                  [ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                   ((ctr_1 * (ctr_1 + 1)) / (2)) -
                   (((-ctr_2 + micro_edges_per_macro_edge) *
                     (-ctr_2 + micro_edges_per_macro_edge + 1) *
                     (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                    (6)) +
                   (((micro_edges_per_macro_edge + 1) *
                     (micro_edges_per_macro_edge + 2) *
                     (micro_edges_per_macro_edge + 3)) /
                    (6)) +
                   1])));
}
}
{{_mm256_storeu_pd(
    &_data_invDiag_
        [ctr_0 + (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
         (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
         (((micro_edges_per_macro_edge + 1) * (micro_edges_per_macro_edge + 2) *
           (micro_edges_per_macro_edge + 3)) /
          (6)) -
         (((-ctr_2 + micro_edges_per_macro_edge + 1) *
           (-ctr_2 + micro_edges_per_macro_edge + 2) *
           (-ctr_2 + micro_edges_per_macro_edge + 3)) /
          (6))],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_11, tmp_moved_constant_11,
                      tmp_moved_constant_11, tmp_moved_constant_11),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
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
                             (6))])));
_mm256_storeu_pd(
    &_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                    ((ctr_1 * (ctr_1 + 1)) / (2)) -
                    (((-ctr_2 + micro_edges_per_macro_edge) *
                      (-ctr_2 + micro_edges_per_macro_edge + 1) *
                      (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                     (6)) +
                    (((micro_edges_per_macro_edge + 1) *
                      (micro_edges_per_macro_edge + 2) *
                      (micro_edges_per_macro_edge + 3)) /
                     (6))],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_12, tmp_moved_constant_12,
                      tmp_moved_constant_12, tmp_moved_constant_12),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
                            ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6))])));
_mm256_storeu_pd(
    &_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                    ((ctr_1 * (ctr_1 + 1)) / (2)) -
                    (((-ctr_2 + micro_edges_per_macro_edge) *
                      (-ctr_2 + micro_edges_per_macro_edge + 1) *
                      (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                     (6)) +
                    (((micro_edges_per_macro_edge + 1) *
                      (micro_edges_per_macro_edge + 2) *
                      (micro_edges_per_macro_edge + 3)) /
                     (6)) +
                    1],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_13, tmp_moved_constant_13,
                      tmp_moved_constant_13, tmp_moved_constant_13),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
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
                            1])));
_mm256_storeu_pd(
    &_data_invDiag_[ctr_0 +
                    (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                    (((-ctr_2 + micro_edges_per_macro_edge) *
                      (-ctr_2 + micro_edges_per_macro_edge + 1) *
                      (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                     (6)) +
                    (((micro_edges_per_macro_edge + 1) *
                      (micro_edges_per_macro_edge + 2) *
                      (micro_edges_per_macro_edge + 3)) /
                     (6))],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_14, tmp_moved_constant_14,
                      tmp_moved_constant_14, tmp_moved_constant_14),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
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
                             (6))])));
}
}
{{_mm256_storeu_pd(
    &_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                    ((ctr_1 * (ctr_1 + 1)) / (2)) +
                    (((micro_edges_per_macro_edge + 1) *
                      (micro_edges_per_macro_edge + 2) *
                      (micro_edges_per_macro_edge + 3)) /
                     (6)) -
                    (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                      (-ctr_2 + micro_edges_per_macro_edge + 2) *
                      (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                     (6)) +
                    1],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_16, tmp_moved_constant_16,
                      tmp_moved_constant_16, tmp_moved_constant_16),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
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
                            1])));
_mm256_storeu_pd(
    &_data_invDiag_[ctr_0 +
                    (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                    (((micro_edges_per_macro_edge + 1) *
                      (micro_edges_per_macro_edge + 2) *
                      (micro_edges_per_macro_edge + 3)) /
                     (6)) -
                    (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                      (-ctr_2 + micro_edges_per_macro_edge + 2) *
                      (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                     (6))],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_17, tmp_moved_constant_17,
                      tmp_moved_constant_17, tmp_moved_constant_17),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
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
                             (6))])));
_mm256_storeu_pd(
    &_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                    ((ctr_1 * (ctr_1 + 1)) / (2)) -
                    (((-ctr_2 + micro_edges_per_macro_edge) *
                      (-ctr_2 + micro_edges_per_macro_edge + 1) *
                      (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                     (6)) +
                    (((micro_edges_per_macro_edge + 1) *
                      (micro_edges_per_macro_edge + 2) *
                      (micro_edges_per_macro_edge + 3)) /
                     (6))],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_18, tmp_moved_constant_18,
                      tmp_moved_constant_18, tmp_moved_constant_18),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
                            ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) -
                            (((-ctr_2 + micro_edges_per_macro_edge) *
                              (-ctr_2 + micro_edges_per_macro_edge + 1) *
                              (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                             (6)) +
                            (((micro_edges_per_macro_edge + 1) *
                              (micro_edges_per_macro_edge + 2) *
                              (micro_edges_per_macro_edge + 3)) /
                             (6))])));
_mm256_storeu_pd(
    &_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                    ((ctr_1 * (ctr_1 + 1)) / (2)) -
                    (((-ctr_2 + micro_edges_per_macro_edge) *
                      (-ctr_2 + micro_edges_per_macro_edge + 1) *
                      (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                     (6)) +
                    (((micro_edges_per_macro_edge + 1) *
                      (micro_edges_per_macro_edge + 2) *
                      (micro_edges_per_macro_edge + 3)) /
                     (6)) +
                    1],
    _mm256_add_pd(
        _mm256_set_pd(tmp_moved_constant_19, tmp_moved_constant_19,
                      tmp_moved_constant_19, tmp_moved_constant_19),
        _mm256_loadu_pd(
            &_data_invDiag_[ctr_0 +
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
                            1])));
}
}
{
  {
    _mm256_storeu_pd(
        &_data_invDiag_[ctr_0 +
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
                         (6))],
        _mm256_add_pd(
            _mm256_set_pd(tmp_moved_constant_21, tmp_moved_constant_21,
                          tmp_moved_constant_21, tmp_moved_constant_21),
            _mm256_loadu_pd(
                &_data_invDiag_[ctr_0 +
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
                                 (6))])));
    _mm256_storeu_pd(
        &_data_invDiag_[ctr_0 +
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
                        1],
        _mm256_add_pd(
            _mm256_set_pd(tmp_moved_constant_22, tmp_moved_constant_22,
                          tmp_moved_constant_22, tmp_moved_constant_22),
            _mm256_loadu_pd(
                &_data_invDiag_[ctr_0 +
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
                                1])));
    _mm256_storeu_pd(
        &_data_invDiag_[ctr_0 +
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
                        1],
        _mm256_add_pd(
            _mm256_set_pd(tmp_moved_constant_23, tmp_moved_constant_23,
                          tmp_moved_constant_23, tmp_moved_constant_23),
            _mm256_loadu_pd(
                &_data_invDiag_
                    [ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) +
                     1])));
    _mm256_storeu_pd(
        &_data_invDiag_[ctr_0 +
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
                         (6))],
        _mm256_add_pd(
            _mm256_set_pd(tmp_moved_constant_24, tmp_moved_constant_24,
                          tmp_moved_constant_24, tmp_moved_constant_24),
            _mm256_loadu_pd(
                &_data_invDiag_[ctr_0 +
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
                                 (6))])));
  }
}
}
for (int64_t ctr_0 =
         (int64_t)((-ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2) / (4)) *
         (4);
     ctr_0 < -ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2; ctr_0 += 1) {
  {{_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
  _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
      _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
      elMatDiag_2 +
      _data_invDiag_[ctr_0 +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6))];
  _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
      _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
{{_data_invDiag_[ctr_0 +
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
      tmp_moved_constant_1 +
      _data_invDiag_[ctr_0 +
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
                     1];
_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
    tmp_moved_constant_2 +
    _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
_data_invDiag_[ctr_0 + (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
               (((-ctr_2 + micro_edges_per_macro_edge) *
                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                (6)) +
               (((micro_edges_per_macro_edge + 1) *
                 (micro_edges_per_macro_edge + 2) *
                 (micro_edges_per_macro_edge + 3)) /
                (6))] =
    tmp_moved_constant_3 +
    _data_invDiag_[ctr_0 +
                   (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                   (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                   (((-ctr_2 + micro_edges_per_macro_edge) *
                     (-ctr_2 + micro_edges_per_macro_edge + 1) *
                     (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                    (6)) +
                   (((micro_edges_per_macro_edge + 1) *
                     (micro_edges_per_macro_edge + 2) *
                     (micro_edges_per_macro_edge + 3)) /
                    (6))];
_data_invDiag_[ctr_0 + (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
    tmp_moved_constant_4 +
    _data_invDiag_[ctr_0 +
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
                   1];
}
}
{{_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
      tmp_moved_constant_6 +
      _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
_data_invDiag_[ctr_0 + (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
               (((micro_edges_per_macro_edge + 1) *
                 (micro_edges_per_macro_edge + 2) *
                 (micro_edges_per_macro_edge + 3)) /
                (6)) -
               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                (6))] =
    tmp_moved_constant_7 +
    _data_invDiag_[ctr_0 +
                   (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                   (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                   (((micro_edges_per_macro_edge + 1) *
                     (micro_edges_per_macro_edge + 2) *
                     (micro_edges_per_macro_edge + 3)) /
                    (6)) -
                   (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                     (-ctr_2 + micro_edges_per_macro_edge + 2) *
                     (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                    (6))];
_data_invDiag_[ctr_0 + (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
    tmp_moved_constant_8 +
    _data_invDiag_[ctr_0 +
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
                   1];
_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
    tmp_moved_constant_9 +
    _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
{{_data_invDiag_[ctr_0 +
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
      tmp_moved_constant_11 +
      _data_invDiag_[ctr_0 +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6))];
_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
               ((ctr_1 * (ctr_1 + 1)) / (2)) -
               (((-ctr_2 + micro_edges_per_macro_edge) *
                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                (6)) +
               (((micro_edges_per_macro_edge + 1) *
                 (micro_edges_per_macro_edge + 2) *
                 (micro_edges_per_macro_edge + 3)) /
                (6))] =
    tmp_moved_constant_12 +
    _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                   ((ctr_1 * (ctr_1 + 1)) / (2)) -
                   (((-ctr_2 + micro_edges_per_macro_edge) *
                     (-ctr_2 + micro_edges_per_macro_edge + 1) *
                     (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                    (6)) +
                   (((micro_edges_per_macro_edge + 1) *
                     (micro_edges_per_macro_edge + 2) *
                     (micro_edges_per_macro_edge + 3)) /
                    (6))];
_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
    tmp_moved_constant_13 +
    _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
_data_invDiag_[ctr_0 + (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
               (((-ctr_2 + micro_edges_per_macro_edge) *
                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                (6)) +
               (((micro_edges_per_macro_edge + 1) *
                 (micro_edges_per_macro_edge + 2) *
                 (micro_edges_per_macro_edge + 3)) /
                (6))] =
    tmp_moved_constant_14 +
    _data_invDiag_[ctr_0 +
                   (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
{{_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
      tmp_moved_constant_16 +
      _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
_data_invDiag_[ctr_0 + (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
               (((micro_edges_per_macro_edge + 1) *
                 (micro_edges_per_macro_edge + 2) *
                 (micro_edges_per_macro_edge + 3)) /
                (6)) -
               (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                 (-ctr_2 + micro_edges_per_macro_edge + 2) *
                 (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                (6))] =
    tmp_moved_constant_17 +
    _data_invDiag_[ctr_0 +
                   (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                   (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                   (((micro_edges_per_macro_edge + 1) *
                     (micro_edges_per_macro_edge + 2) *
                     (micro_edges_per_macro_edge + 3)) /
                    (6)) -
                   (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                     (-ctr_2 + micro_edges_per_macro_edge + 2) *
                     (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                    (6))];
_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
               ((ctr_1 * (ctr_1 + 1)) / (2)) -
               (((-ctr_2 + micro_edges_per_macro_edge) *
                 (-ctr_2 + micro_edges_per_macro_edge + 1) *
                 (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                (6)) +
               (((micro_edges_per_macro_edge + 1) *
                 (micro_edges_per_macro_edge + 2) *
                 (micro_edges_per_macro_edge + 3)) /
                (6))] =
    tmp_moved_constant_18 +
    _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                   ((ctr_1 * (ctr_1 + 1)) / (2)) -
                   (((-ctr_2 + micro_edges_per_macro_edge) *
                     (-ctr_2 + micro_edges_per_macro_edge + 1) *
                     (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                    (6)) +
                   (((micro_edges_per_macro_edge + 1) *
                     (micro_edges_per_macro_edge + 2) *
                     (micro_edges_per_macro_edge + 3)) /
                    (6))];
_data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
    tmp_moved_constant_19 +
    _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
  {
    _data_invDiag_[ctr_0 +
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
        tmp_moved_constant_21 +
        _data_invDiag_[ctr_0 +
                       (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
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
        tmp_moved_constant_22 +
        _data_invDiag_[ctr_0 +
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
                       1];
    _data_invDiag_[ctr_0 + ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
        tmp_moved_constant_23 +
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
                       1];
    _data_invDiag_[ctr_0 +
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
        tmp_moved_constant_24 +
        _data_invDiag_[ctr_0 +
                       (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
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
if (-ctr_1 - ctr_2 + micro_edges_per_macro_edge - 2 >= 0) {
  const walberla::float64 tmp_kernel_op_0 =
      abs_det_jac_affine_WHITE_UP * 0.16666666666666663;
  const walberla::float64 elMatDiag_0 =
      tmp_kernel_op_0 *
      (((-jac_affine_inv_0_0_WHITE_UP - jac_affine_inv_1_0_WHITE_UP -
         jac_affine_inv_2_0_WHITE_UP) *
        (-jac_affine_inv_0_0_WHITE_UP - jac_affine_inv_1_0_WHITE_UP -
         jac_affine_inv_2_0_WHITE_UP)) +
       ((-jac_affine_inv_0_1_WHITE_UP - jac_affine_inv_1_1_WHITE_UP -
         jac_affine_inv_2_1_WHITE_UP) *
        (-jac_affine_inv_0_1_WHITE_UP - jac_affine_inv_1_1_WHITE_UP -
         jac_affine_inv_2_1_WHITE_UP)) +
       ((-jac_affine_inv_0_2_WHITE_UP - jac_affine_inv_1_2_WHITE_UP -
         jac_affine_inv_2_2_WHITE_UP) *
        (-jac_affine_inv_0_2_WHITE_UP - jac_affine_inv_1_2_WHITE_UP -
         jac_affine_inv_2_2_WHITE_UP)));
  const walberla::float64 elMatDiag_1 =
      tmp_kernel_op_0 *
      ((jac_affine_inv_0_0_WHITE_UP * jac_affine_inv_0_0_WHITE_UP) +
       (jac_affine_inv_0_1_WHITE_UP * jac_affine_inv_0_1_WHITE_UP) +
       (jac_affine_inv_0_2_WHITE_UP * jac_affine_inv_0_2_WHITE_UP));
  const walberla::float64 elMatDiag_2 =
      tmp_kernel_op_0 *
      ((jac_affine_inv_1_0_WHITE_UP * jac_affine_inv_1_0_WHITE_UP) +
       (jac_affine_inv_1_1_WHITE_UP * jac_affine_inv_1_1_WHITE_UP) +
       (jac_affine_inv_1_2_WHITE_UP * jac_affine_inv_1_2_WHITE_UP));
  const walberla::float64 elMatDiag_3 =
      tmp_kernel_op_0 *
      ((jac_affine_inv_2_0_WHITE_UP * jac_affine_inv_2_0_WHITE_UP) +
       (jac_affine_inv_2_1_WHITE_UP * jac_affine_inv_2_1_WHITE_UP) +
       (jac_affine_inv_2_2_WHITE_UP * jac_affine_inv_2_2_WHITE_UP));
  {
    {
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          elMatDiag_0 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2) *
                           (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         2];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          elMatDiag_1 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2) *
                           (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          elMatDiag_2 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         2];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          elMatDiag_3 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         2];
    }
  }
  const walberla::float64 tmp_moved_constant_25 =
      abs_det_jac_affine_BLUE_UP * 0.16666666666666663;
  const walberla::float64 tmp_moved_constant_26 =
      tmp_moved_constant_25 *
      (((-jac_affine_inv_0_0_BLUE_UP - jac_affine_inv_1_0_BLUE_UP -
         jac_affine_inv_2_0_BLUE_UP) *
        (-jac_affine_inv_0_0_BLUE_UP - jac_affine_inv_1_0_BLUE_UP -
         jac_affine_inv_2_0_BLUE_UP)) +
       ((-jac_affine_inv_0_1_BLUE_UP - jac_affine_inv_1_1_BLUE_UP -
         jac_affine_inv_2_1_BLUE_UP) *
        (-jac_affine_inv_0_1_BLUE_UP - jac_affine_inv_1_1_BLUE_UP -
         jac_affine_inv_2_1_BLUE_UP)) +
       ((-jac_affine_inv_0_2_BLUE_UP - jac_affine_inv_1_2_BLUE_UP -
         jac_affine_inv_2_2_BLUE_UP) *
        (-jac_affine_inv_0_2_BLUE_UP - jac_affine_inv_1_2_BLUE_UP -
         jac_affine_inv_2_2_BLUE_UP)));
  const walberla::float64 tmp_moved_constant_27 =
      tmp_moved_constant_25 *
      ((jac_affine_inv_0_0_BLUE_UP * jac_affine_inv_0_0_BLUE_UP) +
       (jac_affine_inv_0_1_BLUE_UP * jac_affine_inv_0_1_BLUE_UP) +
       (jac_affine_inv_0_2_BLUE_UP * jac_affine_inv_0_2_BLUE_UP));
  const walberla::float64 tmp_moved_constant_28 =
      tmp_moved_constant_25 *
      ((jac_affine_inv_1_0_BLUE_UP * jac_affine_inv_1_0_BLUE_UP) +
       (jac_affine_inv_1_1_BLUE_UP * jac_affine_inv_1_1_BLUE_UP) +
       (jac_affine_inv_1_2_BLUE_UP * jac_affine_inv_1_2_BLUE_UP));
  const walberla::float64 tmp_moved_constant_29 =
      tmp_moved_constant_25 *
      ((jac_affine_inv_2_0_BLUE_UP * jac_affine_inv_2_0_BLUE_UP) +
       (jac_affine_inv_2_1_BLUE_UP * jac_affine_inv_2_1_BLUE_UP) +
       (jac_affine_inv_2_2_BLUE_UP * jac_affine_inv_2_2_BLUE_UP));
  {
    {
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_26 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2) *
                           (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_27 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         2];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_28 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         1];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_29 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
    }
  }
  const walberla::float64 tmp_moved_constant_30 =
      abs_det_jac_affine_BLUE_DOWN * 0.16666666666666663;
  const walberla::float64 tmp_moved_constant_31 =
      tmp_moved_constant_30 *
      (((-jac_affine_inv_0_0_BLUE_DOWN - jac_affine_inv_1_0_BLUE_DOWN -
         jac_affine_inv_2_0_BLUE_DOWN) *
        (-jac_affine_inv_0_0_BLUE_DOWN - jac_affine_inv_1_0_BLUE_DOWN -
         jac_affine_inv_2_0_BLUE_DOWN)) +
       ((-jac_affine_inv_0_1_BLUE_DOWN - jac_affine_inv_1_1_BLUE_DOWN -
         jac_affine_inv_2_1_BLUE_DOWN) *
        (-jac_affine_inv_0_1_BLUE_DOWN - jac_affine_inv_1_1_BLUE_DOWN -
         jac_affine_inv_2_1_BLUE_DOWN)) +
       ((-jac_affine_inv_0_2_BLUE_DOWN - jac_affine_inv_1_2_BLUE_DOWN -
         jac_affine_inv_2_2_BLUE_DOWN) *
        (-jac_affine_inv_0_2_BLUE_DOWN - jac_affine_inv_1_2_BLUE_DOWN -
         jac_affine_inv_2_2_BLUE_DOWN)));
  const walberla::float64 tmp_moved_constant_32 =
      tmp_moved_constant_30 *
      ((jac_affine_inv_0_0_BLUE_DOWN * jac_affine_inv_0_0_BLUE_DOWN) +
       (jac_affine_inv_0_1_BLUE_DOWN * jac_affine_inv_0_1_BLUE_DOWN) +
       (jac_affine_inv_0_2_BLUE_DOWN * jac_affine_inv_0_2_BLUE_DOWN));
  const walberla::float64 tmp_moved_constant_33 =
      tmp_moved_constant_30 *
      ((jac_affine_inv_1_0_BLUE_DOWN * jac_affine_inv_1_0_BLUE_DOWN) +
       (jac_affine_inv_1_1_BLUE_DOWN * jac_affine_inv_1_1_BLUE_DOWN) +
       (jac_affine_inv_1_2_BLUE_DOWN * jac_affine_inv_1_2_BLUE_DOWN));
  const walberla::float64 tmp_moved_constant_34 =
      tmp_moved_constant_30 *
      ((jac_affine_inv_2_0_BLUE_DOWN * jac_affine_inv_2_0_BLUE_DOWN) +
       (jac_affine_inv_2_1_BLUE_DOWN * jac_affine_inv_2_1_BLUE_DOWN) +
       (jac_affine_inv_2_2_BLUE_DOWN * jac_affine_inv_2_2_BLUE_DOWN));
  {
    {
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_31 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         2];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_32 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         2];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_33 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_34 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         2];
    }
  }
  const walberla::float64 tmp_moved_constant_35 =
      abs_det_jac_affine_GREEN_UP * 0.16666666666666663;
  const walberla::float64 tmp_moved_constant_36 =
      tmp_moved_constant_35 *
      (((-jac_affine_inv_0_0_GREEN_UP - jac_affine_inv_1_0_GREEN_UP -
         jac_affine_inv_2_0_GREEN_UP) *
        (-jac_affine_inv_0_0_GREEN_UP - jac_affine_inv_1_0_GREEN_UP -
         jac_affine_inv_2_0_GREEN_UP)) +
       ((-jac_affine_inv_0_1_GREEN_UP - jac_affine_inv_1_1_GREEN_UP -
         jac_affine_inv_2_1_GREEN_UP) *
        (-jac_affine_inv_0_1_GREEN_UP - jac_affine_inv_1_1_GREEN_UP -
         jac_affine_inv_2_1_GREEN_UP)) +
       ((-jac_affine_inv_0_2_GREEN_UP - jac_affine_inv_1_2_GREEN_UP -
         jac_affine_inv_2_2_GREEN_UP) *
        (-jac_affine_inv_0_2_GREEN_UP - jac_affine_inv_1_2_GREEN_UP -
         jac_affine_inv_2_2_GREEN_UP)));
  const walberla::float64 tmp_moved_constant_37 =
      tmp_moved_constant_35 *
      ((jac_affine_inv_0_0_GREEN_UP * jac_affine_inv_0_0_GREEN_UP) +
       (jac_affine_inv_0_1_GREEN_UP * jac_affine_inv_0_1_GREEN_UP) +
       (jac_affine_inv_0_2_GREEN_UP * jac_affine_inv_0_2_GREEN_UP));
  const walberla::float64 tmp_moved_constant_38 =
      tmp_moved_constant_35 *
      ((jac_affine_inv_1_0_GREEN_UP * jac_affine_inv_1_0_GREEN_UP) +
       (jac_affine_inv_1_1_GREEN_UP * jac_affine_inv_1_1_GREEN_UP) +
       (jac_affine_inv_1_2_GREEN_UP * jac_affine_inv_1_2_GREEN_UP));
  const walberla::float64 tmp_moved_constant_39 =
      tmp_moved_constant_35 *
      ((jac_affine_inv_2_0_GREEN_UP * jac_affine_inv_2_0_GREEN_UP) +
       (jac_affine_inv_2_1_GREEN_UP * jac_affine_inv_2_1_GREEN_UP) +
       (jac_affine_inv_2_2_GREEN_UP * jac_affine_inv_2_2_GREEN_UP));
  {
    {
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_36 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2) *
                           (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_37 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         2];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_38 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         2];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_39 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
    }
  }
  const walberla::float64 tmp_moved_constant_40 =
      abs_det_jac_affine_GREEN_DOWN * 0.16666666666666663;
  const walberla::float64 tmp_moved_constant_41 =
      tmp_moved_constant_40 *
      (((-jac_affine_inv_0_0_GREEN_DOWN - jac_affine_inv_1_0_GREEN_DOWN -
         jac_affine_inv_2_0_GREEN_DOWN) *
        (-jac_affine_inv_0_0_GREEN_DOWN - jac_affine_inv_1_0_GREEN_DOWN -
         jac_affine_inv_2_0_GREEN_DOWN)) +
       ((-jac_affine_inv_0_1_GREEN_DOWN - jac_affine_inv_1_1_GREEN_DOWN -
         jac_affine_inv_2_1_GREEN_DOWN) *
        (-jac_affine_inv_0_1_GREEN_DOWN - jac_affine_inv_1_1_GREEN_DOWN -
         jac_affine_inv_2_1_GREEN_DOWN)) +
       ((-jac_affine_inv_0_2_GREEN_DOWN - jac_affine_inv_1_2_GREEN_DOWN -
         jac_affine_inv_2_2_GREEN_DOWN) *
        (-jac_affine_inv_0_2_GREEN_DOWN - jac_affine_inv_1_2_GREEN_DOWN -
         jac_affine_inv_2_2_GREEN_DOWN)));
  const walberla::float64 tmp_moved_constant_42 =
      tmp_moved_constant_40 *
      ((jac_affine_inv_0_0_GREEN_DOWN * jac_affine_inv_0_0_GREEN_DOWN) +
       (jac_affine_inv_0_1_GREEN_DOWN * jac_affine_inv_0_1_GREEN_DOWN) +
       (jac_affine_inv_0_2_GREEN_DOWN * jac_affine_inv_0_2_GREEN_DOWN));
  const walberla::float64 tmp_moved_constant_43 =
      tmp_moved_constant_40 *
      ((jac_affine_inv_1_0_GREEN_DOWN * jac_affine_inv_1_0_GREEN_DOWN) +
       (jac_affine_inv_1_1_GREEN_DOWN * jac_affine_inv_1_1_GREEN_DOWN) +
       (jac_affine_inv_1_2_GREEN_DOWN * jac_affine_inv_1_2_GREEN_DOWN));
  const walberla::float64 tmp_moved_constant_44 =
      tmp_moved_constant_40 *
      ((jac_affine_inv_2_0_GREEN_DOWN * jac_affine_inv_2_0_GREEN_DOWN) +
       (jac_affine_inv_2_1_GREEN_DOWN * jac_affine_inv_2_1_GREEN_DOWN) +
       (jac_affine_inv_2_2_GREEN_DOWN * jac_affine_inv_2_2_GREEN_DOWN));
  {
    {
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_41 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         2];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_42 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         1];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          tmp_moved_constant_43 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     2] =
          tmp_moved_constant_44 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         2];
    }
  }
}
{
  {
    {
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          elMatDiag_0 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2) *
                           (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6))] =
          elMatDiag_1 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2) *
                           (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                          (6))];
      _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
                     (ctr_1 + 1) * (-ctr_2 + micro_edges_per_macro_edge + 2) -
                     (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     (((-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2) *
                       (-ctr_2 + micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          elMatDiag_2 +
          _data_invDiag_[-ctr_1 - ctr_2 + micro_edges_per_macro_edge +
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
                          (6)) -
                         1];
      _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) - ctr_1 -
                     ctr_2 + micro_edges_per_macro_edge -
                     ((ctr_1 * (ctr_1 + 1)) / (2)) -
                     (((-ctr_2 + micro_edges_per_macro_edge) *
                       (-ctr_2 + micro_edges_per_macro_edge + 1) *
                       (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                      (6)) +
                     (((micro_edges_per_macro_edge + 1) *
                       (micro_edges_per_macro_edge + 2) *
                       (micro_edges_per_macro_edge + 3)) /
                      (6)) -
                     1] =
          elMatDiag_3 +
          _data_invDiag_[ctr_1 * (-ctr_2 + micro_edges_per_macro_edge + 1) -
                         ctr_1 - ctr_2 + micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) -
                         (((-ctr_2 + micro_edges_per_macro_edge) *
                           (-ctr_2 + micro_edges_per_macro_edge + 1) *
                           (-ctr_2 + micro_edges_per_macro_edge + 2)) /
                          (6)) +
                         (((micro_edges_per_macro_edge + 1) *
                           (micro_edges_per_macro_edge + 2) *
                           (micro_edges_per_macro_edge + 3)) /
                          (6)) -
                         1];
    }
  }
}
}
}
}

} // namespace operatorgeneration

} // namespace hyteg
