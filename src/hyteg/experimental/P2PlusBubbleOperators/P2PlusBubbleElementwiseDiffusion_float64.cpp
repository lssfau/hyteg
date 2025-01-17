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

#include "P2PlusBubbleElementwiseDiffusion_float64.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P2PlusBubbleElementwiseDiffusion_float64::
    P2PlusBubbleElementwiseDiffusion_float64(
        const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
        size_t maxLevel)
    : Operator(storage, minLevel, maxLevel) {}

void P2PlusBubbleElementwiseDiffusion_float64::apply(
    const P2PlusBubbleFunction<walberla::float64> &src,
    const P2PlusBubbleFunction<walberla::float64> &dst, uint_t level,
    DoFType flag, UpdateType updateType) const {
  this->startTiming("apply");

  // Make sure that halos are up-to-date
  this->timingTree_->start("pre-communication");
  if (this->storage_->hasGlobalCells()) {
    WALBERLA_ABORT("Not implemented.");
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
    WALBERLA_ABORT("Not implemented.");
  } else {
    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      // get hold of the actual numerical data in the functions
      walberla::float64 *_data_srcVertex =
          face.getData(src.getVertexDoFFunction().getFaceDataID())
              ->getPointer(level);
      walberla::float64 *_data_srcEdge =
          face.getData(src.getEdgeDoFFunction().getFaceDataID())
              ->getPointer(level);
      walberla::float64 *_data_src =
          src.getVolumeDoFFunction().dofMemory(it.first, level);
      walberla::float64 *_data_dstVertex =
          face.getData(dst.getVertexDoFFunction().getFaceDataID())
              ->getPointer(level);
      walberla::float64 *_data_dstEdge =
          face.getData(dst.getEdgeDoFFunction().getFaceDataID())
              ->getPointer(level);
      walberla::float64 *_data_dst =
          dst.getVolumeDoFFunction().dofMemory(it.first, level);

      // Zero out dst halos only
      //
      // This is also necessary when using update type == Add.
      // During additive comm we then skip zeroing the data on the lower-dim
      // primitives.
      for (const auto &idx : vertexdof::macroface::Iterator(level)) {
        if (vertexdof::macroface::isVertexOnBoundary(level, idx)) {
          auto arrayIdx = vertexdof::macroface::index(level, idx.x(), idx.y());
          _data_dstVertex[arrayIdx] =
              walberla::numeric_cast<walberla::float64>(0);
        }
      }
      for (const auto &idx : edgedof::macroface::Iterator(level)) {
        for (const auto &orientation : edgedof::faceLocalEdgeDoFOrientations) {
          if (!edgedof::macroface::isInnerEdgeDoF(level, idx, orientation)) {
            auto arrayIdx =
                edgedof::macroface::index(level, idx.x(), idx.y(), orientation);
            _data_dstEdge[arrayIdx] =
                walberla::numeric_cast<walberla::float64>(0);
          }
        }
      }

      const auto micro_edges_per_macro_edge =
          (int64_t)levelinfo::num_microedges_per_edge(level);
      const auto num_microfaces_per_face =
          (int64_t)levelinfo::num_microfaces_per_face(level);
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

      apply_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(

          _data_dst, _data_dstEdge, _data_dstVertex, _data_src, _data_srcEdge,
          _data_srcVertex, macro_vertex_coord_id_0comp0,
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
    dst.getVertexDoFFunction().communicateAdditively<Face, Edge>(
        level, DoFType::All ^ flag, *storage_, updateType == Replace);
    dst.getVertexDoFFunction().communicateAdditively<Face, Vertex>(
        level, DoFType::All ^ flag, *storage_, updateType == Replace);
    dst.getEdgeDoFFunction().communicateAdditively<Face, Edge>(
        level, DoFType::All ^ flag, *storage_, updateType == Replace);
    this->timingTree_->stop("post-communication");
  }

  this->stopTiming("apply");
}
void P2PlusBubbleElementwiseDiffusion_float64::toMatrix(
    const std::shared_ptr<SparseMatrixProxy> &mat,
    const P2PlusBubbleFunction<idx_t> &src,
    const P2PlusBubbleFunction<idx_t> &dst, uint_t level, DoFType flag) const {
  this->startTiming("toMatrix");

  // We currently ignore the flag provided!
  if (flag != All) {
    WALBERLA_LOG_WARNING_ON_ROOT(
        "Input flag ignored in toMatrix; using flag = All");
  }

  if (storage_->hasGlobalCells()) {
    this->timingTree_->start("pre-communication");

    this->timingTree_->stop("pre-communication");

    WALBERLA_ABORT("Not implemented.");
  } else {
    this->timingTree_->start("pre-communication");

    this->timingTree_->stop("pre-communication");

    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      // get hold of the actual numerical data
      idx_t *_data_srcVertex =
          face.getData(src.getVertexDoFFunction().getFaceDataID())
              ->getPointer(level);
      idx_t *_data_srcEdge =
          face.getData(src.getEdgeDoFFunction().getFaceDataID())
              ->getPointer(level);
      idx_t *_data_src = src.getVolumeDoFFunction().dofMemory(it.first, level);
      idx_t *_data_dstVertex =
          face.getData(dst.getVertexDoFFunction().getFaceDataID())
              ->getPointer(level);
      idx_t *_data_dstEdge =
          face.getData(dst.getEdgeDoFFunction().getFaceDataID())
              ->getPointer(level);
      idx_t *_data_dst = dst.getVolumeDoFFunction().dofMemory(it.first, level);

      const auto micro_edges_per_macro_edge =
          (int64_t)levelinfo::num_microedges_per_edge(level);
      const auto num_microfaces_per_face =
          (int64_t)levelinfo::num_microfaces_per_face(level);
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

      toMatrix_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(

          _data_dst, _data_dstEdge, _data_dstVertex, _data_src, _data_srcEdge,
          _data_srcVertex, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_1comp0,
          macro_vertex_coord_id_1comp1, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, mat, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float);

      this->timingTree_->stop("kernel");
    }
  }
  this->stopTiming("toMatrix");
}
void P2PlusBubbleElementwiseDiffusion_float64::
    computeInverseDiagonalOperatorValues() {
  this->startTiming("computeInverseDiagonalOperatorValues");

  if (invDiag_ == nullptr) {
    invDiag_ = std::make_shared<P2PlusBubbleFunction<walberla::float64>>(
        "inverse diagonal entries", storage_, minLevel_, maxLevel_);
  }

  for (uint_t level = minLevel_; level <= maxLevel_; level++) {
    invDiag_->setToZero(level);

    if (storage_->hasGlobalCells()) {
      this->timingTree_->start("pre-communication");

      this->timingTree_->stop("pre-communication");

      WALBERLA_ABORT("Not implemented.");
      (*invDiag_).invertElementwise(level);
    } else {
      this->timingTree_->start("pre-communication");

      this->timingTree_->stop("pre-communication");

      for (auto &it : storage_->getFaces()) {
        Face &face = *it.second;

        // get hold of the actual numerical data
        walberla::float64 *_data_invDiag_Vertex =
            face.getData((*invDiag_).getVertexDoFFunction().getFaceDataID())
                ->getPointer(level);
        walberla::float64 *_data_invDiag_Edge =
            face.getData((*invDiag_).getEdgeDoFFunction().getFaceDataID())
                ->getPointer(level);
        walberla::float64 *_data_invDiag_ =
            (*invDiag_).getVolumeDoFFunction().dofMemory(it.first, level);

        const auto micro_edges_per_macro_edge =
            (int64_t)levelinfo::num_microedges_per_edge(level);
        const auto num_microfaces_per_face =
            (int64_t)levelinfo::num_microfaces_per_face(level);
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

        computeInverseDiagonalOperatorValues_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(

            _data_invDiag_, _data_invDiag_Edge, _data_invDiag_Vertex,
            macro_vertex_coord_id_0comp0, macro_vertex_coord_id_0comp1,
            macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
            macro_vertex_coord_id_2comp0, macro_vertex_coord_id_2comp1,
            micro_edges_per_macro_edge, micro_edges_per_macro_edge_float);

        this->timingTree_->stop("kernel");
      }

      // Push result to lower-dimensional primitives
      //
      this->timingTree_->start("post-communication");
      // Note: We could avoid communication here by implementing the apply()
      // also for the respective
      //       lower dimensional primitives!
      (*invDiag_).getVertexDoFFunction().communicateAdditively<Face, Edge>(
          level);
      (*invDiag_).getVertexDoFFunction().communicateAdditively<Face, Vertex>(
          level);
      (*invDiag_).getEdgeDoFFunction().communicateAdditively<Face, Edge>(level);
      this->timingTree_->stop("post-communication");
      (*invDiag_).invertElementwise(level);
    }
  }

  this->stopTiming("computeInverseDiagonalOperatorValues");
}
std::shared_ptr<P2PlusBubbleFunction<walberla::float64>>
P2PlusBubbleElementwiseDiffusion_float64::getInverseDiagonalValues() const {
  return invDiag_;
}
void P2PlusBubbleElementwiseDiffusion_float64::
    apply_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(
        walberla::float64 *RESTRICT _data_dst,
        walberla::float64 *RESTRICT _data_dstEdge,
        walberla::float64 *RESTRICT _data_dstVertex,
        walberla::float64 *RESTRICT _data_src,
        walberla::float64 *RESTRICT _data_srcEdge,
        walberla::float64 *RESTRICT _data_srcVertex,
        walberla::float64 macro_vertex_coord_id_0comp0,
        walberla::float64 macro_vertex_coord_id_0comp1,
        walberla::float64 macro_vertex_coord_id_1comp0,
        walberla::float64 macro_vertex_coord_id_1comp1,
        walberla::float64 macro_vertex_coord_id_2comp0,
        walberla::float64 macro_vertex_coord_id_2comp1,
        int64_t micro_edges_per_macro_edge,
        walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    {
      /* FaceType.GRAY */
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
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge;
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
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_0_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_2_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 src_dof_0 =
              _data_srcVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 src_dof_1 =
              _data_srcVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          const walberla::float64 src_dof_2 =
              _data_srcVertex[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          const walberla::float64 src_dof_3 =
              _data_srcEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge + 1)) /
                             (2))];
          const walberla::float64 src_dof_4 =
              _data_srcEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            2 * ((micro_edges_per_macro_edge *
                                  (micro_edges_per_macro_edge + 1)) /
                                 (2))];
          const walberla::float64 src_dof_5 =
              _data_srcEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 src_dof_6 =
              _data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 tmp_kernel_op_0 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_1 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_2 =
              tmp_kernel_op_0 + tmp_kernel_op_1 - 3.0;
          const walberla::float64 tmp_kernel_op_3 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_2 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_5 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_6 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_7 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_8 =
              tmp_kernel_op_6 + tmp_kernel_op_7 - 3.0;
          const walberla::float64 tmp_kernel_op_9 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_8 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_10 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_8 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_11 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_12 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_13 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_12 + tmp_kernel_op_13 - 3.0;
          const walberla::float64 tmp_kernel_op_15 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_14 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_16 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_14 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_17 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_18 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_19 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_20 =
              tmp_kernel_op_18 + tmp_kernel_op_19 - 3.0;
          const walberla::float64 tmp_kernel_op_21 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_20 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_20 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_23 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_24 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_25 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_26 =
              tmp_kernel_op_24 + tmp_kernel_op_25 - 3.0;
          const walberla::float64 tmp_kernel_op_27 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_26 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_28 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_26 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_29 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_30 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_31 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_30 + tmp_kernel_op_31 - 3.0;
          const walberla::float64 tmp_kernel_op_33 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_32 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_34 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_32 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_35 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_36 = tmp_kernel_op_0 - 1.0;
          const walberla::float64 tmp_kernel_op_37 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_38 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_39 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_40 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_41 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_42 = tmp_kernel_op_12 - 1.0;
          const walberla::float64 tmp_kernel_op_43 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_44 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_45 = tmp_kernel_op_18 - 1.0;
          const walberla::float64 tmp_kernel_op_46 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_47 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_48 = tmp_kernel_op_24 - 1.0;
          const walberla::float64 tmp_kernel_op_49 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_50 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_51 = tmp_kernel_op_30 - 1.0;
          const walberla::float64 tmp_kernel_op_52 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_53 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_41 +
                                  tmp_kernel_op_40 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_43 +
                                  tmp_kernel_op_16 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_46 +
                                  tmp_kernel_op_22 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_49 +
                                  tmp_kernel_op_28 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_52 +
                                  tmp_kernel_op_34 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_37 +
                                 tmp_kernel_op_38 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_55 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_56 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_57 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_58 = tmp_kernel_op_7 - 1.0;
          const walberla::float64 tmp_kernel_op_59 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_60 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_61 = tmp_kernel_op_13 - 1.0;
          const walberla::float64 tmp_kernel_op_62 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_63 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_64 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_65 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_66 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_67 = tmp_kernel_op_25 - 1.0;
          const walberla::float64 tmp_kernel_op_68 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_69 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_70 = tmp_kernel_op_31 - 1.0;
          const walberla::float64 tmp_kernel_op_71 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_72 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_60 +
                                  tmp_kernel_op_59 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_62 +
                                  tmp_kernel_op_16 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_65 +
                                  tmp_kernel_op_22 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_68 +
                                  tmp_kernel_op_28 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_71 +
                                  tmp_kernel_op_34 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_56 +
                                 tmp_kernel_op_4 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_74 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_75 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_74 + tmp_kernel_op_75;
          const walberla::float64 tmp_kernel_op_77 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_78 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_77 + tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_81 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_80 + tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_83 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_84 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_83 + tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_86 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_87 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_86 + tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_90 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_91 =
              tmp_kernel_op_89 + tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_92 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_93 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_92 + tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_96 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_95 + tmp_kernel_op_96;
          const walberla::float64 tmp_kernel_op_98 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_99 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_98 + tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_101 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_102 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_103 =
              tmp_kernel_op_101 + tmp_kernel_op_102;
          const walberla::float64 tmp_kernel_op_104 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_105 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_104 + tmp_kernel_op_105;
          const walberla::float64 tmp_kernel_op_107 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_108 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_109 =
              tmp_kernel_op_107 + tmp_kernel_op_108;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_85 +
                                  tmp_kernel_op_82 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_88 +
                                  tmp_kernel_op_16 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_94 +
                                  tmp_kernel_op_22 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_27 +
                                  tmp_kernel_op_103 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_33 +
                                  tmp_kernel_op_109 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_76 +
                                 tmp_kernel_op_4 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_111 =
              -tmp_kernel_op_1 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_112 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_111 - tmp_kernel_op_75;
          const walberla::float64 tmp_kernel_op_113 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_111 - tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_114 =
              -tmp_kernel_op_7 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_115 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_114 - tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_116 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_114 - tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_117 =
              -tmp_kernel_op_13 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_118 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_117 - tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_119 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_117 - tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_120 =
              -tmp_kernel_op_19 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_121 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_120 - tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_122 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_120 - tmp_kernel_op_96;
          const walberla::float64 tmp_kernel_op_123 =
              -tmp_kernel_op_25 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_124 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_123 - tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_125 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_123 - tmp_kernel_op_102;
          const walberla::float64 tmp_kernel_op_126 =
              -tmp_kernel_op_31 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_127 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_126 - tmp_kernel_op_105;
          const walberla::float64 tmp_kernel_op_128 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_126 - tmp_kernel_op_108;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_116 +
                                  tmp_kernel_op_115 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_15 +
                                  tmp_kernel_op_119 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_21 +
                                  tmp_kernel_op_122 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_27 +
                                  tmp_kernel_op_125 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_33 +
                                  tmp_kernel_op_128 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_3 +
                                 tmp_kernel_op_113 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_0 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_131 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_130 - tmp_kernel_op_74;
          const walberla::float64 tmp_kernel_op_132 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_130 - tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_133 =
              -tmp_kernel_op_6 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_134 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_133 - tmp_kernel_op_80;
          const walberla::float64 tmp_kernel_op_135 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_133 - tmp_kernel_op_83;
          const walberla::float64 tmp_kernel_op_136 =
              -tmp_kernel_op_12 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_137 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_136 - tmp_kernel_op_86;
          const walberla::float64 tmp_kernel_op_138 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_136 - tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_139 =
              -tmp_kernel_op_18 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_140 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_139 - tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_141 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_139 - tmp_kernel_op_95;
          const walberla::float64 tmp_kernel_op_142 =
              -tmp_kernel_op_24 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_143 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_142 - tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_144 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_142 - tmp_kernel_op_101;
          const walberla::float64 tmp_kernel_op_145 =
              -tmp_kernel_op_30 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_146 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_145 - tmp_kernel_op_104;
          const walberla::float64 tmp_kernel_op_147 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_145 - tmp_kernel_op_107;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_135 +
                                  tmp_kernel_op_134 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_15 +
                                  tmp_kernel_op_138 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_21 +
                                  tmp_kernel_op_141 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_27 +
                                  tmp_kernel_op_144 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_33 +
                                  tmp_kernel_op_147 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_3 +
                                 tmp_kernel_op_132 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_149 =
              jac_affine_inv_0_0_GRAY * 27.0;
          const int64_t tmp_kernel_op_150 = 0;
          const walberla::float64 tmp_kernel_op_151 =
              jac_affine_inv_1_0_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_152 = 0.15066167873471437;
          const walberla::float64 tmp_kernel_op_153 =
              tmp_kernel_op_149 * ((walberla::float64)(tmp_kernel_op_150)) +
              tmp_kernel_op_151 * tmp_kernel_op_152;
          const walberla::float64 tmp_kernel_op_154 =
              jac_affine_inv_0_1_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_155 =
              jac_affine_inv_1_1_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_152 * tmp_kernel_op_155 +
              tmp_kernel_op_154 * ((walberla::float64)(tmp_kernel_op_150));
          const walberla::float64 tmp_kernel_op_157 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_158 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_149 * tmp_kernel_op_157 +
              tmp_kernel_op_151 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_154 * tmp_kernel_op_157 +
              tmp_kernel_op_155 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_161 = 0.15066167873471437;
          const int64_t tmp_kernel_op_162 = 0;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_149 * tmp_kernel_op_161 +
              tmp_kernel_op_151 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_154 * tmp_kernel_op_161 +
              tmp_kernel_op_155 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_165 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_166 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_149 * tmp_kernel_op_165 +
              tmp_kernel_op_151 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_154 * tmp_kernel_op_165 +
              tmp_kernel_op_155 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_169 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_170 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_149 * tmp_kernel_op_169 +
              tmp_kernel_op_151 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_154 * tmp_kernel_op_169 +
              tmp_kernel_op_155 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_173 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_174 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_149 * tmp_kernel_op_173 +
              tmp_kernel_op_151 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_154 * tmp_kernel_op_173 +
              tmp_kernel_op_155 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_160 +
                                  tmp_kernel_op_159 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_163 +
                                  tmp_kernel_op_16 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_21 +
                                  tmp_kernel_op_168 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_27 +
                                  tmp_kernel_op_172 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_33 +
                                  tmp_kernel_op_176 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_3 +
                                 tmp_kernel_op_156 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_178 =
              (jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY);
          const walberla::float64 tmp_kernel_op_179 =
              (tmp_kernel_op_36 * tmp_kernel_op_36);
          const walberla::float64 tmp_kernel_op_180 =
              (jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY);
          const walberla::float64 tmp_kernel_op_181 =
              (tmp_kernel_op_39 * tmp_kernel_op_39);
          const walberla::float64 tmp_kernel_op_182 =
              (tmp_kernel_op_42 * tmp_kernel_op_42);
          const walberla::float64 tmp_kernel_op_183 =
              (tmp_kernel_op_45 * tmp_kernel_op_45);
          const walberla::float64 tmp_kernel_op_184 =
              (tmp_kernel_op_48 * tmp_kernel_op_48);
          const walberla::float64 tmp_kernel_op_185 =
              (tmp_kernel_op_51 * tmp_kernel_op_51);
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_59 +
                                  tmp_kernel_op_41 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_62 +
                                  tmp_kernel_op_44 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_65 +
                                  tmp_kernel_op_47 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_49 * tmp_kernel_op_68 +
                                  tmp_kernel_op_50 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_52 * tmp_kernel_op_71 +
                                  tmp_kernel_op_53 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_56 +
                                 tmp_kernel_op_38 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_82 +
                                  tmp_kernel_op_41 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_88 +
                                  tmp_kernel_op_44 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_94 +
                                  tmp_kernel_op_47 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_49 +
                                  tmp_kernel_op_103 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_52 +
                                  tmp_kernel_op_109 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_76 +
                                 tmp_kernel_op_38 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_40 +
                                  tmp_kernel_op_116 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_43 +
                                  tmp_kernel_op_119 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_46 +
                                  tmp_kernel_op_122 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_49 +
                                  tmp_kernel_op_125 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_52 +
                                  tmp_kernel_op_128 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_37 +
                                 tmp_kernel_op_113 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_40 +
                                  tmp_kernel_op_135 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_43 +
                                  tmp_kernel_op_138 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_46 +
                                  tmp_kernel_op_141 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_49 +
                                  tmp_kernel_op_144 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_52 +
                                  tmp_kernel_op_147 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_37 +
                                 tmp_kernel_op_132 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_40 +
                                  tmp_kernel_op_160 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_43 +
                                  tmp_kernel_op_164 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_46 +
                                  tmp_kernel_op_168 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_49 +
                                  tmp_kernel_op_172 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_52 +
                                  tmp_kernel_op_176 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_37 +
                                 tmp_kernel_op_156 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_191 =
              (jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY);
          const walberla::float64 tmp_kernel_op_192 =
              (tmp_kernel_op_55 * tmp_kernel_op_55);
          const walberla::float64 tmp_kernel_op_193 =
              (jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY);
          const walberla::float64 tmp_kernel_op_194 =
              (tmp_kernel_op_58 * tmp_kernel_op_58);
          const walberla::float64 tmp_kernel_op_195 =
              (tmp_kernel_op_61 * tmp_kernel_op_61);
          const walberla::float64 tmp_kernel_op_196 =
              (tmp_kernel_op_64 * tmp_kernel_op_64);
          const walberla::float64 tmp_kernel_op_197 =
              (tmp_kernel_op_67 * tmp_kernel_op_67);
          const walberla::float64 tmp_kernel_op_198 =
              (tmp_kernel_op_70 * tmp_kernel_op_70);
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_11 * (tmp_kernel_op_59 * tmp_kernel_op_82 +
                                  tmp_kernel_op_60 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_62 * tmp_kernel_op_88 +
                                  tmp_kernel_op_63 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_65 * tmp_kernel_op_94 +
                                  tmp_kernel_op_66 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_68 +
                                  tmp_kernel_op_103 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_71 +
                                  tmp_kernel_op_109 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_56 * tmp_kernel_op_76 +
                                 tmp_kernel_op_57 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_59 +
                                  tmp_kernel_op_116 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_62 +
                                  tmp_kernel_op_119 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_65 +
                                  tmp_kernel_op_122 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_68 +
                                  tmp_kernel_op_125 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_71 +
                                  tmp_kernel_op_128 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_56 +
                                 tmp_kernel_op_113 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_201 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_59 +
                                  tmp_kernel_op_135 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_62 +
                                  tmp_kernel_op_138 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_65 +
                                  tmp_kernel_op_141 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_68 +
                                  tmp_kernel_op_144 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_71 +
                                  tmp_kernel_op_147 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_56 +
                                 tmp_kernel_op_132 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_59 +
                                  tmp_kernel_op_160 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_62 +
                                  tmp_kernel_op_164 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_65 +
                                  tmp_kernel_op_168 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_68 +
                                  tmp_kernel_op_172 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_71 +
                                  tmp_kernel_op_176 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_56 +
                                 tmp_kernel_op_156 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_203 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_82 +
                                  tmp_kernel_op_116 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_88 +
                                  tmp_kernel_op_119 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_94 +
                                  tmp_kernel_op_122 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_124 +
                                  tmp_kernel_op_103 * tmp_kernel_op_125) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_127 +
                                  tmp_kernel_op_109 * tmp_kernel_op_128) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_76 +
                                 tmp_kernel_op_113 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_134 +
                                  tmp_kernel_op_116 * tmp_kernel_op_135) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_137 +
                                  tmp_kernel_op_119 * tmp_kernel_op_138) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_140 +
                                  tmp_kernel_op_122 * tmp_kernel_op_141) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_143 +
                                  tmp_kernel_op_125 * tmp_kernel_op_144) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_146 +
                                  tmp_kernel_op_128 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_131 +
                                 tmp_kernel_op_113 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_205 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_159 +
                                  tmp_kernel_op_116 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_163 +
                                  tmp_kernel_op_119 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_167 +
                                  tmp_kernel_op_122 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_171 +
                                  tmp_kernel_op_125 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_175 +
                                  tmp_kernel_op_128 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_153 +
                                 tmp_kernel_op_113 * tmp_kernel_op_156);
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_82 +
                                  tmp_kernel_op_135 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_88 +
                                  tmp_kernel_op_138 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_94 +
                                  tmp_kernel_op_141 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_143 +
                                  tmp_kernel_op_103 * tmp_kernel_op_144) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_146 +
                                  tmp_kernel_op_109 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_76 +
                                 tmp_kernel_op_132 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_207 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_82 +
                                  tmp_kernel_op_160 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_88 +
                                  tmp_kernel_op_164 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_94 +
                                  tmp_kernel_op_168 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_171 +
                                  tmp_kernel_op_103 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_175 +
                                  tmp_kernel_op_109 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_76 +
                                 tmp_kernel_op_156 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_159 +
                                  tmp_kernel_op_135 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_163 +
                                  tmp_kernel_op_138 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_167 +
                                  tmp_kernel_op_141 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_171 +
                                  tmp_kernel_op_144 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_175 +
                                  tmp_kernel_op_147 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_153 +
                                 tmp_kernel_op_132 * tmp_kernel_op_156);
          const walberla::float64 elMatVec_0 =
              src_dof_0 *
                  (tmp_kernel_op_11 * ((tmp_kernel_op_10 * tmp_kernel_op_10) +
                                       (tmp_kernel_op_9 * tmp_kernel_op_9)) +
                   tmp_kernel_op_17 * ((tmp_kernel_op_15 * tmp_kernel_op_15) +
                                       (tmp_kernel_op_16 * tmp_kernel_op_16)) +
                   tmp_kernel_op_23 * ((tmp_kernel_op_21 * tmp_kernel_op_21) +
                                       (tmp_kernel_op_22 * tmp_kernel_op_22)) +
                   tmp_kernel_op_29 * ((tmp_kernel_op_27 * tmp_kernel_op_27) +
                                       (tmp_kernel_op_28 * tmp_kernel_op_28)) +
                   tmp_kernel_op_35 * ((tmp_kernel_op_33 * tmp_kernel_op_33) +
                                       (tmp_kernel_op_34 * tmp_kernel_op_34)) +
                   tmp_kernel_op_5 * ((tmp_kernel_op_3 * tmp_kernel_op_3) +
                                      (tmp_kernel_op_4 * tmp_kernel_op_4))) +
              src_dof_1 * tmp_kernel_op_54 + src_dof_2 * tmp_kernel_op_73 +
              src_dof_3 * tmp_kernel_op_129 + src_dof_4 * tmp_kernel_op_110 +
              src_dof_5 * tmp_kernel_op_148 + src_dof_6 * tmp_kernel_op_177;
          const walberla::float64 elMatVec_1 =
              src_dof_0 * tmp_kernel_op_54 +
              src_dof_1 *
                  (tmp_kernel_op_11 * (tmp_kernel_op_178 * tmp_kernel_op_181 +
                                       tmp_kernel_op_180 * tmp_kernel_op_181) +
                   tmp_kernel_op_17 * (tmp_kernel_op_178 * tmp_kernel_op_182 +
                                       tmp_kernel_op_180 * tmp_kernel_op_182) +
                   tmp_kernel_op_23 * (tmp_kernel_op_178 * tmp_kernel_op_183 +
                                       tmp_kernel_op_180 * tmp_kernel_op_183) +
                   tmp_kernel_op_29 * (tmp_kernel_op_178 * tmp_kernel_op_184 +
                                       tmp_kernel_op_180 * tmp_kernel_op_184) +
                   tmp_kernel_op_35 * (tmp_kernel_op_178 * tmp_kernel_op_185 +
                                       tmp_kernel_op_180 * tmp_kernel_op_185) +
                   tmp_kernel_op_5 * (tmp_kernel_op_178 * tmp_kernel_op_179 +
                                      tmp_kernel_op_179 * tmp_kernel_op_180)) +
              src_dof_2 * tmp_kernel_op_186 + src_dof_3 * tmp_kernel_op_188 +
              src_dof_4 * tmp_kernel_op_187 + src_dof_5 * tmp_kernel_op_189 +
              src_dof_6 * tmp_kernel_op_190;
          const walberla::float64 elMatVec_2 =
              src_dof_0 * tmp_kernel_op_73 + src_dof_1 * tmp_kernel_op_186 +
              src_dof_2 *
                  (tmp_kernel_op_11 * (tmp_kernel_op_191 * tmp_kernel_op_194 +
                                       tmp_kernel_op_193 * tmp_kernel_op_194) +
                   tmp_kernel_op_17 * (tmp_kernel_op_191 * tmp_kernel_op_195 +
                                       tmp_kernel_op_193 * tmp_kernel_op_195) +
                   tmp_kernel_op_23 * (tmp_kernel_op_191 * tmp_kernel_op_196 +
                                       tmp_kernel_op_193 * tmp_kernel_op_196) +
                   tmp_kernel_op_29 * (tmp_kernel_op_191 * tmp_kernel_op_197 +
                                       tmp_kernel_op_193 * tmp_kernel_op_197) +
                   tmp_kernel_op_35 * (tmp_kernel_op_191 * tmp_kernel_op_198 +
                                       tmp_kernel_op_193 * tmp_kernel_op_198) +
                   tmp_kernel_op_5 * (tmp_kernel_op_191 * tmp_kernel_op_192 +
                                      tmp_kernel_op_192 * tmp_kernel_op_193)) +
              src_dof_3 * tmp_kernel_op_200 + src_dof_4 * tmp_kernel_op_199 +
              src_dof_5 * tmp_kernel_op_201 + src_dof_6 * tmp_kernel_op_202;
          const walberla::float64 elMatVec_3 =
              src_dof_0 * tmp_kernel_op_129 + src_dof_1 * tmp_kernel_op_188 +
              src_dof_2 * tmp_kernel_op_200 +
              src_dof_3 * (tmp_kernel_op_11 *
                               ((tmp_kernel_op_115 * tmp_kernel_op_115) +
                                (tmp_kernel_op_116 * tmp_kernel_op_116)) +
                           tmp_kernel_op_17 *
                               ((tmp_kernel_op_118 * tmp_kernel_op_118) +
                                (tmp_kernel_op_119 * tmp_kernel_op_119)) +
                           tmp_kernel_op_23 *
                               ((tmp_kernel_op_121 * tmp_kernel_op_121) +
                                (tmp_kernel_op_122 * tmp_kernel_op_122)) +
                           tmp_kernel_op_29 *
                               ((tmp_kernel_op_124 * tmp_kernel_op_124) +
                                (tmp_kernel_op_125 * tmp_kernel_op_125)) +
                           tmp_kernel_op_35 *
                               ((tmp_kernel_op_127 * tmp_kernel_op_127) +
                                (tmp_kernel_op_128 * tmp_kernel_op_128)) +
                           tmp_kernel_op_5 *
                               ((tmp_kernel_op_112 * tmp_kernel_op_112) +
                                (tmp_kernel_op_113 * tmp_kernel_op_113))) +
              src_dof_4 * tmp_kernel_op_203 + src_dof_5 * tmp_kernel_op_204 +
              src_dof_6 * tmp_kernel_op_205;
          const walberla::float64 elMatVec_4 =
              src_dof_0 * tmp_kernel_op_110 + src_dof_1 * tmp_kernel_op_187 +
              src_dof_2 * tmp_kernel_op_199 + src_dof_3 * tmp_kernel_op_203 +
              src_dof_4 *
                  (tmp_kernel_op_11 * ((tmp_kernel_op_82 * tmp_kernel_op_82) +
                                       (tmp_kernel_op_85 * tmp_kernel_op_85)) +
                   tmp_kernel_op_17 * ((tmp_kernel_op_88 * tmp_kernel_op_88) +
                                       (tmp_kernel_op_91 * tmp_kernel_op_91)) +
                   tmp_kernel_op_23 * ((tmp_kernel_op_94 * tmp_kernel_op_94) +
                                       (tmp_kernel_op_97 * tmp_kernel_op_97)) +
                   tmp_kernel_op_29 *
                       ((tmp_kernel_op_100 * tmp_kernel_op_100) +
                        (tmp_kernel_op_103 * tmp_kernel_op_103)) +
                   tmp_kernel_op_35 *
                       ((tmp_kernel_op_106 * tmp_kernel_op_106) +
                        (tmp_kernel_op_109 * tmp_kernel_op_109)) +
                   tmp_kernel_op_5 * ((tmp_kernel_op_76 * tmp_kernel_op_76) +
                                      (tmp_kernel_op_79 * tmp_kernel_op_79))) +
              src_dof_5 * tmp_kernel_op_206 + src_dof_6 * tmp_kernel_op_207;
          const walberla::float64 elMatVec_5 =
              src_dof_0 * tmp_kernel_op_148 + src_dof_1 * tmp_kernel_op_189 +
              src_dof_2 * tmp_kernel_op_201 + src_dof_3 * tmp_kernel_op_204 +
              src_dof_4 * tmp_kernel_op_206 +
              src_dof_5 * (tmp_kernel_op_11 *
                               ((tmp_kernel_op_134 * tmp_kernel_op_134) +
                                (tmp_kernel_op_135 * tmp_kernel_op_135)) +
                           tmp_kernel_op_17 *
                               ((tmp_kernel_op_137 * tmp_kernel_op_137) +
                                (tmp_kernel_op_138 * tmp_kernel_op_138)) +
                           tmp_kernel_op_23 *
                               ((tmp_kernel_op_140 * tmp_kernel_op_140) +
                                (tmp_kernel_op_141 * tmp_kernel_op_141)) +
                           tmp_kernel_op_29 *
                               ((tmp_kernel_op_143 * tmp_kernel_op_143) +
                                (tmp_kernel_op_144 * tmp_kernel_op_144)) +
                           tmp_kernel_op_35 *
                               ((tmp_kernel_op_146 * tmp_kernel_op_146) +
                                (tmp_kernel_op_147 * tmp_kernel_op_147)) +
                           tmp_kernel_op_5 *
                               ((tmp_kernel_op_131 * tmp_kernel_op_131) +
                                (tmp_kernel_op_132 * tmp_kernel_op_132))) +
              src_dof_6 * tmp_kernel_op_208;
          const walberla::float64 elMatVec_6 =
              src_dof_0 * tmp_kernel_op_177 + src_dof_1 * tmp_kernel_op_190 +
              src_dof_2 * tmp_kernel_op_202 + src_dof_3 * tmp_kernel_op_205 +
              src_dof_4 * tmp_kernel_op_207 + src_dof_5 * tmp_kernel_op_208 +
              src_dof_6 *
                  (tmp_kernel_op_11 *
                       ((tmp_kernel_op_159 * tmp_kernel_op_159) +
                        (tmp_kernel_op_160 * tmp_kernel_op_160)) +
                   tmp_kernel_op_17 *
                       ((tmp_kernel_op_163 * tmp_kernel_op_163) +
                        (tmp_kernel_op_164 * tmp_kernel_op_164)) +
                   tmp_kernel_op_23 *
                       ((tmp_kernel_op_167 * tmp_kernel_op_167) +
                        (tmp_kernel_op_168 * tmp_kernel_op_168)) +
                   tmp_kernel_op_29 *
                       ((tmp_kernel_op_171 * tmp_kernel_op_171) +
                        (tmp_kernel_op_172 * tmp_kernel_op_172)) +
                   tmp_kernel_op_35 *
                       ((tmp_kernel_op_175 * tmp_kernel_op_175) +
                        (tmp_kernel_op_176 * tmp_kernel_op_176)) +
                   tmp_kernel_op_5 * ((tmp_kernel_op_153 * tmp_kernel_op_153) +
                                      (tmp_kernel_op_156 * tmp_kernel_op_156)));
          _data_dstVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2))] =
              elMatVec_0 +
              _data_dstVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2))];
          _data_dstVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
              elMatVec_1 +
              _data_dstVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          _data_dstVertex[ctr_0 +
                          (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatVec_2 +
              _data_dstVertex[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1)) /
                         (2))] =
              elMatVec_3 +
              _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge + 1)) /
                             (2))];
          _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge + 1)) /
                             (2))] =
              elMatVec_4 +
              _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            2 * ((micro_edges_per_macro_edge *
                                  (micro_edges_per_macro_edge + 1)) /
                                 (2))];
          _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2))] =
              elMatVec_5 +
              _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2))];
          _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                    ((ctr_1 * (ctr_1 + 1)) / (2))] =
              elMatVec_6 +
              _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2))];
        }
    }
    {
      /* FaceType.BLUE */
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
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1;
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
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_0_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_1_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 src_dof_0 =
              _data_srcVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          const walberla::float64 src_dof_1 =
              _data_srcVertex[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          const walberla::float64 src_dof_2 =
              _data_srcVertex[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
          const walberla::float64 src_dof_3 =
              _data_srcEdge[ctr_0 +
                            (ctr_1 + 1) * (micro_edges_per_macro_edge + 1) -
                            (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          const walberla::float64 src_dof_4 =
              _data_srcEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            2 * ((micro_edges_per_macro_edge *
                                  (micro_edges_per_macro_edge + 1)) /
                                 (2)) +
                            1];
          const walberla::float64 src_dof_5 =
              _data_srcEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge + 1)) /
                             (2))];
          const walberla::float64 src_dof_6 =
              _data_src[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1)) /
                         (2))];
          const walberla::float64 tmp_kernel_op_0 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_1 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_2 =
              tmp_kernel_op_0 + tmp_kernel_op_1 - 3.0;
          const walberla::float64 tmp_kernel_op_3 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_2 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_5 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_6 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_7 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_8 =
              tmp_kernel_op_6 + tmp_kernel_op_7 - 3.0;
          const walberla::float64 tmp_kernel_op_9 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_8 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_10 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_8 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_11 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_12 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_13 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_12 + tmp_kernel_op_13 - 3.0;
          const walberla::float64 tmp_kernel_op_15 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_14 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_16 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_14 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_17 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_18 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_19 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_20 =
              tmp_kernel_op_18 + tmp_kernel_op_19 - 3.0;
          const walberla::float64 tmp_kernel_op_21 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_20 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_20 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_23 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_24 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_25 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_26 =
              tmp_kernel_op_24 + tmp_kernel_op_25 - 3.0;
          const walberla::float64 tmp_kernel_op_27 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_26 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_28 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_26 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_29 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_30 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_31 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_30 + tmp_kernel_op_31 - 3.0;
          const walberla::float64 tmp_kernel_op_33 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_32 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_34 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_32 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_35 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_36 = tmp_kernel_op_0 - 1.0;
          const walberla::float64 tmp_kernel_op_37 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_38 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_39 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_40 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_41 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_42 = tmp_kernel_op_12 - 1.0;
          const walberla::float64 tmp_kernel_op_43 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_44 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_45 = tmp_kernel_op_18 - 1.0;
          const walberla::float64 tmp_kernel_op_46 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_47 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_48 = tmp_kernel_op_24 - 1.0;
          const walberla::float64 tmp_kernel_op_49 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_50 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_51 = tmp_kernel_op_30 - 1.0;
          const walberla::float64 tmp_kernel_op_52 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_53 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_41 +
                                  tmp_kernel_op_40 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_43 +
                                  tmp_kernel_op_16 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_46 +
                                  tmp_kernel_op_22 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_49 +
                                  tmp_kernel_op_28 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_52 +
                                  tmp_kernel_op_34 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_37 +
                                 tmp_kernel_op_38 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_55 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_56 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_57 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_58 = tmp_kernel_op_7 - 1.0;
          const walberla::float64 tmp_kernel_op_59 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_60 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_61 = tmp_kernel_op_13 - 1.0;
          const walberla::float64 tmp_kernel_op_62 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_63 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_64 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_65 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_66 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_67 = tmp_kernel_op_25 - 1.0;
          const walberla::float64 tmp_kernel_op_68 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_69 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_70 = tmp_kernel_op_31 - 1.0;
          const walberla::float64 tmp_kernel_op_71 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_72 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_60 +
                                  tmp_kernel_op_59 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_62 +
                                  tmp_kernel_op_16 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_65 +
                                  tmp_kernel_op_22 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_68 +
                                  tmp_kernel_op_28 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_71 +
                                  tmp_kernel_op_34 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_56 +
                                 tmp_kernel_op_4 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_74 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_75 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_74 + tmp_kernel_op_75;
          const walberla::float64 tmp_kernel_op_77 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_78 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_77 + tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_81 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_80 + tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_83 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_84 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_83 + tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_86 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_87 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_86 + tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_90 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_91 =
              tmp_kernel_op_89 + tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_92 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_93 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_92 + tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_96 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_95 + tmp_kernel_op_96;
          const walberla::float64 tmp_kernel_op_98 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_99 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_98 + tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_101 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_102 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_103 =
              tmp_kernel_op_101 + tmp_kernel_op_102;
          const walberla::float64 tmp_kernel_op_104 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_105 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_104 + tmp_kernel_op_105;
          const walberla::float64 tmp_kernel_op_107 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_108 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_109 =
              tmp_kernel_op_107 + tmp_kernel_op_108;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_85 +
                                  tmp_kernel_op_82 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_88 +
                                  tmp_kernel_op_16 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_94 +
                                  tmp_kernel_op_22 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_27 +
                                  tmp_kernel_op_103 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_33 +
                                  tmp_kernel_op_109 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_76 +
                                 tmp_kernel_op_4 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_111 =
              -tmp_kernel_op_1 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_112 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_111 - tmp_kernel_op_75;
          const walberla::float64 tmp_kernel_op_113 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_111 - tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_114 =
              -tmp_kernel_op_7 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_115 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_114 - tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_116 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_114 - tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_117 =
              -tmp_kernel_op_13 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_118 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_117 - tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_119 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_117 - tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_120 =
              -tmp_kernel_op_19 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_121 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_120 - tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_122 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_120 - tmp_kernel_op_96;
          const walberla::float64 tmp_kernel_op_123 =
              -tmp_kernel_op_25 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_124 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_123 - tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_125 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_123 - tmp_kernel_op_102;
          const walberla::float64 tmp_kernel_op_126 =
              -tmp_kernel_op_31 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_127 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_126 - tmp_kernel_op_105;
          const walberla::float64 tmp_kernel_op_128 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_126 - tmp_kernel_op_108;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_116 +
                                  tmp_kernel_op_115 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_15 +
                                  tmp_kernel_op_119 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_21 +
                                  tmp_kernel_op_122 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_27 +
                                  tmp_kernel_op_125 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_33 +
                                  tmp_kernel_op_128 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_3 +
                                 tmp_kernel_op_113 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_0 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_131 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_130 - tmp_kernel_op_74;
          const walberla::float64 tmp_kernel_op_132 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_130 - tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_133 =
              -tmp_kernel_op_6 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_134 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_133 - tmp_kernel_op_80;
          const walberla::float64 tmp_kernel_op_135 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_133 - tmp_kernel_op_83;
          const walberla::float64 tmp_kernel_op_136 =
              -tmp_kernel_op_12 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_137 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_136 - tmp_kernel_op_86;
          const walberla::float64 tmp_kernel_op_138 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_136 - tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_139 =
              -tmp_kernel_op_18 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_140 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_139 - tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_141 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_139 - tmp_kernel_op_95;
          const walberla::float64 tmp_kernel_op_142 =
              -tmp_kernel_op_24 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_143 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_142 - tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_144 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_142 - tmp_kernel_op_101;
          const walberla::float64 tmp_kernel_op_145 =
              -tmp_kernel_op_30 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_146 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_145 - tmp_kernel_op_104;
          const walberla::float64 tmp_kernel_op_147 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_145 - tmp_kernel_op_107;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_135 +
                                  tmp_kernel_op_134 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_15 +
                                  tmp_kernel_op_138 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_21 +
                                  tmp_kernel_op_141 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_27 +
                                  tmp_kernel_op_144 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_33 +
                                  tmp_kernel_op_147 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_3 +
                                 tmp_kernel_op_132 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_149 =
              jac_affine_inv_0_0_BLUE * 27.0;
          const int64_t tmp_kernel_op_150 = 0;
          const walberla::float64 tmp_kernel_op_151 =
              jac_affine_inv_1_0_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_152 = 0.15066167873471437;
          const walberla::float64 tmp_kernel_op_153 =
              tmp_kernel_op_149 * ((walberla::float64)(tmp_kernel_op_150)) +
              tmp_kernel_op_151 * tmp_kernel_op_152;
          const walberla::float64 tmp_kernel_op_154 =
              jac_affine_inv_0_1_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_155 =
              jac_affine_inv_1_1_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_152 * tmp_kernel_op_155 +
              tmp_kernel_op_154 * ((walberla::float64)(tmp_kernel_op_150));
          const walberla::float64 tmp_kernel_op_157 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_158 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_149 * tmp_kernel_op_157 +
              tmp_kernel_op_151 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_154 * tmp_kernel_op_157 +
              tmp_kernel_op_155 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_161 = 0.15066167873471437;
          const int64_t tmp_kernel_op_162 = 0;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_149 * tmp_kernel_op_161 +
              tmp_kernel_op_151 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_154 * tmp_kernel_op_161 +
              tmp_kernel_op_155 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_165 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_166 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_149 * tmp_kernel_op_165 +
              tmp_kernel_op_151 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_154 * tmp_kernel_op_165 +
              tmp_kernel_op_155 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_169 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_170 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_149 * tmp_kernel_op_169 +
              tmp_kernel_op_151 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_154 * tmp_kernel_op_169 +
              tmp_kernel_op_155 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_173 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_174 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_149 * tmp_kernel_op_173 +
              tmp_kernel_op_151 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_154 * tmp_kernel_op_173 +
              tmp_kernel_op_155 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_160 +
                                  tmp_kernel_op_159 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_163 +
                                  tmp_kernel_op_16 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_21 +
                                  tmp_kernel_op_168 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_27 +
                                  tmp_kernel_op_172 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_33 +
                                  tmp_kernel_op_176 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_3 +
                                 tmp_kernel_op_156 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_178 =
              (jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE);
          const walberla::float64 tmp_kernel_op_179 =
              (tmp_kernel_op_36 * tmp_kernel_op_36);
          const walberla::float64 tmp_kernel_op_180 =
              (jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE);
          const walberla::float64 tmp_kernel_op_181 =
              (tmp_kernel_op_39 * tmp_kernel_op_39);
          const walberla::float64 tmp_kernel_op_182 =
              (tmp_kernel_op_42 * tmp_kernel_op_42);
          const walberla::float64 tmp_kernel_op_183 =
              (tmp_kernel_op_45 * tmp_kernel_op_45);
          const walberla::float64 tmp_kernel_op_184 =
              (tmp_kernel_op_48 * tmp_kernel_op_48);
          const walberla::float64 tmp_kernel_op_185 =
              (tmp_kernel_op_51 * tmp_kernel_op_51);
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_59 +
                                  tmp_kernel_op_41 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_62 +
                                  tmp_kernel_op_44 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_65 +
                                  tmp_kernel_op_47 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_49 * tmp_kernel_op_68 +
                                  tmp_kernel_op_50 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_52 * tmp_kernel_op_71 +
                                  tmp_kernel_op_53 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_56 +
                                 tmp_kernel_op_38 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_82 +
                                  tmp_kernel_op_41 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_88 +
                                  tmp_kernel_op_44 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_94 +
                                  tmp_kernel_op_47 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_49 +
                                  tmp_kernel_op_103 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_52 +
                                  tmp_kernel_op_109 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_76 +
                                 tmp_kernel_op_38 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_40 +
                                  tmp_kernel_op_116 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_43 +
                                  tmp_kernel_op_119 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_46 +
                                  tmp_kernel_op_122 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_49 +
                                  tmp_kernel_op_125 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_52 +
                                  tmp_kernel_op_128 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_37 +
                                 tmp_kernel_op_113 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_40 +
                                  tmp_kernel_op_135 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_43 +
                                  tmp_kernel_op_138 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_46 +
                                  tmp_kernel_op_141 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_49 +
                                  tmp_kernel_op_144 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_52 +
                                  tmp_kernel_op_147 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_37 +
                                 tmp_kernel_op_132 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_40 +
                                  tmp_kernel_op_160 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_43 +
                                  tmp_kernel_op_164 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_46 +
                                  tmp_kernel_op_168 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_49 +
                                  tmp_kernel_op_172 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_52 +
                                  tmp_kernel_op_176 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_37 +
                                 tmp_kernel_op_156 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_191 =
              (jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE);
          const walberla::float64 tmp_kernel_op_192 =
              (tmp_kernel_op_55 * tmp_kernel_op_55);
          const walberla::float64 tmp_kernel_op_193 =
              (jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE);
          const walberla::float64 tmp_kernel_op_194 =
              (tmp_kernel_op_58 * tmp_kernel_op_58);
          const walberla::float64 tmp_kernel_op_195 =
              (tmp_kernel_op_61 * tmp_kernel_op_61);
          const walberla::float64 tmp_kernel_op_196 =
              (tmp_kernel_op_64 * tmp_kernel_op_64);
          const walberla::float64 tmp_kernel_op_197 =
              (tmp_kernel_op_67 * tmp_kernel_op_67);
          const walberla::float64 tmp_kernel_op_198 =
              (tmp_kernel_op_70 * tmp_kernel_op_70);
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_11 * (tmp_kernel_op_59 * tmp_kernel_op_82 +
                                  tmp_kernel_op_60 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_62 * tmp_kernel_op_88 +
                                  tmp_kernel_op_63 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_65 * tmp_kernel_op_94 +
                                  tmp_kernel_op_66 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_68 +
                                  tmp_kernel_op_103 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_71 +
                                  tmp_kernel_op_109 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_56 * tmp_kernel_op_76 +
                                 tmp_kernel_op_57 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_59 +
                                  tmp_kernel_op_116 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_62 +
                                  tmp_kernel_op_119 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_65 +
                                  tmp_kernel_op_122 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_68 +
                                  tmp_kernel_op_125 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_71 +
                                  tmp_kernel_op_128 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_56 +
                                 tmp_kernel_op_113 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_201 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_59 +
                                  tmp_kernel_op_135 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_62 +
                                  tmp_kernel_op_138 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_65 +
                                  tmp_kernel_op_141 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_68 +
                                  tmp_kernel_op_144 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_71 +
                                  tmp_kernel_op_147 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_56 +
                                 tmp_kernel_op_132 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_59 +
                                  tmp_kernel_op_160 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_62 +
                                  tmp_kernel_op_164 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_65 +
                                  tmp_kernel_op_168 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_68 +
                                  tmp_kernel_op_172 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_71 +
                                  tmp_kernel_op_176 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_56 +
                                 tmp_kernel_op_156 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_203 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_82 +
                                  tmp_kernel_op_116 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_88 +
                                  tmp_kernel_op_119 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_94 +
                                  tmp_kernel_op_122 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_124 +
                                  tmp_kernel_op_103 * tmp_kernel_op_125) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_127 +
                                  tmp_kernel_op_109 * tmp_kernel_op_128) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_76 +
                                 tmp_kernel_op_113 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_134 +
                                  tmp_kernel_op_116 * tmp_kernel_op_135) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_137 +
                                  tmp_kernel_op_119 * tmp_kernel_op_138) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_140 +
                                  tmp_kernel_op_122 * tmp_kernel_op_141) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_143 +
                                  tmp_kernel_op_125 * tmp_kernel_op_144) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_146 +
                                  tmp_kernel_op_128 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_131 +
                                 tmp_kernel_op_113 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_205 =
              tmp_kernel_op_11 * (tmp_kernel_op_115 * tmp_kernel_op_159 +
                                  tmp_kernel_op_116 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_118 * tmp_kernel_op_163 +
                                  tmp_kernel_op_119 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_121 * tmp_kernel_op_167 +
                                  tmp_kernel_op_122 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_124 * tmp_kernel_op_171 +
                                  tmp_kernel_op_125 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_127 * tmp_kernel_op_175 +
                                  tmp_kernel_op_128 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_112 * tmp_kernel_op_153 +
                                 tmp_kernel_op_113 * tmp_kernel_op_156);
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_82 +
                                  tmp_kernel_op_135 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_88 +
                                  tmp_kernel_op_138 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_94 +
                                  tmp_kernel_op_141 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_143 +
                                  tmp_kernel_op_103 * tmp_kernel_op_144) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_146 +
                                  tmp_kernel_op_109 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_76 +
                                 tmp_kernel_op_132 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_207 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_82 +
                                  tmp_kernel_op_160 * tmp_kernel_op_85) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_88 +
                                  tmp_kernel_op_164 * tmp_kernel_op_91) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_94 +
                                  tmp_kernel_op_168 * tmp_kernel_op_97) +
              tmp_kernel_op_29 * (tmp_kernel_op_100 * tmp_kernel_op_171 +
                                  tmp_kernel_op_103 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_106 * tmp_kernel_op_175 +
                                  tmp_kernel_op_109 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_76 +
                                 tmp_kernel_op_156 * tmp_kernel_op_79);
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_159 +
                                  tmp_kernel_op_135 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_163 +
                                  tmp_kernel_op_138 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_167 +
                                  tmp_kernel_op_141 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_171 +
                                  tmp_kernel_op_144 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_175 +
                                  tmp_kernel_op_147 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_153 +
                                 tmp_kernel_op_132 * tmp_kernel_op_156);
          const walberla::float64 elMatVec_0 =
              src_dof_0 *
                  (tmp_kernel_op_11 * ((tmp_kernel_op_10 * tmp_kernel_op_10) +
                                       (tmp_kernel_op_9 * tmp_kernel_op_9)) +
                   tmp_kernel_op_17 * ((tmp_kernel_op_15 * tmp_kernel_op_15) +
                                       (tmp_kernel_op_16 * tmp_kernel_op_16)) +
                   tmp_kernel_op_23 * ((tmp_kernel_op_21 * tmp_kernel_op_21) +
                                       (tmp_kernel_op_22 * tmp_kernel_op_22)) +
                   tmp_kernel_op_29 * ((tmp_kernel_op_27 * tmp_kernel_op_27) +
                                       (tmp_kernel_op_28 * tmp_kernel_op_28)) +
                   tmp_kernel_op_35 * ((tmp_kernel_op_33 * tmp_kernel_op_33) +
                                       (tmp_kernel_op_34 * tmp_kernel_op_34)) +
                   tmp_kernel_op_5 * ((tmp_kernel_op_3 * tmp_kernel_op_3) +
                                      (tmp_kernel_op_4 * tmp_kernel_op_4))) +
              src_dof_1 * tmp_kernel_op_54 + src_dof_2 * tmp_kernel_op_73 +
              src_dof_3 * tmp_kernel_op_129 + src_dof_4 * tmp_kernel_op_110 +
              src_dof_5 * tmp_kernel_op_148 + src_dof_6 * tmp_kernel_op_177;
          const walberla::float64 elMatVec_1 =
              src_dof_0 * tmp_kernel_op_54 +
              src_dof_1 *
                  (tmp_kernel_op_11 * (tmp_kernel_op_178 * tmp_kernel_op_181 +
                                       tmp_kernel_op_180 * tmp_kernel_op_181) +
                   tmp_kernel_op_17 * (tmp_kernel_op_178 * tmp_kernel_op_182 +
                                       tmp_kernel_op_180 * tmp_kernel_op_182) +
                   tmp_kernel_op_23 * (tmp_kernel_op_178 * tmp_kernel_op_183 +
                                       tmp_kernel_op_180 * tmp_kernel_op_183) +
                   tmp_kernel_op_29 * (tmp_kernel_op_178 * tmp_kernel_op_184 +
                                       tmp_kernel_op_180 * tmp_kernel_op_184) +
                   tmp_kernel_op_35 * (tmp_kernel_op_178 * tmp_kernel_op_185 +
                                       tmp_kernel_op_180 * tmp_kernel_op_185) +
                   tmp_kernel_op_5 * (tmp_kernel_op_178 * tmp_kernel_op_179 +
                                      tmp_kernel_op_179 * tmp_kernel_op_180)) +
              src_dof_2 * tmp_kernel_op_186 + src_dof_3 * tmp_kernel_op_188 +
              src_dof_4 * tmp_kernel_op_187 + src_dof_5 * tmp_kernel_op_189 +
              src_dof_6 * tmp_kernel_op_190;
          const walberla::float64 elMatVec_2 =
              src_dof_0 * tmp_kernel_op_73 + src_dof_1 * tmp_kernel_op_186 +
              src_dof_2 *
                  (tmp_kernel_op_11 * (tmp_kernel_op_191 * tmp_kernel_op_194 +
                                       tmp_kernel_op_193 * tmp_kernel_op_194) +
                   tmp_kernel_op_17 * (tmp_kernel_op_191 * tmp_kernel_op_195 +
                                       tmp_kernel_op_193 * tmp_kernel_op_195) +
                   tmp_kernel_op_23 * (tmp_kernel_op_191 * tmp_kernel_op_196 +
                                       tmp_kernel_op_193 * tmp_kernel_op_196) +
                   tmp_kernel_op_29 * (tmp_kernel_op_191 * tmp_kernel_op_197 +
                                       tmp_kernel_op_193 * tmp_kernel_op_197) +
                   tmp_kernel_op_35 * (tmp_kernel_op_191 * tmp_kernel_op_198 +
                                       tmp_kernel_op_193 * tmp_kernel_op_198) +
                   tmp_kernel_op_5 * (tmp_kernel_op_191 * tmp_kernel_op_192 +
                                      tmp_kernel_op_192 * tmp_kernel_op_193)) +
              src_dof_3 * tmp_kernel_op_200 + src_dof_4 * tmp_kernel_op_199 +
              src_dof_5 * tmp_kernel_op_201 + src_dof_6 * tmp_kernel_op_202;
          const walberla::float64 elMatVec_3 =
              src_dof_0 * tmp_kernel_op_129 + src_dof_1 * tmp_kernel_op_188 +
              src_dof_2 * tmp_kernel_op_200 +
              src_dof_3 * (tmp_kernel_op_11 *
                               ((tmp_kernel_op_115 * tmp_kernel_op_115) +
                                (tmp_kernel_op_116 * tmp_kernel_op_116)) +
                           tmp_kernel_op_17 *
                               ((tmp_kernel_op_118 * tmp_kernel_op_118) +
                                (tmp_kernel_op_119 * tmp_kernel_op_119)) +
                           tmp_kernel_op_23 *
                               ((tmp_kernel_op_121 * tmp_kernel_op_121) +
                                (tmp_kernel_op_122 * tmp_kernel_op_122)) +
                           tmp_kernel_op_29 *
                               ((tmp_kernel_op_124 * tmp_kernel_op_124) +
                                (tmp_kernel_op_125 * tmp_kernel_op_125)) +
                           tmp_kernel_op_35 *
                               ((tmp_kernel_op_127 * tmp_kernel_op_127) +
                                (tmp_kernel_op_128 * tmp_kernel_op_128)) +
                           tmp_kernel_op_5 *
                               ((tmp_kernel_op_112 * tmp_kernel_op_112) +
                                (tmp_kernel_op_113 * tmp_kernel_op_113))) +
              src_dof_4 * tmp_kernel_op_203 + src_dof_5 * tmp_kernel_op_204 +
              src_dof_6 * tmp_kernel_op_205;
          const walberla::float64 elMatVec_4 =
              src_dof_0 * tmp_kernel_op_110 + src_dof_1 * tmp_kernel_op_187 +
              src_dof_2 * tmp_kernel_op_199 + src_dof_3 * tmp_kernel_op_203 +
              src_dof_4 *
                  (tmp_kernel_op_11 * ((tmp_kernel_op_82 * tmp_kernel_op_82) +
                                       (tmp_kernel_op_85 * tmp_kernel_op_85)) +
                   tmp_kernel_op_17 * ((tmp_kernel_op_88 * tmp_kernel_op_88) +
                                       (tmp_kernel_op_91 * tmp_kernel_op_91)) +
                   tmp_kernel_op_23 * ((tmp_kernel_op_94 * tmp_kernel_op_94) +
                                       (tmp_kernel_op_97 * tmp_kernel_op_97)) +
                   tmp_kernel_op_29 *
                       ((tmp_kernel_op_100 * tmp_kernel_op_100) +
                        (tmp_kernel_op_103 * tmp_kernel_op_103)) +
                   tmp_kernel_op_35 *
                       ((tmp_kernel_op_106 * tmp_kernel_op_106) +
                        (tmp_kernel_op_109 * tmp_kernel_op_109)) +
                   tmp_kernel_op_5 * ((tmp_kernel_op_76 * tmp_kernel_op_76) +
                                      (tmp_kernel_op_79 * tmp_kernel_op_79))) +
              src_dof_5 * tmp_kernel_op_206 + src_dof_6 * tmp_kernel_op_207;
          const walberla::float64 elMatVec_5 =
              src_dof_0 * tmp_kernel_op_148 + src_dof_1 * tmp_kernel_op_189 +
              src_dof_2 * tmp_kernel_op_201 + src_dof_3 * tmp_kernel_op_204 +
              src_dof_4 * tmp_kernel_op_206 +
              src_dof_5 * (tmp_kernel_op_11 *
                               ((tmp_kernel_op_134 * tmp_kernel_op_134) +
                                (tmp_kernel_op_135 * tmp_kernel_op_135)) +
                           tmp_kernel_op_17 *
                               ((tmp_kernel_op_137 * tmp_kernel_op_137) +
                                (tmp_kernel_op_138 * tmp_kernel_op_138)) +
                           tmp_kernel_op_23 *
                               ((tmp_kernel_op_140 * tmp_kernel_op_140) +
                                (tmp_kernel_op_141 * tmp_kernel_op_141)) +
                           tmp_kernel_op_29 *
                               ((tmp_kernel_op_143 * tmp_kernel_op_143) +
                                (tmp_kernel_op_144 * tmp_kernel_op_144)) +
                           tmp_kernel_op_35 *
                               ((tmp_kernel_op_146 * tmp_kernel_op_146) +
                                (tmp_kernel_op_147 * tmp_kernel_op_147)) +
                           tmp_kernel_op_5 *
                               ((tmp_kernel_op_131 * tmp_kernel_op_131) +
                                (tmp_kernel_op_132 * tmp_kernel_op_132))) +
              src_dof_6 * tmp_kernel_op_208;
          const walberla::float64 elMatVec_6 =
              src_dof_0 * tmp_kernel_op_177 + src_dof_1 * tmp_kernel_op_190 +
              src_dof_2 * tmp_kernel_op_202 + src_dof_3 * tmp_kernel_op_205 +
              src_dof_4 * tmp_kernel_op_207 + src_dof_5 * tmp_kernel_op_208 +
              src_dof_6 *
                  (tmp_kernel_op_11 *
                       ((tmp_kernel_op_159 * tmp_kernel_op_159) +
                        (tmp_kernel_op_160 * tmp_kernel_op_160)) +
                   tmp_kernel_op_17 *
                       ((tmp_kernel_op_163 * tmp_kernel_op_163) +
                        (tmp_kernel_op_164 * tmp_kernel_op_164)) +
                   tmp_kernel_op_23 *
                       ((tmp_kernel_op_167 * tmp_kernel_op_167) +
                        (tmp_kernel_op_168 * tmp_kernel_op_168)) +
                   tmp_kernel_op_29 *
                       ((tmp_kernel_op_171 * tmp_kernel_op_171) +
                        (tmp_kernel_op_172 * tmp_kernel_op_172)) +
                   tmp_kernel_op_35 *
                       ((tmp_kernel_op_175 * tmp_kernel_op_175) +
                        (tmp_kernel_op_176 * tmp_kernel_op_176)) +
                   tmp_kernel_op_5 * ((tmp_kernel_op_153 * tmp_kernel_op_153) +
                                      (tmp_kernel_op_156 * tmp_kernel_op_156)));
          _data_dstVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                          ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
              elMatVec_0 +
              _data_dstVertex[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                              ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          _data_dstVertex[ctr_0 +
                          (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatVec_1 +
              _data_dstVertex[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_dstVertex[ctr_0 +
                          (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1] =
              elMatVec_2 +
              _data_dstVertex[ctr_0 +
                              (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                              (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
          _data_dstEdge[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 1) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatVec_3 +
              _data_dstEdge[ctr_0 +
                            (ctr_1 + 1) * (micro_edges_per_macro_edge + 1) -
                            (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        2 * ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge + 1)) /
                             (2)) +
                        1] =
              elMatVec_4 +
              _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            2 * ((micro_edges_per_macro_edge *
                                  (micro_edges_per_macro_edge + 1)) /
                                 (2)) +
                            1];
          _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1)) /
                         (2))] =
              elMatVec_5 +
              _data_dstEdge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                            ((ctr_1 * (ctr_1 + 1)) / (2)) +
                            ((micro_edges_per_macro_edge *
                              (micro_edges_per_macro_edge + 1)) /
                             (2))];
          _data_dst[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                    ((ctr_1 * (ctr_1 + 1)) / (2)) +
                    ((micro_edges_per_macro_edge *
                      (micro_edges_per_macro_edge + 1)) /
                     (2))] =
              elMatVec_6 +
              _data_dst[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) +
                        ((micro_edges_per_macro_edge *
                          (micro_edges_per_macro_edge + 1)) /
                         (2))];
        }
    }
  }
}
void P2PlusBubbleElementwiseDiffusion_float64::
    toMatrix_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(
        idx_t *RESTRICT _data_dst, idx_t *RESTRICT _data_dstEdge,
        idx_t *RESTRICT _data_dstVertex, idx_t *RESTRICT _data_src,
        idx_t *RESTRICT _data_srcEdge, idx_t *RESTRICT _data_srcVertex,
        walberla::float64 macro_vertex_coord_id_0comp0,
        walberla::float64 macro_vertex_coord_id_0comp1,
        walberla::float64 macro_vertex_coord_id_1comp0,
        walberla::float64 macro_vertex_coord_id_1comp1,
        walberla::float64 macro_vertex_coord_id_2comp0,
        walberla::float64 macro_vertex_coord_id_2comp1,
        std::shared_ptr<SparseMatrixProxy> mat,
        int64_t micro_edges_per_macro_edge,
        walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    {
      /* FaceType.GRAY */
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
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge;
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
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_0_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_2_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 tmp_kernel_op_0 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_1 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_2 =
              tmp_kernel_op_0 + tmp_kernel_op_1 - 3.0;
          const walberla::float64 tmp_kernel_op_3 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_2 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_5 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_6 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_7 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_8 =
              tmp_kernel_op_6 + tmp_kernel_op_7 - 3.0;
          const walberla::float64 tmp_kernel_op_9 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_8 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_10 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_8 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_11 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_12 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_13 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_12 + tmp_kernel_op_13 - 3.0;
          const walberla::float64 tmp_kernel_op_15 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_14 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_16 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_14 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_17 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_18 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_19 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_20 =
              tmp_kernel_op_18 + tmp_kernel_op_19 - 3.0;
          const walberla::float64 tmp_kernel_op_21 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_20 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_20 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_23 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_24 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_25 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_26 =
              tmp_kernel_op_24 + tmp_kernel_op_25 - 3.0;
          const walberla::float64 tmp_kernel_op_27 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_26 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_28 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_26 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_29 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_30 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_31 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_30 + tmp_kernel_op_31 - 3.0;
          const walberla::float64 tmp_kernel_op_33 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_32 +
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_34 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_32 +
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_35 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_36 = tmp_kernel_op_0 - 1.0;
          const walberla::float64 tmp_kernel_op_37 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_38 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_39 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_40 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_41 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_42 = tmp_kernel_op_12 - 1.0;
          const walberla::float64 tmp_kernel_op_43 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_44 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_45 = tmp_kernel_op_18 - 1.0;
          const walberla::float64 tmp_kernel_op_46 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_47 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_48 = tmp_kernel_op_24 - 1.0;
          const walberla::float64 tmp_kernel_op_49 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_50 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_51 = tmp_kernel_op_30 - 1.0;
          const walberla::float64 tmp_kernel_op_52 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_53 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_41 +
                                  tmp_kernel_op_40 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_43 +
                                  tmp_kernel_op_16 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_46 +
                                  tmp_kernel_op_22 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_49 +
                                  tmp_kernel_op_28 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_52 +
                                  tmp_kernel_op_34 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_37 +
                                 tmp_kernel_op_38 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_55 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_56 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_57 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_58 = tmp_kernel_op_7 - 1.0;
          const walberla::float64 tmp_kernel_op_59 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_60 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_61 = tmp_kernel_op_13 - 1.0;
          const walberla::float64 tmp_kernel_op_62 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_63 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_64 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_65 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_66 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_67 = tmp_kernel_op_25 - 1.0;
          const walberla::float64 tmp_kernel_op_68 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_69 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_70 = tmp_kernel_op_31 - 1.0;
          const walberla::float64 tmp_kernel_op_71 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_72 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_60 +
                                  tmp_kernel_op_59 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_62 +
                                  tmp_kernel_op_16 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_65 +
                                  tmp_kernel_op_22 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_68 +
                                  tmp_kernel_op_28 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_71 +
                                  tmp_kernel_op_34 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_56 +
                                 tmp_kernel_op_4 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_74 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_75 =
              -tmp_kernel_op_1 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_76 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_75 - tmp_kernel_op_74;
          const walberla::float64 tmp_kernel_op_77 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_78 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_75 - tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_79 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_80 =
              -tmp_kernel_op_7 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_81 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_80 - tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_82 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_83 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_80 - tmp_kernel_op_82;
          const walberla::float64 tmp_kernel_op_84 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_85 =
              -tmp_kernel_op_13 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_86 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_85 - tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_87 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_88 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_85 - tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_90 =
              -tmp_kernel_op_19 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_91 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_90 - tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_92 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_93 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_90 - tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_94 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_95 =
              -tmp_kernel_op_25 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_96 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_95 - tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_97 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_98 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_95 - tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_99 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_100 =
              -tmp_kernel_op_31 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_101 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_100 - tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_102 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_103 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_100 - tmp_kernel_op_102;
          const walberla::float64 tmp_kernel_op_104 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_83 +
                                  tmp_kernel_op_81 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_86 +
                                  tmp_kernel_op_16 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_91 +
                                  tmp_kernel_op_22 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_96 +
                                  tmp_kernel_op_28 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_33 +
                                  tmp_kernel_op_103 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_76 +
                                 tmp_kernel_op_4 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_105 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_105 + tmp_kernel_op_74;
          const walberla::float64 tmp_kernel_op_107 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_108 =
              tmp_kernel_op_107 + tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_109 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_109 + tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_111 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_112 =
              tmp_kernel_op_111 + tmp_kernel_op_82;
          const walberla::float64 tmp_kernel_op_113 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_114 =
              tmp_kernel_op_113 + tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_115 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_116 =
              tmp_kernel_op_115 + tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_117 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_118 =
              tmp_kernel_op_117 + tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_119 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_120 =
              tmp_kernel_op_119 + tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_121 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_122 =
              tmp_kernel_op_121 + tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_123 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_124 =
              tmp_kernel_op_123 + tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_125 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_126 =
              tmp_kernel_op_125 + tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_127 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_128 =
              tmp_kernel_op_102 + tmp_kernel_op_127;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_112 +
                                  tmp_kernel_op_110 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_15 +
                                  tmp_kernel_op_116 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_21 +
                                  tmp_kernel_op_120 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_27 +
                                  tmp_kernel_op_124 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_33 +
                                  tmp_kernel_op_128 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_3 +
                                 tmp_kernel_op_108 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_0 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_131 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_130 - tmp_kernel_op_105;
          const walberla::float64 tmp_kernel_op_132 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_130 - tmp_kernel_op_107;
          const walberla::float64 tmp_kernel_op_133 =
              -tmp_kernel_op_6 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_134 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_133 - tmp_kernel_op_109;
          const walberla::float64 tmp_kernel_op_135 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_133 - tmp_kernel_op_111;
          const walberla::float64 tmp_kernel_op_136 =
              -tmp_kernel_op_12 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_137 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_136 - tmp_kernel_op_113;
          const walberla::float64 tmp_kernel_op_138 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_136 - tmp_kernel_op_115;
          const walberla::float64 tmp_kernel_op_139 =
              -tmp_kernel_op_18 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_140 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_139 - tmp_kernel_op_117;
          const walberla::float64 tmp_kernel_op_141 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_139 - tmp_kernel_op_119;
          const walberla::float64 tmp_kernel_op_142 =
              -tmp_kernel_op_24 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_143 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_142 - tmp_kernel_op_121;
          const walberla::float64 tmp_kernel_op_144 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_142 - tmp_kernel_op_123;
          const walberla::float64 tmp_kernel_op_145 =
              -tmp_kernel_op_30 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_146 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_145 - tmp_kernel_op_125;
          const walberla::float64 tmp_kernel_op_147 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_145 - tmp_kernel_op_127;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_135 +
                                  tmp_kernel_op_134 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_15 +
                                  tmp_kernel_op_138 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_21 +
                                  tmp_kernel_op_141 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_27 +
                                  tmp_kernel_op_144 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_33 +
                                  tmp_kernel_op_147 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_3 +
                                 tmp_kernel_op_132 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_149 =
              jac_affine_inv_0_0_GRAY * 27.0;
          const int64_t tmp_kernel_op_150 = 0;
          const walberla::float64 tmp_kernel_op_151 =
              jac_affine_inv_1_0_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_152 = 0.15066167873471437;
          const walberla::float64 tmp_kernel_op_153 =
              tmp_kernel_op_149 * ((walberla::float64)(tmp_kernel_op_150)) +
              tmp_kernel_op_151 * tmp_kernel_op_152;
          const walberla::float64 tmp_kernel_op_154 =
              jac_affine_inv_0_1_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_155 =
              jac_affine_inv_1_1_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_152 * tmp_kernel_op_155 +
              tmp_kernel_op_154 * ((walberla::float64)(tmp_kernel_op_150));
          const walberla::float64 tmp_kernel_op_157 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_158 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_149 * tmp_kernel_op_157 +
              tmp_kernel_op_151 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_154 * tmp_kernel_op_157 +
              tmp_kernel_op_155 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_161 = 0.15066167873471437;
          const int64_t tmp_kernel_op_162 = 0;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_149 * tmp_kernel_op_161 +
              tmp_kernel_op_151 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_154 * tmp_kernel_op_161 +
              tmp_kernel_op_155 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_165 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_166 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_149 * tmp_kernel_op_165 +
              tmp_kernel_op_151 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_154 * tmp_kernel_op_165 +
              tmp_kernel_op_155 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_169 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_170 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_149 * tmp_kernel_op_169 +
              tmp_kernel_op_151 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_154 * tmp_kernel_op_169 +
              tmp_kernel_op_155 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_173 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_174 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_149 * tmp_kernel_op_173 +
              tmp_kernel_op_151 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_154 * tmp_kernel_op_173 +
              tmp_kernel_op_155 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_160 +
                                  tmp_kernel_op_159 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_163 +
                                  tmp_kernel_op_16 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_21 +
                                  tmp_kernel_op_168 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_27 +
                                  tmp_kernel_op_172 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_33 +
                                  tmp_kernel_op_176 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_3 +
                                 tmp_kernel_op_156 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_178 =
              (jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY);
          const walberla::float64 tmp_kernel_op_179 =
              (tmp_kernel_op_36 * tmp_kernel_op_36);
          const walberla::float64 tmp_kernel_op_180 =
              (jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY);
          const walberla::float64 tmp_kernel_op_181 =
              (tmp_kernel_op_39 * tmp_kernel_op_39);
          const walberla::float64 tmp_kernel_op_182 =
              (tmp_kernel_op_42 * tmp_kernel_op_42);
          const walberla::float64 tmp_kernel_op_183 =
              (tmp_kernel_op_45 * tmp_kernel_op_45);
          const walberla::float64 tmp_kernel_op_184 =
              (tmp_kernel_op_48 * tmp_kernel_op_48);
          const walberla::float64 tmp_kernel_op_185 =
              (tmp_kernel_op_51 * tmp_kernel_op_51);
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_59 +
                                  tmp_kernel_op_41 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_62 +
                                  tmp_kernel_op_44 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_65 +
                                  tmp_kernel_op_47 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_49 * tmp_kernel_op_68 +
                                  tmp_kernel_op_50 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_52 * tmp_kernel_op_71 +
                                  tmp_kernel_op_53 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_56 +
                                 tmp_kernel_op_38 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_81 +
                                  tmp_kernel_op_41 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_86 +
                                  tmp_kernel_op_44 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_91 +
                                  tmp_kernel_op_47 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_49 * tmp_kernel_op_96 +
                                  tmp_kernel_op_50 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_52 +
                                  tmp_kernel_op_103 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_76 +
                                 tmp_kernel_op_38 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_40 +
                                  tmp_kernel_op_112 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_43 +
                                  tmp_kernel_op_116 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_46 +
                                  tmp_kernel_op_120 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_49 +
                                  tmp_kernel_op_124 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_52 +
                                  tmp_kernel_op_128 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_37 +
                                 tmp_kernel_op_108 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_40 +
                                  tmp_kernel_op_135 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_43 +
                                  tmp_kernel_op_138 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_46 +
                                  tmp_kernel_op_141 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_49 +
                                  tmp_kernel_op_144 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_52 +
                                  tmp_kernel_op_147 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_37 +
                                 tmp_kernel_op_132 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_40 +
                                  tmp_kernel_op_160 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_43 +
                                  tmp_kernel_op_164 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_46 +
                                  tmp_kernel_op_168 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_49 +
                                  tmp_kernel_op_172 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_52 +
                                  tmp_kernel_op_176 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_37 +
                                 tmp_kernel_op_156 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_191 =
              (jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY);
          const walberla::float64 tmp_kernel_op_192 =
              (tmp_kernel_op_55 * tmp_kernel_op_55);
          const walberla::float64 tmp_kernel_op_193 =
              (jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY);
          const walberla::float64 tmp_kernel_op_194 =
              (tmp_kernel_op_58 * tmp_kernel_op_58);
          const walberla::float64 tmp_kernel_op_195 =
              (tmp_kernel_op_61 * tmp_kernel_op_61);
          const walberla::float64 tmp_kernel_op_196 =
              (tmp_kernel_op_64 * tmp_kernel_op_64);
          const walberla::float64 tmp_kernel_op_197 =
              (tmp_kernel_op_67 * tmp_kernel_op_67);
          const walberla::float64 tmp_kernel_op_198 =
              (tmp_kernel_op_70 * tmp_kernel_op_70);
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_11 * (tmp_kernel_op_59 * tmp_kernel_op_81 +
                                  tmp_kernel_op_60 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_62 * tmp_kernel_op_86 +
                                  tmp_kernel_op_63 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_65 * tmp_kernel_op_91 +
                                  tmp_kernel_op_66 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_68 * tmp_kernel_op_96 +
                                  tmp_kernel_op_69 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_71 +
                                  tmp_kernel_op_103 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_56 * tmp_kernel_op_76 +
                                 tmp_kernel_op_57 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_59 +
                                  tmp_kernel_op_112 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_62 +
                                  tmp_kernel_op_116 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_65 +
                                  tmp_kernel_op_120 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_68 +
                                  tmp_kernel_op_124 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_71 +
                                  tmp_kernel_op_128 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_56 +
                                 tmp_kernel_op_108 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_201 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_59 +
                                  tmp_kernel_op_135 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_62 +
                                  tmp_kernel_op_138 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_65 +
                                  tmp_kernel_op_141 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_68 +
                                  tmp_kernel_op_144 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_71 +
                                  tmp_kernel_op_147 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_56 +
                                 tmp_kernel_op_132 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_59 +
                                  tmp_kernel_op_160 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_62 +
                                  tmp_kernel_op_164 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_65 +
                                  tmp_kernel_op_168 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_68 +
                                  tmp_kernel_op_172 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_71 +
                                  tmp_kernel_op_176 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_56 +
                                 tmp_kernel_op_156 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_203 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_81 +
                                  tmp_kernel_op_112 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_86 +
                                  tmp_kernel_op_116 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_91 +
                                  tmp_kernel_op_120 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_96 +
                                  tmp_kernel_op_124 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_126 +
                                  tmp_kernel_op_103 * tmp_kernel_op_128) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_76 +
                                 tmp_kernel_op_108 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_81 +
                                  tmp_kernel_op_135 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_86 +
                                  tmp_kernel_op_138 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_91 +
                                  tmp_kernel_op_141 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_96 +
                                  tmp_kernel_op_144 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_146 +
                                  tmp_kernel_op_103 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_76 +
                                 tmp_kernel_op_132 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_205 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_81 +
                                  tmp_kernel_op_160 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_86 +
                                  tmp_kernel_op_164 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_91 +
                                  tmp_kernel_op_168 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_96 +
                                  tmp_kernel_op_172 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_175 +
                                  tmp_kernel_op_103 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_76 +
                                 tmp_kernel_op_156 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_134 +
                                  tmp_kernel_op_112 * tmp_kernel_op_135) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_137 +
                                  tmp_kernel_op_116 * tmp_kernel_op_138) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_140 +
                                  tmp_kernel_op_120 * tmp_kernel_op_141) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_143 +
                                  tmp_kernel_op_124 * tmp_kernel_op_144) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_146 +
                                  tmp_kernel_op_128 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_131 +
                                 tmp_kernel_op_108 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_207 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_159 +
                                  tmp_kernel_op_112 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_163 +
                                  tmp_kernel_op_116 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_167 +
                                  tmp_kernel_op_120 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_171 +
                                  tmp_kernel_op_124 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_175 +
                                  tmp_kernel_op_128 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_153 +
                                 tmp_kernel_op_108 * tmp_kernel_op_156);
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_159 +
                                  tmp_kernel_op_135 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_163 +
                                  tmp_kernel_op_138 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_167 +
                                  tmp_kernel_op_141 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_171 +
                                  tmp_kernel_op_144 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_175 +
                                  tmp_kernel_op_147 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_153 +
                                 tmp_kernel_op_132 * tmp_kernel_op_156);
          const walberla::float64 elMat_0_0 =
              tmp_kernel_op_11 * ((tmp_kernel_op_10 * tmp_kernel_op_10) +
                                  (tmp_kernel_op_9 * tmp_kernel_op_9)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_15 * tmp_kernel_op_15) +
                                  (tmp_kernel_op_16 * tmp_kernel_op_16)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_21 * tmp_kernel_op_21) +
                                  (tmp_kernel_op_22 * tmp_kernel_op_22)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_27 * tmp_kernel_op_27) +
                                  (tmp_kernel_op_28 * tmp_kernel_op_28)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_33 * tmp_kernel_op_33) +
                                  (tmp_kernel_op_34 * tmp_kernel_op_34)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_3 * tmp_kernel_op_3) +
                                 (tmp_kernel_op_4 * tmp_kernel_op_4));
          const walberla::float64 elMat_0_1 = tmp_kernel_op_54;
          const walberla::float64 elMat_0_2 = tmp_kernel_op_73;
          const walberla::float64 elMat_0_3 = tmp_kernel_op_104;
          const walberla::float64 elMat_0_4 = tmp_kernel_op_129;
          const walberla::float64 elMat_0_5 = tmp_kernel_op_148;
          const walberla::float64 elMat_0_6 = tmp_kernel_op_177;
          const walberla::float64 elMat_1_0 = tmp_kernel_op_54;
          const walberla::float64 elMat_1_1 =
              tmp_kernel_op_11 * (tmp_kernel_op_178 * tmp_kernel_op_181 +
                                  tmp_kernel_op_180 * tmp_kernel_op_181) +
              tmp_kernel_op_17 * (tmp_kernel_op_178 * tmp_kernel_op_182 +
                                  tmp_kernel_op_180 * tmp_kernel_op_182) +
              tmp_kernel_op_23 * (tmp_kernel_op_178 * tmp_kernel_op_183 +
                                  tmp_kernel_op_180 * tmp_kernel_op_183) +
              tmp_kernel_op_29 * (tmp_kernel_op_178 * tmp_kernel_op_184 +
                                  tmp_kernel_op_180 * tmp_kernel_op_184) +
              tmp_kernel_op_35 * (tmp_kernel_op_178 * tmp_kernel_op_185 +
                                  tmp_kernel_op_180 * tmp_kernel_op_185) +
              tmp_kernel_op_5 * (tmp_kernel_op_178 * tmp_kernel_op_179 +
                                 tmp_kernel_op_179 * tmp_kernel_op_180);
          const walberla::float64 elMat_1_2 = tmp_kernel_op_186;
          const walberla::float64 elMat_1_3 = tmp_kernel_op_187;
          const walberla::float64 elMat_1_4 = tmp_kernel_op_188;
          const walberla::float64 elMat_1_5 = tmp_kernel_op_189;
          const walberla::float64 elMat_1_6 = tmp_kernel_op_190;
          const walberla::float64 elMat_2_0 = tmp_kernel_op_73;
          const walberla::float64 elMat_2_1 = tmp_kernel_op_186;
          const walberla::float64 elMat_2_2 =
              tmp_kernel_op_11 * (tmp_kernel_op_191 * tmp_kernel_op_194 +
                                  tmp_kernel_op_193 * tmp_kernel_op_194) +
              tmp_kernel_op_17 * (tmp_kernel_op_191 * tmp_kernel_op_195 +
                                  tmp_kernel_op_193 * tmp_kernel_op_195) +
              tmp_kernel_op_23 * (tmp_kernel_op_191 * tmp_kernel_op_196 +
                                  tmp_kernel_op_193 * tmp_kernel_op_196) +
              tmp_kernel_op_29 * (tmp_kernel_op_191 * tmp_kernel_op_197 +
                                  tmp_kernel_op_193 * tmp_kernel_op_197) +
              tmp_kernel_op_35 * (tmp_kernel_op_191 * tmp_kernel_op_198 +
                                  tmp_kernel_op_193 * tmp_kernel_op_198) +
              tmp_kernel_op_5 * (tmp_kernel_op_191 * tmp_kernel_op_192 +
                                 tmp_kernel_op_192 * tmp_kernel_op_193);
          const walberla::float64 elMat_2_3 = tmp_kernel_op_199;
          const walberla::float64 elMat_2_4 = tmp_kernel_op_200;
          const walberla::float64 elMat_2_5 = tmp_kernel_op_201;
          const walberla::float64 elMat_2_6 = tmp_kernel_op_202;
          const walberla::float64 elMat_3_0 = tmp_kernel_op_104;
          const walberla::float64 elMat_3_1 = tmp_kernel_op_187;
          const walberla::float64 elMat_3_2 = tmp_kernel_op_199;
          const walberla::float64 elMat_3_3 =
              tmp_kernel_op_11 * ((tmp_kernel_op_81 * tmp_kernel_op_81) +
                                  (tmp_kernel_op_83 * tmp_kernel_op_83)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_86 * tmp_kernel_op_86) +
                                  (tmp_kernel_op_88 * tmp_kernel_op_88)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_91 * tmp_kernel_op_91) +
                                  (tmp_kernel_op_93 * tmp_kernel_op_93)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_96 * tmp_kernel_op_96) +
                                  (tmp_kernel_op_98 * tmp_kernel_op_98)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_101 * tmp_kernel_op_101) +
                                  (tmp_kernel_op_103 * tmp_kernel_op_103)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_76 * tmp_kernel_op_76) +
                                 (tmp_kernel_op_78 * tmp_kernel_op_78));
          const walberla::float64 elMat_3_4 = tmp_kernel_op_203;
          const walberla::float64 elMat_3_5 = tmp_kernel_op_204;
          const walberla::float64 elMat_3_6 = tmp_kernel_op_205;
          const walberla::float64 elMat_4_0 = tmp_kernel_op_129;
          const walberla::float64 elMat_4_1 = tmp_kernel_op_188;
          const walberla::float64 elMat_4_2 = tmp_kernel_op_200;
          const walberla::float64 elMat_4_3 = tmp_kernel_op_203;
          const walberla::float64 elMat_4_4 =
              tmp_kernel_op_11 * ((tmp_kernel_op_110 * tmp_kernel_op_110) +
                                  (tmp_kernel_op_112 * tmp_kernel_op_112)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_114 * tmp_kernel_op_114) +
                                  (tmp_kernel_op_116 * tmp_kernel_op_116)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_118 * tmp_kernel_op_118) +
                                  (tmp_kernel_op_120 * tmp_kernel_op_120)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_122 * tmp_kernel_op_122) +
                                  (tmp_kernel_op_124 * tmp_kernel_op_124)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_126 * tmp_kernel_op_126) +
                                  (tmp_kernel_op_128 * tmp_kernel_op_128)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_106 * tmp_kernel_op_106) +
                                 (tmp_kernel_op_108 * tmp_kernel_op_108));
          const walberla::float64 elMat_4_5 = tmp_kernel_op_206;
          const walberla::float64 elMat_4_6 = tmp_kernel_op_207;
          const walberla::float64 elMat_5_0 = tmp_kernel_op_148;
          const walberla::float64 elMat_5_1 = tmp_kernel_op_189;
          const walberla::float64 elMat_5_2 = tmp_kernel_op_201;
          const walberla::float64 elMat_5_3 = tmp_kernel_op_204;
          const walberla::float64 elMat_5_4 = tmp_kernel_op_206;
          const walberla::float64 elMat_5_5 =
              tmp_kernel_op_11 * ((tmp_kernel_op_134 * tmp_kernel_op_134) +
                                  (tmp_kernel_op_135 * tmp_kernel_op_135)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_137 * tmp_kernel_op_137) +
                                  (tmp_kernel_op_138 * tmp_kernel_op_138)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_140 * tmp_kernel_op_140) +
                                  (tmp_kernel_op_141 * tmp_kernel_op_141)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_143 * tmp_kernel_op_143) +
                                  (tmp_kernel_op_144 * tmp_kernel_op_144)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_146 * tmp_kernel_op_146) +
                                  (tmp_kernel_op_147 * tmp_kernel_op_147)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_131 * tmp_kernel_op_131) +
                                 (tmp_kernel_op_132 * tmp_kernel_op_132));
          const walberla::float64 elMat_5_6 = tmp_kernel_op_208;
          const walberla::float64 elMat_6_0 = tmp_kernel_op_177;
          const walberla::float64 elMat_6_1 = tmp_kernel_op_190;
          const walberla::float64 elMat_6_2 = tmp_kernel_op_202;
          const walberla::float64 elMat_6_3 = tmp_kernel_op_205;
          const walberla::float64 elMat_6_4 = tmp_kernel_op_207;
          const walberla::float64 elMat_6_5 = tmp_kernel_op_208;
          const walberla::float64 elMat_6_6 =
              tmp_kernel_op_11 * ((tmp_kernel_op_159 * tmp_kernel_op_159) +
                                  (tmp_kernel_op_160 * tmp_kernel_op_160)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_163 * tmp_kernel_op_163) +
                                  (tmp_kernel_op_164 * tmp_kernel_op_164)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_167 * tmp_kernel_op_167) +
                                  (tmp_kernel_op_168 * tmp_kernel_op_168)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_171 * tmp_kernel_op_171) +
                                  (tmp_kernel_op_172 * tmp_kernel_op_172)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_175 * tmp_kernel_op_175) +
                                  (tmp_kernel_op_176 * tmp_kernel_op_176)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_153 * tmp_kernel_op_153) +
                                 (tmp_kernel_op_156 * tmp_kernel_op_156));

          std::vector<uint_t> _data_rowIdx(7);
          std::vector<uint_t> _data_colIdx(7);
          std::vector<real_t> _data_mat(49);

          _data_rowIdx[0] =
              ((uint64_t)(_data_dstVertex[ctr_0 +
                                          ctr_1 *
                                              (micro_edges_per_macro_edge + 2) -
                                          ((ctr_1 * (ctr_1 + 1)) / (2))]));
          _data_rowIdx[1] =
              ((uint64_t)(_data_dstVertex[ctr_0 +
                                          ctr_1 *
                                              (micro_edges_per_macro_edge + 2) -
                                          ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_rowIdx[2] = ((
              uint64_t)(_data_dstVertex[ctr_0 +
                                        (ctr_1 + 1) *
                                            (micro_edges_per_macro_edge + 2) -
                                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_rowIdx[3] = ((
              uint64_t)(_data_dstEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      ((micro_edges_per_macro_edge *
                                        (micro_edges_per_macro_edge + 1)) /
                                       (2))]));
          _data_rowIdx[4] = ((
              uint64_t)(_data_dstEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      2 * ((micro_edges_per_macro_edge *
                                            (micro_edges_per_macro_edge + 1)) /
                                           (2))]));
          _data_rowIdx[5] = ((
              uint64_t)(_data_dstEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2))]));
          _data_rowIdx[6] =
              ((uint64_t)(_data_dst[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 1) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2))]));
          _data_colIdx[0] =
              ((uint64_t)(_data_srcVertex[ctr_0 +
                                          ctr_1 *
                                              (micro_edges_per_macro_edge + 2) -
                                          ((ctr_1 * (ctr_1 + 1)) / (2))]));
          _data_colIdx[1] =
              ((uint64_t)(_data_srcVertex[ctr_0 +
                                          ctr_1 *
                                              (micro_edges_per_macro_edge + 2) -
                                          ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_colIdx[2] = ((
              uint64_t)(_data_srcVertex[ctr_0 +
                                        (ctr_1 + 1) *
                                            (micro_edges_per_macro_edge + 2) -
                                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_colIdx[3] = ((
              uint64_t)(_data_srcEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      ((micro_edges_per_macro_edge *
                                        (micro_edges_per_macro_edge + 1)) /
                                       (2))]));
          _data_colIdx[4] = ((
              uint64_t)(_data_srcEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      2 * ((micro_edges_per_macro_edge *
                                            (micro_edges_per_macro_edge + 1)) /
                                           (2))]));
          _data_colIdx[5] = ((
              uint64_t)(_data_srcEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2))]));
          _data_colIdx[6] =
              ((uint64_t)(_data_src[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 1) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2))]));

          /* Apply basis transformation */

          _data_mat[0] = ((real_t)(elMat_0_0));
          _data_mat[1] = ((real_t)(elMat_0_1));
          _data_mat[2] = ((real_t)(elMat_0_2));
          _data_mat[3] = ((real_t)(elMat_0_3));
          _data_mat[4] = ((real_t)(elMat_0_4));
          _data_mat[5] = ((real_t)(elMat_0_5));
          _data_mat[6] = ((real_t)(elMat_0_6));
          _data_mat[7] = ((real_t)(elMat_1_0));
          _data_mat[8] = ((real_t)(elMat_1_1));
          _data_mat[9] = ((real_t)(elMat_1_2));
          _data_mat[10] = ((real_t)(elMat_1_3));
          _data_mat[11] = ((real_t)(elMat_1_4));
          _data_mat[12] = ((real_t)(elMat_1_5));
          _data_mat[13] = ((real_t)(elMat_1_6));
          _data_mat[14] = ((real_t)(elMat_2_0));
          _data_mat[15] = ((real_t)(elMat_2_1));
          _data_mat[16] = ((real_t)(elMat_2_2));
          _data_mat[17] = ((real_t)(elMat_2_3));
          _data_mat[18] = ((real_t)(elMat_2_4));
          _data_mat[19] = ((real_t)(elMat_2_5));
          _data_mat[20] = ((real_t)(elMat_2_6));
          _data_mat[21] = ((real_t)(elMat_3_0));
          _data_mat[22] = ((real_t)(elMat_3_1));
          _data_mat[23] = ((real_t)(elMat_3_2));
          _data_mat[24] = ((real_t)(elMat_3_3));
          _data_mat[25] = ((real_t)(elMat_3_4));
          _data_mat[26] = ((real_t)(elMat_3_5));
          _data_mat[27] = ((real_t)(elMat_3_6));
          _data_mat[28] = ((real_t)(elMat_4_0));
          _data_mat[29] = ((real_t)(elMat_4_1));
          _data_mat[30] = ((real_t)(elMat_4_2));
          _data_mat[31] = ((real_t)(elMat_4_3));
          _data_mat[32] = ((real_t)(elMat_4_4));
          _data_mat[33] = ((real_t)(elMat_4_5));
          _data_mat[34] = ((real_t)(elMat_4_6));
          _data_mat[35] = ((real_t)(elMat_5_0));
          _data_mat[36] = ((real_t)(elMat_5_1));
          _data_mat[37] = ((real_t)(elMat_5_2));
          _data_mat[38] = ((real_t)(elMat_5_3));
          _data_mat[39] = ((real_t)(elMat_5_4));
          _data_mat[40] = ((real_t)(elMat_5_5));
          _data_mat[41] = ((real_t)(elMat_5_6));
          _data_mat[42] = ((real_t)(elMat_6_0));
          _data_mat[43] = ((real_t)(elMat_6_1));
          _data_mat[44] = ((real_t)(elMat_6_2));
          _data_mat[45] = ((real_t)(elMat_6_3));
          _data_mat[46] = ((real_t)(elMat_6_4));
          _data_mat[47] = ((real_t)(elMat_6_5));
          _data_mat[48] = ((real_t)(elMat_6_6));

          mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
        }
    }
    {
      /* FaceType.BLUE */
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
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1;
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
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_0_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_1_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 tmp_kernel_op_0 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_1 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_2 =
              tmp_kernel_op_0 + tmp_kernel_op_1 - 3.0;
          const walberla::float64 tmp_kernel_op_3 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_2 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_5 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_6 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_7 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_8 =
              tmp_kernel_op_6 + tmp_kernel_op_7 - 3.0;
          const walberla::float64 tmp_kernel_op_9 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_8 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_10 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_8 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_11 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_12 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_13 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_12 + tmp_kernel_op_13 - 3.0;
          const walberla::float64 tmp_kernel_op_15 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_14 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_16 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_14 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_14;
          const walberla::float64 tmp_kernel_op_17 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_18 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_19 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_20 =
              tmp_kernel_op_18 + tmp_kernel_op_19 - 3.0;
          const walberla::float64 tmp_kernel_op_21 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_20 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_20 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_23 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_24 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_25 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_26 =
              tmp_kernel_op_24 + tmp_kernel_op_25 - 3.0;
          const walberla::float64 tmp_kernel_op_27 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_26 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_28 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_26 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_26;
          const walberla::float64 tmp_kernel_op_29 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_30 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_31 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_30 + tmp_kernel_op_31 - 3.0;
          const walberla::float64 tmp_kernel_op_33 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_32 +
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_34 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_32 +
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_35 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_36 = tmp_kernel_op_0 - 1.0;
          const walberla::float64 tmp_kernel_op_37 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_38 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_39 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_40 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_41 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_42 = tmp_kernel_op_12 - 1.0;
          const walberla::float64 tmp_kernel_op_43 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_44 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_45 = tmp_kernel_op_18 - 1.0;
          const walberla::float64 tmp_kernel_op_46 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_47 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_48 = tmp_kernel_op_24 - 1.0;
          const walberla::float64 tmp_kernel_op_49 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_50 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_51 = tmp_kernel_op_30 - 1.0;
          const walberla::float64 tmp_kernel_op_52 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_53 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_41 +
                                  tmp_kernel_op_40 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_43 +
                                  tmp_kernel_op_16 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_46 +
                                  tmp_kernel_op_22 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_49 +
                                  tmp_kernel_op_28 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_52 +
                                  tmp_kernel_op_34 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_37 +
                                 tmp_kernel_op_38 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_55 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_56 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_57 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_55;
          const walberla::float64 tmp_kernel_op_58 = tmp_kernel_op_7 - 1.0;
          const walberla::float64 tmp_kernel_op_59 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_60 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_61 = tmp_kernel_op_13 - 1.0;
          const walberla::float64 tmp_kernel_op_62 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_63 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_64 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_65 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_66 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_67 = tmp_kernel_op_25 - 1.0;
          const walberla::float64 tmp_kernel_op_68 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_69 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_70 = tmp_kernel_op_31 - 1.0;
          const walberla::float64 tmp_kernel_op_71 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_72 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_70;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_60 +
                                  tmp_kernel_op_59 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_62 +
                                  tmp_kernel_op_16 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_65 +
                                  tmp_kernel_op_22 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_68 +
                                  tmp_kernel_op_28 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_33 * tmp_kernel_op_71 +
                                  tmp_kernel_op_34 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_56 +
                                 tmp_kernel_op_4 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_74 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_75 =
              -tmp_kernel_op_1 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_76 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_75 - tmp_kernel_op_74;
          const walberla::float64 tmp_kernel_op_77 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_78 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_75 - tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_79 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_80 =
              -tmp_kernel_op_7 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_81 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_80 - tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_82 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_83 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_80 - tmp_kernel_op_82;
          const walberla::float64 tmp_kernel_op_84 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_85 =
              -tmp_kernel_op_13 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_86 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_85 - tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_87 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_88 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_85 - tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_90 =
              -tmp_kernel_op_19 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_91 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_90 - tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_92 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_93 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_90 - tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_94 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_95 =
              -tmp_kernel_op_25 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_96 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_95 - tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_97 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_98 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_95 - tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_99 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_100 =
              -tmp_kernel_op_31 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_101 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_100 - tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_102 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_103 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_100 - tmp_kernel_op_102;
          const walberla::float64 tmp_kernel_op_104 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_83 +
                                  tmp_kernel_op_81 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_86 +
                                  tmp_kernel_op_16 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_21 * tmp_kernel_op_91 +
                                  tmp_kernel_op_22 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_27 * tmp_kernel_op_96 +
                                  tmp_kernel_op_28 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_33 +
                                  tmp_kernel_op_103 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_3 * tmp_kernel_op_76 +
                                 tmp_kernel_op_4 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_105 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_105 + tmp_kernel_op_74;
          const walberla::float64 tmp_kernel_op_107 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_108 =
              tmp_kernel_op_107 + tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_109 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_109 + tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_111 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_112 =
              tmp_kernel_op_111 + tmp_kernel_op_82;
          const walberla::float64 tmp_kernel_op_113 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_114 =
              tmp_kernel_op_113 + tmp_kernel_op_84;
          const walberla::float64 tmp_kernel_op_115 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_116 =
              tmp_kernel_op_115 + tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_117 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_118 =
              tmp_kernel_op_117 + tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_119 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_19;
          const walberla::float64 tmp_kernel_op_120 =
              tmp_kernel_op_119 + tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_121 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_122 =
              tmp_kernel_op_121 + tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_123 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_124 =
              tmp_kernel_op_123 + tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_125 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_126 =
              tmp_kernel_op_125 + tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_127 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_128 =
              tmp_kernel_op_102 + tmp_kernel_op_127;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_112 +
                                  tmp_kernel_op_110 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_15 +
                                  tmp_kernel_op_116 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_21 +
                                  tmp_kernel_op_120 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_27 +
                                  tmp_kernel_op_124 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_33 +
                                  tmp_kernel_op_128 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_3 +
                                 tmp_kernel_op_108 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_0 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_131 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_130 - tmp_kernel_op_105;
          const walberla::float64 tmp_kernel_op_132 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_130 - tmp_kernel_op_107;
          const walberla::float64 tmp_kernel_op_133 =
              -tmp_kernel_op_6 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_134 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_133 - tmp_kernel_op_109;
          const walberla::float64 tmp_kernel_op_135 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_133 - tmp_kernel_op_111;
          const walberla::float64 tmp_kernel_op_136 =
              -tmp_kernel_op_12 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_137 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_136 - tmp_kernel_op_113;
          const walberla::float64 tmp_kernel_op_138 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_136 - tmp_kernel_op_115;
          const walberla::float64 tmp_kernel_op_139 =
              -tmp_kernel_op_18 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_140 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_139 - tmp_kernel_op_117;
          const walberla::float64 tmp_kernel_op_141 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_139 - tmp_kernel_op_119;
          const walberla::float64 tmp_kernel_op_142 =
              -tmp_kernel_op_24 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_143 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_142 - tmp_kernel_op_121;
          const walberla::float64 tmp_kernel_op_144 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_142 - tmp_kernel_op_123;
          const walberla::float64 tmp_kernel_op_145 =
              -tmp_kernel_op_30 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_146 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_145 - tmp_kernel_op_125;
          const walberla::float64 tmp_kernel_op_147 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_145 - tmp_kernel_op_127;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_135 +
                                  tmp_kernel_op_134 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_15 +
                                  tmp_kernel_op_138 * tmp_kernel_op_16) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_21 +
                                  tmp_kernel_op_141 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_27 +
                                  tmp_kernel_op_144 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_33 +
                                  tmp_kernel_op_147 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_3 +
                                 tmp_kernel_op_132 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_149 =
              jac_affine_inv_0_0_BLUE * 27.0;
          const int64_t tmp_kernel_op_150 = 0;
          const walberla::float64 tmp_kernel_op_151 =
              jac_affine_inv_1_0_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_152 = 0.15066167873471437;
          const walberla::float64 tmp_kernel_op_153 =
              tmp_kernel_op_149 * ((walberla::float64)(tmp_kernel_op_150)) +
              tmp_kernel_op_151 * tmp_kernel_op_152;
          const walberla::float64 tmp_kernel_op_154 =
              jac_affine_inv_0_1_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_155 =
              jac_affine_inv_1_1_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_152 * tmp_kernel_op_155 +
              tmp_kernel_op_154 * ((walberla::float64)(tmp_kernel_op_150));
          const walberla::float64 tmp_kernel_op_157 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_158 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_149 * tmp_kernel_op_157 +
              tmp_kernel_op_151 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_154 * tmp_kernel_op_157 +
              tmp_kernel_op_155 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_161 = 0.15066167873471437;
          const int64_t tmp_kernel_op_162 = 0;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_149 * tmp_kernel_op_161 +
              tmp_kernel_op_151 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_154 * tmp_kernel_op_161 +
              tmp_kernel_op_155 * ((walberla::float64)(tmp_kernel_op_162));
          const walberla::float64 tmp_kernel_op_165 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_166 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_149 * tmp_kernel_op_165 +
              tmp_kernel_op_151 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_154 * tmp_kernel_op_165 +
              tmp_kernel_op_155 * tmp_kernel_op_166;
          const walberla::float64 tmp_kernel_op_169 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_170 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_149 * tmp_kernel_op_169 +
              tmp_kernel_op_151 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_154 * tmp_kernel_op_169 +
              tmp_kernel_op_155 * tmp_kernel_op_170;
          const walberla::float64 tmp_kernel_op_173 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_174 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_149 * tmp_kernel_op_173 +
              tmp_kernel_op_151 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_154 * tmp_kernel_op_173 +
              tmp_kernel_op_155 * tmp_kernel_op_174;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_11 * (tmp_kernel_op_10 * tmp_kernel_op_160 +
                                  tmp_kernel_op_159 * tmp_kernel_op_9) +
              tmp_kernel_op_17 * (tmp_kernel_op_15 * tmp_kernel_op_163 +
                                  tmp_kernel_op_16 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_21 +
                                  tmp_kernel_op_168 * tmp_kernel_op_22) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_27 +
                                  tmp_kernel_op_172 * tmp_kernel_op_28) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_33 +
                                  tmp_kernel_op_176 * tmp_kernel_op_34) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_3 +
                                 tmp_kernel_op_156 * tmp_kernel_op_4);
          const walberla::float64 tmp_kernel_op_178 =
              (jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE);
          const walberla::float64 tmp_kernel_op_179 =
              (tmp_kernel_op_36 * tmp_kernel_op_36);
          const walberla::float64 tmp_kernel_op_180 =
              (jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE);
          const walberla::float64 tmp_kernel_op_181 =
              (tmp_kernel_op_39 * tmp_kernel_op_39);
          const walberla::float64 tmp_kernel_op_182 =
              (tmp_kernel_op_42 * tmp_kernel_op_42);
          const walberla::float64 tmp_kernel_op_183 =
              (tmp_kernel_op_45 * tmp_kernel_op_45);
          const walberla::float64 tmp_kernel_op_184 =
              (tmp_kernel_op_48 * tmp_kernel_op_48);
          const walberla::float64 tmp_kernel_op_185 =
              (tmp_kernel_op_51 * tmp_kernel_op_51);
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_59 +
                                  tmp_kernel_op_41 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_62 +
                                  tmp_kernel_op_44 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_65 +
                                  tmp_kernel_op_47 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_49 * tmp_kernel_op_68 +
                                  tmp_kernel_op_50 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_52 * tmp_kernel_op_71 +
                                  tmp_kernel_op_53 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_56 +
                                 tmp_kernel_op_38 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_11 * (tmp_kernel_op_40 * tmp_kernel_op_81 +
                                  tmp_kernel_op_41 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_43 * tmp_kernel_op_86 +
                                  tmp_kernel_op_44 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_46 * tmp_kernel_op_91 +
                                  tmp_kernel_op_47 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_49 * tmp_kernel_op_96 +
                                  tmp_kernel_op_50 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_52 +
                                  tmp_kernel_op_103 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_37 * tmp_kernel_op_76 +
                                 tmp_kernel_op_38 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_40 +
                                  tmp_kernel_op_112 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_43 +
                                  tmp_kernel_op_116 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_46 +
                                  tmp_kernel_op_120 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_49 +
                                  tmp_kernel_op_124 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_52 +
                                  tmp_kernel_op_128 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_37 +
                                 tmp_kernel_op_108 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_40 +
                                  tmp_kernel_op_135 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_43 +
                                  tmp_kernel_op_138 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_46 +
                                  tmp_kernel_op_141 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_49 +
                                  tmp_kernel_op_144 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_52 +
                                  tmp_kernel_op_147 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_37 +
                                 tmp_kernel_op_132 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_40 +
                                  tmp_kernel_op_160 * tmp_kernel_op_41) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_43 +
                                  tmp_kernel_op_164 * tmp_kernel_op_44) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_46 +
                                  tmp_kernel_op_168 * tmp_kernel_op_47) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_49 +
                                  tmp_kernel_op_172 * tmp_kernel_op_50) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_52 +
                                  tmp_kernel_op_176 * tmp_kernel_op_53) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_37 +
                                 tmp_kernel_op_156 * tmp_kernel_op_38);
          const walberla::float64 tmp_kernel_op_191 =
              (jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE);
          const walberla::float64 tmp_kernel_op_192 =
              (tmp_kernel_op_55 * tmp_kernel_op_55);
          const walberla::float64 tmp_kernel_op_193 =
              (jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE);
          const walberla::float64 tmp_kernel_op_194 =
              (tmp_kernel_op_58 * tmp_kernel_op_58);
          const walberla::float64 tmp_kernel_op_195 =
              (tmp_kernel_op_61 * tmp_kernel_op_61);
          const walberla::float64 tmp_kernel_op_196 =
              (tmp_kernel_op_64 * tmp_kernel_op_64);
          const walberla::float64 tmp_kernel_op_197 =
              (tmp_kernel_op_67 * tmp_kernel_op_67);
          const walberla::float64 tmp_kernel_op_198 =
              (tmp_kernel_op_70 * tmp_kernel_op_70);
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_11 * (tmp_kernel_op_59 * tmp_kernel_op_81 +
                                  tmp_kernel_op_60 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_62 * tmp_kernel_op_86 +
                                  tmp_kernel_op_63 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_65 * tmp_kernel_op_91 +
                                  tmp_kernel_op_66 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_68 * tmp_kernel_op_96 +
                                  tmp_kernel_op_69 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_71 +
                                  tmp_kernel_op_103 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_56 * tmp_kernel_op_76 +
                                 tmp_kernel_op_57 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_59 +
                                  tmp_kernel_op_112 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_62 +
                                  tmp_kernel_op_116 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_65 +
                                  tmp_kernel_op_120 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_68 +
                                  tmp_kernel_op_124 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_71 +
                                  tmp_kernel_op_128 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_56 +
                                 tmp_kernel_op_108 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_201 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_59 +
                                  tmp_kernel_op_135 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_62 +
                                  tmp_kernel_op_138 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_65 +
                                  tmp_kernel_op_141 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_68 +
                                  tmp_kernel_op_144 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_71 +
                                  tmp_kernel_op_147 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_56 +
                                 tmp_kernel_op_132 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_59 +
                                  tmp_kernel_op_160 * tmp_kernel_op_60) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_62 +
                                  tmp_kernel_op_164 * tmp_kernel_op_63) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_65 +
                                  tmp_kernel_op_168 * tmp_kernel_op_66) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_68 +
                                  tmp_kernel_op_172 * tmp_kernel_op_69) +
              tmp_kernel_op_35 * (tmp_kernel_op_175 * tmp_kernel_op_71 +
                                  tmp_kernel_op_176 * tmp_kernel_op_72) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_56 +
                                 tmp_kernel_op_156 * tmp_kernel_op_57);
          const walberla::float64 tmp_kernel_op_203 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_81 +
                                  tmp_kernel_op_112 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_86 +
                                  tmp_kernel_op_116 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_91 +
                                  tmp_kernel_op_120 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_96 +
                                  tmp_kernel_op_124 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_126 +
                                  tmp_kernel_op_103 * tmp_kernel_op_128) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_76 +
                                 tmp_kernel_op_108 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_81 +
                                  tmp_kernel_op_135 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_86 +
                                  tmp_kernel_op_138 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_91 +
                                  tmp_kernel_op_141 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_96 +
                                  tmp_kernel_op_144 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_146 +
                                  tmp_kernel_op_103 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_76 +
                                 tmp_kernel_op_132 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_205 =
              tmp_kernel_op_11 * (tmp_kernel_op_159 * tmp_kernel_op_81 +
                                  tmp_kernel_op_160 * tmp_kernel_op_83) +
              tmp_kernel_op_17 * (tmp_kernel_op_163 * tmp_kernel_op_86 +
                                  tmp_kernel_op_164 * tmp_kernel_op_88) +
              tmp_kernel_op_23 * (tmp_kernel_op_167 * tmp_kernel_op_91 +
                                  tmp_kernel_op_168 * tmp_kernel_op_93) +
              tmp_kernel_op_29 * (tmp_kernel_op_171 * tmp_kernel_op_96 +
                                  tmp_kernel_op_172 * tmp_kernel_op_98) +
              tmp_kernel_op_35 * (tmp_kernel_op_101 * tmp_kernel_op_175 +
                                  tmp_kernel_op_103 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_153 * tmp_kernel_op_76 +
                                 tmp_kernel_op_156 * tmp_kernel_op_78);
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_134 +
                                  tmp_kernel_op_112 * tmp_kernel_op_135) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_137 +
                                  tmp_kernel_op_116 * tmp_kernel_op_138) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_140 +
                                  tmp_kernel_op_120 * tmp_kernel_op_141) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_143 +
                                  tmp_kernel_op_124 * tmp_kernel_op_144) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_146 +
                                  tmp_kernel_op_128 * tmp_kernel_op_147) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_131 +
                                 tmp_kernel_op_108 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_207 =
              tmp_kernel_op_11 * (tmp_kernel_op_110 * tmp_kernel_op_159 +
                                  tmp_kernel_op_112 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_114 * tmp_kernel_op_163 +
                                  tmp_kernel_op_116 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_118 * tmp_kernel_op_167 +
                                  tmp_kernel_op_120 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_122 * tmp_kernel_op_171 +
                                  tmp_kernel_op_124 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_126 * tmp_kernel_op_175 +
                                  tmp_kernel_op_128 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_106 * tmp_kernel_op_153 +
                                 tmp_kernel_op_108 * tmp_kernel_op_156);
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_11 * (tmp_kernel_op_134 * tmp_kernel_op_159 +
                                  tmp_kernel_op_135 * tmp_kernel_op_160) +
              tmp_kernel_op_17 * (tmp_kernel_op_137 * tmp_kernel_op_163 +
                                  tmp_kernel_op_138 * tmp_kernel_op_164) +
              tmp_kernel_op_23 * (tmp_kernel_op_140 * tmp_kernel_op_167 +
                                  tmp_kernel_op_141 * tmp_kernel_op_168) +
              tmp_kernel_op_29 * (tmp_kernel_op_143 * tmp_kernel_op_171 +
                                  tmp_kernel_op_144 * tmp_kernel_op_172) +
              tmp_kernel_op_35 * (tmp_kernel_op_146 * tmp_kernel_op_175 +
                                  tmp_kernel_op_147 * tmp_kernel_op_176) +
              tmp_kernel_op_5 * (tmp_kernel_op_131 * tmp_kernel_op_153 +
                                 tmp_kernel_op_132 * tmp_kernel_op_156);
          const walberla::float64 elMat_0_0 =
              tmp_kernel_op_11 * ((tmp_kernel_op_10 * tmp_kernel_op_10) +
                                  (tmp_kernel_op_9 * tmp_kernel_op_9)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_15 * tmp_kernel_op_15) +
                                  (tmp_kernel_op_16 * tmp_kernel_op_16)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_21 * tmp_kernel_op_21) +
                                  (tmp_kernel_op_22 * tmp_kernel_op_22)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_27 * tmp_kernel_op_27) +
                                  (tmp_kernel_op_28 * tmp_kernel_op_28)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_33 * tmp_kernel_op_33) +
                                  (tmp_kernel_op_34 * tmp_kernel_op_34)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_3 * tmp_kernel_op_3) +
                                 (tmp_kernel_op_4 * tmp_kernel_op_4));
          const walberla::float64 elMat_0_1 = tmp_kernel_op_54;
          const walberla::float64 elMat_0_2 = tmp_kernel_op_73;
          const walberla::float64 elMat_0_3 = tmp_kernel_op_104;
          const walberla::float64 elMat_0_4 = tmp_kernel_op_129;
          const walberla::float64 elMat_0_5 = tmp_kernel_op_148;
          const walberla::float64 elMat_0_6 = tmp_kernel_op_177;
          const walberla::float64 elMat_1_0 = tmp_kernel_op_54;
          const walberla::float64 elMat_1_1 =
              tmp_kernel_op_11 * (tmp_kernel_op_178 * tmp_kernel_op_181 +
                                  tmp_kernel_op_180 * tmp_kernel_op_181) +
              tmp_kernel_op_17 * (tmp_kernel_op_178 * tmp_kernel_op_182 +
                                  tmp_kernel_op_180 * tmp_kernel_op_182) +
              tmp_kernel_op_23 * (tmp_kernel_op_178 * tmp_kernel_op_183 +
                                  tmp_kernel_op_180 * tmp_kernel_op_183) +
              tmp_kernel_op_29 * (tmp_kernel_op_178 * tmp_kernel_op_184 +
                                  tmp_kernel_op_180 * tmp_kernel_op_184) +
              tmp_kernel_op_35 * (tmp_kernel_op_178 * tmp_kernel_op_185 +
                                  tmp_kernel_op_180 * tmp_kernel_op_185) +
              tmp_kernel_op_5 * (tmp_kernel_op_178 * tmp_kernel_op_179 +
                                 tmp_kernel_op_179 * tmp_kernel_op_180);
          const walberla::float64 elMat_1_2 = tmp_kernel_op_186;
          const walberla::float64 elMat_1_3 = tmp_kernel_op_187;
          const walberla::float64 elMat_1_4 = tmp_kernel_op_188;
          const walberla::float64 elMat_1_5 = tmp_kernel_op_189;
          const walberla::float64 elMat_1_6 = tmp_kernel_op_190;
          const walberla::float64 elMat_2_0 = tmp_kernel_op_73;
          const walberla::float64 elMat_2_1 = tmp_kernel_op_186;
          const walberla::float64 elMat_2_2 =
              tmp_kernel_op_11 * (tmp_kernel_op_191 * tmp_kernel_op_194 +
                                  tmp_kernel_op_193 * tmp_kernel_op_194) +
              tmp_kernel_op_17 * (tmp_kernel_op_191 * tmp_kernel_op_195 +
                                  tmp_kernel_op_193 * tmp_kernel_op_195) +
              tmp_kernel_op_23 * (tmp_kernel_op_191 * tmp_kernel_op_196 +
                                  tmp_kernel_op_193 * tmp_kernel_op_196) +
              tmp_kernel_op_29 * (tmp_kernel_op_191 * tmp_kernel_op_197 +
                                  tmp_kernel_op_193 * tmp_kernel_op_197) +
              tmp_kernel_op_35 * (tmp_kernel_op_191 * tmp_kernel_op_198 +
                                  tmp_kernel_op_193 * tmp_kernel_op_198) +
              tmp_kernel_op_5 * (tmp_kernel_op_191 * tmp_kernel_op_192 +
                                 tmp_kernel_op_192 * tmp_kernel_op_193);
          const walberla::float64 elMat_2_3 = tmp_kernel_op_199;
          const walberla::float64 elMat_2_4 = tmp_kernel_op_200;
          const walberla::float64 elMat_2_5 = tmp_kernel_op_201;
          const walberla::float64 elMat_2_6 = tmp_kernel_op_202;
          const walberla::float64 elMat_3_0 = tmp_kernel_op_104;
          const walberla::float64 elMat_3_1 = tmp_kernel_op_187;
          const walberla::float64 elMat_3_2 = tmp_kernel_op_199;
          const walberla::float64 elMat_3_3 =
              tmp_kernel_op_11 * ((tmp_kernel_op_81 * tmp_kernel_op_81) +
                                  (tmp_kernel_op_83 * tmp_kernel_op_83)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_86 * tmp_kernel_op_86) +
                                  (tmp_kernel_op_88 * tmp_kernel_op_88)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_91 * tmp_kernel_op_91) +
                                  (tmp_kernel_op_93 * tmp_kernel_op_93)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_96 * tmp_kernel_op_96) +
                                  (tmp_kernel_op_98 * tmp_kernel_op_98)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_101 * tmp_kernel_op_101) +
                                  (tmp_kernel_op_103 * tmp_kernel_op_103)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_76 * tmp_kernel_op_76) +
                                 (tmp_kernel_op_78 * tmp_kernel_op_78));
          const walberla::float64 elMat_3_4 = tmp_kernel_op_203;
          const walberla::float64 elMat_3_5 = tmp_kernel_op_204;
          const walberla::float64 elMat_3_6 = tmp_kernel_op_205;
          const walberla::float64 elMat_4_0 = tmp_kernel_op_129;
          const walberla::float64 elMat_4_1 = tmp_kernel_op_188;
          const walberla::float64 elMat_4_2 = tmp_kernel_op_200;
          const walberla::float64 elMat_4_3 = tmp_kernel_op_203;
          const walberla::float64 elMat_4_4 =
              tmp_kernel_op_11 * ((tmp_kernel_op_110 * tmp_kernel_op_110) +
                                  (tmp_kernel_op_112 * tmp_kernel_op_112)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_114 * tmp_kernel_op_114) +
                                  (tmp_kernel_op_116 * tmp_kernel_op_116)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_118 * tmp_kernel_op_118) +
                                  (tmp_kernel_op_120 * tmp_kernel_op_120)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_122 * tmp_kernel_op_122) +
                                  (tmp_kernel_op_124 * tmp_kernel_op_124)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_126 * tmp_kernel_op_126) +
                                  (tmp_kernel_op_128 * tmp_kernel_op_128)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_106 * tmp_kernel_op_106) +
                                 (tmp_kernel_op_108 * tmp_kernel_op_108));
          const walberla::float64 elMat_4_5 = tmp_kernel_op_206;
          const walberla::float64 elMat_4_6 = tmp_kernel_op_207;
          const walberla::float64 elMat_5_0 = tmp_kernel_op_148;
          const walberla::float64 elMat_5_1 = tmp_kernel_op_189;
          const walberla::float64 elMat_5_2 = tmp_kernel_op_201;
          const walberla::float64 elMat_5_3 = tmp_kernel_op_204;
          const walberla::float64 elMat_5_4 = tmp_kernel_op_206;
          const walberla::float64 elMat_5_5 =
              tmp_kernel_op_11 * ((tmp_kernel_op_134 * tmp_kernel_op_134) +
                                  (tmp_kernel_op_135 * tmp_kernel_op_135)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_137 * tmp_kernel_op_137) +
                                  (tmp_kernel_op_138 * tmp_kernel_op_138)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_140 * tmp_kernel_op_140) +
                                  (tmp_kernel_op_141 * tmp_kernel_op_141)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_143 * tmp_kernel_op_143) +
                                  (tmp_kernel_op_144 * tmp_kernel_op_144)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_146 * tmp_kernel_op_146) +
                                  (tmp_kernel_op_147 * tmp_kernel_op_147)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_131 * tmp_kernel_op_131) +
                                 (tmp_kernel_op_132 * tmp_kernel_op_132));
          const walberla::float64 elMat_5_6 = tmp_kernel_op_208;
          const walberla::float64 elMat_6_0 = tmp_kernel_op_177;
          const walberla::float64 elMat_6_1 = tmp_kernel_op_190;
          const walberla::float64 elMat_6_2 = tmp_kernel_op_202;
          const walberla::float64 elMat_6_3 = tmp_kernel_op_205;
          const walberla::float64 elMat_6_4 = tmp_kernel_op_207;
          const walberla::float64 elMat_6_5 = tmp_kernel_op_208;
          const walberla::float64 elMat_6_6 =
              tmp_kernel_op_11 * ((tmp_kernel_op_159 * tmp_kernel_op_159) +
                                  (tmp_kernel_op_160 * tmp_kernel_op_160)) +
              tmp_kernel_op_17 * ((tmp_kernel_op_163 * tmp_kernel_op_163) +
                                  (tmp_kernel_op_164 * tmp_kernel_op_164)) +
              tmp_kernel_op_23 * ((tmp_kernel_op_167 * tmp_kernel_op_167) +
                                  (tmp_kernel_op_168 * tmp_kernel_op_168)) +
              tmp_kernel_op_29 * ((tmp_kernel_op_171 * tmp_kernel_op_171) +
                                  (tmp_kernel_op_172 * tmp_kernel_op_172)) +
              tmp_kernel_op_35 * ((tmp_kernel_op_175 * tmp_kernel_op_175) +
                                  (tmp_kernel_op_176 * tmp_kernel_op_176)) +
              tmp_kernel_op_5 * ((tmp_kernel_op_153 * tmp_kernel_op_153) +
                                 (tmp_kernel_op_156 * tmp_kernel_op_156));

          std::vector<uint_t> _data_rowIdx(7);
          std::vector<uint_t> _data_colIdx(7);
          std::vector<real_t> _data_mat(49);

          _data_rowIdx[0] =
              ((uint64_t)(_data_dstVertex[ctr_0 +
                                          ctr_1 *
                                              (micro_edges_per_macro_edge + 2) -
                                          ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_rowIdx[1] = ((
              uint64_t)(_data_dstVertex[ctr_0 +
                                        (ctr_1 + 1) *
                                            (micro_edges_per_macro_edge + 2) -
                                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_rowIdx[2] =
              ((uint64_t)(_data_dstVertex[ctr_0 +
                                          (ctr_1 + 1) *
                                              (micro_edges_per_macro_edge + 2) -
                                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                          1]));
          _data_rowIdx[3] =
              ((uint64_t)(_data_dstEdge[ctr_0 +
                                        (ctr_1 + 1) *
                                            (micro_edges_per_macro_edge + 1) -
                                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_rowIdx[4] = ((
              uint64_t)(_data_dstEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      2 * ((micro_edges_per_macro_edge *
                                            (micro_edges_per_macro_edge + 1)) /
                                           (2)) +
                                      1]));
          _data_rowIdx[5] = ((
              uint64_t)(_data_dstEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      ((micro_edges_per_macro_edge *
                                        (micro_edges_per_macro_edge + 1)) /
                                       (2))]));
          _data_rowIdx[6] =
              ((uint64_t)(_data_dst[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                                    ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                    ((micro_edges_per_macro_edge *
                                      (micro_edges_per_macro_edge + 1)) /
                                     (2))]));
          _data_colIdx[0] =
              ((uint64_t)(_data_srcVertex[ctr_0 +
                                          ctr_1 *
                                              (micro_edges_per_macro_edge + 2) -
                                          ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_colIdx[1] = ((
              uint64_t)(_data_srcVertex[ctr_0 +
                                        (ctr_1 + 1) *
                                            (micro_edges_per_macro_edge + 2) -
                                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_colIdx[2] =
              ((uint64_t)(_data_srcVertex[ctr_0 +
                                          (ctr_1 + 1) *
                                              (micro_edges_per_macro_edge + 2) -
                                          (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) +
                                          1]));
          _data_colIdx[3] =
              ((uint64_t)(_data_srcEdge[ctr_0 +
                                        (ctr_1 + 1) *
                                            (micro_edges_per_macro_edge + 1) -
                                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_colIdx[4] = ((
              uint64_t)(_data_srcEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      2 * ((micro_edges_per_macro_edge *
                                            (micro_edges_per_macro_edge + 1)) /
                                           (2)) +
                                      1]));
          _data_colIdx[5] = ((
              uint64_t)(_data_srcEdge[ctr_0 +
                                      ctr_1 * (micro_edges_per_macro_edge + 1) -
                                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                      ((micro_edges_per_macro_edge *
                                        (micro_edges_per_macro_edge + 1)) /
                                       (2))]));
          _data_colIdx[6] =
              ((uint64_t)(_data_src[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                                    ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                    ((micro_edges_per_macro_edge *
                                      (micro_edges_per_macro_edge + 1)) /
                                     (2))]));

          /* Apply basis transformation */

          _data_mat[0] = ((real_t)(elMat_0_0));
          _data_mat[1] = ((real_t)(elMat_0_1));
          _data_mat[2] = ((real_t)(elMat_0_2));
          _data_mat[3] = ((real_t)(elMat_0_3));
          _data_mat[4] = ((real_t)(elMat_0_4));
          _data_mat[5] = ((real_t)(elMat_0_5));
          _data_mat[6] = ((real_t)(elMat_0_6));
          _data_mat[7] = ((real_t)(elMat_1_0));
          _data_mat[8] = ((real_t)(elMat_1_1));
          _data_mat[9] = ((real_t)(elMat_1_2));
          _data_mat[10] = ((real_t)(elMat_1_3));
          _data_mat[11] = ((real_t)(elMat_1_4));
          _data_mat[12] = ((real_t)(elMat_1_5));
          _data_mat[13] = ((real_t)(elMat_1_6));
          _data_mat[14] = ((real_t)(elMat_2_0));
          _data_mat[15] = ((real_t)(elMat_2_1));
          _data_mat[16] = ((real_t)(elMat_2_2));
          _data_mat[17] = ((real_t)(elMat_2_3));
          _data_mat[18] = ((real_t)(elMat_2_4));
          _data_mat[19] = ((real_t)(elMat_2_5));
          _data_mat[20] = ((real_t)(elMat_2_6));
          _data_mat[21] = ((real_t)(elMat_3_0));
          _data_mat[22] = ((real_t)(elMat_3_1));
          _data_mat[23] = ((real_t)(elMat_3_2));
          _data_mat[24] = ((real_t)(elMat_3_3));
          _data_mat[25] = ((real_t)(elMat_3_4));
          _data_mat[26] = ((real_t)(elMat_3_5));
          _data_mat[27] = ((real_t)(elMat_3_6));
          _data_mat[28] = ((real_t)(elMat_4_0));
          _data_mat[29] = ((real_t)(elMat_4_1));
          _data_mat[30] = ((real_t)(elMat_4_2));
          _data_mat[31] = ((real_t)(elMat_4_3));
          _data_mat[32] = ((real_t)(elMat_4_4));
          _data_mat[33] = ((real_t)(elMat_4_5));
          _data_mat[34] = ((real_t)(elMat_4_6));
          _data_mat[35] = ((real_t)(elMat_5_0));
          _data_mat[36] = ((real_t)(elMat_5_1));
          _data_mat[37] = ((real_t)(elMat_5_2));
          _data_mat[38] = ((real_t)(elMat_5_3));
          _data_mat[39] = ((real_t)(elMat_5_4));
          _data_mat[40] = ((real_t)(elMat_5_5));
          _data_mat[41] = ((real_t)(elMat_5_6));
          _data_mat[42] = ((real_t)(elMat_6_0));
          _data_mat[43] = ((real_t)(elMat_6_1));
          _data_mat[44] = ((real_t)(elMat_6_2));
          _data_mat[45] = ((real_t)(elMat_6_3));
          _data_mat[46] = ((real_t)(elMat_6_4));
          _data_mat[47] = ((real_t)(elMat_6_5));
          _data_mat[48] = ((real_t)(elMat_6_6));

          mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
        }
    }
  }
}
void P2PlusBubbleElementwiseDiffusion_float64::
    computeInverseDiagonalOperatorValues_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(
        walberla::float64 *RESTRICT _data_invDiag_,
        walberla::float64 *RESTRICT _data_invDiag_Edge,
        walberla::float64 *RESTRICT _data_invDiag_Vertex,
        walberla::float64 macro_vertex_coord_id_0comp0,
        walberla::float64 macro_vertex_coord_id_0comp1,
        walberla::float64 macro_vertex_coord_id_1comp0,
        walberla::float64 macro_vertex_coord_id_1comp1,
        walberla::float64 macro_vertex_coord_id_2comp0,
        walberla::float64 macro_vertex_coord_id_2comp1,
        int64_t micro_edges_per_macro_edge,
        walberla::float64 micro_edges_per_macro_edge_float) const {
  {
    {
      /* FaceType.GRAY */
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
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge;
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
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_0_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_2_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 tmp_kernel_op_0 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_1 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_2 =
              tmp_kernel_op_0 + tmp_kernel_op_1 - 3.0;
          const walberla::float64 tmp_kernel_op_3 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_4 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_5 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_6 =
              tmp_kernel_op_4 + tmp_kernel_op_5 - 3.0;
          const walberla::float64 tmp_kernel_op_7 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_8 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_9 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_10 =
              tmp_kernel_op_8 + tmp_kernel_op_9 - 3.0;
          const walberla::float64 tmp_kernel_op_11 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_12 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_13 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_12 + tmp_kernel_op_13 - 3.0;
          const walberla::float64 tmp_kernel_op_15 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_16 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_17 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_18 =
              tmp_kernel_op_16 + tmp_kernel_op_17 - 3.0;
          const walberla::float64 tmp_kernel_op_19 =
              abs_det_jac_affine_GRAY * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_20 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_21 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_20 + tmp_kernel_op_21 - 3.0;
          const walberla::float64 tmp_kernel_op_23 =
              abs_det_jac_affine_GRAY * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_24 =
              (jac_affine_inv_0_0_GRAY * jac_affine_inv_0_0_GRAY);
          const walberla::float64 tmp_kernel_op_25 =
              ((tmp_kernel_op_0 - 1.0) * (tmp_kernel_op_0 - 1.0));
          const walberla::float64 tmp_kernel_op_26 =
              (jac_affine_inv_0_1_GRAY * jac_affine_inv_0_1_GRAY);
          const walberla::float64 tmp_kernel_op_27 =
              ((tmp_kernel_op_4 - 1.0) * (tmp_kernel_op_4 - 1.0));
          const walberla::float64 tmp_kernel_op_28 =
              ((tmp_kernel_op_8 - 1.0) * (tmp_kernel_op_8 - 1.0));
          const walberla::float64 tmp_kernel_op_29 =
              ((tmp_kernel_op_12 - 1.0) * (tmp_kernel_op_12 - 1.0));
          const walberla::float64 tmp_kernel_op_30 =
              ((tmp_kernel_op_16 - 1.0) * (tmp_kernel_op_16 - 1.0));
          const walberla::float64 tmp_kernel_op_31 =
              ((tmp_kernel_op_20 - 1.0) * (tmp_kernel_op_20 - 1.0));
          const walberla::float64 tmp_kernel_op_32 =
              (jac_affine_inv_1_0_GRAY * jac_affine_inv_1_0_GRAY);
          const walberla::float64 tmp_kernel_op_33 =
              ((tmp_kernel_op_1 - 1.0) * (tmp_kernel_op_1 - 1.0));
          const walberla::float64 tmp_kernel_op_34 =
              (jac_affine_inv_1_1_GRAY * jac_affine_inv_1_1_GRAY);
          const walberla::float64 tmp_kernel_op_35 =
              ((tmp_kernel_op_5 - 1.0) * (tmp_kernel_op_5 - 1.0));
          const walberla::float64 tmp_kernel_op_36 =
              ((tmp_kernel_op_9 - 1.0) * (tmp_kernel_op_9 - 1.0));
          const walberla::float64 tmp_kernel_op_37 =
              ((tmp_kernel_op_13 - 1.0) * (tmp_kernel_op_13 - 1.0));
          const walberla::float64 tmp_kernel_op_38 =
              ((tmp_kernel_op_17 - 1.0) * (tmp_kernel_op_17 - 1.0));
          const walberla::float64 tmp_kernel_op_39 =
              ((tmp_kernel_op_21 - 1.0) * (tmp_kernel_op_21 - 1.0));
          const walberla::float64 tmp_kernel_op_40 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_41 =
              -tmp_kernel_op_1 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_42 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_43 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_44 =
              -tmp_kernel_op_5 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_45 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_46 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_47 =
              -tmp_kernel_op_9 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_48 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_49 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_50 =
              -tmp_kernel_op_13 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_51 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_52 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_16;
          const walberla::float64 tmp_kernel_op_53 =
              -tmp_kernel_op_17 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_54 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_16;
          const walberla::float64 tmp_kernel_op_55 =
              jac_affine_inv_1_0_GRAY * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_56 =
              -tmp_kernel_op_21 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_57 =
              jac_affine_inv_1_1_GRAY * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_58 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_59 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_60 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_61 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_62 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_63 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_64 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_65 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_66 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_67 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_68 =
              jac_affine_inv_0_0_GRAY * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_69 =
              jac_affine_inv_0_1_GRAY * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_70 =
              -tmp_kernel_op_0 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_71 =
              -tmp_kernel_op_4 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_72 =
              -tmp_kernel_op_8 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_73 =
              -tmp_kernel_op_12 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_74 =
              -tmp_kernel_op_16 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_75 =
              -tmp_kernel_op_20 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_76 =
              jac_affine_inv_0_0_GRAY * 27.0;
          const int64_t tmp_kernel_op_77 = 0;
          const walberla::float64 tmp_kernel_op_78 =
              jac_affine_inv_1_0_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_79 = 0.15066167873471437;
          const walberla::float64 tmp_kernel_op_80 =
              jac_affine_inv_0_1_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_81 =
              jac_affine_inv_1_1_GRAY * 27.0;
          const walberla::float64 tmp_kernel_op_82 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_83 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_84 = 0.15066167873471437;
          const int64_t tmp_kernel_op_85 = 0;
          const walberla::float64 tmp_kernel_op_86 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_87 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_88 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_89 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_90 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_91 = 0.066417604867409386;
          const walberla::float64 elMatDiag_0 =
              tmp_kernel_op_11 *
                  (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_10 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_10) *
                    (jac_affine_inv_0_0_GRAY * tmp_kernel_op_10 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_10)) +
                   ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_10 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_10) *
                    (jac_affine_inv_0_1_GRAY * tmp_kernel_op_10 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_10))) +
              tmp_kernel_op_15 *
                  (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_14 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_14) *
                    (jac_affine_inv_0_0_GRAY * tmp_kernel_op_14 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_14)) +
                   ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_14 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_14) *
                    (jac_affine_inv_0_1_GRAY * tmp_kernel_op_14 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_14))) +
              tmp_kernel_op_19 *
                  (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_18 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_18) *
                    (jac_affine_inv_0_0_GRAY * tmp_kernel_op_18 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_18)) +
                   ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_18 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_18) *
                    (jac_affine_inv_0_1_GRAY * tmp_kernel_op_18 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_18))) +
              tmp_kernel_op_23 *
                  (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_22 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_22) *
                    (jac_affine_inv_0_0_GRAY * tmp_kernel_op_22 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_22)) +
                   ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_22 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_22) *
                    (jac_affine_inv_0_1_GRAY * tmp_kernel_op_22 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_22))) +
              tmp_kernel_op_3 *
                  (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_2 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_2) *
                    (jac_affine_inv_0_0_GRAY * tmp_kernel_op_2 +
                     jac_affine_inv_1_0_GRAY * tmp_kernel_op_2)) +
                   ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_2) *
                    (jac_affine_inv_0_1_GRAY * tmp_kernel_op_2 +
                     jac_affine_inv_1_1_GRAY * tmp_kernel_op_2))) +
              tmp_kernel_op_7 * (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_6 +
                                   jac_affine_inv_1_0_GRAY * tmp_kernel_op_6) *
                                  (jac_affine_inv_0_0_GRAY * tmp_kernel_op_6 +
                                   jac_affine_inv_1_0_GRAY * tmp_kernel_op_6)) +
                                 ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_6 +
                                   jac_affine_inv_1_1_GRAY * tmp_kernel_op_6) *
                                  (jac_affine_inv_0_1_GRAY * tmp_kernel_op_6 +
                                   jac_affine_inv_1_1_GRAY * tmp_kernel_op_6)));
          const walberla::float64 elMatDiag_1 =
              tmp_kernel_op_11 * (tmp_kernel_op_24 * tmp_kernel_op_28 +
                                  tmp_kernel_op_26 * tmp_kernel_op_28) +
              tmp_kernel_op_15 * (tmp_kernel_op_24 * tmp_kernel_op_29 +
                                  tmp_kernel_op_26 * tmp_kernel_op_29) +
              tmp_kernel_op_19 * (tmp_kernel_op_24 * tmp_kernel_op_30 +
                                  tmp_kernel_op_26 * tmp_kernel_op_30) +
              tmp_kernel_op_23 * (tmp_kernel_op_24 * tmp_kernel_op_31 +
                                  tmp_kernel_op_26 * tmp_kernel_op_31) +
              tmp_kernel_op_3 * (tmp_kernel_op_24 * tmp_kernel_op_25 +
                                 tmp_kernel_op_25 * tmp_kernel_op_26) +
              tmp_kernel_op_7 * (tmp_kernel_op_24 * tmp_kernel_op_27 +
                                 tmp_kernel_op_26 * tmp_kernel_op_27);
          const walberla::float64 elMatDiag_2 =
              tmp_kernel_op_11 * (tmp_kernel_op_32 * tmp_kernel_op_36 +
                                  tmp_kernel_op_34 * tmp_kernel_op_36) +
              tmp_kernel_op_15 * (tmp_kernel_op_32 * tmp_kernel_op_37 +
                                  tmp_kernel_op_34 * tmp_kernel_op_37) +
              tmp_kernel_op_19 * (tmp_kernel_op_32 * tmp_kernel_op_38 +
                                  tmp_kernel_op_34 * tmp_kernel_op_38) +
              tmp_kernel_op_23 * (tmp_kernel_op_32 * tmp_kernel_op_39 +
                                  tmp_kernel_op_34 * tmp_kernel_op_39) +
              tmp_kernel_op_3 * (tmp_kernel_op_32 * tmp_kernel_op_33 +
                                 tmp_kernel_op_33 * tmp_kernel_op_34) +
              tmp_kernel_op_7 * (tmp_kernel_op_32 * tmp_kernel_op_35 +
                                 tmp_kernel_op_34 * tmp_kernel_op_35);
          const walberla::float64 elMatDiag_3 =
              tmp_kernel_op_11 * (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_47 -
                                    tmp_kernel_op_46) *
                                   (jac_affine_inv_0_0_GRAY * tmp_kernel_op_47 -
                                    tmp_kernel_op_46)) +
                                  ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_47 -
                                    tmp_kernel_op_48) *
                                   (jac_affine_inv_0_1_GRAY * tmp_kernel_op_47 -
                                    tmp_kernel_op_48))) +
              tmp_kernel_op_15 * (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_50 -
                                    tmp_kernel_op_49) *
                                   (jac_affine_inv_0_0_GRAY * tmp_kernel_op_50 -
                                    tmp_kernel_op_49)) +
                                  ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_50 -
                                    tmp_kernel_op_51) *
                                   (jac_affine_inv_0_1_GRAY * tmp_kernel_op_50 -
                                    tmp_kernel_op_51))) +
              tmp_kernel_op_19 * (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_53 -
                                    tmp_kernel_op_52) *
                                   (jac_affine_inv_0_0_GRAY * tmp_kernel_op_53 -
                                    tmp_kernel_op_52)) +
                                  ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_53 -
                                    tmp_kernel_op_54) *
                                   (jac_affine_inv_0_1_GRAY * tmp_kernel_op_53 -
                                    tmp_kernel_op_54))) +
              tmp_kernel_op_23 * (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_56 -
                                    tmp_kernel_op_55) *
                                   (jac_affine_inv_0_0_GRAY * tmp_kernel_op_56 -
                                    tmp_kernel_op_55)) +
                                  ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_56 -
                                    tmp_kernel_op_57) *
                                   (jac_affine_inv_0_1_GRAY * tmp_kernel_op_56 -
                                    tmp_kernel_op_57))) +
              tmp_kernel_op_3 * (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_41 -
                                   tmp_kernel_op_40) *
                                  (jac_affine_inv_0_0_GRAY * tmp_kernel_op_41 -
                                   tmp_kernel_op_40)) +
                                 ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_41 -
                                   tmp_kernel_op_42) *
                                  (jac_affine_inv_0_1_GRAY * tmp_kernel_op_41 -
                                   tmp_kernel_op_42))) +
              tmp_kernel_op_7 * (((jac_affine_inv_0_0_GRAY * tmp_kernel_op_44 -
                                   tmp_kernel_op_43) *
                                  (jac_affine_inv_0_0_GRAY * tmp_kernel_op_44 -
                                   tmp_kernel_op_43)) +
                                 ((jac_affine_inv_0_1_GRAY * tmp_kernel_op_44 -
                                   tmp_kernel_op_45) *
                                  (jac_affine_inv_0_1_GRAY * tmp_kernel_op_44 -
                                   tmp_kernel_op_45)));
          const walberla::float64 elMatDiag_4 =
              tmp_kernel_op_11 * (((tmp_kernel_op_46 + tmp_kernel_op_62) *
                                   (tmp_kernel_op_46 + tmp_kernel_op_62)) +
                                  ((tmp_kernel_op_48 + tmp_kernel_op_63) *
                                   (tmp_kernel_op_48 + tmp_kernel_op_63))) +
              tmp_kernel_op_15 * (((tmp_kernel_op_49 + tmp_kernel_op_64) *
                                   (tmp_kernel_op_49 + tmp_kernel_op_64)) +
                                  ((tmp_kernel_op_51 + tmp_kernel_op_65) *
                                   (tmp_kernel_op_51 + tmp_kernel_op_65))) +
              tmp_kernel_op_19 * (((tmp_kernel_op_52 + tmp_kernel_op_66) *
                                   (tmp_kernel_op_52 + tmp_kernel_op_66)) +
                                  ((tmp_kernel_op_54 + tmp_kernel_op_67) *
                                   (tmp_kernel_op_54 + tmp_kernel_op_67))) +
              tmp_kernel_op_23 * (((tmp_kernel_op_55 + tmp_kernel_op_68) *
                                   (tmp_kernel_op_55 + tmp_kernel_op_68)) +
                                  ((tmp_kernel_op_57 + tmp_kernel_op_69) *
                                   (tmp_kernel_op_57 + tmp_kernel_op_69))) +
              tmp_kernel_op_3 * (((tmp_kernel_op_40 + tmp_kernel_op_58) *
                                  (tmp_kernel_op_40 + tmp_kernel_op_58)) +
                                 ((tmp_kernel_op_42 + tmp_kernel_op_59) *
                                  (tmp_kernel_op_42 + tmp_kernel_op_59))) +
              tmp_kernel_op_7 * (((tmp_kernel_op_43 + tmp_kernel_op_60) *
                                  (tmp_kernel_op_43 + tmp_kernel_op_60)) +
                                 ((tmp_kernel_op_45 + tmp_kernel_op_61) *
                                  (tmp_kernel_op_45 + tmp_kernel_op_61)));
          const walberla::float64 elMatDiag_5 =
              tmp_kernel_op_11 * (((jac_affine_inv_1_0_GRAY * tmp_kernel_op_72 -
                                    tmp_kernel_op_62) *
                                   (jac_affine_inv_1_0_GRAY * tmp_kernel_op_72 -
                                    tmp_kernel_op_62)) +
                                  ((jac_affine_inv_1_1_GRAY * tmp_kernel_op_72 -
                                    tmp_kernel_op_63) *
                                   (jac_affine_inv_1_1_GRAY * tmp_kernel_op_72 -
                                    tmp_kernel_op_63))) +
              tmp_kernel_op_15 * (((jac_affine_inv_1_0_GRAY * tmp_kernel_op_73 -
                                    tmp_kernel_op_64) *
                                   (jac_affine_inv_1_0_GRAY * tmp_kernel_op_73 -
                                    tmp_kernel_op_64)) +
                                  ((jac_affine_inv_1_1_GRAY * tmp_kernel_op_73 -
                                    tmp_kernel_op_65) *
                                   (jac_affine_inv_1_1_GRAY * tmp_kernel_op_73 -
                                    tmp_kernel_op_65))) +
              tmp_kernel_op_19 * (((jac_affine_inv_1_0_GRAY * tmp_kernel_op_74 -
                                    tmp_kernel_op_66) *
                                   (jac_affine_inv_1_0_GRAY * tmp_kernel_op_74 -
                                    tmp_kernel_op_66)) +
                                  ((jac_affine_inv_1_1_GRAY * tmp_kernel_op_74 -
                                    tmp_kernel_op_67) *
                                   (jac_affine_inv_1_1_GRAY * tmp_kernel_op_74 -
                                    tmp_kernel_op_67))) +
              tmp_kernel_op_23 * (((jac_affine_inv_1_0_GRAY * tmp_kernel_op_75 -
                                    tmp_kernel_op_68) *
                                   (jac_affine_inv_1_0_GRAY * tmp_kernel_op_75 -
                                    tmp_kernel_op_68)) +
                                  ((jac_affine_inv_1_1_GRAY * tmp_kernel_op_75 -
                                    tmp_kernel_op_69) *
                                   (jac_affine_inv_1_1_GRAY * tmp_kernel_op_75 -
                                    tmp_kernel_op_69))) +
              tmp_kernel_op_3 * (((jac_affine_inv_1_0_GRAY * tmp_kernel_op_70 -
                                   tmp_kernel_op_58) *
                                  (jac_affine_inv_1_0_GRAY * tmp_kernel_op_70 -
                                   tmp_kernel_op_58)) +
                                 ((jac_affine_inv_1_1_GRAY * tmp_kernel_op_70 -
                                   tmp_kernel_op_59) *
                                  (jac_affine_inv_1_1_GRAY * tmp_kernel_op_70 -
                                   tmp_kernel_op_59))) +
              tmp_kernel_op_7 * (((jac_affine_inv_1_0_GRAY * tmp_kernel_op_71 -
                                   tmp_kernel_op_60) *
                                  (jac_affine_inv_1_0_GRAY * tmp_kernel_op_71 -
                                   tmp_kernel_op_60)) +
                                 ((jac_affine_inv_1_1_GRAY * tmp_kernel_op_71 -
                                   tmp_kernel_op_61) *
                                  (jac_affine_inv_1_1_GRAY * tmp_kernel_op_71 -
                                   tmp_kernel_op_61)));
          const walberla::float64 elMatDiag_6 =
              tmp_kernel_op_11 *
                  (((tmp_kernel_op_76 * tmp_kernel_op_84 +
                     tmp_kernel_op_78 *
                         ((walberla::float64)(tmp_kernel_op_85))) *
                    (tmp_kernel_op_76 * tmp_kernel_op_84 +
                     tmp_kernel_op_78 *
                         ((walberla::float64)(tmp_kernel_op_85)))) +
                   ((tmp_kernel_op_80 * tmp_kernel_op_84 +
                     tmp_kernel_op_81 *
                         ((walberla::float64)(tmp_kernel_op_85))) *
                    (tmp_kernel_op_80 * tmp_kernel_op_84 +
                     tmp_kernel_op_81 *
                         ((walberla::float64)(tmp_kernel_op_85))))) +
              tmp_kernel_op_15 * (((tmp_kernel_op_76 * tmp_kernel_op_86 +
                                    tmp_kernel_op_78 * tmp_kernel_op_87) *
                                   (tmp_kernel_op_76 * tmp_kernel_op_86 +
                                    tmp_kernel_op_78 * tmp_kernel_op_87)) +
                                  ((tmp_kernel_op_80 * tmp_kernel_op_86 +
                                    tmp_kernel_op_81 * tmp_kernel_op_87) *
                                   (tmp_kernel_op_80 * tmp_kernel_op_86 +
                                    tmp_kernel_op_81 * tmp_kernel_op_87))) +
              tmp_kernel_op_19 * (((tmp_kernel_op_76 * tmp_kernel_op_88 +
                                    tmp_kernel_op_78 * tmp_kernel_op_89) *
                                   (tmp_kernel_op_76 * tmp_kernel_op_88 +
                                    tmp_kernel_op_78 * tmp_kernel_op_89)) +
                                  ((tmp_kernel_op_80 * tmp_kernel_op_88 +
                                    tmp_kernel_op_81 * tmp_kernel_op_89) *
                                   (tmp_kernel_op_80 * tmp_kernel_op_88 +
                                    tmp_kernel_op_81 * tmp_kernel_op_89))) +
              tmp_kernel_op_23 * (((tmp_kernel_op_76 * tmp_kernel_op_90 +
                                    tmp_kernel_op_78 * tmp_kernel_op_91) *
                                   (tmp_kernel_op_76 * tmp_kernel_op_90 +
                                    tmp_kernel_op_78 * tmp_kernel_op_91)) +
                                  ((tmp_kernel_op_80 * tmp_kernel_op_90 +
                                    tmp_kernel_op_81 * tmp_kernel_op_91) *
                                   (tmp_kernel_op_80 * tmp_kernel_op_90 +
                                    tmp_kernel_op_81 * tmp_kernel_op_91))) +
              tmp_kernel_op_3 *
                  (((tmp_kernel_op_76 *
                         ((walberla::float64)(tmp_kernel_op_77)) +
                     tmp_kernel_op_78 * tmp_kernel_op_79) *
                    (tmp_kernel_op_76 *
                         ((walberla::float64)(tmp_kernel_op_77)) +
                     tmp_kernel_op_78 * tmp_kernel_op_79)) +
                   ((tmp_kernel_op_79 * tmp_kernel_op_81 +
                     tmp_kernel_op_80 *
                         ((walberla::float64)(tmp_kernel_op_77))) *
                    (tmp_kernel_op_79 * tmp_kernel_op_81 +
                     tmp_kernel_op_80 *
                         ((walberla::float64)(tmp_kernel_op_77))))) +
              tmp_kernel_op_7 * (((tmp_kernel_op_76 * tmp_kernel_op_82 +
                                   tmp_kernel_op_78 * tmp_kernel_op_83) *
                                  (tmp_kernel_op_76 * tmp_kernel_op_82 +
                                   tmp_kernel_op_78 * tmp_kernel_op_83)) +
                                 ((tmp_kernel_op_80 * tmp_kernel_op_82 +
                                   tmp_kernel_op_81 * tmp_kernel_op_83) *
                                  (tmp_kernel_op_80 * tmp_kernel_op_82 +
                                   tmp_kernel_op_81 * tmp_kernel_op_83)));
          _data_invDiag_Vertex[ctr_0 +
                               ctr_1 * (micro_edges_per_macro_edge + 2) -
                               ((ctr_1 * (ctr_1 + 1)) / (2))] =
              elMatDiag_0 +
              _data_invDiag_Vertex[ctr_0 +
                                   ctr_1 * (micro_edges_per_macro_edge + 2) -
                                   ((ctr_1 * (ctr_1 + 1)) / (2))];
          _data_invDiag_Vertex[ctr_0 +
                               ctr_1 * (micro_edges_per_macro_edge + 2) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
              elMatDiag_1 +
              _data_invDiag_Vertex[ctr_0 +
                                   ctr_1 * (micro_edges_per_macro_edge + 2) -
                                   ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          _data_invDiag_Vertex[ctr_0 +
                               (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatDiag_2 +
              _data_invDiag_Vertex[ctr_0 +
                                   (ctr_1 + 1) *
                                       (micro_edges_per_macro_edge + 2) -
                                   (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_invDiag_Edge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) +
                             ((micro_edges_per_macro_edge *
                               (micro_edges_per_macro_edge + 1)) /
                              (2))] =
              elMatDiag_3 +
              _data_invDiag_Edge[ctr_0 +
                                 ctr_1 * (micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 ((micro_edges_per_macro_edge *
                                   (micro_edges_per_macro_edge + 1)) /
                                  (2))];
          _data_invDiag_Edge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) +
                             2 * ((micro_edges_per_macro_edge *
                                   (micro_edges_per_macro_edge + 1)) /
                                  (2))] =
              elMatDiag_4 +
              _data_invDiag_Edge[ctr_0 +
                                 ctr_1 * (micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 2 * ((micro_edges_per_macro_edge *
                                       (micro_edges_per_macro_edge + 1)) /
                                      (2))];
          _data_invDiag_Edge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                             ((ctr_1 * (ctr_1 + 1)) / (2))] =
              elMatDiag_5 +
              _data_invDiag_Edge[ctr_0 +
                                 ctr_1 * (micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2))];
          _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                         ((ctr_1 * (ctr_1 + 1)) / (2))] =
              elMatDiag_6 +
              _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                             ((ctr_1 * (ctr_1 + 1)) / (2))];
        }
    }
    {
      /* FaceType.BLUE */
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
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1;
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
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_0_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1));
          const walberla::float64 p_affine_1_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_1_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_0 =
              macro_vertex_coord_id_0comp0 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_1comp0) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp0 +
                   macro_vertex_coord_id_2comp0) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 p_affine_2_1 =
              macro_vertex_coord_id_0comp1 +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_1comp1) *
                  1.0 * ((walberla::float64)(ctr_0 + 1)) +
              1.0 / (micro_edges_per_macro_edge_float) *
                  (-macro_vertex_coord_id_0comp1 +
                   macro_vertex_coord_id_2comp1) *
                  1.0 * ((walberla::float64)(ctr_1 + 1));
          const walberla::float64 tmp_kernel_op_0 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_1 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_2 =
              tmp_kernel_op_0 + tmp_kernel_op_1 - 3.0;
          const walberla::float64 tmp_kernel_op_3 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_4 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_5 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_6 =
              tmp_kernel_op_4 + tmp_kernel_op_5 - 3.0;
          const walberla::float64 tmp_kernel_op_7 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_8 = 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_9 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_10 =
              tmp_kernel_op_8 + tmp_kernel_op_9 - 3.0;
          const walberla::float64 tmp_kernel_op_11 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_12 = 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_13 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_12 + tmp_kernel_op_13 - 3.0;
          const walberla::float64 tmp_kernel_op_15 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_16 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_17 = 1.7837939636638596;
          const walberla::float64 tmp_kernel_op_18 =
              tmp_kernel_op_16 + tmp_kernel_op_17 - 3.0;
          const walberla::float64 tmp_kernel_op_19 =
              abs_det_jac_affine_BLUE * 0.11169079483900581;
          const walberla::float64 tmp_kernel_op_20 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_21 = 0.36630485403908286;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_20 + tmp_kernel_op_21 - 3.0;
          const walberla::float64 tmp_kernel_op_23 =
              abs_det_jac_affine_BLUE * 0.054975871827660949;
          const walberla::float64 tmp_kernel_op_24 =
              (jac_affine_inv_0_0_BLUE * jac_affine_inv_0_0_BLUE);
          const walberla::float64 tmp_kernel_op_25 =
              ((tmp_kernel_op_0 - 1.0) * (tmp_kernel_op_0 - 1.0));
          const walberla::float64 tmp_kernel_op_26 =
              (jac_affine_inv_0_1_BLUE * jac_affine_inv_0_1_BLUE);
          const walberla::float64 tmp_kernel_op_27 =
              ((tmp_kernel_op_4 - 1.0) * (tmp_kernel_op_4 - 1.0));
          const walberla::float64 tmp_kernel_op_28 =
              ((tmp_kernel_op_8 - 1.0) * (tmp_kernel_op_8 - 1.0));
          const walberla::float64 tmp_kernel_op_29 =
              ((tmp_kernel_op_12 - 1.0) * (tmp_kernel_op_12 - 1.0));
          const walberla::float64 tmp_kernel_op_30 =
              ((tmp_kernel_op_16 - 1.0) * (tmp_kernel_op_16 - 1.0));
          const walberla::float64 tmp_kernel_op_31 =
              ((tmp_kernel_op_20 - 1.0) * (tmp_kernel_op_20 - 1.0));
          const walberla::float64 tmp_kernel_op_32 =
              (jac_affine_inv_1_0_BLUE * jac_affine_inv_1_0_BLUE);
          const walberla::float64 tmp_kernel_op_33 =
              ((tmp_kernel_op_1 - 1.0) * (tmp_kernel_op_1 - 1.0));
          const walberla::float64 tmp_kernel_op_34 =
              (jac_affine_inv_1_1_BLUE * jac_affine_inv_1_1_BLUE);
          const walberla::float64 tmp_kernel_op_35 =
              ((tmp_kernel_op_5 - 1.0) * (tmp_kernel_op_5 - 1.0));
          const walberla::float64 tmp_kernel_op_36 =
              ((tmp_kernel_op_9 - 1.0) * (tmp_kernel_op_9 - 1.0));
          const walberla::float64 tmp_kernel_op_37 =
              ((tmp_kernel_op_13 - 1.0) * (tmp_kernel_op_13 - 1.0));
          const walberla::float64 tmp_kernel_op_38 =
              ((tmp_kernel_op_17 - 1.0) * (tmp_kernel_op_17 - 1.0));
          const walberla::float64 tmp_kernel_op_39 =
              ((tmp_kernel_op_21 - 1.0) * (tmp_kernel_op_21 - 1.0));
          const walberla::float64 tmp_kernel_op_40 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_41 =
              -tmp_kernel_op_1 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_42 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_0;
          const walberla::float64 tmp_kernel_op_43 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_44 =
              -tmp_kernel_op_5 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_45 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_46 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_47 =
              -tmp_kernel_op_9 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_48 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_49 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_50 =
              -tmp_kernel_op_13 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_51 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_52 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_16;
          const walberla::float64 tmp_kernel_op_53 =
              -tmp_kernel_op_17 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_54 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_16;
          const walberla::float64 tmp_kernel_op_55 =
              jac_affine_inv_1_0_BLUE * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_56 =
              -tmp_kernel_op_21 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_57 =
              jac_affine_inv_1_1_BLUE * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_58 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_59 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_1;
          const walberla::float64 tmp_kernel_op_60 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_61 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_62 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_63 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_64 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_65 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_66 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_67 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_68 =
              jac_affine_inv_0_0_BLUE * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_69 =
              jac_affine_inv_0_1_BLUE * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_70 =
              -tmp_kernel_op_0 + 3.1351758546554382;
          const walberla::float64 tmp_kernel_op_71 =
              -tmp_kernel_op_4 - 2.5347805838436681;
          const walberla::float64 tmp_kernel_op_72 =
              -tmp_kernel_op_8 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_73 =
              -tmp_kernel_op_12 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_74 =
              -tmp_kernel_op_16 + 0.43241207267228088;
          const walberla::float64 tmp_kernel_op_75 =
              -tmp_kernel_op_20 + 3.2673902919218341;
          const walberla::float64 tmp_kernel_op_76 =
              jac_affine_inv_0_0_BLUE * 27.0;
          const int64_t tmp_kernel_op_77 = 0;
          const walberla::float64 tmp_kernel_op_78 =
              jac_affine_inv_1_0_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_79 = 0.15066167873471437;
          const walberla::float64 tmp_kernel_op_80 =
              jac_affine_inv_0_1_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_81 =
              jac_affine_inv_1_1_BLUE * 27.0;
          const walberla::float64 tmp_kernel_op_82 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_83 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_84 = 0.15066167873471437;
          const int64_t tmp_kernel_op_85 = 0;
          const walberla::float64 tmp_kernel_op_86 = -0.066417604867409372;
          const walberla::float64 tmp_kernel_op_87 = 4.5344149156604147e-17;
          const walberla::float64 tmp_kernel_op_88 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_89 = -0.15066167873471437;
          const walberla::float64 tmp_kernel_op_90 = 0.066417604867409386;
          const walberla::float64 tmp_kernel_op_91 = 0.066417604867409386;
          const walberla::float64 elMatDiag_0 =
              tmp_kernel_op_11 *
                  (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_10 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_10) *
                    (jac_affine_inv_0_0_BLUE * tmp_kernel_op_10 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_10)) +
                   ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_10 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_10) *
                    (jac_affine_inv_0_1_BLUE * tmp_kernel_op_10 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_10))) +
              tmp_kernel_op_15 *
                  (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_14 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_14) *
                    (jac_affine_inv_0_0_BLUE * tmp_kernel_op_14 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_14)) +
                   ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_14 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_14) *
                    (jac_affine_inv_0_1_BLUE * tmp_kernel_op_14 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_14))) +
              tmp_kernel_op_19 *
                  (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_18 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_18) *
                    (jac_affine_inv_0_0_BLUE * tmp_kernel_op_18 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_18)) +
                   ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_18 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_18) *
                    (jac_affine_inv_0_1_BLUE * tmp_kernel_op_18 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_18))) +
              tmp_kernel_op_23 *
                  (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_22 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_22) *
                    (jac_affine_inv_0_0_BLUE * tmp_kernel_op_22 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_22)) +
                   ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_22 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_22) *
                    (jac_affine_inv_0_1_BLUE * tmp_kernel_op_22 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_22))) +
              tmp_kernel_op_3 *
                  (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_2 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_2) *
                    (jac_affine_inv_0_0_BLUE * tmp_kernel_op_2 +
                     jac_affine_inv_1_0_BLUE * tmp_kernel_op_2)) +
                   ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_2) *
                    (jac_affine_inv_0_1_BLUE * tmp_kernel_op_2 +
                     jac_affine_inv_1_1_BLUE * tmp_kernel_op_2))) +
              tmp_kernel_op_7 * (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_6 +
                                   jac_affine_inv_1_0_BLUE * tmp_kernel_op_6) *
                                  (jac_affine_inv_0_0_BLUE * tmp_kernel_op_6 +
                                   jac_affine_inv_1_0_BLUE * tmp_kernel_op_6)) +
                                 ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_6 +
                                   jac_affine_inv_1_1_BLUE * tmp_kernel_op_6) *
                                  (jac_affine_inv_0_1_BLUE * tmp_kernel_op_6 +
                                   jac_affine_inv_1_1_BLUE * tmp_kernel_op_6)));
          const walberla::float64 elMatDiag_1 =
              tmp_kernel_op_11 * (tmp_kernel_op_24 * tmp_kernel_op_28 +
                                  tmp_kernel_op_26 * tmp_kernel_op_28) +
              tmp_kernel_op_15 * (tmp_kernel_op_24 * tmp_kernel_op_29 +
                                  tmp_kernel_op_26 * tmp_kernel_op_29) +
              tmp_kernel_op_19 * (tmp_kernel_op_24 * tmp_kernel_op_30 +
                                  tmp_kernel_op_26 * tmp_kernel_op_30) +
              tmp_kernel_op_23 * (tmp_kernel_op_24 * tmp_kernel_op_31 +
                                  tmp_kernel_op_26 * tmp_kernel_op_31) +
              tmp_kernel_op_3 * (tmp_kernel_op_24 * tmp_kernel_op_25 +
                                 tmp_kernel_op_25 * tmp_kernel_op_26) +
              tmp_kernel_op_7 * (tmp_kernel_op_24 * tmp_kernel_op_27 +
                                 tmp_kernel_op_26 * tmp_kernel_op_27);
          const walberla::float64 elMatDiag_2 =
              tmp_kernel_op_11 * (tmp_kernel_op_32 * tmp_kernel_op_36 +
                                  tmp_kernel_op_34 * tmp_kernel_op_36) +
              tmp_kernel_op_15 * (tmp_kernel_op_32 * tmp_kernel_op_37 +
                                  tmp_kernel_op_34 * tmp_kernel_op_37) +
              tmp_kernel_op_19 * (tmp_kernel_op_32 * tmp_kernel_op_38 +
                                  tmp_kernel_op_34 * tmp_kernel_op_38) +
              tmp_kernel_op_23 * (tmp_kernel_op_32 * tmp_kernel_op_39 +
                                  tmp_kernel_op_34 * tmp_kernel_op_39) +
              tmp_kernel_op_3 * (tmp_kernel_op_32 * tmp_kernel_op_33 +
                                 tmp_kernel_op_33 * tmp_kernel_op_34) +
              tmp_kernel_op_7 * (tmp_kernel_op_32 * tmp_kernel_op_35 +
                                 tmp_kernel_op_34 * tmp_kernel_op_35);
          const walberla::float64 elMatDiag_3 =
              tmp_kernel_op_11 * (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_47 -
                                    tmp_kernel_op_46) *
                                   (jac_affine_inv_0_0_BLUE * tmp_kernel_op_47 -
                                    tmp_kernel_op_46)) +
                                  ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_47 -
                                    tmp_kernel_op_48) *
                                   (jac_affine_inv_0_1_BLUE * tmp_kernel_op_47 -
                                    tmp_kernel_op_48))) +
              tmp_kernel_op_15 * (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_50 -
                                    tmp_kernel_op_49) *
                                   (jac_affine_inv_0_0_BLUE * tmp_kernel_op_50 -
                                    tmp_kernel_op_49)) +
                                  ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_50 -
                                    tmp_kernel_op_51) *
                                   (jac_affine_inv_0_1_BLUE * tmp_kernel_op_50 -
                                    tmp_kernel_op_51))) +
              tmp_kernel_op_19 * (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_53 -
                                    tmp_kernel_op_52) *
                                   (jac_affine_inv_0_0_BLUE * tmp_kernel_op_53 -
                                    tmp_kernel_op_52)) +
                                  ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_53 -
                                    tmp_kernel_op_54) *
                                   (jac_affine_inv_0_1_BLUE * tmp_kernel_op_53 -
                                    tmp_kernel_op_54))) +
              tmp_kernel_op_23 * (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_56 -
                                    tmp_kernel_op_55) *
                                   (jac_affine_inv_0_0_BLUE * tmp_kernel_op_56 -
                                    tmp_kernel_op_55)) +
                                  ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_56 -
                                    tmp_kernel_op_57) *
                                   (jac_affine_inv_0_1_BLUE * tmp_kernel_op_56 -
                                    tmp_kernel_op_57))) +
              tmp_kernel_op_3 * (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_41 -
                                   tmp_kernel_op_40) *
                                  (jac_affine_inv_0_0_BLUE * tmp_kernel_op_41 -
                                   tmp_kernel_op_40)) +
                                 ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_41 -
                                   tmp_kernel_op_42) *
                                  (jac_affine_inv_0_1_BLUE * tmp_kernel_op_41 -
                                   tmp_kernel_op_42))) +
              tmp_kernel_op_7 * (((jac_affine_inv_0_0_BLUE * tmp_kernel_op_44 -
                                   tmp_kernel_op_43) *
                                  (jac_affine_inv_0_0_BLUE * tmp_kernel_op_44 -
                                   tmp_kernel_op_43)) +
                                 ((jac_affine_inv_0_1_BLUE * tmp_kernel_op_44 -
                                   tmp_kernel_op_45) *
                                  (jac_affine_inv_0_1_BLUE * tmp_kernel_op_44 -
                                   tmp_kernel_op_45)));
          const walberla::float64 elMatDiag_4 =
              tmp_kernel_op_11 * (((tmp_kernel_op_46 + tmp_kernel_op_62) *
                                   (tmp_kernel_op_46 + tmp_kernel_op_62)) +
                                  ((tmp_kernel_op_48 + tmp_kernel_op_63) *
                                   (tmp_kernel_op_48 + tmp_kernel_op_63))) +
              tmp_kernel_op_15 * (((tmp_kernel_op_49 + tmp_kernel_op_64) *
                                   (tmp_kernel_op_49 + tmp_kernel_op_64)) +
                                  ((tmp_kernel_op_51 + tmp_kernel_op_65) *
                                   (tmp_kernel_op_51 + tmp_kernel_op_65))) +
              tmp_kernel_op_19 * (((tmp_kernel_op_52 + tmp_kernel_op_66) *
                                   (tmp_kernel_op_52 + tmp_kernel_op_66)) +
                                  ((tmp_kernel_op_54 + tmp_kernel_op_67) *
                                   (tmp_kernel_op_54 + tmp_kernel_op_67))) +
              tmp_kernel_op_23 * (((tmp_kernel_op_55 + tmp_kernel_op_68) *
                                   (tmp_kernel_op_55 + tmp_kernel_op_68)) +
                                  ((tmp_kernel_op_57 + tmp_kernel_op_69) *
                                   (tmp_kernel_op_57 + tmp_kernel_op_69))) +
              tmp_kernel_op_3 * (((tmp_kernel_op_40 + tmp_kernel_op_58) *
                                  (tmp_kernel_op_40 + tmp_kernel_op_58)) +
                                 ((tmp_kernel_op_42 + tmp_kernel_op_59) *
                                  (tmp_kernel_op_42 + tmp_kernel_op_59))) +
              tmp_kernel_op_7 * (((tmp_kernel_op_43 + tmp_kernel_op_60) *
                                  (tmp_kernel_op_43 + tmp_kernel_op_60)) +
                                 ((tmp_kernel_op_45 + tmp_kernel_op_61) *
                                  (tmp_kernel_op_45 + tmp_kernel_op_61)));
          const walberla::float64 elMatDiag_5 =
              tmp_kernel_op_11 * (((jac_affine_inv_1_0_BLUE * tmp_kernel_op_72 -
                                    tmp_kernel_op_62) *
                                   (jac_affine_inv_1_0_BLUE * tmp_kernel_op_72 -
                                    tmp_kernel_op_62)) +
                                  ((jac_affine_inv_1_1_BLUE * tmp_kernel_op_72 -
                                    tmp_kernel_op_63) *
                                   (jac_affine_inv_1_1_BLUE * tmp_kernel_op_72 -
                                    tmp_kernel_op_63))) +
              tmp_kernel_op_15 * (((jac_affine_inv_1_0_BLUE * tmp_kernel_op_73 -
                                    tmp_kernel_op_64) *
                                   (jac_affine_inv_1_0_BLUE * tmp_kernel_op_73 -
                                    tmp_kernel_op_64)) +
                                  ((jac_affine_inv_1_1_BLUE * tmp_kernel_op_73 -
                                    tmp_kernel_op_65) *
                                   (jac_affine_inv_1_1_BLUE * tmp_kernel_op_73 -
                                    tmp_kernel_op_65))) +
              tmp_kernel_op_19 * (((jac_affine_inv_1_0_BLUE * tmp_kernel_op_74 -
                                    tmp_kernel_op_66) *
                                   (jac_affine_inv_1_0_BLUE * tmp_kernel_op_74 -
                                    tmp_kernel_op_66)) +
                                  ((jac_affine_inv_1_1_BLUE * tmp_kernel_op_74 -
                                    tmp_kernel_op_67) *
                                   (jac_affine_inv_1_1_BLUE * tmp_kernel_op_74 -
                                    tmp_kernel_op_67))) +
              tmp_kernel_op_23 * (((jac_affine_inv_1_0_BLUE * tmp_kernel_op_75 -
                                    tmp_kernel_op_68) *
                                   (jac_affine_inv_1_0_BLUE * tmp_kernel_op_75 -
                                    tmp_kernel_op_68)) +
                                  ((jac_affine_inv_1_1_BLUE * tmp_kernel_op_75 -
                                    tmp_kernel_op_69) *
                                   (jac_affine_inv_1_1_BLUE * tmp_kernel_op_75 -
                                    tmp_kernel_op_69))) +
              tmp_kernel_op_3 * (((jac_affine_inv_1_0_BLUE * tmp_kernel_op_70 -
                                   tmp_kernel_op_58) *
                                  (jac_affine_inv_1_0_BLUE * tmp_kernel_op_70 -
                                   tmp_kernel_op_58)) +
                                 ((jac_affine_inv_1_1_BLUE * tmp_kernel_op_70 -
                                   tmp_kernel_op_59) *
                                  (jac_affine_inv_1_1_BLUE * tmp_kernel_op_70 -
                                   tmp_kernel_op_59))) +
              tmp_kernel_op_7 * (((jac_affine_inv_1_0_BLUE * tmp_kernel_op_71 -
                                   tmp_kernel_op_60) *
                                  (jac_affine_inv_1_0_BLUE * tmp_kernel_op_71 -
                                   tmp_kernel_op_60)) +
                                 ((jac_affine_inv_1_1_BLUE * tmp_kernel_op_71 -
                                   tmp_kernel_op_61) *
                                  (jac_affine_inv_1_1_BLUE * tmp_kernel_op_71 -
                                   tmp_kernel_op_61)));
          const walberla::float64 elMatDiag_6 =
              tmp_kernel_op_11 *
                  (((tmp_kernel_op_76 * tmp_kernel_op_84 +
                     tmp_kernel_op_78 *
                         ((walberla::float64)(tmp_kernel_op_85))) *
                    (tmp_kernel_op_76 * tmp_kernel_op_84 +
                     tmp_kernel_op_78 *
                         ((walberla::float64)(tmp_kernel_op_85)))) +
                   ((tmp_kernel_op_80 * tmp_kernel_op_84 +
                     tmp_kernel_op_81 *
                         ((walberla::float64)(tmp_kernel_op_85))) *
                    (tmp_kernel_op_80 * tmp_kernel_op_84 +
                     tmp_kernel_op_81 *
                         ((walberla::float64)(tmp_kernel_op_85))))) +
              tmp_kernel_op_15 * (((tmp_kernel_op_76 * tmp_kernel_op_86 +
                                    tmp_kernel_op_78 * tmp_kernel_op_87) *
                                   (tmp_kernel_op_76 * tmp_kernel_op_86 +
                                    tmp_kernel_op_78 * tmp_kernel_op_87)) +
                                  ((tmp_kernel_op_80 * tmp_kernel_op_86 +
                                    tmp_kernel_op_81 * tmp_kernel_op_87) *
                                   (tmp_kernel_op_80 * tmp_kernel_op_86 +
                                    tmp_kernel_op_81 * tmp_kernel_op_87))) +
              tmp_kernel_op_19 * (((tmp_kernel_op_76 * tmp_kernel_op_88 +
                                    tmp_kernel_op_78 * tmp_kernel_op_89) *
                                   (tmp_kernel_op_76 * tmp_kernel_op_88 +
                                    tmp_kernel_op_78 * tmp_kernel_op_89)) +
                                  ((tmp_kernel_op_80 * tmp_kernel_op_88 +
                                    tmp_kernel_op_81 * tmp_kernel_op_89) *
                                   (tmp_kernel_op_80 * tmp_kernel_op_88 +
                                    tmp_kernel_op_81 * tmp_kernel_op_89))) +
              tmp_kernel_op_23 * (((tmp_kernel_op_76 * tmp_kernel_op_90 +
                                    tmp_kernel_op_78 * tmp_kernel_op_91) *
                                   (tmp_kernel_op_76 * tmp_kernel_op_90 +
                                    tmp_kernel_op_78 * tmp_kernel_op_91)) +
                                  ((tmp_kernel_op_80 * tmp_kernel_op_90 +
                                    tmp_kernel_op_81 * tmp_kernel_op_91) *
                                   (tmp_kernel_op_80 * tmp_kernel_op_90 +
                                    tmp_kernel_op_81 * tmp_kernel_op_91))) +
              tmp_kernel_op_3 *
                  (((tmp_kernel_op_76 *
                         ((walberla::float64)(tmp_kernel_op_77)) +
                     tmp_kernel_op_78 * tmp_kernel_op_79) *
                    (tmp_kernel_op_76 *
                         ((walberla::float64)(tmp_kernel_op_77)) +
                     tmp_kernel_op_78 * tmp_kernel_op_79)) +
                   ((tmp_kernel_op_79 * tmp_kernel_op_81 +
                     tmp_kernel_op_80 *
                         ((walberla::float64)(tmp_kernel_op_77))) *
                    (tmp_kernel_op_79 * tmp_kernel_op_81 +
                     tmp_kernel_op_80 *
                         ((walberla::float64)(tmp_kernel_op_77))))) +
              tmp_kernel_op_7 * (((tmp_kernel_op_76 * tmp_kernel_op_82 +
                                   tmp_kernel_op_78 * tmp_kernel_op_83) *
                                  (tmp_kernel_op_76 * tmp_kernel_op_82 +
                                   tmp_kernel_op_78 * tmp_kernel_op_83)) +
                                 ((tmp_kernel_op_80 * tmp_kernel_op_82 +
                                   tmp_kernel_op_81 * tmp_kernel_op_83) *
                                  (tmp_kernel_op_80 * tmp_kernel_op_82 +
                                   tmp_kernel_op_81 * tmp_kernel_op_83)));
          _data_invDiag_Vertex[ctr_0 +
                               ctr_1 * (micro_edges_per_macro_edge + 2) -
                               ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
              elMatDiag_0 +
              _data_invDiag_Vertex[ctr_0 +
                                   ctr_1 * (micro_edges_per_macro_edge + 2) -
                                   ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          _data_invDiag_Vertex[ctr_0 +
                               (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatDiag_1 +
              _data_invDiag_Vertex[ctr_0 +
                                   (ctr_1 + 1) *
                                       (micro_edges_per_macro_edge + 2) -
                                   (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_invDiag_Vertex[ctr_0 +
                               (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                               (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1] =
              elMatDiag_2 +
              _data_invDiag_Vertex[ctr_0 +
                                   (ctr_1 + 1) *
                                       (micro_edges_per_macro_edge + 2) -
                                   (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
          _data_invDiag_Edge[ctr_0 +
                             (ctr_1 + 1) * (micro_edges_per_macro_edge + 1) -
                             (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatDiag_3 +
              _data_invDiag_Edge[ctr_0 +
                                 (ctr_1 + 1) *
                                     (micro_edges_per_macro_edge + 1) -
                                 (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_invDiag_Edge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) +
                             2 * ((micro_edges_per_macro_edge *
                                   (micro_edges_per_macro_edge + 1)) /
                                  (2)) +
                             1] =
              elMatDiag_4 +
              _data_invDiag_Edge[ctr_0 +
                                 ctr_1 * (micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 2 * ((micro_edges_per_macro_edge *
                                       (micro_edges_per_macro_edge + 1)) /
                                      (2)) +
                                 1];
          _data_invDiag_Edge[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) +
                             ((micro_edges_per_macro_edge *
                               (micro_edges_per_macro_edge + 1)) /
                              (2))] =
              elMatDiag_5 +
              _data_invDiag_Edge[ctr_0 +
                                 ctr_1 * (micro_edges_per_macro_edge + 1) -
                                 ((ctr_1 * (ctr_1 + 1)) / (2)) +
                                 ((micro_edges_per_macro_edge *
                                   (micro_edges_per_macro_edge + 1)) /
                                  (2))];
          _data_invDiag_[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) +
                         ((micro_edges_per_macro_edge *
                           (micro_edges_per_macro_edge + 1)) /
                          (2))] =
              elMatDiag_6 +
              _data_invDiag_[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) +
                             ((micro_edges_per_macro_edge *
                               (micro_edges_per_macro_edge + 1)) /
                              (2))];
        }
    }
  }
}

} // namespace operatorgeneration

} // namespace hyteg
