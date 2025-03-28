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

#include "P2PlusBubbleElementwiseMass_AffineMap2D_float64.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P2PlusBubbleElementwiseMass_AffineMap2D_float64::
    P2PlusBubbleElementwiseMass_AffineMap2D_float64(
        const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
        size_t maxLevel)
    : Operator(storage, minLevel, maxLevel) {}

void P2PlusBubbleElementwiseMass_AffineMap2D_float64::apply(
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
      // const auto num_microfaces_per_face =
      //     (int64_t)levelinfo::num_microfaces_per_face(level);
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
      WALBERLA_CHECK_NOT_NULLPTR(
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap()),
          "This operator requires the AffineMap2D to be registered as "
          "GeometryMap on every macro-face.")
      real_t bMat_00 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(0, 0);
      real_t bMat_01 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(0, 1);
      real_t bMat_10 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(1, 0);
      real_t bMat_11 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(1, 1);
      // real_t bVec_0 =
      //     std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
      //         ->getVector()[0];
      // real_t bVec_1 =
      //     std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
      //         ->getVector()[1];

      this->timingTree_->start("kernel");

      apply_P2PlusBubbleElementwiseMass_AffineMap2D_float64_macro_2D(

          _data_dst, _data_dstEdge, _data_dstVertex, _data_src, _data_srcEdge,
          _data_srcVertex, bMat_00, bMat_01, bMat_10, bMat_11,
          macro_vertex_coord_id_0comp0, macro_vertex_coord_id_0comp1,
          macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
          macro_vertex_coord_id_2comp0, macro_vertex_coord_id_2comp1,
          micro_edges_per_macro_edge, micro_edges_per_macro_edge_float);

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
void P2PlusBubbleElementwiseMass_AffineMap2D_float64::toMatrix(
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
      // const auto num_microfaces_per_face =
      //     (int64_t)levelinfo::num_microfaces_per_face(level);
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
      WALBERLA_CHECK_NOT_NULLPTR(
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap()),
          "This operator requires the AffineMap2D to be registered as "
          "GeometryMap on every macro-face.")
      real_t bMat_00 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(0, 0);
      real_t bMat_01 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(0, 1);
      real_t bMat_10 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(1, 0);
      real_t bMat_11 =
          std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
              ->getMatrix()(1, 1);
      // real_t bVec_0 =
      //     std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
      //         ->getVector()[0];
      // real_t bVec_1 =
      //     std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
      //         ->getVector()[1];

      this->timingTree_->start("kernel");

      toMatrix_P2PlusBubbleElementwiseMass_AffineMap2D_float64_macro_2D(

          _data_dst, _data_dstEdge, _data_dstVertex, _data_src, _data_srcEdge,
          _data_srcVertex, bMat_00, bMat_01, bMat_10, bMat_11,
          macro_vertex_coord_id_0comp0, macro_vertex_coord_id_0comp1,
          macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
          macro_vertex_coord_id_2comp0, macro_vertex_coord_id_2comp1, mat,
          micro_edges_per_macro_edge, micro_edges_per_macro_edge_float);

      this->timingTree_->stop("kernel");
    }
  }
  this->stopTiming("toMatrix");
}
void P2PlusBubbleElementwiseMass_AffineMap2D_float64::
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
        // const auto num_microfaces_per_face =
        //     (int64_t)levelinfo::num_microfaces_per_face(level);
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
        WALBERLA_CHECK_NOT_NULLPTR(
            std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap()),
            "This operator requires the AffineMap2D to be registered as "
            "GeometryMap on every macro-face.")
        real_t bMat_00 =
            std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
                ->getMatrix()(0, 0);
        real_t bMat_01 =
            std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
                ->getMatrix()(0, 1);
        real_t bMat_10 =
            std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
                ->getMatrix()(1, 0);
        real_t bMat_11 =
            std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
                ->getMatrix()(1, 1);
        // real_t bVec_0 =
        //     std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
        //         ->getVector()[0];
        // real_t bVec_1 =
        //     std::dynamic_pointer_cast<AffineMap2D>(face.getGeometryMap())
        //         ->getVector()[1];

        this->timingTree_->start("kernel");

        computeInverseDiagonalOperatorValues_P2PlusBubbleElementwiseMass_AffineMap2D_float64_macro_2D(

            _data_invDiag_, _data_invDiag_Edge, _data_invDiag_Vertex, bMat_00,
            bMat_01, bMat_10, bMat_11, macro_vertex_coord_id_0comp0,
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
P2PlusBubbleElementwiseMass_AffineMap2D_float64::getInverseDiagonalValues()
    const {
  return invDiag_;
}
void P2PlusBubbleElementwiseMass_AffineMap2D_float64::
    apply_P2PlusBubbleElementwiseMass_AffineMap2D_float64_macro_2D(
        walberla::float64 *RESTRICT _data_dst,
        walberla::float64 *RESTRICT _data_dstEdge,
        walberla::float64 *RESTRICT _data_dstVertex,
        walberla::float64 *RESTRICT _data_src,
        walberla::float64 *RESTRICT _data_srcEdge,
        walberla::float64 *RESTRICT _data_srcVertex, walberla::float64 bMat_00,
        walberla::float64 bMat_01, walberla::float64 bMat_10,
        walberla::float64 bMat_11,
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
      const walberla::float64 abs_det_jac_affine_GRAY =
          abs(jac_affine_0_0_GRAY * jac_affine_1_1_GRAY -
              jac_affine_0_1_GRAY * jac_affine_1_0_GRAY);
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge;
             ctr_0 += 1) {
#if 0
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
#endif
          const walberla::float64 jac_blending_1_1_q_9 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_9 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_9 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_9 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_9 =
              jac_blending_0_0_q_9 * jac_blending_1_1_q_9 -
              jac_blending_0_1_q_9 * jac_blending_1_0_q_9;
          const walberla::float64 jac_blending_1_1_q_8 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_8 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_8 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_8 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_8 =
              jac_blending_0_0_q_8 * jac_blending_1_1_q_8 -
              jac_blending_0_1_q_8 * jac_blending_1_0_q_8;
          const walberla::float64 jac_blending_1_1_q_7 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_7 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_7 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_7 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_7 =
              jac_blending_0_0_q_7 * jac_blending_1_1_q_7 -
              jac_blending_0_1_q_7 * jac_blending_1_0_q_7;
          const walberla::float64 jac_blending_1_1_q_6 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_6 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_6 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_6 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_6 =
              jac_blending_0_0_q_6 * jac_blending_1_1_q_6 -
              jac_blending_0_1_q_6 * jac_blending_1_0_q_6;
          const walberla::float64 jac_blending_1_1_q_5 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_5 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_5 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_5 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_5 =
              jac_blending_0_0_q_5 * jac_blending_1_1_q_5 -
              jac_blending_0_1_q_5 * jac_blending_1_0_q_5;
          const walberla::float64 jac_blending_1_1_q_4 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_4 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_4 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_4 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_4 =
              jac_blending_0_0_q_4 * jac_blending_1_1_q_4 -
              jac_blending_0_1_q_4 * jac_blending_1_0_q_4;
          const walberla::float64 jac_blending_1_1_q_3 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_3 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_3 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_3 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_2 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_2 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_2 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_1 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_1 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_1 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_0 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_0 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_0 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
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
          const walberla::float64 tmp_kernel_op_0 = -0.035167988576825245;
          const walberla::float64 tmp_kernel_op_1 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_0 *
                                                    0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_2 = 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_3 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_4 = -1.0175839942884126;
          const walberla::float64 tmp_kernel_op_5 =
              -tmp_kernel_op_3 - tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_6 =
              tmp_kernel_op_2 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_7 =
              tmp_kernel_op_1 * tmp_kernel_op_6 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_8 = 0.25875522422404362;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_1 *
                                                    0.052397656566423402;
          const walberla::float64 tmp_kernel_op_10 = 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_11 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_12 = -0.87062238788797819;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_15 =
              tmp_kernel_op_14 * tmp_kernel_op_9 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_16 = 0.30730286739632839;
          const walberla::float64 tmp_kernel_op_17 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_2 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_18 = 0.40172877323475986;
          const walberla::float64 tmp_kernel_op_19 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_20 = -0.84634856630183575;
          const walberla::float64 tmp_kernel_op_21 =
              -tmp_kernel_op_19 - tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_18 * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_17 * tmp_kernel_op_22 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_24 = 1.1801004774776607;
          const walberla::float64 tmp_kernel_op_25 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_3 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_26 = 0.64933214716985033;
          const walberla::float64 tmp_kernel_op_27 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_28 = -0.40994976126116967;
          const walberla::float64 tmp_kernel_op_29 =
              -tmp_kernel_op_27 - tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_26 * tmp_kernel_op_29;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_25 * tmp_kernel_op_30 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_32 = 5.1547196302936094;
          const walberla::float64 tmp_kernel_op_33 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_4 *
                                                     0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_34 = -0.0087919971442063094;
          const walberla::float64 tmp_kernel_op_35 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_36 = 1.5773598151468047;
          const walberla::float64 tmp_kernel_op_37 =
              -tmp_kernel_op_35 - tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_34 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_39 =
              tmp_kernel_op_33 * tmp_kernel_op_38 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_40 = 3.1529168382674562;
          const walberla::float64 tmp_kernel_op_41 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_5 *
                                                     0.052397656566423402;
          const walberla::float64 tmp_kernel_op_42 = 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_43 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_44 = 0.57645841913372808;
          const walberla::float64 tmp_kernel_op_45 =
              -tmp_kernel_op_43 - tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_46 =
              tmp_kernel_op_42 * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_47 =
              tmp_kernel_op_41 * tmp_kernel_op_46 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_48 = 1.6069150929390392;
          const walberla::float64 tmp_kernel_op_49 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_6 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_50 = 0.076825716849082126;
          const walberla::float64 tmp_kernel_op_51 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_52 = -0.19654245353048039;
          const walberla::float64 tmp_kernel_op_53 =
              -tmp_kernel_op_51 - tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_50 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_49 * tmp_kernel_op_54 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_56 = 2.5973285886794009;
          const walberla::float64 tmp_kernel_op_57 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_7 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_58 = 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_59 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_60 = 0.29866429433970043;
          const walberla::float64 tmp_kernel_op_61 =
              -tmp_kernel_op_59 - tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_58 * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_63 =
              tmp_kernel_op_57 * tmp_kernel_op_62 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_64 = 0.3413326446812428;
          const walberla::float64 tmp_kernel_op_65 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_8 *
                                                     0.045496761795224737;
          const walberla::float64 tmp_kernel_op_66 = 0.085333161170310645;
          const walberla::float64 tmp_kernel_op_67 = 1.6586673553187572;
          const walberla::float64 tmp_kernel_op_68 = -0.8293336776593786;
          const walberla::float64 tmp_kernel_op_69 =
              -tmp_kernel_op_67 - tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_66 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_65 * tmp_kernel_op_70 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_72 = 1.4147621866798579;
          const walberla::float64 tmp_kernel_op_73 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_9 *
                                                     0.10566414783403316;
          const walberla::float64 tmp_kernel_op_74 = 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_75 = 0.58523781332014213;
          const walberla::float64 tmp_kernel_op_76 = -0.29261890666007107;
          const walberla::float64 tmp_kernel_op_77 =
              -tmp_kernel_op_75 - tmp_kernel_op_76;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_74 * tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_73 * tmp_kernel_op_78 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_0 * tmp_kernel_op_7 +
              tmp_kernel_op_15 * tmp_kernel_op_8 +
              tmp_kernel_op_16 * tmp_kernel_op_23 +
              tmp_kernel_op_24 * tmp_kernel_op_31 +
              tmp_kernel_op_32 * tmp_kernel_op_39 +
              tmp_kernel_op_40 * tmp_kernel_op_47 +
              tmp_kernel_op_48 * tmp_kernel_op_55 +
              tmp_kernel_op_56 * tmp_kernel_op_63 +
              tmp_kernel_op_64 * tmp_kernel_op_71 +
              tmp_kernel_op_72 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_81 =
              tmp_kernel_op_1 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_6 * tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_83 =
              tmp_kernel_op_9 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_14 * tmp_kernel_op_83;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_17 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_22 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_87 =
              tmp_kernel_op_25 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_30 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 =
              tmp_kernel_op_33 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_38 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_91 =
              tmp_kernel_op_41 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_46 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_93 =
              tmp_kernel_op_49 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_54 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              tmp_kernel_op_57 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_62 * tmp_kernel_op_95;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_65 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_70 * tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_99 =
              tmp_kernel_op_73 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_78 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_100 * tmp_kernel_op_76 +
              tmp_kernel_op_12 * tmp_kernel_op_84 +
              tmp_kernel_op_20 * tmp_kernel_op_86 +
              tmp_kernel_op_28 * tmp_kernel_op_88 +
              tmp_kernel_op_36 * tmp_kernel_op_90 +
              tmp_kernel_op_4 * tmp_kernel_op_82 +
              tmp_kernel_op_44 * tmp_kernel_op_92 +
              tmp_kernel_op_52 * tmp_kernel_op_94 +
              tmp_kernel_op_60 * tmp_kernel_op_96 +
              tmp_kernel_op_68 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_102 = tmp_kernel_op_3 - 1.0;
          const walberla::float64 tmp_kernel_op_103 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_104 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_105 = tmp_kernel_op_27 - 1.0;
          const walberla::float64 tmp_kernel_op_106 = tmp_kernel_op_35 - 1.0;
          const walberla::float64 tmp_kernel_op_107 = tmp_kernel_op_43 - 1.0;
          const walberla::float64 tmp_kernel_op_108 = tmp_kernel_op_51 - 1.0;
          const walberla::float64 tmp_kernel_op_109 = tmp_kernel_op_59 - 1.0;
          const walberla::float64 tmp_kernel_op_110 = tmp_kernel_op_67 - 1.0;
          const walberla::float64 tmp_kernel_op_111 = tmp_kernel_op_75 - 1.0;
          const walberla::float64 tmp_kernel_op_112 =
              tmp_kernel_op_102 * tmp_kernel_op_7 +
              tmp_kernel_op_103 * tmp_kernel_op_15 +
              tmp_kernel_op_104 * tmp_kernel_op_23 +
              tmp_kernel_op_105 * tmp_kernel_op_31 +
              tmp_kernel_op_106 * tmp_kernel_op_39 +
              tmp_kernel_op_107 * tmp_kernel_op_47 +
              tmp_kernel_op_108 * tmp_kernel_op_55 +
              tmp_kernel_op_109 * tmp_kernel_op_63 +
              tmp_kernel_op_110 * tmp_kernel_op_71 +
              tmp_kernel_op_111 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_113 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_114 =
              -tmp_kernel_op_0 - tmp_kernel_op_113 + 4.0;
          const walberla::float64 tmp_kernel_op_115 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_116 =
              -tmp_kernel_op_115 - tmp_kernel_op_8 + 4.0;
          const walberla::float64 tmp_kernel_op_117 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_118 =
              -tmp_kernel_op_117 - tmp_kernel_op_16 + 4.0;
          const walberla::float64 tmp_kernel_op_119 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_120 =
              -tmp_kernel_op_119 - tmp_kernel_op_24 + 4.0;
          const walberla::float64 tmp_kernel_op_121 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_122 =
              -tmp_kernel_op_121 - tmp_kernel_op_32 + 4.0;
          const walberla::float64 tmp_kernel_op_123 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_124 =
              -tmp_kernel_op_123 - tmp_kernel_op_40 + 4.0;
          const walberla::float64 tmp_kernel_op_125 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_126 =
              -tmp_kernel_op_125 - tmp_kernel_op_48 + 4.0;
          const walberla::float64 tmp_kernel_op_127 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_128 =
              -tmp_kernel_op_127 - tmp_kernel_op_56 + 4.0;
          const walberla::float64 tmp_kernel_op_129 = 3.3173347106375144;
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_129 - tmp_kernel_op_64 + 4.0;
          const walberla::float64 tmp_kernel_op_131 = 1.1704756266402843;
          const walberla::float64 tmp_kernel_op_132 =
              -tmp_kernel_op_131 - tmp_kernel_op_72 + 4.0;
          const walberla::float64 tmp_kernel_op_133 =
              tmp_kernel_op_114 * tmp_kernel_op_7 +
              tmp_kernel_op_116 * tmp_kernel_op_15 +
              tmp_kernel_op_118 * tmp_kernel_op_23 +
              tmp_kernel_op_120 * tmp_kernel_op_31 +
              tmp_kernel_op_122 * tmp_kernel_op_39 +
              tmp_kernel_op_124 * tmp_kernel_op_47 +
              tmp_kernel_op_126 * tmp_kernel_op_55 +
              tmp_kernel_op_128 * tmp_kernel_op_63 +
              tmp_kernel_op_130 * tmp_kernel_op_71 +
              tmp_kernel_op_132 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_134 =
              tmp_kernel_op_100 * tmp_kernel_op_132 +
              tmp_kernel_op_114 * tmp_kernel_op_82 +
              tmp_kernel_op_116 * tmp_kernel_op_84 +
              tmp_kernel_op_118 * tmp_kernel_op_86 +
              tmp_kernel_op_120 * tmp_kernel_op_88 +
              tmp_kernel_op_122 * tmp_kernel_op_90 +
              tmp_kernel_op_124 * tmp_kernel_op_92 +
              tmp_kernel_op_126 * tmp_kernel_op_94 +
              tmp_kernel_op_128 * tmp_kernel_op_96 +
              tmp_kernel_op_130 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_135 = 34.794357504481866;
          const walberla::float64 tmp_kernel_op_136 =
              tmp_kernel_op_81 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_137 = 21.28218865830533;
          const walberla::float64 tmp_kernel_op_138 =
              tmp_kernel_op_83 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_139 = 10.846676877338515;
          const walberla::float64 tmp_kernel_op_140 =
              tmp_kernel_op_85 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_141 = 17.531967973585957;
          const walberla::float64 tmp_kernel_op_142 =
              tmp_kernel_op_87 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_143 = -0.23738392289357257;
          const walberla::float64 tmp_kernel_op_144 =
              tmp_kernel_op_89 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_145 = 1.7465977635122938;
          const walberla::float64 tmp_kernel_op_146 =
              tmp_kernel_op_91 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_147 = 2.0742943549252164;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_93 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_149 = 7.9656782229742085;
          const walberla::float64 tmp_kernel_op_150 =
              tmp_kernel_op_95 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_151 = 2.3039953515983882;
          const walberla::float64 tmp_kernel_op_152 =
              tmp_kernel_op_97 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_153 = 9.5496447600890395;
          const walberla::float64 tmp_kernel_op_154 =
              tmp_kernel_op_99 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_155 =
              tmp_kernel_op_135 * tmp_kernel_op_136 * tmp_kernel_op_6 +
              tmp_kernel_op_137 * tmp_kernel_op_138 * tmp_kernel_op_14 +
              tmp_kernel_op_139 * tmp_kernel_op_140 * tmp_kernel_op_22 +
              tmp_kernel_op_141 * tmp_kernel_op_142 * tmp_kernel_op_30 +
              tmp_kernel_op_143 * tmp_kernel_op_144 * tmp_kernel_op_38 +
              tmp_kernel_op_145 * tmp_kernel_op_146 * tmp_kernel_op_46 +
              tmp_kernel_op_147 * tmp_kernel_op_148 * tmp_kernel_op_54 +
              tmp_kernel_op_149 * tmp_kernel_op_150 * tmp_kernel_op_62 +
              tmp_kernel_op_151 * tmp_kernel_op_152 * tmp_kernel_op_70 +
              tmp_kernel_op_153 * tmp_kernel_op_154 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_1 * 7.7299213783731935e-5;
          const walberla::float64 tmp_kernel_op_157 =
              tmp_kernel_op_9 * 0.0041846416289521935;
          const walberla::float64 tmp_kernel_op_158 =
              tmp_kernel_op_17 * 0.0059021907693753367;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_25 * 0.087039821058937664;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_33 * 1.6606959041833929;
          const walberla::float64 tmp_kernel_op_161 =
              tmp_kernel_op_41 * 0.62130528681440322;
          const walberla::float64 tmp_kernel_op_162 =
              tmp_kernel_op_49 * 0.16138600724470506;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_57 * 0.42163223734820804;
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_65 * 0.0072817483953182219;
          const walberla::float64 tmp_kernel_op_165 =
              tmp_kernel_op_73 * 0.12509700280369832;
          const walberla::float64 tmp_kernel_op_166 =
              tmp_kernel_op_156 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_12 * tmp_kernel_op_157;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_158 * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_169 =
              tmp_kernel_op_159 * tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_170 =
              tmp_kernel_op_160 * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_161 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_162 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_173 =
              tmp_kernel_op_163 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_174 =
              tmp_kernel_op_164 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_165 * tmp_kernel_op_76;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_113 * tmp_kernel_op_166 +
              tmp_kernel_op_115 * tmp_kernel_op_167 +
              tmp_kernel_op_117 * tmp_kernel_op_168 +
              tmp_kernel_op_119 * tmp_kernel_op_169 +
              tmp_kernel_op_121 * tmp_kernel_op_170 +
              tmp_kernel_op_123 * tmp_kernel_op_171 +
              tmp_kernel_op_125 * tmp_kernel_op_172 +
              tmp_kernel_op_127 * tmp_kernel_op_173 +
              tmp_kernel_op_129 * tmp_kernel_op_174 +
              tmp_kernel_op_131 * tmp_kernel_op_175;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_136 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_178 =
              tmp_kernel_op_12 * tmp_kernel_op_138;
          const walberla::float64 tmp_kernel_op_179 =
              tmp_kernel_op_140 * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_180 =
              tmp_kernel_op_142 * tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_181 =
              tmp_kernel_op_144 * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_182 =
              tmp_kernel_op_146 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_183 =
              tmp_kernel_op_148 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_184 =
              tmp_kernel_op_150 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_185 =
              tmp_kernel_op_152 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_154 * tmp_kernel_op_76;
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_102 * tmp_kernel_op_177 +
              tmp_kernel_op_103 * tmp_kernel_op_178 +
              tmp_kernel_op_104 * tmp_kernel_op_179 +
              tmp_kernel_op_105 * tmp_kernel_op_180 +
              tmp_kernel_op_106 * tmp_kernel_op_181 +
              tmp_kernel_op_107 * tmp_kernel_op_182 +
              tmp_kernel_op_108 * tmp_kernel_op_183 +
              tmp_kernel_op_109 * tmp_kernel_op_184 +
              tmp_kernel_op_110 * tmp_kernel_op_185 +
              tmp_kernel_op_111 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_114 * tmp_kernel_op_177 +
              tmp_kernel_op_116 * tmp_kernel_op_178 +
              tmp_kernel_op_118 * tmp_kernel_op_179 +
              tmp_kernel_op_120 * tmp_kernel_op_180 +
              tmp_kernel_op_122 * tmp_kernel_op_181 +
              tmp_kernel_op_124 * tmp_kernel_op_182 +
              tmp_kernel_op_126 * tmp_kernel_op_183 +
              tmp_kernel_op_128 * tmp_kernel_op_184 +
              tmp_kernel_op_130 * tmp_kernel_op_185 +
              tmp_kernel_op_132 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_114 * tmp_kernel_op_166 +
              tmp_kernel_op_116 * tmp_kernel_op_167 +
              tmp_kernel_op_118 * tmp_kernel_op_168 +
              tmp_kernel_op_120 * tmp_kernel_op_169 +
              tmp_kernel_op_122 * tmp_kernel_op_170 +
              tmp_kernel_op_124 * tmp_kernel_op_171 +
              tmp_kernel_op_126 * tmp_kernel_op_172 +
              tmp_kernel_op_128 * tmp_kernel_op_173 +
              tmp_kernel_op_130 * tmp_kernel_op_174 +
              tmp_kernel_op_132 * tmp_kernel_op_175;
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_135 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_191 =
              tmp_kernel_op_137 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_192 =
              tmp_kernel_op_139 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_193 =
              tmp_kernel_op_141 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_194 =
              tmp_kernel_op_143 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_195 =
              tmp_kernel_op_145 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_196 =
              tmp_kernel_op_147 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_197 =
              tmp_kernel_op_149 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_198 =
              tmp_kernel_op_151 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_153 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_166 * tmp_kernel_op_190 +
              tmp_kernel_op_167 * tmp_kernel_op_191 +
              tmp_kernel_op_168 * tmp_kernel_op_192 +
              tmp_kernel_op_169 * tmp_kernel_op_193 +
              tmp_kernel_op_170 * tmp_kernel_op_194 +
              tmp_kernel_op_171 * tmp_kernel_op_195 +
              tmp_kernel_op_172 * tmp_kernel_op_196 +
              tmp_kernel_op_173 * tmp_kernel_op_197 +
              tmp_kernel_op_174 * tmp_kernel_op_198 +
              tmp_kernel_op_175 * tmp_kernel_op_199;
          const walberla::float64 tmp_kernel_op_201 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_1 * tmp_kernel_op_201;
          const walberla::float64 tmp_kernel_op_203 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_203 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_205 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_17 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_207 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_207 * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_209 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_210 =
              tmp_kernel_op_209 * tmp_kernel_op_33;
          const walberla::float64 tmp_kernel_op_211 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_212 =
              tmp_kernel_op_211 * tmp_kernel_op_41;
          const walberla::float64 tmp_kernel_op_213 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_214 =
              tmp_kernel_op_213 * tmp_kernel_op_49;
          const walberla::float64 tmp_kernel_op_215 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_216 =
              tmp_kernel_op_215 * tmp_kernel_op_57;
          const walberla::float64 tmp_kernel_op_217 = 0.68779434890003011;
          const walberla::float64 tmp_kernel_op_218 =
              tmp_kernel_op_217 * tmp_kernel_op_65;
          const walberla::float64 tmp_kernel_op_219 = 0.085625824534935377;
          const walberla::float64 tmp_kernel_op_220 =
              tmp_kernel_op_219 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_221 =
              tmp_kernel_op_102 * tmp_kernel_op_202;
          const walberla::float64 tmp_kernel_op_222 =
              tmp_kernel_op_103 * tmp_kernel_op_204;
          const walberla::float64 tmp_kernel_op_223 =
              tmp_kernel_op_104 * tmp_kernel_op_206;
          const walberla::float64 tmp_kernel_op_224 =
              tmp_kernel_op_105 * tmp_kernel_op_208;
          const walberla::float64 tmp_kernel_op_225 =
              tmp_kernel_op_106 * tmp_kernel_op_210;
          const walberla::float64 tmp_kernel_op_226 =
              tmp_kernel_op_107 * tmp_kernel_op_212;
          const walberla::float64 tmp_kernel_op_227 =
              tmp_kernel_op_108 * tmp_kernel_op_214;
          const walberla::float64 tmp_kernel_op_228 =
              tmp_kernel_op_109 * tmp_kernel_op_216;
          const walberla::float64 tmp_kernel_op_229 =
              tmp_kernel_op_110 * tmp_kernel_op_218;
          const walberla::float64 tmp_kernel_op_230 =
              tmp_kernel_op_111 * tmp_kernel_op_220;
          const walberla::float64 tmp_kernel_op_231 =
              tmp_kernel_op_0 * tmp_kernel_op_221 +
              tmp_kernel_op_16 * tmp_kernel_op_223 +
              tmp_kernel_op_222 * tmp_kernel_op_8 +
              tmp_kernel_op_224 * tmp_kernel_op_24 +
              tmp_kernel_op_225 * tmp_kernel_op_32 +
              tmp_kernel_op_226 * tmp_kernel_op_40 +
              tmp_kernel_op_227 * tmp_kernel_op_48 +
              tmp_kernel_op_228 * tmp_kernel_op_56 +
              tmp_kernel_op_229 * tmp_kernel_op_64 +
              tmp_kernel_op_230 * tmp_kernel_op_72;
          const walberla::float64 tmp_kernel_op_232 =
              tmp_kernel_op_102 * tmp_kernel_op_114 * tmp_kernel_op_136 +
              tmp_kernel_op_103 * tmp_kernel_op_116 * tmp_kernel_op_138 +
              tmp_kernel_op_104 * tmp_kernel_op_118 * tmp_kernel_op_140 +
              tmp_kernel_op_105 * tmp_kernel_op_120 * tmp_kernel_op_142 +
              tmp_kernel_op_106 * tmp_kernel_op_122 * tmp_kernel_op_144 +
              tmp_kernel_op_107 * tmp_kernel_op_124 * tmp_kernel_op_146 +
              tmp_kernel_op_108 * tmp_kernel_op_126 * tmp_kernel_op_148 +
              tmp_kernel_op_109 * tmp_kernel_op_128 * tmp_kernel_op_150 +
              tmp_kernel_op_110 * tmp_kernel_op_130 * tmp_kernel_op_152 +
              tmp_kernel_op_111 * tmp_kernel_op_132 * tmp_kernel_op_154;
          const walberla::float64 tmp_kernel_op_233 =
              tmp_kernel_op_114 * tmp_kernel_op_221 +
              tmp_kernel_op_116 * tmp_kernel_op_222 +
              tmp_kernel_op_118 * tmp_kernel_op_223 +
              tmp_kernel_op_120 * tmp_kernel_op_224 +
              tmp_kernel_op_122 * tmp_kernel_op_225 +
              tmp_kernel_op_124 * tmp_kernel_op_226 +
              tmp_kernel_op_126 * tmp_kernel_op_227 +
              tmp_kernel_op_128 * tmp_kernel_op_228 +
              tmp_kernel_op_130 * tmp_kernel_op_229 +
              tmp_kernel_op_132 * tmp_kernel_op_230;
          const walberla::float64 tmp_kernel_op_234 =
              tmp_kernel_op_135 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_235 =
              tmp_kernel_op_137 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_236 =
              tmp_kernel_op_139 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_237 =
              tmp_kernel_op_141 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_238 =
              tmp_kernel_op_143 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_239 =
              tmp_kernel_op_145 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_240 =
              tmp_kernel_op_147 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_241 =
              tmp_kernel_op_149 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_242 =
              tmp_kernel_op_151 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_243 =
              tmp_kernel_op_153 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_244 =
              tmp_kernel_op_221 * tmp_kernel_op_234 +
              tmp_kernel_op_222 * tmp_kernel_op_235 +
              tmp_kernel_op_223 * tmp_kernel_op_236 +
              tmp_kernel_op_224 * tmp_kernel_op_237 +
              tmp_kernel_op_225 * tmp_kernel_op_238 +
              tmp_kernel_op_226 * tmp_kernel_op_239 +
              tmp_kernel_op_227 * tmp_kernel_op_240 +
              tmp_kernel_op_228 * tmp_kernel_op_241 +
              tmp_kernel_op_229 * tmp_kernel_op_242 +
              tmp_kernel_op_230 * tmp_kernel_op_243;
          const walberla::float64 tmp_kernel_op_245 =
              tmp_kernel_op_156 * tmp_kernel_op_201;
          const walberla::float64 tmp_kernel_op_246 =
              tmp_kernel_op_157 * tmp_kernel_op_203;
          const walberla::float64 tmp_kernel_op_247 =
              tmp_kernel_op_158 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_248 =
              tmp_kernel_op_159 * tmp_kernel_op_207;
          const walberla::float64 tmp_kernel_op_249 =
              tmp_kernel_op_160 * tmp_kernel_op_209;
          const walberla::float64 tmp_kernel_op_250 =
              tmp_kernel_op_161 * tmp_kernel_op_211;
          const walberla::float64 tmp_kernel_op_251 =
              tmp_kernel_op_162 * tmp_kernel_op_213;
          const walberla::float64 tmp_kernel_op_252 =
              tmp_kernel_op_163 * tmp_kernel_op_215;
          const walberla::float64 tmp_kernel_op_253 =
              tmp_kernel_op_164 * tmp_kernel_op_217;
          const walberla::float64 tmp_kernel_op_254 =
              tmp_kernel_op_165 * tmp_kernel_op_219;
          const walberla::float64 tmp_kernel_op_255 =
              tmp_kernel_op_114 * tmp_kernel_op_202;
          const walberla::float64 tmp_kernel_op_256 =
              tmp_kernel_op_116 * tmp_kernel_op_204;
          const walberla::float64 tmp_kernel_op_257 =
              tmp_kernel_op_118 * tmp_kernel_op_206;
          const walberla::float64 tmp_kernel_op_258 =
              tmp_kernel_op_120 * tmp_kernel_op_208;
          const walberla::float64 tmp_kernel_op_259 =
              tmp_kernel_op_122 * tmp_kernel_op_210;
          const walberla::float64 tmp_kernel_op_260 =
              tmp_kernel_op_124 * tmp_kernel_op_212;
          const walberla::float64 tmp_kernel_op_261 =
              tmp_kernel_op_126 * tmp_kernel_op_214;
          const walberla::float64 tmp_kernel_op_262 =
              tmp_kernel_op_128 * tmp_kernel_op_216;
          const walberla::float64 tmp_kernel_op_263 =
              tmp_kernel_op_130 * tmp_kernel_op_218;
          const walberla::float64 tmp_kernel_op_264 =
              tmp_kernel_op_132 * tmp_kernel_op_220;
          const walberla::float64 tmp_kernel_op_265 =
              tmp_kernel_op_0 * tmp_kernel_op_255 +
              tmp_kernel_op_16 * tmp_kernel_op_257 +
              tmp_kernel_op_24 * tmp_kernel_op_258 +
              tmp_kernel_op_256 * tmp_kernel_op_8 +
              tmp_kernel_op_259 * tmp_kernel_op_32 +
              tmp_kernel_op_260 * tmp_kernel_op_40 +
              tmp_kernel_op_261 * tmp_kernel_op_48 +
              tmp_kernel_op_262 * tmp_kernel_op_56 +
              tmp_kernel_op_263 * tmp_kernel_op_64 +
              tmp_kernel_op_264 * tmp_kernel_op_72;
          const walberla::float64 tmp_kernel_op_266 =
              tmp_kernel_op_114 * tmp_kernel_op_156;
          const walberla::float64 tmp_kernel_op_267 =
              tmp_kernel_op_116 * tmp_kernel_op_157;
          const walberla::float64 tmp_kernel_op_268 =
              tmp_kernel_op_118 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_269 =
              tmp_kernel_op_120 * tmp_kernel_op_159;
          const walberla::float64 tmp_kernel_op_270 =
              tmp_kernel_op_122 * tmp_kernel_op_160;
          const walberla::float64 tmp_kernel_op_271 =
              tmp_kernel_op_124 * tmp_kernel_op_161;
          const walberla::float64 tmp_kernel_op_272 =
              tmp_kernel_op_126 * tmp_kernel_op_162;
          const walberla::float64 tmp_kernel_op_273 =
              tmp_kernel_op_128 * tmp_kernel_op_163;
          const walberla::float64 tmp_kernel_op_274 =
              tmp_kernel_op_130 * tmp_kernel_op_164;
          const walberla::float64 tmp_kernel_op_275 =
              tmp_kernel_op_132 * tmp_kernel_op_165;
          const walberla::float64 tmp_kernel_op_276 =
              tmp_kernel_op_113 * tmp_kernel_op_266 +
              tmp_kernel_op_115 * tmp_kernel_op_267 +
              tmp_kernel_op_117 * tmp_kernel_op_268 +
              tmp_kernel_op_119 * tmp_kernel_op_269 +
              tmp_kernel_op_121 * tmp_kernel_op_270 +
              tmp_kernel_op_123 * tmp_kernel_op_271 +
              tmp_kernel_op_125 * tmp_kernel_op_272 +
              tmp_kernel_op_127 * tmp_kernel_op_273 +
              tmp_kernel_op_129 * tmp_kernel_op_274 +
              tmp_kernel_op_131 * tmp_kernel_op_275;
          const walberla::float64 tmp_kernel_op_277 =
              tmp_kernel_op_135 * tmp_kernel_op_245 * 4.0 +
              tmp_kernel_op_137 * tmp_kernel_op_246 * 4.0 +
              tmp_kernel_op_139 * tmp_kernel_op_247 * 4.0 +
              tmp_kernel_op_141 * tmp_kernel_op_248 * 4.0 +
              tmp_kernel_op_143 * tmp_kernel_op_249 * 4.0 +
              tmp_kernel_op_145 * tmp_kernel_op_250 * 4.0 +
              tmp_kernel_op_147 * tmp_kernel_op_251 * 4.0 +
              tmp_kernel_op_149 * tmp_kernel_op_252 * 4.0 +
              tmp_kernel_op_151 * tmp_kernel_op_253 * 4.0 +
              tmp_kernel_op_153 * tmp_kernel_op_254 * 4.0;
          const walberla::float64 tmp_kernel_op_278 =
              (tmp_kernel_op_114 * tmp_kernel_op_114);
          const walberla::float64 tmp_kernel_op_279 =
              (tmp_kernel_op_116 * tmp_kernel_op_116);
          const walberla::float64 tmp_kernel_op_280 =
              (tmp_kernel_op_118 * tmp_kernel_op_118);
          const walberla::float64 tmp_kernel_op_281 =
              (tmp_kernel_op_120 * tmp_kernel_op_120);
          const walberla::float64 tmp_kernel_op_282 =
              (tmp_kernel_op_122 * tmp_kernel_op_122);
          const walberla::float64 tmp_kernel_op_283 =
              (tmp_kernel_op_124 * tmp_kernel_op_124);
          const walberla::float64 tmp_kernel_op_284 =
              (tmp_kernel_op_126 * tmp_kernel_op_126);
          const walberla::float64 tmp_kernel_op_285 =
              (tmp_kernel_op_128 * tmp_kernel_op_128);
          const walberla::float64 tmp_kernel_op_286 =
              (tmp_kernel_op_130 * tmp_kernel_op_130);
          const walberla::float64 tmp_kernel_op_287 =
              (tmp_kernel_op_132 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_288 =
              tmp_kernel_op_136 * tmp_kernel_op_278 +
              tmp_kernel_op_138 * tmp_kernel_op_279 +
              tmp_kernel_op_140 * tmp_kernel_op_280 +
              tmp_kernel_op_142 * tmp_kernel_op_281 +
              tmp_kernel_op_144 * tmp_kernel_op_282 +
              tmp_kernel_op_146 * tmp_kernel_op_283 +
              tmp_kernel_op_148 * tmp_kernel_op_284 +
              tmp_kernel_op_150 * tmp_kernel_op_285 +
              tmp_kernel_op_152 * tmp_kernel_op_286 +
              tmp_kernel_op_154 * tmp_kernel_op_287;
          const walberla::float64 tmp_kernel_op_289 =
              tmp_kernel_op_234 * tmp_kernel_op_255 +
              tmp_kernel_op_235 * tmp_kernel_op_256 +
              tmp_kernel_op_236 * tmp_kernel_op_257 +
              tmp_kernel_op_237 * tmp_kernel_op_258 +
              tmp_kernel_op_238 * tmp_kernel_op_259 +
              tmp_kernel_op_239 * tmp_kernel_op_260 +
              tmp_kernel_op_240 * tmp_kernel_op_261 +
              tmp_kernel_op_241 * tmp_kernel_op_262 +
              tmp_kernel_op_242 * tmp_kernel_op_263 +
              tmp_kernel_op_243 * tmp_kernel_op_264;
          const walberla::float64 tmp_kernel_op_290 =
              tmp_kernel_op_190 * tmp_kernel_op_266 +
              tmp_kernel_op_191 * tmp_kernel_op_267 +
              tmp_kernel_op_192 * tmp_kernel_op_268 +
              tmp_kernel_op_193 * tmp_kernel_op_269 +
              tmp_kernel_op_194 * tmp_kernel_op_270 +
              tmp_kernel_op_195 * tmp_kernel_op_271 +
              tmp_kernel_op_196 * tmp_kernel_op_272 +
              tmp_kernel_op_197 * tmp_kernel_op_273 +
              tmp_kernel_op_198 * tmp_kernel_op_274 +
              tmp_kernel_op_199 * tmp_kernel_op_275;
          const walberla::float64 elMatVec_0 =
              src_dof_0 *
                  (tmp_kernel_op_1 * (tmp_kernel_op_2 * tmp_kernel_op_2) *
                       (tmp_kernel_op_5 * tmp_kernel_op_5) +
                   (tmp_kernel_op_10 * tmp_kernel_op_10) *
                       (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_9 +
                   tmp_kernel_op_17 * (tmp_kernel_op_18 * tmp_kernel_op_18) *
                       (tmp_kernel_op_21 * tmp_kernel_op_21) +
                   tmp_kernel_op_25 * (tmp_kernel_op_26 * tmp_kernel_op_26) *
                       (tmp_kernel_op_29 * tmp_kernel_op_29) +
                   tmp_kernel_op_33 * (tmp_kernel_op_34 * tmp_kernel_op_34) *
                       (tmp_kernel_op_37 * tmp_kernel_op_37) +
                   tmp_kernel_op_41 * (tmp_kernel_op_42 * tmp_kernel_op_42) *
                       (tmp_kernel_op_45 * tmp_kernel_op_45) +
                   tmp_kernel_op_49 * (tmp_kernel_op_50 * tmp_kernel_op_50) *
                       (tmp_kernel_op_53 * tmp_kernel_op_53) +
                   tmp_kernel_op_57 * (tmp_kernel_op_58 * tmp_kernel_op_58) *
                       (tmp_kernel_op_61 * tmp_kernel_op_61) +
                   tmp_kernel_op_65 * (tmp_kernel_op_66 * tmp_kernel_op_66) *
                       (tmp_kernel_op_69 * tmp_kernel_op_69) +
                   tmp_kernel_op_73 * (tmp_kernel_op_74 * tmp_kernel_op_74) *
                       (tmp_kernel_op_77 * tmp_kernel_op_77)) +
              src_dof_1 * tmp_kernel_op_101 + src_dof_2 * tmp_kernel_op_112 +
              src_dof_3 * tmp_kernel_op_80 + src_dof_4 * tmp_kernel_op_133 +
              src_dof_5 * tmp_kernel_op_134 + src_dof_6 * tmp_kernel_op_155;
          const walberla::float64 elMatVec_1 =
              src_dof_0 * tmp_kernel_op_101 +
              src_dof_1 *
                  ((tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_157 +
                   tmp_kernel_op_156 * (tmp_kernel_op_4 * tmp_kernel_op_4) +
                   tmp_kernel_op_158 * (tmp_kernel_op_20 * tmp_kernel_op_20) +
                   tmp_kernel_op_159 * (tmp_kernel_op_28 * tmp_kernel_op_28) +
                   tmp_kernel_op_160 * (tmp_kernel_op_36 * tmp_kernel_op_36) +
                   tmp_kernel_op_161 * (tmp_kernel_op_44 * tmp_kernel_op_44) +
                   tmp_kernel_op_162 * (tmp_kernel_op_52 * tmp_kernel_op_52) +
                   tmp_kernel_op_163 * (tmp_kernel_op_60 * tmp_kernel_op_60) +
                   tmp_kernel_op_164 * (tmp_kernel_op_68 * tmp_kernel_op_68) +
                   tmp_kernel_op_165 * (tmp_kernel_op_76 * tmp_kernel_op_76)) +
              src_dof_2 * tmp_kernel_op_187 + src_dof_3 * tmp_kernel_op_176 +
              src_dof_4 * tmp_kernel_op_188 + src_dof_5 * tmp_kernel_op_189 +
              src_dof_6 * tmp_kernel_op_200;
          const walberla::float64 elMatVec_2 =
              src_dof_0 * tmp_kernel_op_112 + src_dof_1 * tmp_kernel_op_187 +
              src_dof_2 *
                  ((tmp_kernel_op_102 * tmp_kernel_op_102) * tmp_kernel_op_202 +
                   (tmp_kernel_op_103 * tmp_kernel_op_103) * tmp_kernel_op_204 +
                   (tmp_kernel_op_104 * tmp_kernel_op_104) * tmp_kernel_op_206 +
                   (tmp_kernel_op_105 * tmp_kernel_op_105) * tmp_kernel_op_208 +
                   (tmp_kernel_op_106 * tmp_kernel_op_106) * tmp_kernel_op_210 +
                   (tmp_kernel_op_107 * tmp_kernel_op_107) * tmp_kernel_op_212 +
                   (tmp_kernel_op_108 * tmp_kernel_op_108) * tmp_kernel_op_214 +
                   (tmp_kernel_op_109 * tmp_kernel_op_109) * tmp_kernel_op_216 +
                   (tmp_kernel_op_110 * tmp_kernel_op_110) * tmp_kernel_op_218 +
                   (tmp_kernel_op_111 * tmp_kernel_op_111) *
                       tmp_kernel_op_220) +
              src_dof_3 * tmp_kernel_op_231 + src_dof_4 * tmp_kernel_op_233 +
              src_dof_5 * tmp_kernel_op_232 + src_dof_6 * tmp_kernel_op_244;
          const walberla::float64 elMatVec_3 =
              src_dof_0 * tmp_kernel_op_80 + src_dof_1 * tmp_kernel_op_176 +
              src_dof_2 * tmp_kernel_op_231 +
              src_dof_3 *
                  (tmp_kernel_op_245 * 16.0 + tmp_kernel_op_246 * 16.0 +
                   tmp_kernel_op_247 * 16.0 + tmp_kernel_op_248 * 16.0 +
                   tmp_kernel_op_249 * 16.0 + tmp_kernel_op_250 * 16.0 +
                   tmp_kernel_op_251 * 16.0 + tmp_kernel_op_252 * 16.0 +
                   tmp_kernel_op_253 * 16.0 + tmp_kernel_op_254 * 16.0) +
              src_dof_4 * tmp_kernel_op_265 + src_dof_5 * tmp_kernel_op_276 +
              src_dof_6 * tmp_kernel_op_277;
          const walberla::float64 elMatVec_4 =
              src_dof_0 * tmp_kernel_op_133 + src_dof_1 * tmp_kernel_op_188 +
              src_dof_2 * tmp_kernel_op_233 + src_dof_3 * tmp_kernel_op_265 +
              src_dof_4 * (tmp_kernel_op_202 * tmp_kernel_op_278 +
                           tmp_kernel_op_204 * tmp_kernel_op_279 +
                           tmp_kernel_op_206 * tmp_kernel_op_280 +
                           tmp_kernel_op_208 * tmp_kernel_op_281 +
                           tmp_kernel_op_210 * tmp_kernel_op_282 +
                           tmp_kernel_op_212 * tmp_kernel_op_283 +
                           tmp_kernel_op_214 * tmp_kernel_op_284 +
                           tmp_kernel_op_216 * tmp_kernel_op_285 +
                           tmp_kernel_op_218 * tmp_kernel_op_286 +
                           tmp_kernel_op_220 * tmp_kernel_op_287) +
              src_dof_5 * tmp_kernel_op_288 + src_dof_6 * tmp_kernel_op_289;
          const walberla::float64 elMatVec_5 =
              src_dof_0 * tmp_kernel_op_134 + src_dof_1 * tmp_kernel_op_189 +
              src_dof_2 * tmp_kernel_op_232 + src_dof_3 * tmp_kernel_op_276 +
              src_dof_4 * tmp_kernel_op_288 +
              src_dof_5 * (tmp_kernel_op_156 * tmp_kernel_op_278 +
                           tmp_kernel_op_157 * tmp_kernel_op_279 +
                           tmp_kernel_op_158 * tmp_kernel_op_280 +
                           tmp_kernel_op_159 * tmp_kernel_op_281 +
                           tmp_kernel_op_160 * tmp_kernel_op_282 +
                           tmp_kernel_op_161 * tmp_kernel_op_283 +
                           tmp_kernel_op_162 * tmp_kernel_op_284 +
                           tmp_kernel_op_163 * tmp_kernel_op_285 +
                           tmp_kernel_op_164 * tmp_kernel_op_286 +
                           tmp_kernel_op_165 * tmp_kernel_op_287) +
              src_dof_6 * tmp_kernel_op_290;
          const walberla::float64 elMatVec_6 =
              src_dof_0 * tmp_kernel_op_155 + src_dof_1 * tmp_kernel_op_200 +
              src_dof_2 * tmp_kernel_op_244 + src_dof_3 * tmp_kernel_op_277 +
              src_dof_4 * tmp_kernel_op_289 + src_dof_5 * tmp_kernel_op_290 +
              src_dof_6 *
                  ((tmp_kernel_op_135 * tmp_kernel_op_135) * tmp_kernel_op_245 +
                   (tmp_kernel_op_137 * tmp_kernel_op_137) * tmp_kernel_op_246 +
                   (tmp_kernel_op_139 * tmp_kernel_op_139) * tmp_kernel_op_247 +
                   (tmp_kernel_op_141 * tmp_kernel_op_141) * tmp_kernel_op_248 +
                   (tmp_kernel_op_143 * tmp_kernel_op_143) * tmp_kernel_op_249 +
                   (tmp_kernel_op_145 * tmp_kernel_op_145) * tmp_kernel_op_250 +
                   (tmp_kernel_op_147 * tmp_kernel_op_147) * tmp_kernel_op_251 +
                   (tmp_kernel_op_149 * tmp_kernel_op_149) * tmp_kernel_op_252 +
                   (tmp_kernel_op_151 * tmp_kernel_op_151) * tmp_kernel_op_253 +
                   (tmp_kernel_op_153 * tmp_kernel_op_153) * tmp_kernel_op_254);
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
      const walberla::float64 abs_det_jac_affine_BLUE =
          abs(jac_affine_0_0_BLUE * jac_affine_1_1_BLUE -
              jac_affine_0_1_BLUE * jac_affine_1_0_BLUE);
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1;
             ctr_0 += 1) {
#if 0
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
#endif
          const walberla::float64 jac_blending_1_1_q_9 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_9 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_9 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_9 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_9 =
              jac_blending_0_0_q_9 * jac_blending_1_1_q_9 -
              jac_blending_0_1_q_9 * jac_blending_1_0_q_9;
          const walberla::float64 jac_blending_1_1_q_8 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_8 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_8 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_8 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_8 =
              jac_blending_0_0_q_8 * jac_blending_1_1_q_8 -
              jac_blending_0_1_q_8 * jac_blending_1_0_q_8;
          const walberla::float64 jac_blending_1_1_q_7 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_7 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_7 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_7 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_7 =
              jac_blending_0_0_q_7 * jac_blending_1_1_q_7 -
              jac_blending_0_1_q_7 * jac_blending_1_0_q_7;
          const walberla::float64 jac_blending_1_1_q_6 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_6 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_6 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_6 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_6 =
              jac_blending_0_0_q_6 * jac_blending_1_1_q_6 -
              jac_blending_0_1_q_6 * jac_blending_1_0_q_6;
          const walberla::float64 jac_blending_1_1_q_5 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_5 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_5 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_5 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_5 =
              jac_blending_0_0_q_5 * jac_blending_1_1_q_5 -
              jac_blending_0_1_q_5 * jac_blending_1_0_q_5;
          const walberla::float64 jac_blending_1_1_q_4 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_4 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_4 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_4 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_4 =
              jac_blending_0_0_q_4 * jac_blending_1_1_q_4 -
              jac_blending_0_1_q_4 * jac_blending_1_0_q_4;
          const walberla::float64 jac_blending_1_1_q_3 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_3 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_3 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_3 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_2 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_2 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_2 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_1 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_1 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_1 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_0 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_0 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_0 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
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
          const walberla::float64 tmp_kernel_op_0 = -0.035167988576825245;
          const walberla::float64 tmp_kernel_op_1 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_0 *
                                                    0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_2 = 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_3 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_4 = -1.0175839942884126;
          const walberla::float64 tmp_kernel_op_5 =
              -tmp_kernel_op_3 - tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_6 =
              tmp_kernel_op_2 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_7 =
              tmp_kernel_op_1 * tmp_kernel_op_6 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_8 = 0.25875522422404362;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_1 *
                                                    0.052397656566423402;
          const walberla::float64 tmp_kernel_op_10 = 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_11 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_12 = -0.87062238788797819;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_15 =
              tmp_kernel_op_14 * tmp_kernel_op_9 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_16 = 0.30730286739632839;
          const walberla::float64 tmp_kernel_op_17 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_2 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_18 = 0.40172877323475986;
          const walberla::float64 tmp_kernel_op_19 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_20 = -0.84634856630183575;
          const walberla::float64 tmp_kernel_op_21 =
              -tmp_kernel_op_19 - tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_18 * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_17 * tmp_kernel_op_22 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_24 = 1.1801004774776607;
          const walberla::float64 tmp_kernel_op_25 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_3 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_26 = 0.64933214716985033;
          const walberla::float64 tmp_kernel_op_27 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_28 = -0.40994976126116967;
          const walberla::float64 tmp_kernel_op_29 =
              -tmp_kernel_op_27 - tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_26 * tmp_kernel_op_29;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_25 * tmp_kernel_op_30 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_32 = 5.1547196302936094;
          const walberla::float64 tmp_kernel_op_33 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_4 *
                                                     0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_34 = -0.0087919971442063094;
          const walberla::float64 tmp_kernel_op_35 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_36 = 1.5773598151468047;
          const walberla::float64 tmp_kernel_op_37 =
              -tmp_kernel_op_35 - tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_34 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_39 =
              tmp_kernel_op_33 * tmp_kernel_op_38 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_40 = 3.1529168382674562;
          const walberla::float64 tmp_kernel_op_41 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_5 *
                                                     0.052397656566423402;
          const walberla::float64 tmp_kernel_op_42 = 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_43 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_44 = 0.57645841913372808;
          const walberla::float64 tmp_kernel_op_45 =
              -tmp_kernel_op_43 - tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_46 =
              tmp_kernel_op_42 * tmp_kernel_op_45;
          const walberla::float64 tmp_kernel_op_47 =
              tmp_kernel_op_41 * tmp_kernel_op_46 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_48 = 1.6069150929390392;
          const walberla::float64 tmp_kernel_op_49 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_6 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_50 = 0.076825716849082126;
          const walberla::float64 tmp_kernel_op_51 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_52 = -0.19654245353048039;
          const walberla::float64 tmp_kernel_op_53 =
              -tmp_kernel_op_51 - tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_50 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_49 * tmp_kernel_op_54 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_56 = 2.5973285886794009;
          const walberla::float64 tmp_kernel_op_57 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_7 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_58 = 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_59 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_60 = 0.29866429433970043;
          const walberla::float64 tmp_kernel_op_61 =
              -tmp_kernel_op_59 - tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_58 * tmp_kernel_op_61;
          const walberla::float64 tmp_kernel_op_63 =
              tmp_kernel_op_57 * tmp_kernel_op_62 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_64 = 0.3413326446812428;
          const walberla::float64 tmp_kernel_op_65 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_8 *
                                                     0.045496761795224737;
          const walberla::float64 tmp_kernel_op_66 = 0.085333161170310645;
          const walberla::float64 tmp_kernel_op_67 = 1.6586673553187572;
          const walberla::float64 tmp_kernel_op_68 = -0.8293336776593786;
          const walberla::float64 tmp_kernel_op_69 =
              -tmp_kernel_op_67 - tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_66 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_65 * tmp_kernel_op_70 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_72 = 1.4147621866798579;
          const walberla::float64 tmp_kernel_op_73 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_9 *
                                                     0.10566414783403316;
          const walberla::float64 tmp_kernel_op_74 = 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_75 = 0.58523781332014213;
          const walberla::float64 tmp_kernel_op_76 = -0.29261890666007107;
          const walberla::float64 tmp_kernel_op_77 =
              -tmp_kernel_op_75 - tmp_kernel_op_76;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_74 * tmp_kernel_op_77;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_73 * tmp_kernel_op_78 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_0 * tmp_kernel_op_7 +
              tmp_kernel_op_15 * tmp_kernel_op_8 +
              tmp_kernel_op_16 * tmp_kernel_op_23 +
              tmp_kernel_op_24 * tmp_kernel_op_31 +
              tmp_kernel_op_32 * tmp_kernel_op_39 +
              tmp_kernel_op_40 * tmp_kernel_op_47 +
              tmp_kernel_op_48 * tmp_kernel_op_55 +
              tmp_kernel_op_56 * tmp_kernel_op_63 +
              tmp_kernel_op_64 * tmp_kernel_op_71 +
              tmp_kernel_op_72 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_81 =
              tmp_kernel_op_1 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_6 * tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_83 =
              tmp_kernel_op_9 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_14 * tmp_kernel_op_83;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_17 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_22 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_87 =
              tmp_kernel_op_25 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_30 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 =
              tmp_kernel_op_33 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_38 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_91 =
              tmp_kernel_op_41 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_46 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_93 =
              tmp_kernel_op_49 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_54 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              tmp_kernel_op_57 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_62 * tmp_kernel_op_95;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_65 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_70 * tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_99 =
              tmp_kernel_op_73 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_78 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_100 * tmp_kernel_op_76 +
              tmp_kernel_op_12 * tmp_kernel_op_84 +
              tmp_kernel_op_20 * tmp_kernel_op_86 +
              tmp_kernel_op_28 * tmp_kernel_op_88 +
              tmp_kernel_op_36 * tmp_kernel_op_90 +
              tmp_kernel_op_4 * tmp_kernel_op_82 +
              tmp_kernel_op_44 * tmp_kernel_op_92 +
              tmp_kernel_op_52 * tmp_kernel_op_94 +
              tmp_kernel_op_60 * tmp_kernel_op_96 +
              tmp_kernel_op_68 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_102 = tmp_kernel_op_3 - 1.0;
          const walberla::float64 tmp_kernel_op_103 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_104 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_105 = tmp_kernel_op_27 - 1.0;
          const walberla::float64 tmp_kernel_op_106 = tmp_kernel_op_35 - 1.0;
          const walberla::float64 tmp_kernel_op_107 = tmp_kernel_op_43 - 1.0;
          const walberla::float64 tmp_kernel_op_108 = tmp_kernel_op_51 - 1.0;
          const walberla::float64 tmp_kernel_op_109 = tmp_kernel_op_59 - 1.0;
          const walberla::float64 tmp_kernel_op_110 = tmp_kernel_op_67 - 1.0;
          const walberla::float64 tmp_kernel_op_111 = tmp_kernel_op_75 - 1.0;
          const walberla::float64 tmp_kernel_op_112 =
              tmp_kernel_op_102 * tmp_kernel_op_7 +
              tmp_kernel_op_103 * tmp_kernel_op_15 +
              tmp_kernel_op_104 * tmp_kernel_op_23 +
              tmp_kernel_op_105 * tmp_kernel_op_31 +
              tmp_kernel_op_106 * tmp_kernel_op_39 +
              tmp_kernel_op_107 * tmp_kernel_op_47 +
              tmp_kernel_op_108 * tmp_kernel_op_55 +
              tmp_kernel_op_109 * tmp_kernel_op_63 +
              tmp_kernel_op_110 * tmp_kernel_op_71 +
              tmp_kernel_op_111 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_113 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_114 =
              -tmp_kernel_op_0 - tmp_kernel_op_113 + 4.0;
          const walberla::float64 tmp_kernel_op_115 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_116 =
              -tmp_kernel_op_115 - tmp_kernel_op_8 + 4.0;
          const walberla::float64 tmp_kernel_op_117 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_118 =
              -tmp_kernel_op_117 - tmp_kernel_op_16 + 4.0;
          const walberla::float64 tmp_kernel_op_119 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_120 =
              -tmp_kernel_op_119 - tmp_kernel_op_24 + 4.0;
          const walberla::float64 tmp_kernel_op_121 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_122 =
              -tmp_kernel_op_121 - tmp_kernel_op_32 + 4.0;
          const walberla::float64 tmp_kernel_op_123 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_124 =
              -tmp_kernel_op_123 - tmp_kernel_op_40 + 4.0;
          const walberla::float64 tmp_kernel_op_125 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_126 =
              -tmp_kernel_op_125 - tmp_kernel_op_48 + 4.0;
          const walberla::float64 tmp_kernel_op_127 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_128 =
              -tmp_kernel_op_127 - tmp_kernel_op_56 + 4.0;
          const walberla::float64 tmp_kernel_op_129 = 3.3173347106375144;
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_129 - tmp_kernel_op_64 + 4.0;
          const walberla::float64 tmp_kernel_op_131 = 1.1704756266402843;
          const walberla::float64 tmp_kernel_op_132 =
              -tmp_kernel_op_131 - tmp_kernel_op_72 + 4.0;
          const walberla::float64 tmp_kernel_op_133 =
              tmp_kernel_op_114 * tmp_kernel_op_7 +
              tmp_kernel_op_116 * tmp_kernel_op_15 +
              tmp_kernel_op_118 * tmp_kernel_op_23 +
              tmp_kernel_op_120 * tmp_kernel_op_31 +
              tmp_kernel_op_122 * tmp_kernel_op_39 +
              tmp_kernel_op_124 * tmp_kernel_op_47 +
              tmp_kernel_op_126 * tmp_kernel_op_55 +
              tmp_kernel_op_128 * tmp_kernel_op_63 +
              tmp_kernel_op_130 * tmp_kernel_op_71 +
              tmp_kernel_op_132 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_134 =
              tmp_kernel_op_100 * tmp_kernel_op_132 +
              tmp_kernel_op_114 * tmp_kernel_op_82 +
              tmp_kernel_op_116 * tmp_kernel_op_84 +
              tmp_kernel_op_118 * tmp_kernel_op_86 +
              tmp_kernel_op_120 * tmp_kernel_op_88 +
              tmp_kernel_op_122 * tmp_kernel_op_90 +
              tmp_kernel_op_124 * tmp_kernel_op_92 +
              tmp_kernel_op_126 * tmp_kernel_op_94 +
              tmp_kernel_op_128 * tmp_kernel_op_96 +
              tmp_kernel_op_130 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_135 = 34.794357504481866;
          const walberla::float64 tmp_kernel_op_136 =
              tmp_kernel_op_81 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_137 = 21.28218865830533;
          const walberla::float64 tmp_kernel_op_138 =
              tmp_kernel_op_83 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_139 = 10.846676877338515;
          const walberla::float64 tmp_kernel_op_140 =
              tmp_kernel_op_85 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_141 = 17.531967973585957;
          const walberla::float64 tmp_kernel_op_142 =
              tmp_kernel_op_87 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_143 = -0.23738392289357257;
          const walberla::float64 tmp_kernel_op_144 =
              tmp_kernel_op_89 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_145 = 1.7465977635122938;
          const walberla::float64 tmp_kernel_op_146 =
              tmp_kernel_op_91 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_147 = 2.0742943549252164;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_93 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_149 = 7.9656782229742085;
          const walberla::float64 tmp_kernel_op_150 =
              tmp_kernel_op_95 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_151 = 2.3039953515983882;
          const walberla::float64 tmp_kernel_op_152 =
              tmp_kernel_op_97 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_153 = 9.5496447600890395;
          const walberla::float64 tmp_kernel_op_154 =
              tmp_kernel_op_99 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_155 =
              tmp_kernel_op_135 * tmp_kernel_op_136 * tmp_kernel_op_6 +
              tmp_kernel_op_137 * tmp_kernel_op_138 * tmp_kernel_op_14 +
              tmp_kernel_op_139 * tmp_kernel_op_140 * tmp_kernel_op_22 +
              tmp_kernel_op_141 * tmp_kernel_op_142 * tmp_kernel_op_30 +
              tmp_kernel_op_143 * tmp_kernel_op_144 * tmp_kernel_op_38 +
              tmp_kernel_op_145 * tmp_kernel_op_146 * tmp_kernel_op_46 +
              tmp_kernel_op_147 * tmp_kernel_op_148 * tmp_kernel_op_54 +
              tmp_kernel_op_149 * tmp_kernel_op_150 * tmp_kernel_op_62 +
              tmp_kernel_op_151 * tmp_kernel_op_152 * tmp_kernel_op_70 +
              tmp_kernel_op_153 * tmp_kernel_op_154 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_1 * 7.7299213783731935e-5;
          const walberla::float64 tmp_kernel_op_157 =
              tmp_kernel_op_9 * 0.0041846416289521935;
          const walberla::float64 tmp_kernel_op_158 =
              tmp_kernel_op_17 * 0.0059021907693753367;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_25 * 0.087039821058937664;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_33 * 1.6606959041833929;
          const walberla::float64 tmp_kernel_op_161 =
              tmp_kernel_op_41 * 0.62130528681440322;
          const walberla::float64 tmp_kernel_op_162 =
              tmp_kernel_op_49 * 0.16138600724470506;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_57 * 0.42163223734820804;
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_65 * 0.0072817483953182219;
          const walberla::float64 tmp_kernel_op_165 =
              tmp_kernel_op_73 * 0.12509700280369832;
          const walberla::float64 tmp_kernel_op_166 =
              tmp_kernel_op_156 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_12 * tmp_kernel_op_157;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_158 * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_169 =
              tmp_kernel_op_159 * tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_170 =
              tmp_kernel_op_160 * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_161 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_162 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_173 =
              tmp_kernel_op_163 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_174 =
              tmp_kernel_op_164 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_165 * tmp_kernel_op_76;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_113 * tmp_kernel_op_166 +
              tmp_kernel_op_115 * tmp_kernel_op_167 +
              tmp_kernel_op_117 * tmp_kernel_op_168 +
              tmp_kernel_op_119 * tmp_kernel_op_169 +
              tmp_kernel_op_121 * tmp_kernel_op_170 +
              tmp_kernel_op_123 * tmp_kernel_op_171 +
              tmp_kernel_op_125 * tmp_kernel_op_172 +
              tmp_kernel_op_127 * tmp_kernel_op_173 +
              tmp_kernel_op_129 * tmp_kernel_op_174 +
              tmp_kernel_op_131 * tmp_kernel_op_175;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_136 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_178 =
              tmp_kernel_op_12 * tmp_kernel_op_138;
          const walberla::float64 tmp_kernel_op_179 =
              tmp_kernel_op_140 * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_180 =
              tmp_kernel_op_142 * tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_181 =
              tmp_kernel_op_144 * tmp_kernel_op_36;
          const walberla::float64 tmp_kernel_op_182 =
              tmp_kernel_op_146 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_183 =
              tmp_kernel_op_148 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_184 =
              tmp_kernel_op_150 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_185 =
              tmp_kernel_op_152 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_154 * tmp_kernel_op_76;
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_102 * tmp_kernel_op_177 +
              tmp_kernel_op_103 * tmp_kernel_op_178 +
              tmp_kernel_op_104 * tmp_kernel_op_179 +
              tmp_kernel_op_105 * tmp_kernel_op_180 +
              tmp_kernel_op_106 * tmp_kernel_op_181 +
              tmp_kernel_op_107 * tmp_kernel_op_182 +
              tmp_kernel_op_108 * tmp_kernel_op_183 +
              tmp_kernel_op_109 * tmp_kernel_op_184 +
              tmp_kernel_op_110 * tmp_kernel_op_185 +
              tmp_kernel_op_111 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_114 * tmp_kernel_op_177 +
              tmp_kernel_op_116 * tmp_kernel_op_178 +
              tmp_kernel_op_118 * tmp_kernel_op_179 +
              tmp_kernel_op_120 * tmp_kernel_op_180 +
              tmp_kernel_op_122 * tmp_kernel_op_181 +
              tmp_kernel_op_124 * tmp_kernel_op_182 +
              tmp_kernel_op_126 * tmp_kernel_op_183 +
              tmp_kernel_op_128 * tmp_kernel_op_184 +
              tmp_kernel_op_130 * tmp_kernel_op_185 +
              tmp_kernel_op_132 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_114 * tmp_kernel_op_166 +
              tmp_kernel_op_116 * tmp_kernel_op_167 +
              tmp_kernel_op_118 * tmp_kernel_op_168 +
              tmp_kernel_op_120 * tmp_kernel_op_169 +
              tmp_kernel_op_122 * tmp_kernel_op_170 +
              tmp_kernel_op_124 * tmp_kernel_op_171 +
              tmp_kernel_op_126 * tmp_kernel_op_172 +
              tmp_kernel_op_128 * tmp_kernel_op_173 +
              tmp_kernel_op_130 * tmp_kernel_op_174 +
              tmp_kernel_op_132 * tmp_kernel_op_175;
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_135 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_191 =
              tmp_kernel_op_137 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_192 =
              tmp_kernel_op_139 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_193 =
              tmp_kernel_op_141 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_194 =
              tmp_kernel_op_143 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_195 =
              tmp_kernel_op_145 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_196 =
              tmp_kernel_op_147 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_197 =
              tmp_kernel_op_149 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_198 =
              tmp_kernel_op_151 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_153 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_166 * tmp_kernel_op_190 +
              tmp_kernel_op_167 * tmp_kernel_op_191 +
              tmp_kernel_op_168 * tmp_kernel_op_192 +
              tmp_kernel_op_169 * tmp_kernel_op_193 +
              tmp_kernel_op_170 * tmp_kernel_op_194 +
              tmp_kernel_op_171 * tmp_kernel_op_195 +
              tmp_kernel_op_172 * tmp_kernel_op_196 +
              tmp_kernel_op_173 * tmp_kernel_op_197 +
              tmp_kernel_op_174 * tmp_kernel_op_198 +
              tmp_kernel_op_175 * tmp_kernel_op_199;
          const walberla::float64 tmp_kernel_op_201 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_1 * tmp_kernel_op_201;
          const walberla::float64 tmp_kernel_op_203 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_203 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_205 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_17 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_207 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_207 * tmp_kernel_op_25;
          const walberla::float64 tmp_kernel_op_209 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_210 =
              tmp_kernel_op_209 * tmp_kernel_op_33;
          const walberla::float64 tmp_kernel_op_211 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_212 =
              tmp_kernel_op_211 * tmp_kernel_op_41;
          const walberla::float64 tmp_kernel_op_213 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_214 =
              tmp_kernel_op_213 * tmp_kernel_op_49;
          const walberla::float64 tmp_kernel_op_215 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_216 =
              tmp_kernel_op_215 * tmp_kernel_op_57;
          const walberla::float64 tmp_kernel_op_217 = 0.68779434890003011;
          const walberla::float64 tmp_kernel_op_218 =
              tmp_kernel_op_217 * tmp_kernel_op_65;
          const walberla::float64 tmp_kernel_op_219 = 0.085625824534935377;
          const walberla::float64 tmp_kernel_op_220 =
              tmp_kernel_op_219 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_221 =
              tmp_kernel_op_102 * tmp_kernel_op_202;
          const walberla::float64 tmp_kernel_op_222 =
              tmp_kernel_op_103 * tmp_kernel_op_204;
          const walberla::float64 tmp_kernel_op_223 =
              tmp_kernel_op_104 * tmp_kernel_op_206;
          const walberla::float64 tmp_kernel_op_224 =
              tmp_kernel_op_105 * tmp_kernel_op_208;
          const walberla::float64 tmp_kernel_op_225 =
              tmp_kernel_op_106 * tmp_kernel_op_210;
          const walberla::float64 tmp_kernel_op_226 =
              tmp_kernel_op_107 * tmp_kernel_op_212;
          const walberla::float64 tmp_kernel_op_227 =
              tmp_kernel_op_108 * tmp_kernel_op_214;
          const walberla::float64 tmp_kernel_op_228 =
              tmp_kernel_op_109 * tmp_kernel_op_216;
          const walberla::float64 tmp_kernel_op_229 =
              tmp_kernel_op_110 * tmp_kernel_op_218;
          const walberla::float64 tmp_kernel_op_230 =
              tmp_kernel_op_111 * tmp_kernel_op_220;
          const walberla::float64 tmp_kernel_op_231 =
              tmp_kernel_op_0 * tmp_kernel_op_221 +
              tmp_kernel_op_16 * tmp_kernel_op_223 +
              tmp_kernel_op_222 * tmp_kernel_op_8 +
              tmp_kernel_op_224 * tmp_kernel_op_24 +
              tmp_kernel_op_225 * tmp_kernel_op_32 +
              tmp_kernel_op_226 * tmp_kernel_op_40 +
              tmp_kernel_op_227 * tmp_kernel_op_48 +
              tmp_kernel_op_228 * tmp_kernel_op_56 +
              tmp_kernel_op_229 * tmp_kernel_op_64 +
              tmp_kernel_op_230 * tmp_kernel_op_72;
          const walberla::float64 tmp_kernel_op_232 =
              tmp_kernel_op_102 * tmp_kernel_op_114 * tmp_kernel_op_136 +
              tmp_kernel_op_103 * tmp_kernel_op_116 * tmp_kernel_op_138 +
              tmp_kernel_op_104 * tmp_kernel_op_118 * tmp_kernel_op_140 +
              tmp_kernel_op_105 * tmp_kernel_op_120 * tmp_kernel_op_142 +
              tmp_kernel_op_106 * tmp_kernel_op_122 * tmp_kernel_op_144 +
              tmp_kernel_op_107 * tmp_kernel_op_124 * tmp_kernel_op_146 +
              tmp_kernel_op_108 * tmp_kernel_op_126 * tmp_kernel_op_148 +
              tmp_kernel_op_109 * tmp_kernel_op_128 * tmp_kernel_op_150 +
              tmp_kernel_op_110 * tmp_kernel_op_130 * tmp_kernel_op_152 +
              tmp_kernel_op_111 * tmp_kernel_op_132 * tmp_kernel_op_154;
          const walberla::float64 tmp_kernel_op_233 =
              tmp_kernel_op_114 * tmp_kernel_op_221 +
              tmp_kernel_op_116 * tmp_kernel_op_222 +
              tmp_kernel_op_118 * tmp_kernel_op_223 +
              tmp_kernel_op_120 * tmp_kernel_op_224 +
              tmp_kernel_op_122 * tmp_kernel_op_225 +
              tmp_kernel_op_124 * tmp_kernel_op_226 +
              tmp_kernel_op_126 * tmp_kernel_op_227 +
              tmp_kernel_op_128 * tmp_kernel_op_228 +
              tmp_kernel_op_130 * tmp_kernel_op_229 +
              tmp_kernel_op_132 * tmp_kernel_op_230;
          const walberla::float64 tmp_kernel_op_234 =
              tmp_kernel_op_135 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_235 =
              tmp_kernel_op_137 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_236 =
              tmp_kernel_op_139 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_237 =
              tmp_kernel_op_141 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_238 =
              tmp_kernel_op_143 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_239 =
              tmp_kernel_op_145 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_240 =
              tmp_kernel_op_147 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_241 =
              tmp_kernel_op_149 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_242 =
              tmp_kernel_op_151 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_243 =
              tmp_kernel_op_153 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_244 =
              tmp_kernel_op_221 * tmp_kernel_op_234 +
              tmp_kernel_op_222 * tmp_kernel_op_235 +
              tmp_kernel_op_223 * tmp_kernel_op_236 +
              tmp_kernel_op_224 * tmp_kernel_op_237 +
              tmp_kernel_op_225 * tmp_kernel_op_238 +
              tmp_kernel_op_226 * tmp_kernel_op_239 +
              tmp_kernel_op_227 * tmp_kernel_op_240 +
              tmp_kernel_op_228 * tmp_kernel_op_241 +
              tmp_kernel_op_229 * tmp_kernel_op_242 +
              tmp_kernel_op_230 * tmp_kernel_op_243;
          const walberla::float64 tmp_kernel_op_245 =
              tmp_kernel_op_156 * tmp_kernel_op_201;
          const walberla::float64 tmp_kernel_op_246 =
              tmp_kernel_op_157 * tmp_kernel_op_203;
          const walberla::float64 tmp_kernel_op_247 =
              tmp_kernel_op_158 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_248 =
              tmp_kernel_op_159 * tmp_kernel_op_207;
          const walberla::float64 tmp_kernel_op_249 =
              tmp_kernel_op_160 * tmp_kernel_op_209;
          const walberla::float64 tmp_kernel_op_250 =
              tmp_kernel_op_161 * tmp_kernel_op_211;
          const walberla::float64 tmp_kernel_op_251 =
              tmp_kernel_op_162 * tmp_kernel_op_213;
          const walberla::float64 tmp_kernel_op_252 =
              tmp_kernel_op_163 * tmp_kernel_op_215;
          const walberla::float64 tmp_kernel_op_253 =
              tmp_kernel_op_164 * tmp_kernel_op_217;
          const walberla::float64 tmp_kernel_op_254 =
              tmp_kernel_op_165 * tmp_kernel_op_219;
          const walberla::float64 tmp_kernel_op_255 =
              tmp_kernel_op_114 * tmp_kernel_op_202;
          const walberla::float64 tmp_kernel_op_256 =
              tmp_kernel_op_116 * tmp_kernel_op_204;
          const walberla::float64 tmp_kernel_op_257 =
              tmp_kernel_op_118 * tmp_kernel_op_206;
          const walberla::float64 tmp_kernel_op_258 =
              tmp_kernel_op_120 * tmp_kernel_op_208;
          const walberla::float64 tmp_kernel_op_259 =
              tmp_kernel_op_122 * tmp_kernel_op_210;
          const walberla::float64 tmp_kernel_op_260 =
              tmp_kernel_op_124 * tmp_kernel_op_212;
          const walberla::float64 tmp_kernel_op_261 =
              tmp_kernel_op_126 * tmp_kernel_op_214;
          const walberla::float64 tmp_kernel_op_262 =
              tmp_kernel_op_128 * tmp_kernel_op_216;
          const walberla::float64 tmp_kernel_op_263 =
              tmp_kernel_op_130 * tmp_kernel_op_218;
          const walberla::float64 tmp_kernel_op_264 =
              tmp_kernel_op_132 * tmp_kernel_op_220;
          const walberla::float64 tmp_kernel_op_265 =
              tmp_kernel_op_0 * tmp_kernel_op_255 +
              tmp_kernel_op_16 * tmp_kernel_op_257 +
              tmp_kernel_op_24 * tmp_kernel_op_258 +
              tmp_kernel_op_256 * tmp_kernel_op_8 +
              tmp_kernel_op_259 * tmp_kernel_op_32 +
              tmp_kernel_op_260 * tmp_kernel_op_40 +
              tmp_kernel_op_261 * tmp_kernel_op_48 +
              tmp_kernel_op_262 * tmp_kernel_op_56 +
              tmp_kernel_op_263 * tmp_kernel_op_64 +
              tmp_kernel_op_264 * tmp_kernel_op_72;
          const walberla::float64 tmp_kernel_op_266 =
              tmp_kernel_op_114 * tmp_kernel_op_156;
          const walberla::float64 tmp_kernel_op_267 =
              tmp_kernel_op_116 * tmp_kernel_op_157;
          const walberla::float64 tmp_kernel_op_268 =
              tmp_kernel_op_118 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_269 =
              tmp_kernel_op_120 * tmp_kernel_op_159;
          const walberla::float64 tmp_kernel_op_270 =
              tmp_kernel_op_122 * tmp_kernel_op_160;
          const walberla::float64 tmp_kernel_op_271 =
              tmp_kernel_op_124 * tmp_kernel_op_161;
          const walberla::float64 tmp_kernel_op_272 =
              tmp_kernel_op_126 * tmp_kernel_op_162;
          const walberla::float64 tmp_kernel_op_273 =
              tmp_kernel_op_128 * tmp_kernel_op_163;
          const walberla::float64 tmp_kernel_op_274 =
              tmp_kernel_op_130 * tmp_kernel_op_164;
          const walberla::float64 tmp_kernel_op_275 =
              tmp_kernel_op_132 * tmp_kernel_op_165;
          const walberla::float64 tmp_kernel_op_276 =
              tmp_kernel_op_113 * tmp_kernel_op_266 +
              tmp_kernel_op_115 * tmp_kernel_op_267 +
              tmp_kernel_op_117 * tmp_kernel_op_268 +
              tmp_kernel_op_119 * tmp_kernel_op_269 +
              tmp_kernel_op_121 * tmp_kernel_op_270 +
              tmp_kernel_op_123 * tmp_kernel_op_271 +
              tmp_kernel_op_125 * tmp_kernel_op_272 +
              tmp_kernel_op_127 * tmp_kernel_op_273 +
              tmp_kernel_op_129 * tmp_kernel_op_274 +
              tmp_kernel_op_131 * tmp_kernel_op_275;
          const walberla::float64 tmp_kernel_op_277 =
              tmp_kernel_op_135 * tmp_kernel_op_245 * 4.0 +
              tmp_kernel_op_137 * tmp_kernel_op_246 * 4.0 +
              tmp_kernel_op_139 * tmp_kernel_op_247 * 4.0 +
              tmp_kernel_op_141 * tmp_kernel_op_248 * 4.0 +
              tmp_kernel_op_143 * tmp_kernel_op_249 * 4.0 +
              tmp_kernel_op_145 * tmp_kernel_op_250 * 4.0 +
              tmp_kernel_op_147 * tmp_kernel_op_251 * 4.0 +
              tmp_kernel_op_149 * tmp_kernel_op_252 * 4.0 +
              tmp_kernel_op_151 * tmp_kernel_op_253 * 4.0 +
              tmp_kernel_op_153 * tmp_kernel_op_254 * 4.0;
          const walberla::float64 tmp_kernel_op_278 =
              (tmp_kernel_op_114 * tmp_kernel_op_114);
          const walberla::float64 tmp_kernel_op_279 =
              (tmp_kernel_op_116 * tmp_kernel_op_116);
          const walberla::float64 tmp_kernel_op_280 =
              (tmp_kernel_op_118 * tmp_kernel_op_118);
          const walberla::float64 tmp_kernel_op_281 =
              (tmp_kernel_op_120 * tmp_kernel_op_120);
          const walberla::float64 tmp_kernel_op_282 =
              (tmp_kernel_op_122 * tmp_kernel_op_122);
          const walberla::float64 tmp_kernel_op_283 =
              (tmp_kernel_op_124 * tmp_kernel_op_124);
          const walberla::float64 tmp_kernel_op_284 =
              (tmp_kernel_op_126 * tmp_kernel_op_126);
          const walberla::float64 tmp_kernel_op_285 =
              (tmp_kernel_op_128 * tmp_kernel_op_128);
          const walberla::float64 tmp_kernel_op_286 =
              (tmp_kernel_op_130 * tmp_kernel_op_130);
          const walberla::float64 tmp_kernel_op_287 =
              (tmp_kernel_op_132 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_288 =
              tmp_kernel_op_136 * tmp_kernel_op_278 +
              tmp_kernel_op_138 * tmp_kernel_op_279 +
              tmp_kernel_op_140 * tmp_kernel_op_280 +
              tmp_kernel_op_142 * tmp_kernel_op_281 +
              tmp_kernel_op_144 * tmp_kernel_op_282 +
              tmp_kernel_op_146 * tmp_kernel_op_283 +
              tmp_kernel_op_148 * tmp_kernel_op_284 +
              tmp_kernel_op_150 * tmp_kernel_op_285 +
              tmp_kernel_op_152 * tmp_kernel_op_286 +
              tmp_kernel_op_154 * tmp_kernel_op_287;
          const walberla::float64 tmp_kernel_op_289 =
              tmp_kernel_op_234 * tmp_kernel_op_255 +
              tmp_kernel_op_235 * tmp_kernel_op_256 +
              tmp_kernel_op_236 * tmp_kernel_op_257 +
              tmp_kernel_op_237 * tmp_kernel_op_258 +
              tmp_kernel_op_238 * tmp_kernel_op_259 +
              tmp_kernel_op_239 * tmp_kernel_op_260 +
              tmp_kernel_op_240 * tmp_kernel_op_261 +
              tmp_kernel_op_241 * tmp_kernel_op_262 +
              tmp_kernel_op_242 * tmp_kernel_op_263 +
              tmp_kernel_op_243 * tmp_kernel_op_264;
          const walberla::float64 tmp_kernel_op_290 =
              tmp_kernel_op_190 * tmp_kernel_op_266 +
              tmp_kernel_op_191 * tmp_kernel_op_267 +
              tmp_kernel_op_192 * tmp_kernel_op_268 +
              tmp_kernel_op_193 * tmp_kernel_op_269 +
              tmp_kernel_op_194 * tmp_kernel_op_270 +
              tmp_kernel_op_195 * tmp_kernel_op_271 +
              tmp_kernel_op_196 * tmp_kernel_op_272 +
              tmp_kernel_op_197 * tmp_kernel_op_273 +
              tmp_kernel_op_198 * tmp_kernel_op_274 +
              tmp_kernel_op_199 * tmp_kernel_op_275;
          const walberla::float64 elMatVec_0 =
              src_dof_0 *
                  (tmp_kernel_op_1 * (tmp_kernel_op_2 * tmp_kernel_op_2) *
                       (tmp_kernel_op_5 * tmp_kernel_op_5) +
                   (tmp_kernel_op_10 * tmp_kernel_op_10) *
                       (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_9 +
                   tmp_kernel_op_17 * (tmp_kernel_op_18 * tmp_kernel_op_18) *
                       (tmp_kernel_op_21 * tmp_kernel_op_21) +
                   tmp_kernel_op_25 * (tmp_kernel_op_26 * tmp_kernel_op_26) *
                       (tmp_kernel_op_29 * tmp_kernel_op_29) +
                   tmp_kernel_op_33 * (tmp_kernel_op_34 * tmp_kernel_op_34) *
                       (tmp_kernel_op_37 * tmp_kernel_op_37) +
                   tmp_kernel_op_41 * (tmp_kernel_op_42 * tmp_kernel_op_42) *
                       (tmp_kernel_op_45 * tmp_kernel_op_45) +
                   tmp_kernel_op_49 * (tmp_kernel_op_50 * tmp_kernel_op_50) *
                       (tmp_kernel_op_53 * tmp_kernel_op_53) +
                   tmp_kernel_op_57 * (tmp_kernel_op_58 * tmp_kernel_op_58) *
                       (tmp_kernel_op_61 * tmp_kernel_op_61) +
                   tmp_kernel_op_65 * (tmp_kernel_op_66 * tmp_kernel_op_66) *
                       (tmp_kernel_op_69 * tmp_kernel_op_69) +
                   tmp_kernel_op_73 * (tmp_kernel_op_74 * tmp_kernel_op_74) *
                       (tmp_kernel_op_77 * tmp_kernel_op_77)) +
              src_dof_1 * tmp_kernel_op_101 + src_dof_2 * tmp_kernel_op_112 +
              src_dof_3 * tmp_kernel_op_80 + src_dof_4 * tmp_kernel_op_133 +
              src_dof_5 * tmp_kernel_op_134 + src_dof_6 * tmp_kernel_op_155;
          const walberla::float64 elMatVec_1 =
              src_dof_0 * tmp_kernel_op_101 +
              src_dof_1 *
                  ((tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_157 +
                   tmp_kernel_op_156 * (tmp_kernel_op_4 * tmp_kernel_op_4) +
                   tmp_kernel_op_158 * (tmp_kernel_op_20 * tmp_kernel_op_20) +
                   tmp_kernel_op_159 * (tmp_kernel_op_28 * tmp_kernel_op_28) +
                   tmp_kernel_op_160 * (tmp_kernel_op_36 * tmp_kernel_op_36) +
                   tmp_kernel_op_161 * (tmp_kernel_op_44 * tmp_kernel_op_44) +
                   tmp_kernel_op_162 * (tmp_kernel_op_52 * tmp_kernel_op_52) +
                   tmp_kernel_op_163 * (tmp_kernel_op_60 * tmp_kernel_op_60) +
                   tmp_kernel_op_164 * (tmp_kernel_op_68 * tmp_kernel_op_68) +
                   tmp_kernel_op_165 * (tmp_kernel_op_76 * tmp_kernel_op_76)) +
              src_dof_2 * tmp_kernel_op_187 + src_dof_3 * tmp_kernel_op_176 +
              src_dof_4 * tmp_kernel_op_188 + src_dof_5 * tmp_kernel_op_189 +
              src_dof_6 * tmp_kernel_op_200;
          const walberla::float64 elMatVec_2 =
              src_dof_0 * tmp_kernel_op_112 + src_dof_1 * tmp_kernel_op_187 +
              src_dof_2 *
                  ((tmp_kernel_op_102 * tmp_kernel_op_102) * tmp_kernel_op_202 +
                   (tmp_kernel_op_103 * tmp_kernel_op_103) * tmp_kernel_op_204 +
                   (tmp_kernel_op_104 * tmp_kernel_op_104) * tmp_kernel_op_206 +
                   (tmp_kernel_op_105 * tmp_kernel_op_105) * tmp_kernel_op_208 +
                   (tmp_kernel_op_106 * tmp_kernel_op_106) * tmp_kernel_op_210 +
                   (tmp_kernel_op_107 * tmp_kernel_op_107) * tmp_kernel_op_212 +
                   (tmp_kernel_op_108 * tmp_kernel_op_108) * tmp_kernel_op_214 +
                   (tmp_kernel_op_109 * tmp_kernel_op_109) * tmp_kernel_op_216 +
                   (tmp_kernel_op_110 * tmp_kernel_op_110) * tmp_kernel_op_218 +
                   (tmp_kernel_op_111 * tmp_kernel_op_111) *
                       tmp_kernel_op_220) +
              src_dof_3 * tmp_kernel_op_231 + src_dof_4 * tmp_kernel_op_233 +
              src_dof_5 * tmp_kernel_op_232 + src_dof_6 * tmp_kernel_op_244;
          const walberla::float64 elMatVec_3 =
              src_dof_0 * tmp_kernel_op_80 + src_dof_1 * tmp_kernel_op_176 +
              src_dof_2 * tmp_kernel_op_231 +
              src_dof_3 *
                  (tmp_kernel_op_245 * 16.0 + tmp_kernel_op_246 * 16.0 +
                   tmp_kernel_op_247 * 16.0 + tmp_kernel_op_248 * 16.0 +
                   tmp_kernel_op_249 * 16.0 + tmp_kernel_op_250 * 16.0 +
                   tmp_kernel_op_251 * 16.0 + tmp_kernel_op_252 * 16.0 +
                   tmp_kernel_op_253 * 16.0 + tmp_kernel_op_254 * 16.0) +
              src_dof_4 * tmp_kernel_op_265 + src_dof_5 * tmp_kernel_op_276 +
              src_dof_6 * tmp_kernel_op_277;
          const walberla::float64 elMatVec_4 =
              src_dof_0 * tmp_kernel_op_133 + src_dof_1 * tmp_kernel_op_188 +
              src_dof_2 * tmp_kernel_op_233 + src_dof_3 * tmp_kernel_op_265 +
              src_dof_4 * (tmp_kernel_op_202 * tmp_kernel_op_278 +
                           tmp_kernel_op_204 * tmp_kernel_op_279 +
                           tmp_kernel_op_206 * tmp_kernel_op_280 +
                           tmp_kernel_op_208 * tmp_kernel_op_281 +
                           tmp_kernel_op_210 * tmp_kernel_op_282 +
                           tmp_kernel_op_212 * tmp_kernel_op_283 +
                           tmp_kernel_op_214 * tmp_kernel_op_284 +
                           tmp_kernel_op_216 * tmp_kernel_op_285 +
                           tmp_kernel_op_218 * tmp_kernel_op_286 +
                           tmp_kernel_op_220 * tmp_kernel_op_287) +
              src_dof_5 * tmp_kernel_op_288 + src_dof_6 * tmp_kernel_op_289;
          const walberla::float64 elMatVec_5 =
              src_dof_0 * tmp_kernel_op_134 + src_dof_1 * tmp_kernel_op_189 +
              src_dof_2 * tmp_kernel_op_232 + src_dof_3 * tmp_kernel_op_276 +
              src_dof_4 * tmp_kernel_op_288 +
              src_dof_5 * (tmp_kernel_op_156 * tmp_kernel_op_278 +
                           tmp_kernel_op_157 * tmp_kernel_op_279 +
                           tmp_kernel_op_158 * tmp_kernel_op_280 +
                           tmp_kernel_op_159 * tmp_kernel_op_281 +
                           tmp_kernel_op_160 * tmp_kernel_op_282 +
                           tmp_kernel_op_161 * tmp_kernel_op_283 +
                           tmp_kernel_op_162 * tmp_kernel_op_284 +
                           tmp_kernel_op_163 * tmp_kernel_op_285 +
                           tmp_kernel_op_164 * tmp_kernel_op_286 +
                           tmp_kernel_op_165 * tmp_kernel_op_287) +
              src_dof_6 * tmp_kernel_op_290;
          const walberla::float64 elMatVec_6 =
              src_dof_0 * tmp_kernel_op_155 + src_dof_1 * tmp_kernel_op_200 +
              src_dof_2 * tmp_kernel_op_244 + src_dof_3 * tmp_kernel_op_277 +
              src_dof_4 * tmp_kernel_op_289 + src_dof_5 * tmp_kernel_op_290 +
              src_dof_6 *
                  ((tmp_kernel_op_135 * tmp_kernel_op_135) * tmp_kernel_op_245 +
                   (tmp_kernel_op_137 * tmp_kernel_op_137) * tmp_kernel_op_246 +
                   (tmp_kernel_op_139 * tmp_kernel_op_139) * tmp_kernel_op_247 +
                   (tmp_kernel_op_141 * tmp_kernel_op_141) * tmp_kernel_op_248 +
                   (tmp_kernel_op_143 * tmp_kernel_op_143) * tmp_kernel_op_249 +
                   (tmp_kernel_op_145 * tmp_kernel_op_145) * tmp_kernel_op_250 +
                   (tmp_kernel_op_147 * tmp_kernel_op_147) * tmp_kernel_op_251 +
                   (tmp_kernel_op_149 * tmp_kernel_op_149) * tmp_kernel_op_252 +
                   (tmp_kernel_op_151 * tmp_kernel_op_151) * tmp_kernel_op_253 +
                   (tmp_kernel_op_153 * tmp_kernel_op_153) * tmp_kernel_op_254);
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
void P2PlusBubbleElementwiseMass_AffineMap2D_float64::
    toMatrix_P2PlusBubbleElementwiseMass_AffineMap2D_float64_macro_2D(
        idx_t *RESTRICT _data_dst, idx_t *RESTRICT _data_dstEdge,
        idx_t *RESTRICT _data_dstVertex, idx_t *RESTRICT _data_src,
        idx_t *RESTRICT _data_srcEdge, idx_t *RESTRICT _data_srcVertex,
        walberla::float64 bMat_00, walberla::float64 bMat_01,
        walberla::float64 bMat_10, walberla::float64 bMat_11,
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
      const walberla::float64 abs_det_jac_affine_GRAY =
          abs(jac_affine_0_0_GRAY * jac_affine_1_1_GRAY -
              jac_affine_0_1_GRAY * jac_affine_1_0_GRAY);
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge;
             ctr_0 += 1) {
#if 0
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
#endif
          const walberla::float64 jac_blending_1_1_q_9 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_9 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_9 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_9 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_9 =
              jac_blending_0_0_q_9 * jac_blending_1_1_q_9 -
              jac_blending_0_1_q_9 * jac_blending_1_0_q_9;
          const walberla::float64 jac_blending_1_1_q_8 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_8 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_8 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_8 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_8 =
              jac_blending_0_0_q_8 * jac_blending_1_1_q_8 -
              jac_blending_0_1_q_8 * jac_blending_1_0_q_8;
          const walberla::float64 jac_blending_1_1_q_7 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_7 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_7 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_7 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_7 =
              jac_blending_0_0_q_7 * jac_blending_1_1_q_7 -
              jac_blending_0_1_q_7 * jac_blending_1_0_q_7;
          const walberla::float64 jac_blending_1_1_q_6 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_6 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_6 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_6 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_6 =
              jac_blending_0_0_q_6 * jac_blending_1_1_q_6 -
              jac_blending_0_1_q_6 * jac_blending_1_0_q_6;
          const walberla::float64 jac_blending_1_1_q_5 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_5 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_5 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_5 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_5 =
              jac_blending_0_0_q_5 * jac_blending_1_1_q_5 -
              jac_blending_0_1_q_5 * jac_blending_1_0_q_5;
          const walberla::float64 jac_blending_1_1_q_4 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_4 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_4 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_4 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_4 =
              jac_blending_0_0_q_4 * jac_blending_1_1_q_4 -
              jac_blending_0_1_q_4 * jac_blending_1_0_q_4;
          const walberla::float64 jac_blending_1_1_q_3 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_3 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_3 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_3 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_2 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_2 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_2 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_1 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_1 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_1 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_0 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_0 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_0 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_1 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_2 = -1.0175839942884126;
          const walberla::float64 tmp_kernel_op_3 =
              -tmp_kernel_op_1 - tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_0 *
                                                    0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_5 = 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_6 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_7 = -0.87062238788797819;
          const walberla::float64 tmp_kernel_op_8 =
              -tmp_kernel_op_6 - tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_1 *
                                                    0.052397656566423402;
          const walberla::float64 tmp_kernel_op_10 = 0.40172877323475986;
          const walberla::float64 tmp_kernel_op_11 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_12 = -0.84634856630183575;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_2 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_15 = 0.64933214716985033;
          const walberla::float64 tmp_kernel_op_16 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_17 = -0.40994976126116967;
          const walberla::float64 tmp_kernel_op_18 =
              -tmp_kernel_op_16 - tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_19 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_3 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_20 = -0.0087919971442063094;
          const walberla::float64 tmp_kernel_op_21 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_22 = 1.5773598151468047;
          const walberla::float64 tmp_kernel_op_23 =
              -tmp_kernel_op_21 - tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_24 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_4 *
                                                     0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_25 = 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_26 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_27 = 0.57645841913372808;
          const walberla::float64 tmp_kernel_op_28 =
              -tmp_kernel_op_26 - tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_29 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_5 *
                                                     0.052397656566423402;
          const walberla::float64 tmp_kernel_op_30 = 0.076825716849082126;
          const walberla::float64 tmp_kernel_op_31 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_32 = -0.19654245353048039;
          const walberla::float64 tmp_kernel_op_33 =
              -tmp_kernel_op_31 - tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_34 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_6 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_35 = 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_36 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_37 = 0.29866429433970043;
          const walberla::float64 tmp_kernel_op_38 =
              -tmp_kernel_op_36 - tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_39 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_7 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_40 = 0.085333161170310645;
          const walberla::float64 tmp_kernel_op_41 = 1.6586673553187572;
          const walberla::float64 tmp_kernel_op_42 = -0.8293336776593786;
          const walberla::float64 tmp_kernel_op_43 =
              -tmp_kernel_op_41 - tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_44 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_8 *
                                                     0.045496761795224737;
          const walberla::float64 tmp_kernel_op_45 = 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_46 = 0.58523781332014213;
          const walberla::float64 tmp_kernel_op_47 = -0.29261890666007107;
          const walberla::float64 tmp_kernel_op_48 =
              -tmp_kernel_op_46 - tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_49 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_9 *
                                                     0.10566414783403316;
          const walberla::float64 tmp_kernel_op_50 =
              tmp_kernel_op_4 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_51 =
              tmp_kernel_op_0 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_52 =
              tmp_kernel_op_50 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_53 =
              tmp_kernel_op_9 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_5 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_53 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_56 =
              tmp_kernel_op_14 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_57 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_58 =
              tmp_kernel_op_56 * tmp_kernel_op_57;
          const walberla::float64 tmp_kernel_op_59 =
              tmp_kernel_op_19 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_15 * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_61 =
              tmp_kernel_op_59 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_24 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_63 =
              tmp_kernel_op_20 * tmp_kernel_op_23;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_62 * tmp_kernel_op_63;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_29 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_25 * tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_65 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_34 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_30 * tmp_kernel_op_33;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_68 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_39 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_72 =
              tmp_kernel_op_35 * tmp_kernel_op_38;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_71 * tmp_kernel_op_72;
          const walberla::float64 tmp_kernel_op_74 =
              tmp_kernel_op_44 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_75 =
              tmp_kernel_op_40 * tmp_kernel_op_43;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_74 * tmp_kernel_op_75;
          const walberla::float64 tmp_kernel_op_77 =
              tmp_kernel_op_49 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_45 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_77 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_12 * tmp_kernel_op_58 +
              tmp_kernel_op_17 * tmp_kernel_op_61 +
              tmp_kernel_op_2 * tmp_kernel_op_52 +
              tmp_kernel_op_22 * tmp_kernel_op_64 +
              tmp_kernel_op_27 * tmp_kernel_op_67 +
              tmp_kernel_op_32 * tmp_kernel_op_70 +
              tmp_kernel_op_37 * tmp_kernel_op_73 +
              tmp_kernel_op_42 * tmp_kernel_op_76 +
              tmp_kernel_op_47 * tmp_kernel_op_79 +
              tmp_kernel_op_55 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_81 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_4 * tmp_kernel_op_51 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_83 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_54 * tmp_kernel_op_9 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_85 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_14 * tmp_kernel_op_57 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_87 = tmp_kernel_op_16 - 1.0;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_19 * tmp_kernel_op_60 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_89 = tmp_kernel_op_21 - 1.0;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_24 * tmp_kernel_op_63 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_91 = tmp_kernel_op_26 - 1.0;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_29 * tmp_kernel_op_66 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_93 = tmp_kernel_op_31 - 1.0;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_34 * tmp_kernel_op_69 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_95 = tmp_kernel_op_36 - 1.0;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_39 * tmp_kernel_op_72 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_97 = tmp_kernel_op_41 - 1.0;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_44 * tmp_kernel_op_75 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_99 = tmp_kernel_op_46 - 1.0;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_49 * tmp_kernel_op_78 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_100 * tmp_kernel_op_99 +
              tmp_kernel_op_81 * tmp_kernel_op_82 +
              tmp_kernel_op_83 * tmp_kernel_op_84 +
              tmp_kernel_op_85 * tmp_kernel_op_86 +
              tmp_kernel_op_87 * tmp_kernel_op_88 +
              tmp_kernel_op_89 * tmp_kernel_op_90 +
              tmp_kernel_op_91 * tmp_kernel_op_92 +
              tmp_kernel_op_93 * tmp_kernel_op_94 +
              tmp_kernel_op_95 * tmp_kernel_op_96 +
              tmp_kernel_op_97 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_102 = -0.035167988576825245;
          const walberla::float64 tmp_kernel_op_103 = 0.25875522422404362;
          const walberla::float64 tmp_kernel_op_104 = 0.30730286739632839;
          const walberla::float64 tmp_kernel_op_105 = 1.1801004774776607;
          const walberla::float64 tmp_kernel_op_106 = 5.1547196302936094;
          const walberla::float64 tmp_kernel_op_107 = 3.1529168382674562;
          const walberla::float64 tmp_kernel_op_108 = 1.6069150929390392;
          const walberla::float64 tmp_kernel_op_109 = 2.5973285886794009;
          const walberla::float64 tmp_kernel_op_110 = 0.3413326446812428;
          const walberla::float64 tmp_kernel_op_111 = 1.4147621866798579;
          const walberla::float64 tmp_kernel_op_112 =
              tmp_kernel_op_100 * tmp_kernel_op_111 +
              tmp_kernel_op_102 * tmp_kernel_op_82 +
              tmp_kernel_op_103 * tmp_kernel_op_84 +
              tmp_kernel_op_104 * tmp_kernel_op_86 +
              tmp_kernel_op_105 * tmp_kernel_op_88 +
              tmp_kernel_op_106 * tmp_kernel_op_90 +
              tmp_kernel_op_107 * tmp_kernel_op_92 +
              tmp_kernel_op_108 * tmp_kernel_op_94 +
              tmp_kernel_op_109 * tmp_kernel_op_96 +
              tmp_kernel_op_110 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_113 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_114 =
              -tmp_kernel_op_102 - tmp_kernel_op_113 + 4.0;
          const walberla::float64 tmp_kernel_op_115 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_116 =
              -tmp_kernel_op_103 - tmp_kernel_op_115 + 4.0;
          const walberla::float64 tmp_kernel_op_117 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_118 =
              -tmp_kernel_op_104 - tmp_kernel_op_117 + 4.0;
          const walberla::float64 tmp_kernel_op_119 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_120 =
              -tmp_kernel_op_105 - tmp_kernel_op_119 + 4.0;
          const walberla::float64 tmp_kernel_op_121 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_122 =
              -tmp_kernel_op_106 - tmp_kernel_op_121 + 4.0;
          const walberla::float64 tmp_kernel_op_123 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_124 =
              -tmp_kernel_op_107 - tmp_kernel_op_123 + 4.0;
          const walberla::float64 tmp_kernel_op_125 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_126 =
              -tmp_kernel_op_108 - tmp_kernel_op_125 + 4.0;
          const walberla::float64 tmp_kernel_op_127 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_128 =
              -tmp_kernel_op_109 - tmp_kernel_op_127 + 4.0;
          const walberla::float64 tmp_kernel_op_129 = 3.3173347106375144;
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_110 - tmp_kernel_op_129 + 4.0;
          const walberla::float64 tmp_kernel_op_131 = 1.1704756266402843;
          const walberla::float64 tmp_kernel_op_132 =
              -tmp_kernel_op_111 - tmp_kernel_op_131 + 4.0;
          const walberla::float64 tmp_kernel_op_133 =
              tmp_kernel_op_100 * tmp_kernel_op_132 +
              tmp_kernel_op_114 * tmp_kernel_op_82 +
              tmp_kernel_op_116 * tmp_kernel_op_84 +
              tmp_kernel_op_118 * tmp_kernel_op_86 +
              tmp_kernel_op_120 * tmp_kernel_op_88 +
              tmp_kernel_op_122 * tmp_kernel_op_90 +
              tmp_kernel_op_124 * tmp_kernel_op_92 +
              tmp_kernel_op_126 * tmp_kernel_op_94 +
              tmp_kernel_op_128 * tmp_kernel_op_96 +
              tmp_kernel_op_130 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_134 =
              tmp_kernel_op_114 * tmp_kernel_op_52 +
              tmp_kernel_op_116 * tmp_kernel_op_55 +
              tmp_kernel_op_118 * tmp_kernel_op_58 +
              tmp_kernel_op_120 * tmp_kernel_op_61 +
              tmp_kernel_op_122 * tmp_kernel_op_64 +
              tmp_kernel_op_124 * tmp_kernel_op_67 +
              tmp_kernel_op_126 * tmp_kernel_op_70 +
              tmp_kernel_op_128 * tmp_kernel_op_73 +
              tmp_kernel_op_130 * tmp_kernel_op_76 +
              tmp_kernel_op_132 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_135 = 34.794357504481866;
          const walberla::float64 tmp_kernel_op_136 =
              tmp_kernel_op_50 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_137 = 21.28218865830533;
          const walberla::float64 tmp_kernel_op_138 =
              tmp_kernel_op_53 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_139 = 10.846676877338515;
          const walberla::float64 tmp_kernel_op_140 =
              tmp_kernel_op_56 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_141 = 17.531967973585957;
          const walberla::float64 tmp_kernel_op_142 =
              tmp_kernel_op_59 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_143 = -0.23738392289357257;
          const walberla::float64 tmp_kernel_op_144 =
              tmp_kernel_op_62 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_145 = 1.7465977635122938;
          const walberla::float64 tmp_kernel_op_146 =
              tmp_kernel_op_65 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_147 = 2.0742943549252164;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_68 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_149 = 7.9656782229742085;
          const walberla::float64 tmp_kernel_op_150 =
              tmp_kernel_op_71 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_151 = 2.3039953515983882;
          const walberla::float64 tmp_kernel_op_152 =
              tmp_kernel_op_74 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_153 = 9.5496447600890395;
          const walberla::float64 tmp_kernel_op_154 =
              tmp_kernel_op_77 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_155 =
              tmp_kernel_op_135 * tmp_kernel_op_136 * tmp_kernel_op_51 +
              tmp_kernel_op_137 * tmp_kernel_op_138 * tmp_kernel_op_54 +
              tmp_kernel_op_139 * tmp_kernel_op_140 * tmp_kernel_op_57 +
              tmp_kernel_op_141 * tmp_kernel_op_142 * tmp_kernel_op_60 +
              tmp_kernel_op_143 * tmp_kernel_op_144 * tmp_kernel_op_63 +
              tmp_kernel_op_145 * tmp_kernel_op_146 * tmp_kernel_op_66 +
              tmp_kernel_op_147 * tmp_kernel_op_148 * tmp_kernel_op_69 +
              tmp_kernel_op_149 * tmp_kernel_op_150 * tmp_kernel_op_72 +
              tmp_kernel_op_151 * tmp_kernel_op_152 * tmp_kernel_op_75 +
              tmp_kernel_op_153 * tmp_kernel_op_154 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_4 * 7.7299213783731935e-5;
          const walberla::float64 tmp_kernel_op_157 =
              tmp_kernel_op_9 * 0.0041846416289521935;
          const walberla::float64 tmp_kernel_op_158 =
              tmp_kernel_op_14 * 0.0059021907693753367;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_19 * 0.087039821058937664;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_24 * 1.6606959041833929;
          const walberla::float64 tmp_kernel_op_161 =
              tmp_kernel_op_29 * 0.62130528681440322;
          const walberla::float64 tmp_kernel_op_162 =
              tmp_kernel_op_34 * 0.16138600724470506;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_39 * 0.42163223734820804;
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_44 * 0.0072817483953182219;
          const walberla::float64 tmp_kernel_op_165 =
              tmp_kernel_op_49 * 0.12509700280369832;
          const walberla::float64 tmp_kernel_op_166 =
              tmp_kernel_op_136 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_138 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_12 * tmp_kernel_op_140;
          const walberla::float64 tmp_kernel_op_169 =
              tmp_kernel_op_142 * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_170 =
              tmp_kernel_op_144 * tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_146 * tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_148 * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_173 =
              tmp_kernel_op_150 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_174 =
              tmp_kernel_op_152 * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_154 * tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_166 * tmp_kernel_op_81 +
              tmp_kernel_op_167 * tmp_kernel_op_83 +
              tmp_kernel_op_168 * tmp_kernel_op_85 +
              tmp_kernel_op_169 * tmp_kernel_op_87 +
              tmp_kernel_op_170 * tmp_kernel_op_89 +
              tmp_kernel_op_171 * tmp_kernel_op_91 +
              tmp_kernel_op_172 * tmp_kernel_op_93 +
              tmp_kernel_op_173 * tmp_kernel_op_95 +
              tmp_kernel_op_174 * tmp_kernel_op_97 +
              tmp_kernel_op_175 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_156 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_178 =
              tmp_kernel_op_157 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_179 =
              tmp_kernel_op_12 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_180 =
              tmp_kernel_op_159 * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_181 =
              tmp_kernel_op_160 * tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_182 =
              tmp_kernel_op_161 * tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_183 =
              tmp_kernel_op_162 * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_184 =
              tmp_kernel_op_163 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_185 =
              tmp_kernel_op_164 * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_165 * tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_113 * tmp_kernel_op_177 +
              tmp_kernel_op_115 * tmp_kernel_op_178 +
              tmp_kernel_op_117 * tmp_kernel_op_179 +
              tmp_kernel_op_119 * tmp_kernel_op_180 +
              tmp_kernel_op_121 * tmp_kernel_op_181 +
              tmp_kernel_op_123 * tmp_kernel_op_182 +
              tmp_kernel_op_125 * tmp_kernel_op_183 +
              tmp_kernel_op_127 * tmp_kernel_op_184 +
              tmp_kernel_op_129 * tmp_kernel_op_185 +
              tmp_kernel_op_131 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_114 * tmp_kernel_op_166 +
              tmp_kernel_op_116 * tmp_kernel_op_167 +
              tmp_kernel_op_118 * tmp_kernel_op_168 +
              tmp_kernel_op_120 * tmp_kernel_op_169 +
              tmp_kernel_op_122 * tmp_kernel_op_170 +
              tmp_kernel_op_124 * tmp_kernel_op_171 +
              tmp_kernel_op_126 * tmp_kernel_op_172 +
              tmp_kernel_op_128 * tmp_kernel_op_173 +
              tmp_kernel_op_130 * tmp_kernel_op_174 +
              tmp_kernel_op_132 * tmp_kernel_op_175;
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_114 * tmp_kernel_op_177 +
              tmp_kernel_op_116 * tmp_kernel_op_178 +
              tmp_kernel_op_118 * tmp_kernel_op_179 +
              tmp_kernel_op_120 * tmp_kernel_op_180 +
              tmp_kernel_op_122 * tmp_kernel_op_181 +
              tmp_kernel_op_124 * tmp_kernel_op_182 +
              tmp_kernel_op_126 * tmp_kernel_op_183 +
              tmp_kernel_op_128 * tmp_kernel_op_184 +
              tmp_kernel_op_130 * tmp_kernel_op_185 +
              tmp_kernel_op_132 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_135 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_191 =
              tmp_kernel_op_137 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_192 =
              tmp_kernel_op_139 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_193 =
              tmp_kernel_op_141 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_194 =
              tmp_kernel_op_143 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_195 =
              tmp_kernel_op_145 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_196 =
              tmp_kernel_op_147 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_197 =
              tmp_kernel_op_149 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_198 =
              tmp_kernel_op_151 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_153 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_177 * tmp_kernel_op_190 +
              tmp_kernel_op_178 * tmp_kernel_op_191 +
              tmp_kernel_op_179 * tmp_kernel_op_192 +
              tmp_kernel_op_180 * tmp_kernel_op_193 +
              tmp_kernel_op_181 * tmp_kernel_op_194 +
              tmp_kernel_op_182 * tmp_kernel_op_195 +
              tmp_kernel_op_183 * tmp_kernel_op_196 +
              tmp_kernel_op_184 * tmp_kernel_op_197 +
              tmp_kernel_op_185 * tmp_kernel_op_198 +
              tmp_kernel_op_186 * tmp_kernel_op_199;
          const walberla::float64 tmp_kernel_op_201 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_201 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_203 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_203 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_205 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_14 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_207 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_19 * tmp_kernel_op_207;
          const walberla::float64 tmp_kernel_op_209 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_210 =
              tmp_kernel_op_209 * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_211 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_212 =
              tmp_kernel_op_211 * tmp_kernel_op_29;
          const walberla::float64 tmp_kernel_op_213 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_214 =
              tmp_kernel_op_213 * tmp_kernel_op_34;
          const walberla::float64 tmp_kernel_op_215 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_216 =
              tmp_kernel_op_215 * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_217 = 0.68779434890003011;
          const walberla::float64 tmp_kernel_op_218 =
              tmp_kernel_op_217 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_219 = 0.085625824534935377;
          const walberla::float64 tmp_kernel_op_220 =
              tmp_kernel_op_219 * tmp_kernel_op_49;
          const walberla::float64 tmp_kernel_op_221 =
              tmp_kernel_op_202 * tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_222 =
              tmp_kernel_op_204 * tmp_kernel_op_83;
          const walberla::float64 tmp_kernel_op_223 =
              tmp_kernel_op_206 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_224 =
              tmp_kernel_op_208 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_225 =
              tmp_kernel_op_210 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_226 =
              tmp_kernel_op_212 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_227 =
              tmp_kernel_op_214 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_228 =
              tmp_kernel_op_216 * tmp_kernel_op_95;
          const walberla::float64 tmp_kernel_op_229 =
              tmp_kernel_op_218 * tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_230 =
              tmp_kernel_op_220 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_231 =
              tmp_kernel_op_102 * tmp_kernel_op_221 +
              tmp_kernel_op_103 * tmp_kernel_op_222 +
              tmp_kernel_op_104 * tmp_kernel_op_223 +
              tmp_kernel_op_105 * tmp_kernel_op_224 +
              tmp_kernel_op_106 * tmp_kernel_op_225 +
              tmp_kernel_op_107 * tmp_kernel_op_226 +
              tmp_kernel_op_108 * tmp_kernel_op_227 +
              tmp_kernel_op_109 * tmp_kernel_op_228 +
              tmp_kernel_op_110 * tmp_kernel_op_229 +
              tmp_kernel_op_111 * tmp_kernel_op_230;
          const walberla::float64 tmp_kernel_op_232 =
              tmp_kernel_op_114 * tmp_kernel_op_221 +
              tmp_kernel_op_116 * tmp_kernel_op_222 +
              tmp_kernel_op_118 * tmp_kernel_op_223 +
              tmp_kernel_op_120 * tmp_kernel_op_224 +
              tmp_kernel_op_122 * tmp_kernel_op_225 +
              tmp_kernel_op_124 * tmp_kernel_op_226 +
              tmp_kernel_op_126 * tmp_kernel_op_227 +
              tmp_kernel_op_128 * tmp_kernel_op_228 +
              tmp_kernel_op_130 * tmp_kernel_op_229 +
              tmp_kernel_op_132 * tmp_kernel_op_230;
          const walberla::float64 tmp_kernel_op_233 =
              tmp_kernel_op_114 * tmp_kernel_op_136 * tmp_kernel_op_81 +
              tmp_kernel_op_116 * tmp_kernel_op_138 * tmp_kernel_op_83 +
              tmp_kernel_op_118 * tmp_kernel_op_140 * tmp_kernel_op_85 +
              tmp_kernel_op_120 * tmp_kernel_op_142 * tmp_kernel_op_87 +
              tmp_kernel_op_122 * tmp_kernel_op_144 * tmp_kernel_op_89 +
              tmp_kernel_op_124 * tmp_kernel_op_146 * tmp_kernel_op_91 +
              tmp_kernel_op_126 * tmp_kernel_op_148 * tmp_kernel_op_93 +
              tmp_kernel_op_128 * tmp_kernel_op_150 * tmp_kernel_op_95 +
              tmp_kernel_op_130 * tmp_kernel_op_152 * tmp_kernel_op_97 +
              tmp_kernel_op_132 * tmp_kernel_op_154 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_234 =
              tmp_kernel_op_135 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_235 =
              tmp_kernel_op_137 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_236 =
              tmp_kernel_op_139 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_237 =
              tmp_kernel_op_141 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_238 =
              tmp_kernel_op_143 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_239 =
              tmp_kernel_op_145 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_240 =
              tmp_kernel_op_147 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_241 =
              tmp_kernel_op_149 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_242 =
              tmp_kernel_op_151 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_243 =
              tmp_kernel_op_153 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_244 =
              tmp_kernel_op_221 * tmp_kernel_op_234 +
              tmp_kernel_op_222 * tmp_kernel_op_235 +
              tmp_kernel_op_223 * tmp_kernel_op_236 +
              tmp_kernel_op_224 * tmp_kernel_op_237 +
              tmp_kernel_op_225 * tmp_kernel_op_238 +
              tmp_kernel_op_226 * tmp_kernel_op_239 +
              tmp_kernel_op_227 * tmp_kernel_op_240 +
              tmp_kernel_op_228 * tmp_kernel_op_241 +
              tmp_kernel_op_229 * tmp_kernel_op_242 +
              tmp_kernel_op_230 * tmp_kernel_op_243;
          const walberla::float64 tmp_kernel_op_245 =
              tmp_kernel_op_156 * tmp_kernel_op_201;
          const walberla::float64 tmp_kernel_op_246 =
              tmp_kernel_op_157 * tmp_kernel_op_203;
          const walberla::float64 tmp_kernel_op_247 =
              tmp_kernel_op_158 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_248 =
              tmp_kernel_op_159 * tmp_kernel_op_207;
          const walberla::float64 tmp_kernel_op_249 =
              tmp_kernel_op_160 * tmp_kernel_op_209;
          const walberla::float64 tmp_kernel_op_250 =
              tmp_kernel_op_161 * tmp_kernel_op_211;
          const walberla::float64 tmp_kernel_op_251 =
              tmp_kernel_op_162 * tmp_kernel_op_213;
          const walberla::float64 tmp_kernel_op_252 =
              tmp_kernel_op_163 * tmp_kernel_op_215;
          const walberla::float64 tmp_kernel_op_253 =
              tmp_kernel_op_164 * tmp_kernel_op_217;
          const walberla::float64 tmp_kernel_op_254 =
              tmp_kernel_op_165 * tmp_kernel_op_219;
          const walberla::float64 tmp_kernel_op_255 =
              tmp_kernel_op_114 * tmp_kernel_op_202;
          const walberla::float64 tmp_kernel_op_256 =
              tmp_kernel_op_116 * tmp_kernel_op_204;
          const walberla::float64 tmp_kernel_op_257 =
              tmp_kernel_op_118 * tmp_kernel_op_206;
          const walberla::float64 tmp_kernel_op_258 =
              tmp_kernel_op_120 * tmp_kernel_op_208;
          const walberla::float64 tmp_kernel_op_259 =
              tmp_kernel_op_122 * tmp_kernel_op_210;
          const walberla::float64 tmp_kernel_op_260 =
              tmp_kernel_op_124 * tmp_kernel_op_212;
          const walberla::float64 tmp_kernel_op_261 =
              tmp_kernel_op_126 * tmp_kernel_op_214;
          const walberla::float64 tmp_kernel_op_262 =
              tmp_kernel_op_128 * tmp_kernel_op_216;
          const walberla::float64 tmp_kernel_op_263 =
              tmp_kernel_op_130 * tmp_kernel_op_218;
          const walberla::float64 tmp_kernel_op_264 =
              tmp_kernel_op_132 * tmp_kernel_op_220;
          const walberla::float64 tmp_kernel_op_265 =
              tmp_kernel_op_102 * tmp_kernel_op_255 +
              tmp_kernel_op_103 * tmp_kernel_op_256 +
              tmp_kernel_op_104 * tmp_kernel_op_257 +
              tmp_kernel_op_105 * tmp_kernel_op_258 +
              tmp_kernel_op_106 * tmp_kernel_op_259 +
              tmp_kernel_op_107 * tmp_kernel_op_260 +
              tmp_kernel_op_108 * tmp_kernel_op_261 +
              tmp_kernel_op_109 * tmp_kernel_op_262 +
              tmp_kernel_op_110 * tmp_kernel_op_263 +
              tmp_kernel_op_111 * tmp_kernel_op_264;
          const walberla::float64 tmp_kernel_op_266 =
              tmp_kernel_op_114 * tmp_kernel_op_156;
          const walberla::float64 tmp_kernel_op_267 =
              tmp_kernel_op_116 * tmp_kernel_op_157;
          const walberla::float64 tmp_kernel_op_268 =
              tmp_kernel_op_118 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_269 =
              tmp_kernel_op_120 * tmp_kernel_op_159;
          const walberla::float64 tmp_kernel_op_270 =
              tmp_kernel_op_122 * tmp_kernel_op_160;
          const walberla::float64 tmp_kernel_op_271 =
              tmp_kernel_op_124 * tmp_kernel_op_161;
          const walberla::float64 tmp_kernel_op_272 =
              tmp_kernel_op_126 * tmp_kernel_op_162;
          const walberla::float64 tmp_kernel_op_273 =
              tmp_kernel_op_128 * tmp_kernel_op_163;
          const walberla::float64 tmp_kernel_op_274 =
              tmp_kernel_op_130 * tmp_kernel_op_164;
          const walberla::float64 tmp_kernel_op_275 =
              tmp_kernel_op_132 * tmp_kernel_op_165;
          const walberla::float64 tmp_kernel_op_276 =
              tmp_kernel_op_113 * tmp_kernel_op_266 +
              tmp_kernel_op_115 * tmp_kernel_op_267 +
              tmp_kernel_op_117 * tmp_kernel_op_268 +
              tmp_kernel_op_119 * tmp_kernel_op_269 +
              tmp_kernel_op_121 * tmp_kernel_op_270 +
              tmp_kernel_op_123 * tmp_kernel_op_271 +
              tmp_kernel_op_125 * tmp_kernel_op_272 +
              tmp_kernel_op_127 * tmp_kernel_op_273 +
              tmp_kernel_op_129 * tmp_kernel_op_274 +
              tmp_kernel_op_131 * tmp_kernel_op_275;
          const walberla::float64 tmp_kernel_op_277 =
              tmp_kernel_op_135 * tmp_kernel_op_245 * 4.0 +
              tmp_kernel_op_137 * tmp_kernel_op_246 * 4.0 +
              tmp_kernel_op_139 * tmp_kernel_op_247 * 4.0 +
              tmp_kernel_op_141 * tmp_kernel_op_248 * 4.0 +
              tmp_kernel_op_143 * tmp_kernel_op_249 * 4.0 +
              tmp_kernel_op_145 * tmp_kernel_op_250 * 4.0 +
              tmp_kernel_op_147 * tmp_kernel_op_251 * 4.0 +
              tmp_kernel_op_149 * tmp_kernel_op_252 * 4.0 +
              tmp_kernel_op_151 * tmp_kernel_op_253 * 4.0 +
              tmp_kernel_op_153 * tmp_kernel_op_254 * 4.0;
          const walberla::float64 tmp_kernel_op_278 =
              (tmp_kernel_op_114 * tmp_kernel_op_114);
          const walberla::float64 tmp_kernel_op_279 =
              (tmp_kernel_op_116 * tmp_kernel_op_116);
          const walberla::float64 tmp_kernel_op_280 =
              (tmp_kernel_op_118 * tmp_kernel_op_118);
          const walberla::float64 tmp_kernel_op_281 =
              (tmp_kernel_op_120 * tmp_kernel_op_120);
          const walberla::float64 tmp_kernel_op_282 =
              (tmp_kernel_op_122 * tmp_kernel_op_122);
          const walberla::float64 tmp_kernel_op_283 =
              (tmp_kernel_op_124 * tmp_kernel_op_124);
          const walberla::float64 tmp_kernel_op_284 =
              (tmp_kernel_op_126 * tmp_kernel_op_126);
          const walberla::float64 tmp_kernel_op_285 =
              (tmp_kernel_op_128 * tmp_kernel_op_128);
          const walberla::float64 tmp_kernel_op_286 =
              (tmp_kernel_op_130 * tmp_kernel_op_130);
          const walberla::float64 tmp_kernel_op_287 =
              (tmp_kernel_op_132 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_288 =
              tmp_kernel_op_136 * tmp_kernel_op_278 +
              tmp_kernel_op_138 * tmp_kernel_op_279 +
              tmp_kernel_op_140 * tmp_kernel_op_280 +
              tmp_kernel_op_142 * tmp_kernel_op_281 +
              tmp_kernel_op_144 * tmp_kernel_op_282 +
              tmp_kernel_op_146 * tmp_kernel_op_283 +
              tmp_kernel_op_148 * tmp_kernel_op_284 +
              tmp_kernel_op_150 * tmp_kernel_op_285 +
              tmp_kernel_op_152 * tmp_kernel_op_286 +
              tmp_kernel_op_154 * tmp_kernel_op_287;
          const walberla::float64 tmp_kernel_op_289 =
              tmp_kernel_op_234 * tmp_kernel_op_255 +
              tmp_kernel_op_235 * tmp_kernel_op_256 +
              tmp_kernel_op_236 * tmp_kernel_op_257 +
              tmp_kernel_op_237 * tmp_kernel_op_258 +
              tmp_kernel_op_238 * tmp_kernel_op_259 +
              tmp_kernel_op_239 * tmp_kernel_op_260 +
              tmp_kernel_op_240 * tmp_kernel_op_261 +
              tmp_kernel_op_241 * tmp_kernel_op_262 +
              tmp_kernel_op_242 * tmp_kernel_op_263 +
              tmp_kernel_op_243 * tmp_kernel_op_264;
          const walberla::float64 tmp_kernel_op_290 =
              tmp_kernel_op_190 * tmp_kernel_op_266 +
              tmp_kernel_op_191 * tmp_kernel_op_267 +
              tmp_kernel_op_192 * tmp_kernel_op_268 +
              tmp_kernel_op_193 * tmp_kernel_op_269 +
              tmp_kernel_op_194 * tmp_kernel_op_270 +
              tmp_kernel_op_195 * tmp_kernel_op_271 +
              tmp_kernel_op_196 * tmp_kernel_op_272 +
              tmp_kernel_op_197 * tmp_kernel_op_273 +
              tmp_kernel_op_198 * tmp_kernel_op_274 +
              tmp_kernel_op_199 * tmp_kernel_op_275;
          const walberla::float64 elMat_0_0 =
              (tmp_kernel_op_0 * tmp_kernel_op_0) *
                  (tmp_kernel_op_3 * tmp_kernel_op_3) * tmp_kernel_op_4 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) *
                  (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_14 +
              (tmp_kernel_op_15 * tmp_kernel_op_15) *
                  (tmp_kernel_op_18 * tmp_kernel_op_18) * tmp_kernel_op_19 +
              (tmp_kernel_op_20 * tmp_kernel_op_20) *
                  (tmp_kernel_op_23 * tmp_kernel_op_23) * tmp_kernel_op_24 +
              (tmp_kernel_op_25 * tmp_kernel_op_25) *
                  (tmp_kernel_op_28 * tmp_kernel_op_28) * tmp_kernel_op_29 +
              (tmp_kernel_op_30 * tmp_kernel_op_30) *
                  (tmp_kernel_op_33 * tmp_kernel_op_33) * tmp_kernel_op_34 +
              (tmp_kernel_op_35 * tmp_kernel_op_35) *
                  (tmp_kernel_op_38 * tmp_kernel_op_38) * tmp_kernel_op_39 +
              (tmp_kernel_op_40 * tmp_kernel_op_40) *
                  (tmp_kernel_op_43 * tmp_kernel_op_43) * tmp_kernel_op_44 +
              (tmp_kernel_op_45 * tmp_kernel_op_45) *
                  (tmp_kernel_op_48 * tmp_kernel_op_48) * tmp_kernel_op_49 +
              (tmp_kernel_op_5 * tmp_kernel_op_5) *
                  (tmp_kernel_op_8 * tmp_kernel_op_8) * tmp_kernel_op_9;
          const walberla::float64 elMat_0_1 = tmp_kernel_op_80;
          const walberla::float64 elMat_0_2 = tmp_kernel_op_101;
          const walberla::float64 elMat_0_3 = tmp_kernel_op_112;
          const walberla::float64 elMat_0_4 = tmp_kernel_op_133;
          const walberla::float64 elMat_0_5 = tmp_kernel_op_134;
          const walberla::float64 elMat_0_6 = tmp_kernel_op_155;
          const walberla::float64 elMat_1_0 = tmp_kernel_op_80;
          const walberla::float64 elMat_1_1 =
              (tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_158 +
              tmp_kernel_op_156 * (tmp_kernel_op_2 * tmp_kernel_op_2) +
              tmp_kernel_op_157 * (tmp_kernel_op_7 * tmp_kernel_op_7) +
              tmp_kernel_op_159 * (tmp_kernel_op_17 * tmp_kernel_op_17) +
              tmp_kernel_op_160 * (tmp_kernel_op_22 * tmp_kernel_op_22) +
              tmp_kernel_op_161 * (tmp_kernel_op_27 * tmp_kernel_op_27) +
              tmp_kernel_op_162 * (tmp_kernel_op_32 * tmp_kernel_op_32) +
              tmp_kernel_op_163 * (tmp_kernel_op_37 * tmp_kernel_op_37) +
              tmp_kernel_op_164 * (tmp_kernel_op_42 * tmp_kernel_op_42) +
              tmp_kernel_op_165 * (tmp_kernel_op_47 * tmp_kernel_op_47);
          const walberla::float64 elMat_1_2 = tmp_kernel_op_176;
          const walberla::float64 elMat_1_3 = tmp_kernel_op_187;
          const walberla::float64 elMat_1_4 = tmp_kernel_op_188;
          const walberla::float64 elMat_1_5 = tmp_kernel_op_189;
          const walberla::float64 elMat_1_6 = tmp_kernel_op_200;
          const walberla::float64 elMat_2_0 = tmp_kernel_op_101;
          const walberla::float64 elMat_2_1 = tmp_kernel_op_176;
          const walberla::float64 elMat_2_2 =
              tmp_kernel_op_202 * (tmp_kernel_op_81 * tmp_kernel_op_81) +
              tmp_kernel_op_204 * (tmp_kernel_op_83 * tmp_kernel_op_83) +
              tmp_kernel_op_206 * (tmp_kernel_op_85 * tmp_kernel_op_85) +
              tmp_kernel_op_208 * (tmp_kernel_op_87 * tmp_kernel_op_87) +
              tmp_kernel_op_210 * (tmp_kernel_op_89 * tmp_kernel_op_89) +
              tmp_kernel_op_212 * (tmp_kernel_op_91 * tmp_kernel_op_91) +
              tmp_kernel_op_214 * (tmp_kernel_op_93 * tmp_kernel_op_93) +
              tmp_kernel_op_216 * (tmp_kernel_op_95 * tmp_kernel_op_95) +
              tmp_kernel_op_218 * (tmp_kernel_op_97 * tmp_kernel_op_97) +
              tmp_kernel_op_220 * (tmp_kernel_op_99 * tmp_kernel_op_99);
          const walberla::float64 elMat_2_3 = tmp_kernel_op_231;
          const walberla::float64 elMat_2_4 = tmp_kernel_op_232;
          const walberla::float64 elMat_2_5 = tmp_kernel_op_233;
          const walberla::float64 elMat_2_6 = tmp_kernel_op_244;
          const walberla::float64 elMat_3_0 = tmp_kernel_op_112;
          const walberla::float64 elMat_3_1 = tmp_kernel_op_187;
          const walberla::float64 elMat_3_2 = tmp_kernel_op_231;
          const walberla::float64 elMat_3_3 =
              tmp_kernel_op_245 * 16.0 + tmp_kernel_op_246 * 16.0 +
              tmp_kernel_op_247 * 16.0 + tmp_kernel_op_248 * 16.0 +
              tmp_kernel_op_249 * 16.0 + tmp_kernel_op_250 * 16.0 +
              tmp_kernel_op_251 * 16.0 + tmp_kernel_op_252 * 16.0 +
              tmp_kernel_op_253 * 16.0 + tmp_kernel_op_254 * 16.0;
          const walberla::float64 elMat_3_4 = tmp_kernel_op_265;
          const walberla::float64 elMat_3_5 = tmp_kernel_op_276;
          const walberla::float64 elMat_3_6 = tmp_kernel_op_277;
          const walberla::float64 elMat_4_0 = tmp_kernel_op_133;
          const walberla::float64 elMat_4_1 = tmp_kernel_op_188;
          const walberla::float64 elMat_4_2 = tmp_kernel_op_232;
          const walberla::float64 elMat_4_3 = tmp_kernel_op_265;
          const walberla::float64 elMat_4_4 =
              tmp_kernel_op_202 * tmp_kernel_op_278 +
              tmp_kernel_op_204 * tmp_kernel_op_279 +
              tmp_kernel_op_206 * tmp_kernel_op_280 +
              tmp_kernel_op_208 * tmp_kernel_op_281 +
              tmp_kernel_op_210 * tmp_kernel_op_282 +
              tmp_kernel_op_212 * tmp_kernel_op_283 +
              tmp_kernel_op_214 * tmp_kernel_op_284 +
              tmp_kernel_op_216 * tmp_kernel_op_285 +
              tmp_kernel_op_218 * tmp_kernel_op_286 +
              tmp_kernel_op_220 * tmp_kernel_op_287;
          const walberla::float64 elMat_4_5 = tmp_kernel_op_288;
          const walberla::float64 elMat_4_6 = tmp_kernel_op_289;
          const walberla::float64 elMat_5_0 = tmp_kernel_op_134;
          const walberla::float64 elMat_5_1 = tmp_kernel_op_189;
          const walberla::float64 elMat_5_2 = tmp_kernel_op_233;
          const walberla::float64 elMat_5_3 = tmp_kernel_op_276;
          const walberla::float64 elMat_5_4 = tmp_kernel_op_288;
          const walberla::float64 elMat_5_5 =
              tmp_kernel_op_156 * tmp_kernel_op_278 +
              tmp_kernel_op_157 * tmp_kernel_op_279 +
              tmp_kernel_op_158 * tmp_kernel_op_280 +
              tmp_kernel_op_159 * tmp_kernel_op_281 +
              tmp_kernel_op_160 * tmp_kernel_op_282 +
              tmp_kernel_op_161 * tmp_kernel_op_283 +
              tmp_kernel_op_162 * tmp_kernel_op_284 +
              tmp_kernel_op_163 * tmp_kernel_op_285 +
              tmp_kernel_op_164 * tmp_kernel_op_286 +
              tmp_kernel_op_165 * tmp_kernel_op_287;
          const walberla::float64 elMat_5_6 = tmp_kernel_op_290;
          const walberla::float64 elMat_6_0 = tmp_kernel_op_155;
          const walberla::float64 elMat_6_1 = tmp_kernel_op_200;
          const walberla::float64 elMat_6_2 = tmp_kernel_op_244;
          const walberla::float64 elMat_6_3 = tmp_kernel_op_277;
          const walberla::float64 elMat_6_4 = tmp_kernel_op_289;
          const walberla::float64 elMat_6_5 = tmp_kernel_op_290;
          const walberla::float64 elMat_6_6 =
              (tmp_kernel_op_135 * tmp_kernel_op_135) * tmp_kernel_op_245 +
              (tmp_kernel_op_137 * tmp_kernel_op_137) * tmp_kernel_op_246 +
              (tmp_kernel_op_139 * tmp_kernel_op_139) * tmp_kernel_op_247 +
              (tmp_kernel_op_141 * tmp_kernel_op_141) * tmp_kernel_op_248 +
              (tmp_kernel_op_143 * tmp_kernel_op_143) * tmp_kernel_op_249 +
              (tmp_kernel_op_145 * tmp_kernel_op_145) * tmp_kernel_op_250 +
              (tmp_kernel_op_147 * tmp_kernel_op_147) * tmp_kernel_op_251 +
              (tmp_kernel_op_149 * tmp_kernel_op_149) * tmp_kernel_op_252 +
              (tmp_kernel_op_151 * tmp_kernel_op_151) * tmp_kernel_op_253 +
              (tmp_kernel_op_153 * tmp_kernel_op_153) * tmp_kernel_op_254;

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
      const walberla::float64 abs_det_jac_affine_BLUE =
          abs(jac_affine_0_0_BLUE * jac_affine_1_1_BLUE -
              jac_affine_0_1_BLUE * jac_affine_1_0_BLUE);
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1;
             ctr_0 += 1) {
#if 0
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

#endif
          const walberla::float64 jac_blending_1_1_q_9 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_9 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_9 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_9 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_9 =
              jac_blending_0_0_q_9 * jac_blending_1_1_q_9 -
              jac_blending_0_1_q_9 * jac_blending_1_0_q_9;
          const walberla::float64 jac_blending_1_1_q_8 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_8 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_8 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_8 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_8 =
              jac_blending_0_0_q_8 * jac_blending_1_1_q_8 -
              jac_blending_0_1_q_8 * jac_blending_1_0_q_8;
          const walberla::float64 jac_blending_1_1_q_7 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_7 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_7 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_7 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_7 =
              jac_blending_0_0_q_7 * jac_blending_1_1_q_7 -
              jac_blending_0_1_q_7 * jac_blending_1_0_q_7;
          const walberla::float64 jac_blending_1_1_q_6 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_6 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_6 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_6 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_6 =
              jac_blending_0_0_q_6 * jac_blending_1_1_q_6 -
              jac_blending_0_1_q_6 * jac_blending_1_0_q_6;
          const walberla::float64 jac_blending_1_1_q_5 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_5 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_5 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_5 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_5 =
              jac_blending_0_0_q_5 * jac_blending_1_1_q_5 -
              jac_blending_0_1_q_5 * jac_blending_1_0_q_5;
          const walberla::float64 jac_blending_1_1_q_4 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_4 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_4 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_4 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_4 =
              jac_blending_0_0_q_4 * jac_blending_1_1_q_4 -
              jac_blending_0_1_q_4 * jac_blending_1_0_q_4;
          const walberla::float64 jac_blending_1_1_q_3 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_3 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_3 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_3 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_2 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_2 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_2 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_1 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_1 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_1 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_0 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_0 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_0 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_1 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_2 = -1.0175839942884126;
          const walberla::float64 tmp_kernel_op_3 =
              -tmp_kernel_op_1 - tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_0 *
                                                    0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_5 = 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_6 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_7 = -0.87062238788797819;
          const walberla::float64 tmp_kernel_op_8 =
              -tmp_kernel_op_6 - tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_1 *
                                                    0.052397656566423402;
          const walberla::float64 tmp_kernel_op_10 = 0.40172877323475986;
          const walberla::float64 tmp_kernel_op_11 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_12 = -0.84634856630183575;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_2 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_15 = 0.64933214716985033;
          const walberla::float64 tmp_kernel_op_16 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_17 = -0.40994976126116967;
          const walberla::float64 tmp_kernel_op_18 =
              -tmp_kernel_op_16 - tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_19 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_3 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_20 = -0.0087919971442063094;
          const walberla::float64 tmp_kernel_op_21 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_22 = 1.5773598151468047;
          const walberla::float64 tmp_kernel_op_23 =
              -tmp_kernel_op_21 - tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_24 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_4 *
                                                     0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_25 = 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_26 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_27 = 0.57645841913372808;
          const walberla::float64 tmp_kernel_op_28 =
              -tmp_kernel_op_26 - tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_29 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_5 *
                                                     0.052397656566423402;
          const walberla::float64 tmp_kernel_op_30 = 0.076825716849082126;
          const walberla::float64 tmp_kernel_op_31 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_32 = -0.19654245353048039;
          const walberla::float64 tmp_kernel_op_33 =
              -tmp_kernel_op_31 - tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_34 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_6 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_35 = 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_36 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_37 = 0.29866429433970043;
          const walberla::float64 tmp_kernel_op_38 =
              -tmp_kernel_op_36 - tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_39 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_7 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_40 = 0.085333161170310645;
          const walberla::float64 tmp_kernel_op_41 = 1.6586673553187572;
          const walberla::float64 tmp_kernel_op_42 = -0.8293336776593786;
          const walberla::float64 tmp_kernel_op_43 =
              -tmp_kernel_op_41 - tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_44 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_8 *
                                                     0.045496761795224737;
          const walberla::float64 tmp_kernel_op_45 = 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_46 = 0.58523781332014213;
          const walberla::float64 tmp_kernel_op_47 = -0.29261890666007107;
          const walberla::float64 tmp_kernel_op_48 =
              -tmp_kernel_op_46 - tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_49 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_9 *
                                                     0.10566414783403316;
          const walberla::float64 tmp_kernel_op_50 =
              tmp_kernel_op_4 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_51 =
              tmp_kernel_op_0 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_52 =
              tmp_kernel_op_50 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_53 =
              tmp_kernel_op_9 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_5 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_53 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_56 =
              tmp_kernel_op_14 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_57 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_58 =
              tmp_kernel_op_56 * tmp_kernel_op_57;
          const walberla::float64 tmp_kernel_op_59 =
              tmp_kernel_op_19 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_15 * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_61 =
              tmp_kernel_op_59 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_24 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_63 =
              tmp_kernel_op_20 * tmp_kernel_op_23;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_62 * tmp_kernel_op_63;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_29 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_25 * tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_65 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_34 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_30 * tmp_kernel_op_33;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_68 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_39 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_72 =
              tmp_kernel_op_35 * tmp_kernel_op_38;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_71 * tmp_kernel_op_72;
          const walberla::float64 tmp_kernel_op_74 =
              tmp_kernel_op_44 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_75 =
              tmp_kernel_op_40 * tmp_kernel_op_43;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_74 * tmp_kernel_op_75;
          const walberla::float64 tmp_kernel_op_77 =
              tmp_kernel_op_49 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_45 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_77 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_12 * tmp_kernel_op_58 +
              tmp_kernel_op_17 * tmp_kernel_op_61 +
              tmp_kernel_op_2 * tmp_kernel_op_52 +
              tmp_kernel_op_22 * tmp_kernel_op_64 +
              tmp_kernel_op_27 * tmp_kernel_op_67 +
              tmp_kernel_op_32 * tmp_kernel_op_70 +
              tmp_kernel_op_37 * tmp_kernel_op_73 +
              tmp_kernel_op_42 * tmp_kernel_op_76 +
              tmp_kernel_op_47 * tmp_kernel_op_79 +
              tmp_kernel_op_55 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_81 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_4 * tmp_kernel_op_51 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_83 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_54 * tmp_kernel_op_9 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_85 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_14 * tmp_kernel_op_57 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_87 = tmp_kernel_op_16 - 1.0;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_19 * tmp_kernel_op_60 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_89 = tmp_kernel_op_21 - 1.0;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_24 * tmp_kernel_op_63 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_91 = tmp_kernel_op_26 - 1.0;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_29 * tmp_kernel_op_66 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_93 = tmp_kernel_op_31 - 1.0;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_34 * tmp_kernel_op_69 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_95 = tmp_kernel_op_36 - 1.0;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_39 * tmp_kernel_op_72 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_97 = tmp_kernel_op_41 - 1.0;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_44 * tmp_kernel_op_75 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_99 = tmp_kernel_op_46 - 1.0;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_49 * tmp_kernel_op_78 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_100 * tmp_kernel_op_99 +
              tmp_kernel_op_81 * tmp_kernel_op_82 +
              tmp_kernel_op_83 * tmp_kernel_op_84 +
              tmp_kernel_op_85 * tmp_kernel_op_86 +
              tmp_kernel_op_87 * tmp_kernel_op_88 +
              tmp_kernel_op_89 * tmp_kernel_op_90 +
              tmp_kernel_op_91 * tmp_kernel_op_92 +
              tmp_kernel_op_93 * tmp_kernel_op_94 +
              tmp_kernel_op_95 * tmp_kernel_op_96 +
              tmp_kernel_op_97 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_102 = -0.035167988576825245;
          const walberla::float64 tmp_kernel_op_103 = 0.25875522422404362;
          const walberla::float64 tmp_kernel_op_104 = 0.30730286739632839;
          const walberla::float64 tmp_kernel_op_105 = 1.1801004774776607;
          const walberla::float64 tmp_kernel_op_106 = 5.1547196302936094;
          const walberla::float64 tmp_kernel_op_107 = 3.1529168382674562;
          const walberla::float64 tmp_kernel_op_108 = 1.6069150929390392;
          const walberla::float64 tmp_kernel_op_109 = 2.5973285886794009;
          const walberla::float64 tmp_kernel_op_110 = 0.3413326446812428;
          const walberla::float64 tmp_kernel_op_111 = 1.4147621866798579;
          const walberla::float64 tmp_kernel_op_112 =
              tmp_kernel_op_100 * tmp_kernel_op_111 +
              tmp_kernel_op_102 * tmp_kernel_op_82 +
              tmp_kernel_op_103 * tmp_kernel_op_84 +
              tmp_kernel_op_104 * tmp_kernel_op_86 +
              tmp_kernel_op_105 * tmp_kernel_op_88 +
              tmp_kernel_op_106 * tmp_kernel_op_90 +
              tmp_kernel_op_107 * tmp_kernel_op_92 +
              tmp_kernel_op_108 * tmp_kernel_op_94 +
              tmp_kernel_op_109 * tmp_kernel_op_96 +
              tmp_kernel_op_110 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_113 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_114 =
              -tmp_kernel_op_102 - tmp_kernel_op_113 + 4.0;
          const walberla::float64 tmp_kernel_op_115 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_116 =
              -tmp_kernel_op_103 - tmp_kernel_op_115 + 4.0;
          const walberla::float64 tmp_kernel_op_117 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_118 =
              -tmp_kernel_op_104 - tmp_kernel_op_117 + 4.0;
          const walberla::float64 tmp_kernel_op_119 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_120 =
              -tmp_kernel_op_105 - tmp_kernel_op_119 + 4.0;
          const walberla::float64 tmp_kernel_op_121 = -1.1195516417167841;
          const walberla::float64 tmp_kernel_op_122 =
              -tmp_kernel_op_106 - tmp_kernel_op_121 + 4.0;
          const walberla::float64 tmp_kernel_op_123 = 0.58832793750850021;
          const walberla::float64 tmp_kernel_op_124 =
              -tmp_kernel_op_107 - tmp_kernel_op_123 + 4.0;
          const walberla::float64 tmp_kernel_op_125 = 2.085782039664632;
          const walberla::float64 tmp_kernel_op_126 =
              -tmp_kernel_op_108 - tmp_kernel_op_125 + 4.0;
          const walberla::float64 tmp_kernel_op_127 = 0.22257093384293847;
          const walberla::float64 tmp_kernel_op_128 =
              -tmp_kernel_op_109 - tmp_kernel_op_127 + 4.0;
          const walberla::float64 tmp_kernel_op_129 = 3.3173347106375144;
          const walberla::float64 tmp_kernel_op_130 =
              -tmp_kernel_op_110 - tmp_kernel_op_129 + 4.0;
          const walberla::float64 tmp_kernel_op_131 = 1.1704756266402843;
          const walberla::float64 tmp_kernel_op_132 =
              -tmp_kernel_op_111 - tmp_kernel_op_131 + 4.0;
          const walberla::float64 tmp_kernel_op_133 =
              tmp_kernel_op_100 * tmp_kernel_op_132 +
              tmp_kernel_op_114 * tmp_kernel_op_82 +
              tmp_kernel_op_116 * tmp_kernel_op_84 +
              tmp_kernel_op_118 * tmp_kernel_op_86 +
              tmp_kernel_op_120 * tmp_kernel_op_88 +
              tmp_kernel_op_122 * tmp_kernel_op_90 +
              tmp_kernel_op_124 * tmp_kernel_op_92 +
              tmp_kernel_op_126 * tmp_kernel_op_94 +
              tmp_kernel_op_128 * tmp_kernel_op_96 +
              tmp_kernel_op_130 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_134 =
              tmp_kernel_op_114 * tmp_kernel_op_52 +
              tmp_kernel_op_116 * tmp_kernel_op_55 +
              tmp_kernel_op_118 * tmp_kernel_op_58 +
              tmp_kernel_op_120 * tmp_kernel_op_61 +
              tmp_kernel_op_122 * tmp_kernel_op_64 +
              tmp_kernel_op_124 * tmp_kernel_op_67 +
              tmp_kernel_op_126 * tmp_kernel_op_70 +
              tmp_kernel_op_128 * tmp_kernel_op_73 +
              tmp_kernel_op_130 * tmp_kernel_op_76 +
              tmp_kernel_op_132 * tmp_kernel_op_79;
          const walberla::float64 tmp_kernel_op_135 = 34.794357504481866;
          const walberla::float64 tmp_kernel_op_136 =
              tmp_kernel_op_50 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_137 = 21.28218865830533;
          const walberla::float64 tmp_kernel_op_138 =
              tmp_kernel_op_53 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_139 = 10.846676877338515;
          const walberla::float64 tmp_kernel_op_140 =
              tmp_kernel_op_56 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_141 = 17.531967973585957;
          const walberla::float64 tmp_kernel_op_142 =
              tmp_kernel_op_59 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_143 = -0.23738392289357257;
          const walberla::float64 tmp_kernel_op_144 =
              tmp_kernel_op_62 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_145 = 1.7465977635122938;
          const walberla::float64 tmp_kernel_op_146 =
              tmp_kernel_op_65 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_147 = 2.0742943549252164;
          const walberla::float64 tmp_kernel_op_148 =
              tmp_kernel_op_68 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_149 = 7.9656782229742085;
          const walberla::float64 tmp_kernel_op_150 =
              tmp_kernel_op_71 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_151 = 2.3039953515983882;
          const walberla::float64 tmp_kernel_op_152 =
              tmp_kernel_op_74 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_153 = 9.5496447600890395;
          const walberla::float64 tmp_kernel_op_154 =
              tmp_kernel_op_77 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_155 =
              tmp_kernel_op_135 * tmp_kernel_op_136 * tmp_kernel_op_51 +
              tmp_kernel_op_137 * tmp_kernel_op_138 * tmp_kernel_op_54 +
              tmp_kernel_op_139 * tmp_kernel_op_140 * tmp_kernel_op_57 +
              tmp_kernel_op_141 * tmp_kernel_op_142 * tmp_kernel_op_60 +
              tmp_kernel_op_143 * tmp_kernel_op_144 * tmp_kernel_op_63 +
              tmp_kernel_op_145 * tmp_kernel_op_146 * tmp_kernel_op_66 +
              tmp_kernel_op_147 * tmp_kernel_op_148 * tmp_kernel_op_69 +
              tmp_kernel_op_149 * tmp_kernel_op_150 * tmp_kernel_op_72 +
              tmp_kernel_op_151 * tmp_kernel_op_152 * tmp_kernel_op_75 +
              tmp_kernel_op_153 * tmp_kernel_op_154 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_156 =
              tmp_kernel_op_4 * 7.7299213783731935e-5;
          const walberla::float64 tmp_kernel_op_157 =
              tmp_kernel_op_9 * 0.0041846416289521935;
          const walberla::float64 tmp_kernel_op_158 =
              tmp_kernel_op_14 * 0.0059021907693753367;
          const walberla::float64 tmp_kernel_op_159 =
              tmp_kernel_op_19 * 0.087039821058937664;
          const walberla::float64 tmp_kernel_op_160 =
              tmp_kernel_op_24 * 1.6606959041833929;
          const walberla::float64 tmp_kernel_op_161 =
              tmp_kernel_op_29 * 0.62130528681440322;
          const walberla::float64 tmp_kernel_op_162 =
              tmp_kernel_op_34 * 0.16138600724470506;
          const walberla::float64 tmp_kernel_op_163 =
              tmp_kernel_op_39 * 0.42163223734820804;
          const walberla::float64 tmp_kernel_op_164 =
              tmp_kernel_op_44 * 0.0072817483953182219;
          const walberla::float64 tmp_kernel_op_165 =
              tmp_kernel_op_49 * 0.12509700280369832;
          const walberla::float64 tmp_kernel_op_166 =
              tmp_kernel_op_136 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_167 =
              tmp_kernel_op_138 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_168 =
              tmp_kernel_op_12 * tmp_kernel_op_140;
          const walberla::float64 tmp_kernel_op_169 =
              tmp_kernel_op_142 * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_170 =
              tmp_kernel_op_144 * tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_171 =
              tmp_kernel_op_146 * tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_172 =
              tmp_kernel_op_148 * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_173 =
              tmp_kernel_op_150 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_174 =
              tmp_kernel_op_152 * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_175 =
              tmp_kernel_op_154 * tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_176 =
              tmp_kernel_op_166 * tmp_kernel_op_81 +
              tmp_kernel_op_167 * tmp_kernel_op_83 +
              tmp_kernel_op_168 * tmp_kernel_op_85 +
              tmp_kernel_op_169 * tmp_kernel_op_87 +
              tmp_kernel_op_170 * tmp_kernel_op_89 +
              tmp_kernel_op_171 * tmp_kernel_op_91 +
              tmp_kernel_op_172 * tmp_kernel_op_93 +
              tmp_kernel_op_173 * tmp_kernel_op_95 +
              tmp_kernel_op_174 * tmp_kernel_op_97 +
              tmp_kernel_op_175 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_177 =
              tmp_kernel_op_156 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_178 =
              tmp_kernel_op_157 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_179 =
              tmp_kernel_op_12 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_180 =
              tmp_kernel_op_159 * tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_181 =
              tmp_kernel_op_160 * tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_182 =
              tmp_kernel_op_161 * tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_183 =
              tmp_kernel_op_162 * tmp_kernel_op_32;
          const walberla::float64 tmp_kernel_op_184 =
              tmp_kernel_op_163 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_185 =
              tmp_kernel_op_164 * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_186 =
              tmp_kernel_op_165 * tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_187 =
              tmp_kernel_op_113 * tmp_kernel_op_177 +
              tmp_kernel_op_115 * tmp_kernel_op_178 +
              tmp_kernel_op_117 * tmp_kernel_op_179 +
              tmp_kernel_op_119 * tmp_kernel_op_180 +
              tmp_kernel_op_121 * tmp_kernel_op_181 +
              tmp_kernel_op_123 * tmp_kernel_op_182 +
              tmp_kernel_op_125 * tmp_kernel_op_183 +
              tmp_kernel_op_127 * tmp_kernel_op_184 +
              tmp_kernel_op_129 * tmp_kernel_op_185 +
              tmp_kernel_op_131 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_188 =
              tmp_kernel_op_114 * tmp_kernel_op_166 +
              tmp_kernel_op_116 * tmp_kernel_op_167 +
              tmp_kernel_op_118 * tmp_kernel_op_168 +
              tmp_kernel_op_120 * tmp_kernel_op_169 +
              tmp_kernel_op_122 * tmp_kernel_op_170 +
              tmp_kernel_op_124 * tmp_kernel_op_171 +
              tmp_kernel_op_126 * tmp_kernel_op_172 +
              tmp_kernel_op_128 * tmp_kernel_op_173 +
              tmp_kernel_op_130 * tmp_kernel_op_174 +
              tmp_kernel_op_132 * tmp_kernel_op_175;
          const walberla::float64 tmp_kernel_op_189 =
              tmp_kernel_op_114 * tmp_kernel_op_177 +
              tmp_kernel_op_116 * tmp_kernel_op_178 +
              tmp_kernel_op_118 * tmp_kernel_op_179 +
              tmp_kernel_op_120 * tmp_kernel_op_180 +
              tmp_kernel_op_122 * tmp_kernel_op_181 +
              tmp_kernel_op_124 * tmp_kernel_op_182 +
              tmp_kernel_op_126 * tmp_kernel_op_183 +
              tmp_kernel_op_128 * tmp_kernel_op_184 +
              tmp_kernel_op_130 * tmp_kernel_op_185 +
              tmp_kernel_op_132 * tmp_kernel_op_186;
          const walberla::float64 tmp_kernel_op_190 =
              tmp_kernel_op_135 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_191 =
              tmp_kernel_op_137 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_192 =
              tmp_kernel_op_139 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_193 =
              tmp_kernel_op_141 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_194 =
              tmp_kernel_op_143 * -0.27988791042919603;
          const walberla::float64 tmp_kernel_op_195 =
              tmp_kernel_op_145 * 0.14708198437712505;
          const walberla::float64 tmp_kernel_op_196 =
              tmp_kernel_op_147 * 0.52144550991615801;
          const walberla::float64 tmp_kernel_op_197 =
              tmp_kernel_op_149 * 0.055642733460734617;
          const walberla::float64 tmp_kernel_op_198 =
              tmp_kernel_op_151 * 0.8293336776593786;
          const walberla::float64 tmp_kernel_op_199 =
              tmp_kernel_op_153 * 0.29261890666007107;
          const walberla::float64 tmp_kernel_op_200 =
              tmp_kernel_op_177 * tmp_kernel_op_190 +
              tmp_kernel_op_178 * tmp_kernel_op_191 +
              tmp_kernel_op_179 * tmp_kernel_op_192 +
              tmp_kernel_op_180 * tmp_kernel_op_193 +
              tmp_kernel_op_181 * tmp_kernel_op_194 +
              tmp_kernel_op_182 * tmp_kernel_op_195 +
              tmp_kernel_op_183 * tmp_kernel_op_196 +
              tmp_kernel_op_184 * tmp_kernel_op_197 +
              tmp_kernel_op_185 * tmp_kernel_op_198 +
              tmp_kernel_op_186 * tmp_kernel_op_199;
          const walberla::float64 tmp_kernel_op_201 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_202 =
              tmp_kernel_op_201 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_203 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_204 =
              tmp_kernel_op_203 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_205 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_206 =
              tmp_kernel_op_14 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_207 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_208 =
              tmp_kernel_op_19 * tmp_kernel_op_207;
          const walberla::float64 tmp_kernel_op_209 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_210 =
              tmp_kernel_op_209 * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_211 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_212 =
              tmp_kernel_op_211 * tmp_kernel_op_29;
          const walberla::float64 tmp_kernel_op_213 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_214 =
              tmp_kernel_op_213 * tmp_kernel_op_34;
          const walberla::float64 tmp_kernel_op_215 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_216 =
              tmp_kernel_op_215 * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_217 = 0.68779434890003011;
          const walberla::float64 tmp_kernel_op_218 =
              tmp_kernel_op_217 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_219 = 0.085625824534935377;
          const walberla::float64 tmp_kernel_op_220 =
              tmp_kernel_op_219 * tmp_kernel_op_49;
          const walberla::float64 tmp_kernel_op_221 =
              tmp_kernel_op_202 * tmp_kernel_op_81;
          const walberla::float64 tmp_kernel_op_222 =
              tmp_kernel_op_204 * tmp_kernel_op_83;
          const walberla::float64 tmp_kernel_op_223 =
              tmp_kernel_op_206 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_224 =
              tmp_kernel_op_208 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_225 =
              tmp_kernel_op_210 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_226 =
              tmp_kernel_op_212 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_227 =
              tmp_kernel_op_214 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_228 =
              tmp_kernel_op_216 * tmp_kernel_op_95;
          const walberla::float64 tmp_kernel_op_229 =
              tmp_kernel_op_218 * tmp_kernel_op_97;
          const walberla::float64 tmp_kernel_op_230 =
              tmp_kernel_op_220 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_231 =
              tmp_kernel_op_102 * tmp_kernel_op_221 +
              tmp_kernel_op_103 * tmp_kernel_op_222 +
              tmp_kernel_op_104 * tmp_kernel_op_223 +
              tmp_kernel_op_105 * tmp_kernel_op_224 +
              tmp_kernel_op_106 * tmp_kernel_op_225 +
              tmp_kernel_op_107 * tmp_kernel_op_226 +
              tmp_kernel_op_108 * tmp_kernel_op_227 +
              tmp_kernel_op_109 * tmp_kernel_op_228 +
              tmp_kernel_op_110 * tmp_kernel_op_229 +
              tmp_kernel_op_111 * tmp_kernel_op_230;
          const walberla::float64 tmp_kernel_op_232 =
              tmp_kernel_op_114 * tmp_kernel_op_221 +
              tmp_kernel_op_116 * tmp_kernel_op_222 +
              tmp_kernel_op_118 * tmp_kernel_op_223 +
              tmp_kernel_op_120 * tmp_kernel_op_224 +
              tmp_kernel_op_122 * tmp_kernel_op_225 +
              tmp_kernel_op_124 * tmp_kernel_op_226 +
              tmp_kernel_op_126 * tmp_kernel_op_227 +
              tmp_kernel_op_128 * tmp_kernel_op_228 +
              tmp_kernel_op_130 * tmp_kernel_op_229 +
              tmp_kernel_op_132 * tmp_kernel_op_230;
          const walberla::float64 tmp_kernel_op_233 =
              tmp_kernel_op_114 * tmp_kernel_op_136 * tmp_kernel_op_81 +
              tmp_kernel_op_116 * tmp_kernel_op_138 * tmp_kernel_op_83 +
              tmp_kernel_op_118 * tmp_kernel_op_140 * tmp_kernel_op_85 +
              tmp_kernel_op_120 * tmp_kernel_op_142 * tmp_kernel_op_87 +
              tmp_kernel_op_122 * tmp_kernel_op_144 * tmp_kernel_op_89 +
              tmp_kernel_op_124 * tmp_kernel_op_146 * tmp_kernel_op_91 +
              tmp_kernel_op_126 * tmp_kernel_op_148 * tmp_kernel_op_93 +
              tmp_kernel_op_128 * tmp_kernel_op_150 * tmp_kernel_op_95 +
              tmp_kernel_op_130 * tmp_kernel_op_152 * tmp_kernel_op_97 +
              tmp_kernel_op_132 * tmp_kernel_op_154 * tmp_kernel_op_99;
          const walberla::float64 tmp_kernel_op_234 =
              tmp_kernel_op_135 * -0.0087919971442063111;
          const walberla::float64 tmp_kernel_op_235 =
              tmp_kernel_op_137 * 0.064688806056010906;
          const walberla::float64 tmp_kernel_op_236 =
              tmp_kernel_op_139 * 0.076825716849082099;
          const walberla::float64 tmp_kernel_op_237 =
              tmp_kernel_op_141 * 0.29502511936941517;
          const walberla::float64 tmp_kernel_op_238 =
              tmp_kernel_op_143 * 1.2886799075734023;
          const walberla::float64 tmp_kernel_op_239 =
              tmp_kernel_op_145 * 0.78822920956686404;
          const walberla::float64 tmp_kernel_op_240 =
              tmp_kernel_op_147 * 0.40172877323475981;
          const walberla::float64 tmp_kernel_op_241 =
              tmp_kernel_op_149 * 0.64933214716985022;
          const walberla::float64 tmp_kernel_op_242 =
              tmp_kernel_op_151 * 0.085333161170310701;
          const walberla::float64 tmp_kernel_op_243 =
              tmp_kernel_op_153 * 0.35369054666996447;
          const walberla::float64 tmp_kernel_op_244 =
              tmp_kernel_op_221 * tmp_kernel_op_234 +
              tmp_kernel_op_222 * tmp_kernel_op_235 +
              tmp_kernel_op_223 * tmp_kernel_op_236 +
              tmp_kernel_op_224 * tmp_kernel_op_237 +
              tmp_kernel_op_225 * tmp_kernel_op_238 +
              tmp_kernel_op_226 * tmp_kernel_op_239 +
              tmp_kernel_op_227 * tmp_kernel_op_240 +
              tmp_kernel_op_228 * tmp_kernel_op_241 +
              tmp_kernel_op_229 * tmp_kernel_op_242 +
              tmp_kernel_op_230 * tmp_kernel_op_243;
          const walberla::float64 tmp_kernel_op_245 =
              tmp_kernel_op_156 * tmp_kernel_op_201;
          const walberla::float64 tmp_kernel_op_246 =
              tmp_kernel_op_157 * tmp_kernel_op_203;
          const walberla::float64 tmp_kernel_op_247 =
              tmp_kernel_op_158 * tmp_kernel_op_205;
          const walberla::float64 tmp_kernel_op_248 =
              tmp_kernel_op_159 * tmp_kernel_op_207;
          const walberla::float64 tmp_kernel_op_249 =
              tmp_kernel_op_160 * tmp_kernel_op_209;
          const walberla::float64 tmp_kernel_op_250 =
              tmp_kernel_op_161 * tmp_kernel_op_211;
          const walberla::float64 tmp_kernel_op_251 =
              tmp_kernel_op_162 * tmp_kernel_op_213;
          const walberla::float64 tmp_kernel_op_252 =
              tmp_kernel_op_163 * tmp_kernel_op_215;
          const walberla::float64 tmp_kernel_op_253 =
              tmp_kernel_op_164 * tmp_kernel_op_217;
          const walberla::float64 tmp_kernel_op_254 =
              tmp_kernel_op_165 * tmp_kernel_op_219;
          const walberla::float64 tmp_kernel_op_255 =
              tmp_kernel_op_114 * tmp_kernel_op_202;
          const walberla::float64 tmp_kernel_op_256 =
              tmp_kernel_op_116 * tmp_kernel_op_204;
          const walberla::float64 tmp_kernel_op_257 =
              tmp_kernel_op_118 * tmp_kernel_op_206;
          const walberla::float64 tmp_kernel_op_258 =
              tmp_kernel_op_120 * tmp_kernel_op_208;
          const walberla::float64 tmp_kernel_op_259 =
              tmp_kernel_op_122 * tmp_kernel_op_210;
          const walberla::float64 tmp_kernel_op_260 =
              tmp_kernel_op_124 * tmp_kernel_op_212;
          const walberla::float64 tmp_kernel_op_261 =
              tmp_kernel_op_126 * tmp_kernel_op_214;
          const walberla::float64 tmp_kernel_op_262 =
              tmp_kernel_op_128 * tmp_kernel_op_216;
          const walberla::float64 tmp_kernel_op_263 =
              tmp_kernel_op_130 * tmp_kernel_op_218;
          const walberla::float64 tmp_kernel_op_264 =
              tmp_kernel_op_132 * tmp_kernel_op_220;
          const walberla::float64 tmp_kernel_op_265 =
              tmp_kernel_op_102 * tmp_kernel_op_255 +
              tmp_kernel_op_103 * tmp_kernel_op_256 +
              tmp_kernel_op_104 * tmp_kernel_op_257 +
              tmp_kernel_op_105 * tmp_kernel_op_258 +
              tmp_kernel_op_106 * tmp_kernel_op_259 +
              tmp_kernel_op_107 * tmp_kernel_op_260 +
              tmp_kernel_op_108 * tmp_kernel_op_261 +
              tmp_kernel_op_109 * tmp_kernel_op_262 +
              tmp_kernel_op_110 * tmp_kernel_op_263 +
              tmp_kernel_op_111 * tmp_kernel_op_264;
          const walberla::float64 tmp_kernel_op_266 =
              tmp_kernel_op_114 * tmp_kernel_op_156;
          const walberla::float64 tmp_kernel_op_267 =
              tmp_kernel_op_116 * tmp_kernel_op_157;
          const walberla::float64 tmp_kernel_op_268 =
              tmp_kernel_op_118 * tmp_kernel_op_158;
          const walberla::float64 tmp_kernel_op_269 =
              tmp_kernel_op_120 * tmp_kernel_op_159;
          const walberla::float64 tmp_kernel_op_270 =
              tmp_kernel_op_122 * tmp_kernel_op_160;
          const walberla::float64 tmp_kernel_op_271 =
              tmp_kernel_op_124 * tmp_kernel_op_161;
          const walberla::float64 tmp_kernel_op_272 =
              tmp_kernel_op_126 * tmp_kernel_op_162;
          const walberla::float64 tmp_kernel_op_273 =
              tmp_kernel_op_128 * tmp_kernel_op_163;
          const walberla::float64 tmp_kernel_op_274 =
              tmp_kernel_op_130 * tmp_kernel_op_164;
          const walberla::float64 tmp_kernel_op_275 =
              tmp_kernel_op_132 * tmp_kernel_op_165;
          const walberla::float64 tmp_kernel_op_276 =
              tmp_kernel_op_113 * tmp_kernel_op_266 +
              tmp_kernel_op_115 * tmp_kernel_op_267 +
              tmp_kernel_op_117 * tmp_kernel_op_268 +
              tmp_kernel_op_119 * tmp_kernel_op_269 +
              tmp_kernel_op_121 * tmp_kernel_op_270 +
              tmp_kernel_op_123 * tmp_kernel_op_271 +
              tmp_kernel_op_125 * tmp_kernel_op_272 +
              tmp_kernel_op_127 * tmp_kernel_op_273 +
              tmp_kernel_op_129 * tmp_kernel_op_274 +
              tmp_kernel_op_131 * tmp_kernel_op_275;
          const walberla::float64 tmp_kernel_op_277 =
              tmp_kernel_op_135 * tmp_kernel_op_245 * 4.0 +
              tmp_kernel_op_137 * tmp_kernel_op_246 * 4.0 +
              tmp_kernel_op_139 * tmp_kernel_op_247 * 4.0 +
              tmp_kernel_op_141 * tmp_kernel_op_248 * 4.0 +
              tmp_kernel_op_143 * tmp_kernel_op_249 * 4.0 +
              tmp_kernel_op_145 * tmp_kernel_op_250 * 4.0 +
              tmp_kernel_op_147 * tmp_kernel_op_251 * 4.0 +
              tmp_kernel_op_149 * tmp_kernel_op_252 * 4.0 +
              tmp_kernel_op_151 * tmp_kernel_op_253 * 4.0 +
              tmp_kernel_op_153 * tmp_kernel_op_254 * 4.0;
          const walberla::float64 tmp_kernel_op_278 =
              (tmp_kernel_op_114 * tmp_kernel_op_114);
          const walberla::float64 tmp_kernel_op_279 =
              (tmp_kernel_op_116 * tmp_kernel_op_116);
          const walberla::float64 tmp_kernel_op_280 =
              (tmp_kernel_op_118 * tmp_kernel_op_118);
          const walberla::float64 tmp_kernel_op_281 =
              (tmp_kernel_op_120 * tmp_kernel_op_120);
          const walberla::float64 tmp_kernel_op_282 =
              (tmp_kernel_op_122 * tmp_kernel_op_122);
          const walberla::float64 tmp_kernel_op_283 =
              (tmp_kernel_op_124 * tmp_kernel_op_124);
          const walberla::float64 tmp_kernel_op_284 =
              (tmp_kernel_op_126 * tmp_kernel_op_126);
          const walberla::float64 tmp_kernel_op_285 =
              (tmp_kernel_op_128 * tmp_kernel_op_128);
          const walberla::float64 tmp_kernel_op_286 =
              (tmp_kernel_op_130 * tmp_kernel_op_130);
          const walberla::float64 tmp_kernel_op_287 =
              (tmp_kernel_op_132 * tmp_kernel_op_132);
          const walberla::float64 tmp_kernel_op_288 =
              tmp_kernel_op_136 * tmp_kernel_op_278 +
              tmp_kernel_op_138 * tmp_kernel_op_279 +
              tmp_kernel_op_140 * tmp_kernel_op_280 +
              tmp_kernel_op_142 * tmp_kernel_op_281 +
              tmp_kernel_op_144 * tmp_kernel_op_282 +
              tmp_kernel_op_146 * tmp_kernel_op_283 +
              tmp_kernel_op_148 * tmp_kernel_op_284 +
              tmp_kernel_op_150 * tmp_kernel_op_285 +
              tmp_kernel_op_152 * tmp_kernel_op_286 +
              tmp_kernel_op_154 * tmp_kernel_op_287;
          const walberla::float64 tmp_kernel_op_289 =
              tmp_kernel_op_234 * tmp_kernel_op_255 +
              tmp_kernel_op_235 * tmp_kernel_op_256 +
              tmp_kernel_op_236 * tmp_kernel_op_257 +
              tmp_kernel_op_237 * tmp_kernel_op_258 +
              tmp_kernel_op_238 * tmp_kernel_op_259 +
              tmp_kernel_op_239 * tmp_kernel_op_260 +
              tmp_kernel_op_240 * tmp_kernel_op_261 +
              tmp_kernel_op_241 * tmp_kernel_op_262 +
              tmp_kernel_op_242 * tmp_kernel_op_263 +
              tmp_kernel_op_243 * tmp_kernel_op_264;
          const walberla::float64 tmp_kernel_op_290 =
              tmp_kernel_op_190 * tmp_kernel_op_266 +
              tmp_kernel_op_191 * tmp_kernel_op_267 +
              tmp_kernel_op_192 * tmp_kernel_op_268 +
              tmp_kernel_op_193 * tmp_kernel_op_269 +
              tmp_kernel_op_194 * tmp_kernel_op_270 +
              tmp_kernel_op_195 * tmp_kernel_op_271 +
              tmp_kernel_op_196 * tmp_kernel_op_272 +
              tmp_kernel_op_197 * tmp_kernel_op_273 +
              tmp_kernel_op_198 * tmp_kernel_op_274 +
              tmp_kernel_op_199 * tmp_kernel_op_275;
          const walberla::float64 elMat_0_0 =
              (tmp_kernel_op_0 * tmp_kernel_op_0) *
                  (tmp_kernel_op_3 * tmp_kernel_op_3) * tmp_kernel_op_4 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) *
                  (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_14 +
              (tmp_kernel_op_15 * tmp_kernel_op_15) *
                  (tmp_kernel_op_18 * tmp_kernel_op_18) * tmp_kernel_op_19 +
              (tmp_kernel_op_20 * tmp_kernel_op_20) *
                  (tmp_kernel_op_23 * tmp_kernel_op_23) * tmp_kernel_op_24 +
              (tmp_kernel_op_25 * tmp_kernel_op_25) *
                  (tmp_kernel_op_28 * tmp_kernel_op_28) * tmp_kernel_op_29 +
              (tmp_kernel_op_30 * tmp_kernel_op_30) *
                  (tmp_kernel_op_33 * tmp_kernel_op_33) * tmp_kernel_op_34 +
              (tmp_kernel_op_35 * tmp_kernel_op_35) *
                  (tmp_kernel_op_38 * tmp_kernel_op_38) * tmp_kernel_op_39 +
              (tmp_kernel_op_40 * tmp_kernel_op_40) *
                  (tmp_kernel_op_43 * tmp_kernel_op_43) * tmp_kernel_op_44 +
              (tmp_kernel_op_45 * tmp_kernel_op_45) *
                  (tmp_kernel_op_48 * tmp_kernel_op_48) * tmp_kernel_op_49 +
              (tmp_kernel_op_5 * tmp_kernel_op_5) *
                  (tmp_kernel_op_8 * tmp_kernel_op_8) * tmp_kernel_op_9;
          const walberla::float64 elMat_0_1 = tmp_kernel_op_80;
          const walberla::float64 elMat_0_2 = tmp_kernel_op_101;
          const walberla::float64 elMat_0_3 = tmp_kernel_op_112;
          const walberla::float64 elMat_0_4 = tmp_kernel_op_133;
          const walberla::float64 elMat_0_5 = tmp_kernel_op_134;
          const walberla::float64 elMat_0_6 = tmp_kernel_op_155;
          const walberla::float64 elMat_1_0 = tmp_kernel_op_80;
          const walberla::float64 elMat_1_1 =
              (tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_158 +
              tmp_kernel_op_156 * (tmp_kernel_op_2 * tmp_kernel_op_2) +
              tmp_kernel_op_157 * (tmp_kernel_op_7 * tmp_kernel_op_7) +
              tmp_kernel_op_159 * (tmp_kernel_op_17 * tmp_kernel_op_17) +
              tmp_kernel_op_160 * (tmp_kernel_op_22 * tmp_kernel_op_22) +
              tmp_kernel_op_161 * (tmp_kernel_op_27 * tmp_kernel_op_27) +
              tmp_kernel_op_162 * (tmp_kernel_op_32 * tmp_kernel_op_32) +
              tmp_kernel_op_163 * (tmp_kernel_op_37 * tmp_kernel_op_37) +
              tmp_kernel_op_164 * (tmp_kernel_op_42 * tmp_kernel_op_42) +
              tmp_kernel_op_165 * (tmp_kernel_op_47 * tmp_kernel_op_47);
          const walberla::float64 elMat_1_2 = tmp_kernel_op_176;
          const walberla::float64 elMat_1_3 = tmp_kernel_op_187;
          const walberla::float64 elMat_1_4 = tmp_kernel_op_188;
          const walberla::float64 elMat_1_5 = tmp_kernel_op_189;
          const walberla::float64 elMat_1_6 = tmp_kernel_op_200;
          const walberla::float64 elMat_2_0 = tmp_kernel_op_101;
          const walberla::float64 elMat_2_1 = tmp_kernel_op_176;
          const walberla::float64 elMat_2_2 =
              tmp_kernel_op_202 * (tmp_kernel_op_81 * tmp_kernel_op_81) +
              tmp_kernel_op_204 * (tmp_kernel_op_83 * tmp_kernel_op_83) +
              tmp_kernel_op_206 * (tmp_kernel_op_85 * tmp_kernel_op_85) +
              tmp_kernel_op_208 * (tmp_kernel_op_87 * tmp_kernel_op_87) +
              tmp_kernel_op_210 * (tmp_kernel_op_89 * tmp_kernel_op_89) +
              tmp_kernel_op_212 * (tmp_kernel_op_91 * tmp_kernel_op_91) +
              tmp_kernel_op_214 * (tmp_kernel_op_93 * tmp_kernel_op_93) +
              tmp_kernel_op_216 * (tmp_kernel_op_95 * tmp_kernel_op_95) +
              tmp_kernel_op_218 * (tmp_kernel_op_97 * tmp_kernel_op_97) +
              tmp_kernel_op_220 * (tmp_kernel_op_99 * tmp_kernel_op_99);
          const walberla::float64 elMat_2_3 = tmp_kernel_op_231;
          const walberla::float64 elMat_2_4 = tmp_kernel_op_232;
          const walberla::float64 elMat_2_5 = tmp_kernel_op_233;
          const walberla::float64 elMat_2_6 = tmp_kernel_op_244;
          const walberla::float64 elMat_3_0 = tmp_kernel_op_112;
          const walberla::float64 elMat_3_1 = tmp_kernel_op_187;
          const walberla::float64 elMat_3_2 = tmp_kernel_op_231;
          const walberla::float64 elMat_3_3 =
              tmp_kernel_op_245 * 16.0 + tmp_kernel_op_246 * 16.0 +
              tmp_kernel_op_247 * 16.0 + tmp_kernel_op_248 * 16.0 +
              tmp_kernel_op_249 * 16.0 + tmp_kernel_op_250 * 16.0 +
              tmp_kernel_op_251 * 16.0 + tmp_kernel_op_252 * 16.0 +
              tmp_kernel_op_253 * 16.0 + tmp_kernel_op_254 * 16.0;
          const walberla::float64 elMat_3_4 = tmp_kernel_op_265;
          const walberla::float64 elMat_3_5 = tmp_kernel_op_276;
          const walberla::float64 elMat_3_6 = tmp_kernel_op_277;
          const walberla::float64 elMat_4_0 = tmp_kernel_op_133;
          const walberla::float64 elMat_4_1 = tmp_kernel_op_188;
          const walberla::float64 elMat_4_2 = tmp_kernel_op_232;
          const walberla::float64 elMat_4_3 = tmp_kernel_op_265;
          const walberla::float64 elMat_4_4 =
              tmp_kernel_op_202 * tmp_kernel_op_278 +
              tmp_kernel_op_204 * tmp_kernel_op_279 +
              tmp_kernel_op_206 * tmp_kernel_op_280 +
              tmp_kernel_op_208 * tmp_kernel_op_281 +
              tmp_kernel_op_210 * tmp_kernel_op_282 +
              tmp_kernel_op_212 * tmp_kernel_op_283 +
              tmp_kernel_op_214 * tmp_kernel_op_284 +
              tmp_kernel_op_216 * tmp_kernel_op_285 +
              tmp_kernel_op_218 * tmp_kernel_op_286 +
              tmp_kernel_op_220 * tmp_kernel_op_287;
          const walberla::float64 elMat_4_5 = tmp_kernel_op_288;
          const walberla::float64 elMat_4_6 = tmp_kernel_op_289;
          const walberla::float64 elMat_5_0 = tmp_kernel_op_134;
          const walberla::float64 elMat_5_1 = tmp_kernel_op_189;
          const walberla::float64 elMat_5_2 = tmp_kernel_op_233;
          const walberla::float64 elMat_5_3 = tmp_kernel_op_276;
          const walberla::float64 elMat_5_4 = tmp_kernel_op_288;
          const walberla::float64 elMat_5_5 =
              tmp_kernel_op_156 * tmp_kernel_op_278 +
              tmp_kernel_op_157 * tmp_kernel_op_279 +
              tmp_kernel_op_158 * tmp_kernel_op_280 +
              tmp_kernel_op_159 * tmp_kernel_op_281 +
              tmp_kernel_op_160 * tmp_kernel_op_282 +
              tmp_kernel_op_161 * tmp_kernel_op_283 +
              tmp_kernel_op_162 * tmp_kernel_op_284 +
              tmp_kernel_op_163 * tmp_kernel_op_285 +
              tmp_kernel_op_164 * tmp_kernel_op_286 +
              tmp_kernel_op_165 * tmp_kernel_op_287;
          const walberla::float64 elMat_5_6 = tmp_kernel_op_290;
          const walberla::float64 elMat_6_0 = tmp_kernel_op_155;
          const walberla::float64 elMat_6_1 = tmp_kernel_op_200;
          const walberla::float64 elMat_6_2 = tmp_kernel_op_244;
          const walberla::float64 elMat_6_3 = tmp_kernel_op_277;
          const walberla::float64 elMat_6_4 = tmp_kernel_op_289;
          const walberla::float64 elMat_6_5 = tmp_kernel_op_290;
          const walberla::float64 elMat_6_6 =
              (tmp_kernel_op_135 * tmp_kernel_op_135) * tmp_kernel_op_245 +
              (tmp_kernel_op_137 * tmp_kernel_op_137) * tmp_kernel_op_246 +
              (tmp_kernel_op_139 * tmp_kernel_op_139) * tmp_kernel_op_247 +
              (tmp_kernel_op_141 * tmp_kernel_op_141) * tmp_kernel_op_248 +
              (tmp_kernel_op_143 * tmp_kernel_op_143) * tmp_kernel_op_249 +
              (tmp_kernel_op_145 * tmp_kernel_op_145) * tmp_kernel_op_250 +
              (tmp_kernel_op_147 * tmp_kernel_op_147) * tmp_kernel_op_251 +
              (tmp_kernel_op_149 * tmp_kernel_op_149) * tmp_kernel_op_252 +
              (tmp_kernel_op_151 * tmp_kernel_op_151) * tmp_kernel_op_253 +
              (tmp_kernel_op_153 * tmp_kernel_op_153) * tmp_kernel_op_254;

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
void P2PlusBubbleElementwiseMass_AffineMap2D_float64::
    computeInverseDiagonalOperatorValues_P2PlusBubbleElementwiseMass_AffineMap2D_float64_macro_2D(
        walberla::float64 *RESTRICT _data_invDiag_,
        walberla::float64 *RESTRICT _data_invDiag_Edge,
        walberla::float64 *RESTRICT _data_invDiag_Vertex,
        walberla::float64 bMat_00, walberla::float64 bMat_01,
        walberla::float64 bMat_10, walberla::float64 bMat_11,
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
      const walberla::float64 abs_det_jac_affine_GRAY =
          abs(jac_affine_0_0_GRAY * jac_affine_1_1_GRAY -
              jac_affine_0_1_GRAY * jac_affine_1_0_GRAY);
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge;
             ctr_0 += 1) {
#if 0
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
#endif
          const walberla::float64 jac_blending_1_1_q_9 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_9 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_9 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_9 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_9 =
              jac_blending_0_0_q_9 * jac_blending_1_1_q_9 -
              jac_blending_0_1_q_9 * jac_blending_1_0_q_9;
          const walberla::float64 jac_blending_1_1_q_8 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_8 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_8 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_8 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_8 =
              jac_blending_0_0_q_8 * jac_blending_1_1_q_8 -
              jac_blending_0_1_q_8 * jac_blending_1_0_q_8;
          const walberla::float64 jac_blending_1_1_q_7 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_7 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_7 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_7 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_7 =
              jac_blending_0_0_q_7 * jac_blending_1_1_q_7 -
              jac_blending_0_1_q_7 * jac_blending_1_0_q_7;
          const walberla::float64 jac_blending_1_1_q_6 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_6 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_6 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_6 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_6 =
              jac_blending_0_0_q_6 * jac_blending_1_1_q_6 -
              jac_blending_0_1_q_6 * jac_blending_1_0_q_6;
          const walberla::float64 jac_blending_1_1_q_5 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_5 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_5 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_5 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_5 =
              jac_blending_0_0_q_5 * jac_blending_1_1_q_5 -
              jac_blending_0_1_q_5 * jac_blending_1_0_q_5;
          const walberla::float64 jac_blending_1_1_q_4 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_4 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_4 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_4 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_4 =
              jac_blending_0_0_q_4 * jac_blending_1_1_q_4 -
              jac_blending_0_1_q_4 * jac_blending_1_0_q_4;
          const walberla::float64 jac_blending_1_1_q_3 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_3 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_3 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_3 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_2 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_2 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_2 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_1 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_1 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_1 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_0 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_0 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_0 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_1 = -1.0175839942884126;
          const walberla::float64 tmp_kernel_op_2 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_0 *
                                                    0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_3 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_4 = -0.87062238788797819;
          const walberla::float64 tmp_kernel_op_5 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_1 *
                                                    0.052397656566423402;
          const walberla::float64 tmp_kernel_op_6 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_7 = -0.84634856630183575;
          const walberla::float64 tmp_kernel_op_8 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_2 *
                                                    0.074275554650521658;
          const walberla::float64 tmp_kernel_op_9 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_10 = -0.40994976126116967;
          const walberla::float64 tmp_kernel_op_11 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_3 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_12 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_13 = 1.5773598151468047;
          const walberla::float64 tmp_kernel_op_14 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_4 *
                                                     0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_15 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_16 = 0.57645841913372808;
          const walberla::float64 tmp_kernel_op_17 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_5 *
                                                     0.052397656566423402;
          const walberla::float64 tmp_kernel_op_18 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_19 = -0.19654245353048039;
          const walberla::float64 tmp_kernel_op_20 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_6 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_21 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_22 = 0.29866429433970043;
          const walberla::float64 tmp_kernel_op_23 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_7 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_24 = 1.6586673553187572;
          const walberla::float64 tmp_kernel_op_25 = -0.8293336776593786;
          const walberla::float64 tmp_kernel_op_26 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_8 *
                                                     0.045496761795224737;
          const walberla::float64 tmp_kernel_op_27 = 0.58523781332014213;
          const walberla::float64 tmp_kernel_op_28 = -0.29261890666007107;
          const walberla::float64 tmp_kernel_op_29 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_9 *
                                                     0.10566414783403316;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_2 * 7.7299213783731935e-5;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_5 * 0.0041846416289521935;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_8 * 0.0059021907693753367;
          const walberla::float64 tmp_kernel_op_33 =
              tmp_kernel_op_11 * 0.087039821058937664;
          const walberla::float64 tmp_kernel_op_34 =
              tmp_kernel_op_14 * 1.6606959041833929;
          const walberla::float64 tmp_kernel_op_35 =
              tmp_kernel_op_17 * 0.62130528681440322;
          const walberla::float64 tmp_kernel_op_36 =
              tmp_kernel_op_20 * 0.16138600724470506;
          const walberla::float64 tmp_kernel_op_37 =
              tmp_kernel_op_23 * 0.42163223734820804;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_26 * 0.0072817483953182219;
          const walberla::float64 tmp_kernel_op_39 =
              tmp_kernel_op_29 * 0.12509700280369832;
          const walberla::float64 tmp_kernel_op_40 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_41 =
              tmp_kernel_op_2 * tmp_kernel_op_40;
          const walberla::float64 tmp_kernel_op_42 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_43 =
              tmp_kernel_op_42 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_44 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_45 =
              tmp_kernel_op_44 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_46 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_47 =
              tmp_kernel_op_11 * tmp_kernel_op_46;
          const walberla::float64 tmp_kernel_op_48 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_49 =
              tmp_kernel_op_14 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_50 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_51 =
              tmp_kernel_op_17 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_52 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_53 =
              tmp_kernel_op_20 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_54 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_23 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_56 = 0.68779434890003011;
          const walberla::float64 tmp_kernel_op_57 =
              tmp_kernel_op_26 * tmp_kernel_op_56;
          const walberla::float64 tmp_kernel_op_58 = 0.085625824534935377;
          const walberla::float64 tmp_kernel_op_59 =
              tmp_kernel_op_29 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_30 * tmp_kernel_op_40;
          const walberla::float64 tmp_kernel_op_61 =
              tmp_kernel_op_31 * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_32 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_63 =
              tmp_kernel_op_33 * tmp_kernel_op_46;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_34 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_35 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_36 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_37 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_38 * tmp_kernel_op_56;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_39 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_70 = 26.571134466934286;
          const walberla::float64 tmp_kernel_op_71 = 9.9408845890304516;
          const walberla::float64 tmp_kernel_op_72 = 2.5821761159152818;
          const walberla::float64 tmp_kernel_op_73 = 6.7461157975713304;
          const walberla::float64 tmp_kernel_op_74 = 0.0012367874205397103;
          const walberla::float64 tmp_kernel_op_75 = 0.066954266063235096;
          const walberla::float64 tmp_kernel_op_76 = 0.094435052310005457;
          const walberla::float64 tmp_kernel_op_77 = 1.3926371369430026;
          const walberla::float64 tmp_kernel_op_78 = 0.1165079743250914;
          const walberla::float64 tmp_kernel_op_79 = 2.001552044859173;
          const walberla::float64 elMatDiag_0 =
              tmp_kernel_op_11 *
                  ((-tmp_kernel_op_10 - tmp_kernel_op_9) *
                   (-tmp_kernel_op_10 - tmp_kernel_op_9)) *
                  0.42163223734820815 +
              tmp_kernel_op_14 *
                  ((-tmp_kernel_op_12 - tmp_kernel_op_13) *
                   (-tmp_kernel_op_12 - tmp_kernel_op_13)) *
                  7.7299213783731894e-5 +
              tmp_kernel_op_17 *
                  ((-tmp_kernel_op_15 - tmp_kernel_op_16) *
                   (-tmp_kernel_op_15 - tmp_kernel_op_16)) *
                  0.0041846416289521935 +
              tmp_kernel_op_2 *
                  ((-tmp_kernel_op_0 - tmp_kernel_op_1) *
                   (-tmp_kernel_op_0 - tmp_kernel_op_1)) *
                  1.6606959041833929 +
              tmp_kernel_op_20 *
                  ((-tmp_kernel_op_18 - tmp_kernel_op_19) *
                   (-tmp_kernel_op_18 - tmp_kernel_op_19)) *
                  0.0059021907693753411 +
              tmp_kernel_op_23 *
                  ((-tmp_kernel_op_21 - tmp_kernel_op_22) *
                   (-tmp_kernel_op_21 - tmp_kernel_op_22)) *
                  0.087039821058937664 +
              tmp_kernel_op_26 *
                  ((-tmp_kernel_op_24 - tmp_kernel_op_25) *
                   (-tmp_kernel_op_24 - tmp_kernel_op_25)) *
                  0.0072817483953182124 +
              tmp_kernel_op_29 *
                  ((-tmp_kernel_op_27 - tmp_kernel_op_28) *
                   (-tmp_kernel_op_27 - tmp_kernel_op_28)) *
                  0.12509700280369832 +
              tmp_kernel_op_5 *
                  ((-tmp_kernel_op_3 - tmp_kernel_op_4) *
                   (-tmp_kernel_op_3 - tmp_kernel_op_4)) *
                  0.62130528681440322 +
              tmp_kernel_op_8 *
                  ((-tmp_kernel_op_6 - tmp_kernel_op_7) *
                   (-tmp_kernel_op_6 - tmp_kernel_op_7)) *
                  0.16138600724470512;
          const walberla::float64 elMatDiag_1 =
              (tmp_kernel_op_1 * tmp_kernel_op_1) * tmp_kernel_op_30 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) * tmp_kernel_op_33 +
              (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_34 +
              (tmp_kernel_op_16 * tmp_kernel_op_16) * tmp_kernel_op_35 +
              (tmp_kernel_op_19 * tmp_kernel_op_19) * tmp_kernel_op_36 +
              (tmp_kernel_op_22 * tmp_kernel_op_22) * tmp_kernel_op_37 +
              (tmp_kernel_op_25 * tmp_kernel_op_25) * tmp_kernel_op_38 +
              (tmp_kernel_op_28 * tmp_kernel_op_28) * tmp_kernel_op_39 +
              tmp_kernel_op_31 * (tmp_kernel_op_4 * tmp_kernel_op_4) +
              tmp_kernel_op_32 * (tmp_kernel_op_7 * tmp_kernel_op_7);
          const walberla::float64 elMatDiag_2 =
              tmp_kernel_op_41 *
                  ((tmp_kernel_op_0 - 1.0) * (tmp_kernel_op_0 - 1.0)) +
              tmp_kernel_op_43 *
                  ((tmp_kernel_op_3 - 1.0) * (tmp_kernel_op_3 - 1.0)) +
              tmp_kernel_op_45 *
                  ((tmp_kernel_op_6 - 1.0) * (tmp_kernel_op_6 - 1.0)) +
              tmp_kernel_op_47 *
                  ((tmp_kernel_op_9 - 1.0) * (tmp_kernel_op_9 - 1.0)) +
              tmp_kernel_op_49 *
                  ((tmp_kernel_op_12 - 1.0) * (tmp_kernel_op_12 - 1.0)) +
              tmp_kernel_op_51 *
                  ((tmp_kernel_op_15 - 1.0) * (tmp_kernel_op_15 - 1.0)) +
              tmp_kernel_op_53 *
                  ((tmp_kernel_op_18 - 1.0) * (tmp_kernel_op_18 - 1.0)) +
              tmp_kernel_op_55 *
                  ((tmp_kernel_op_21 - 1.0) * (tmp_kernel_op_21 - 1.0)) +
              tmp_kernel_op_57 *
                  ((tmp_kernel_op_24 - 1.0) * (tmp_kernel_op_24 - 1.0)) +
              tmp_kernel_op_59 *
                  ((tmp_kernel_op_27 - 1.0) * (tmp_kernel_op_27 - 1.0));
          const walberla::float64 elMatDiag_3 =
              tmp_kernel_op_60 * 16.0 + tmp_kernel_op_61 * 16.0 +
              tmp_kernel_op_62 * 16.0 + tmp_kernel_op_63 * 16.0 +
              tmp_kernel_op_64 * 16.0 + tmp_kernel_op_65 * 16.0 +
              tmp_kernel_op_66 * 16.0 + tmp_kernel_op_67 * 16.0 +
              tmp_kernel_op_68 * 16.0 + tmp_kernel_op_69 * 16.0;
          const walberla::float64 elMatDiag_4 =
              tmp_kernel_op_41 * tmp_kernel_op_70 +
              tmp_kernel_op_43 * tmp_kernel_op_71 +
              tmp_kernel_op_45 * tmp_kernel_op_72 +
              tmp_kernel_op_47 * tmp_kernel_op_73 +
              tmp_kernel_op_49 * tmp_kernel_op_74 +
              tmp_kernel_op_51 * tmp_kernel_op_75 +
              tmp_kernel_op_53 * tmp_kernel_op_76 +
              tmp_kernel_op_55 * tmp_kernel_op_77 +
              tmp_kernel_op_57 * tmp_kernel_op_78 +
              tmp_kernel_op_59 * tmp_kernel_op_79;
          const walberla::float64 elMatDiag_5 =
              tmp_kernel_op_30 * tmp_kernel_op_70 +
              tmp_kernel_op_31 * tmp_kernel_op_71 +
              tmp_kernel_op_32 * tmp_kernel_op_72 +
              tmp_kernel_op_33 * tmp_kernel_op_73 +
              tmp_kernel_op_34 * tmp_kernel_op_74 +
              tmp_kernel_op_35 * tmp_kernel_op_75 +
              tmp_kernel_op_36 * tmp_kernel_op_76 +
              tmp_kernel_op_37 * tmp_kernel_op_77 +
              tmp_kernel_op_38 * tmp_kernel_op_78 +
              tmp_kernel_op_39 * tmp_kernel_op_79;
          const walberla::float64 elMatDiag_6 =
              tmp_kernel_op_60 * 1210.6473141496936 +
              tmp_kernel_op_61 * 452.93155408770002 +
              tmp_kernel_op_62 * 117.65039928139001 +
              tmp_kernel_op_63 * 307.36990102684365 +
              tmp_kernel_op_64 * 0.056351126848341607 +
              tmp_kernel_op_65 * 3.0506037475061465 +
              tmp_kernel_op_66 * 4.3026970708746193 +
              tmp_kernel_op_67 * 63.452029551965545 +
              tmp_kernel_op_68 * 5.3083945801869801 +
              tmp_kernel_op_69 * 91.195715043896044;
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
      const walberla::float64 abs_det_jac_affine_BLUE =
          abs(jac_affine_0_0_BLUE * jac_affine_1_1_BLUE -
              jac_affine_0_1_BLUE * jac_affine_1_0_BLUE);
      for (int64_t ctr_1 = 0; ctr_1 < micro_edges_per_macro_edge; ctr_1 += 1)
        for (int64_t ctr_0 = 0; ctr_0 < -ctr_1 + micro_edges_per_macro_edge - 1;
             ctr_0 += 1) {
#if 0
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
#endif
          const walberla::float64 jac_blending_1_1_q_9 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_9 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_9 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_9 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_9 =
              jac_blending_0_0_q_9 * jac_blending_1_1_q_9 -
              jac_blending_0_1_q_9 * jac_blending_1_0_q_9;
          const walberla::float64 jac_blending_1_1_q_8 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_8 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_8 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_8 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_8 =
              jac_blending_0_0_q_8 * jac_blending_1_1_q_8 -
              jac_blending_0_1_q_8 * jac_blending_1_0_q_8;
          const walberla::float64 jac_blending_1_1_q_7 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_7 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_7 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_7 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_7 =
              jac_blending_0_0_q_7 * jac_blending_1_1_q_7 -
              jac_blending_0_1_q_7 * jac_blending_1_0_q_7;
          const walberla::float64 jac_blending_1_1_q_6 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_6 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_6 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_6 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_6 =
              jac_blending_0_0_q_6 * jac_blending_1_1_q_6 -
              jac_blending_0_1_q_6 * jac_blending_1_0_q_6;
          const walberla::float64 jac_blending_1_1_q_5 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_5 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_5 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_5 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_5 =
              jac_blending_0_0_q_5 * jac_blending_1_1_q_5 -
              jac_blending_0_1_q_5 * jac_blending_1_0_q_5;
          const walberla::float64 jac_blending_1_1_q_4 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_4 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_4 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_4 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_4 =
              jac_blending_0_0_q_4 * jac_blending_1_1_q_4 -
              jac_blending_0_1_q_4 * jac_blending_1_0_q_4;
          const walberla::float64 jac_blending_1_1_q_3 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_3 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_3 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_3 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_2 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_2 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_2 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_1 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_1 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_1 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 = bMat_11;
          const walberla::float64 jac_blending_1_0_q_0 = bMat_10;
          const walberla::float64 jac_blending_0_1_q_0 = bMat_01;
          const walberla::float64 jac_blending_0_0_q_0 = bMat_00;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_1 = -1.0175839942884126;
          const walberla::float64 tmp_kernel_op_2 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_0 *
                                                    0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_3 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_4 = -0.87062238788797819;
          const walberla::float64 tmp_kernel_op_5 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_1 *
                                                    0.052397656566423402;
          const walberla::float64 tmp_kernel_op_6 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_7 = -0.84634856630183575;
          const walberla::float64 tmp_kernel_op_8 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_2 *
                                                    0.074275554650521658;
          const walberla::float64 tmp_kernel_op_9 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_10 = -0.40994976126116967;
          const walberla::float64 tmp_kernel_op_11 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_3 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_12 = -0.55977582085839206;
          const walberla::float64 tmp_kernel_op_13 = 1.5773598151468047;
          const walberla::float64 tmp_kernel_op_14 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_4 *
                                                     0.00025771437964227723;
          const walberla::float64 tmp_kernel_op_15 = 0.2941639687542501;
          const walberla::float64 tmp_kernel_op_16 = 0.57645841913372808;
          const walberla::float64 tmp_kernel_op_17 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_5 *
                                                     0.052397656566423402;
          const walberla::float64 tmp_kernel_op_18 = 1.042891019832316;
          const walberla::float64 tmp_kernel_op_19 = -0.19654245353048039;
          const walberla::float64 tmp_kernel_op_20 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_6 *
                                                     0.074275554650521658;
          const walberla::float64 tmp_kernel_op_21 = 0.11128546692146923;
          const walberla::float64 tmp_kernel_op_22 = 0.29866429433970043;
          const walberla::float64 tmp_kernel_op_23 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_7 *
                                                     0.047488619588783712;
          const walberla::float64 tmp_kernel_op_24 = 1.6586673553187572;
          const walberla::float64 tmp_kernel_op_25 = -0.8293336776593786;
          const walberla::float64 tmp_kernel_op_26 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_8 *
                                                     0.045496761795224737;
          const walberla::float64 tmp_kernel_op_27 = 0.58523781332014213;
          const walberla::float64 tmp_kernel_op_28 = -0.29261890666007107;
          const walberla::float64 tmp_kernel_op_29 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_9 *
                                                     0.10566414783403316;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_2 * 7.7299213783731935e-5;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_5 * 0.0041846416289521935;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_8 * 0.0059021907693753367;
          const walberla::float64 tmp_kernel_op_33 =
              tmp_kernel_op_11 * 0.087039821058937664;
          const walberla::float64 tmp_kernel_op_34 =
              tmp_kernel_op_14 * 1.6606959041833929;
          const walberla::float64 tmp_kernel_op_35 =
              tmp_kernel_op_17 * 0.62130528681440322;
          const walberla::float64 tmp_kernel_op_36 =
              tmp_kernel_op_20 * 0.16138600724470506;
          const walberla::float64 tmp_kernel_op_37 =
              tmp_kernel_op_23 * 0.42163223734820804;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_26 * 0.0072817483953182219;
          const walberla::float64 tmp_kernel_op_39 =
              tmp_kernel_op_29 * 0.12509700280369832;
          const walberla::float64 tmp_kernel_op_40 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_41 =
              tmp_kernel_op_2 * tmp_kernel_op_40;
          const walberla::float64 tmp_kernel_op_42 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_43 =
              tmp_kernel_op_42 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_44 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_45 =
              tmp_kernel_op_44 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_46 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_47 =
              tmp_kernel_op_11 * tmp_kernel_op_46;
          const walberla::float64 tmp_kernel_op_48 = 0.078337242404421664;
          const walberla::float64 tmp_kernel_op_49 =
              tmp_kernel_op_14 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_50 = 0.021633110128312857;
          const walberla::float64 tmp_kernel_op_51 =
              tmp_kernel_op_17 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_52 = 0.27190541981172206;
          const walberla::float64 tmp_kernel_op_53 =
              tmp_kernel_op_20 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_54 = 0.0030961137869823557;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_23 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_56 = 0.68779434890003011;
          const walberla::float64 tmp_kernel_op_57 =
              tmp_kernel_op_26 * tmp_kernel_op_56;
          const walberla::float64 tmp_kernel_op_58 = 0.085625824534935377;
          const walberla::float64 tmp_kernel_op_59 =
              tmp_kernel_op_29 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_30 * tmp_kernel_op_40;
          const walberla::float64 tmp_kernel_op_61 =
              tmp_kernel_op_31 * tmp_kernel_op_42;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_32 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_63 =
              tmp_kernel_op_33 * tmp_kernel_op_46;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_34 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_35 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_36 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_37 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_38 * tmp_kernel_op_56;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_39 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_70 = 26.571134466934286;
          const walberla::float64 tmp_kernel_op_71 = 9.9408845890304516;
          const walberla::float64 tmp_kernel_op_72 = 2.5821761159152818;
          const walberla::float64 tmp_kernel_op_73 = 6.7461157975713304;
          const walberla::float64 tmp_kernel_op_74 = 0.0012367874205397103;
          const walberla::float64 tmp_kernel_op_75 = 0.066954266063235096;
          const walberla::float64 tmp_kernel_op_76 = 0.094435052310005457;
          const walberla::float64 tmp_kernel_op_77 = 1.3926371369430026;
          const walberla::float64 tmp_kernel_op_78 = 0.1165079743250914;
          const walberla::float64 tmp_kernel_op_79 = 2.001552044859173;
          const walberla::float64 elMatDiag_0 =
              tmp_kernel_op_11 *
                  ((-tmp_kernel_op_10 - tmp_kernel_op_9) *
                   (-tmp_kernel_op_10 - tmp_kernel_op_9)) *
                  0.42163223734820815 +
              tmp_kernel_op_14 *
                  ((-tmp_kernel_op_12 - tmp_kernel_op_13) *
                   (-tmp_kernel_op_12 - tmp_kernel_op_13)) *
                  7.7299213783731894e-5 +
              tmp_kernel_op_17 *
                  ((-tmp_kernel_op_15 - tmp_kernel_op_16) *
                   (-tmp_kernel_op_15 - tmp_kernel_op_16)) *
                  0.0041846416289521935 +
              tmp_kernel_op_2 *
                  ((-tmp_kernel_op_0 - tmp_kernel_op_1) *
                   (-tmp_kernel_op_0 - tmp_kernel_op_1)) *
                  1.6606959041833929 +
              tmp_kernel_op_20 *
                  ((-tmp_kernel_op_18 - tmp_kernel_op_19) *
                   (-tmp_kernel_op_18 - tmp_kernel_op_19)) *
                  0.0059021907693753411 +
              tmp_kernel_op_23 *
                  ((-tmp_kernel_op_21 - tmp_kernel_op_22) *
                   (-tmp_kernel_op_21 - tmp_kernel_op_22)) *
                  0.087039821058937664 +
              tmp_kernel_op_26 *
                  ((-tmp_kernel_op_24 - tmp_kernel_op_25) *
                   (-tmp_kernel_op_24 - tmp_kernel_op_25)) *
                  0.0072817483953182124 +
              tmp_kernel_op_29 *
                  ((-tmp_kernel_op_27 - tmp_kernel_op_28) *
                   (-tmp_kernel_op_27 - tmp_kernel_op_28)) *
                  0.12509700280369832 +
              tmp_kernel_op_5 *
                  ((-tmp_kernel_op_3 - tmp_kernel_op_4) *
                   (-tmp_kernel_op_3 - tmp_kernel_op_4)) *
                  0.62130528681440322 +
              tmp_kernel_op_8 *
                  ((-tmp_kernel_op_6 - tmp_kernel_op_7) *
                   (-tmp_kernel_op_6 - tmp_kernel_op_7)) *
                  0.16138600724470512;
          const walberla::float64 elMatDiag_1 =
              (tmp_kernel_op_1 * tmp_kernel_op_1) * tmp_kernel_op_30 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) * tmp_kernel_op_33 +
              (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_34 +
              (tmp_kernel_op_16 * tmp_kernel_op_16) * tmp_kernel_op_35 +
              (tmp_kernel_op_19 * tmp_kernel_op_19) * tmp_kernel_op_36 +
              (tmp_kernel_op_22 * tmp_kernel_op_22) * tmp_kernel_op_37 +
              (tmp_kernel_op_25 * tmp_kernel_op_25) * tmp_kernel_op_38 +
              (tmp_kernel_op_28 * tmp_kernel_op_28) * tmp_kernel_op_39 +
              tmp_kernel_op_31 * (tmp_kernel_op_4 * tmp_kernel_op_4) +
              tmp_kernel_op_32 * (tmp_kernel_op_7 * tmp_kernel_op_7);
          const walberla::float64 elMatDiag_2 =
              tmp_kernel_op_41 *
                  ((tmp_kernel_op_0 - 1.0) * (tmp_kernel_op_0 - 1.0)) +
              tmp_kernel_op_43 *
                  ((tmp_kernel_op_3 - 1.0) * (tmp_kernel_op_3 - 1.0)) +
              tmp_kernel_op_45 *
                  ((tmp_kernel_op_6 - 1.0) * (tmp_kernel_op_6 - 1.0)) +
              tmp_kernel_op_47 *
                  ((tmp_kernel_op_9 - 1.0) * (tmp_kernel_op_9 - 1.0)) +
              tmp_kernel_op_49 *
                  ((tmp_kernel_op_12 - 1.0) * (tmp_kernel_op_12 - 1.0)) +
              tmp_kernel_op_51 *
                  ((tmp_kernel_op_15 - 1.0) * (tmp_kernel_op_15 - 1.0)) +
              tmp_kernel_op_53 *
                  ((tmp_kernel_op_18 - 1.0) * (tmp_kernel_op_18 - 1.0)) +
              tmp_kernel_op_55 *
                  ((tmp_kernel_op_21 - 1.0) * (tmp_kernel_op_21 - 1.0)) +
              tmp_kernel_op_57 *
                  ((tmp_kernel_op_24 - 1.0) * (tmp_kernel_op_24 - 1.0)) +
              tmp_kernel_op_59 *
                  ((tmp_kernel_op_27 - 1.0) * (tmp_kernel_op_27 - 1.0));
          const walberla::float64 elMatDiag_3 =
              tmp_kernel_op_60 * 16.0 + tmp_kernel_op_61 * 16.0 +
              tmp_kernel_op_62 * 16.0 + tmp_kernel_op_63 * 16.0 +
              tmp_kernel_op_64 * 16.0 + tmp_kernel_op_65 * 16.0 +
              tmp_kernel_op_66 * 16.0 + tmp_kernel_op_67 * 16.0 +
              tmp_kernel_op_68 * 16.0 + tmp_kernel_op_69 * 16.0;
          const walberla::float64 elMatDiag_4 =
              tmp_kernel_op_41 * tmp_kernel_op_70 +
              tmp_kernel_op_43 * tmp_kernel_op_71 +
              tmp_kernel_op_45 * tmp_kernel_op_72 +
              tmp_kernel_op_47 * tmp_kernel_op_73 +
              tmp_kernel_op_49 * tmp_kernel_op_74 +
              tmp_kernel_op_51 * tmp_kernel_op_75 +
              tmp_kernel_op_53 * tmp_kernel_op_76 +
              tmp_kernel_op_55 * tmp_kernel_op_77 +
              tmp_kernel_op_57 * tmp_kernel_op_78 +
              tmp_kernel_op_59 * tmp_kernel_op_79;
          const walberla::float64 elMatDiag_5 =
              tmp_kernel_op_30 * tmp_kernel_op_70 +
              tmp_kernel_op_31 * tmp_kernel_op_71 +
              tmp_kernel_op_32 * tmp_kernel_op_72 +
              tmp_kernel_op_33 * tmp_kernel_op_73 +
              tmp_kernel_op_34 * tmp_kernel_op_74 +
              tmp_kernel_op_35 * tmp_kernel_op_75 +
              tmp_kernel_op_36 * tmp_kernel_op_76 +
              tmp_kernel_op_37 * tmp_kernel_op_77 +
              tmp_kernel_op_38 * tmp_kernel_op_78 +
              tmp_kernel_op_39 * tmp_kernel_op_79;
          const walberla::float64 elMatDiag_6 =
              tmp_kernel_op_60 * 1210.6473141496936 +
              tmp_kernel_op_61 * 452.93155408770002 +
              tmp_kernel_op_62 * 117.65039928139001 +
              tmp_kernel_op_63 * 307.36990102684365 +
              tmp_kernel_op_64 * 0.056351126848341607 +
              tmp_kernel_op_65 * 3.0506037475061465 +
              tmp_kernel_op_66 * 4.3026970708746193 +
              tmp_kernel_op_67 * 63.452029551965545 +
              tmp_kernel_op_68 * 5.3083945801869801 +
              tmp_kernel_op_69 * 91.195715043896044;
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
