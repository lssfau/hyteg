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

#include "P2PlusBubbleElementwiseMass_AnnulusMap_float64.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P2PlusBubbleElementwiseMass_AnnulusMap_float64::
    P2PlusBubbleElementwiseMass_AnnulusMap_float64(
        const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
        size_t maxLevel)
    : Operator(storage, minLevel, maxLevel) {}

void P2PlusBubbleElementwiseMass_AnnulusMap_float64::apply(
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
      WALBERLA_CHECK_NOT_NULLPTR(
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap()),
          "This operator requires the AnnulusMap to be registered as "
          "GeometryMap on every macro-cell.")
      real_t radRefVertex =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->radRefVertex();
      real_t radRayVertex =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->radRayVertex();
      real_t refVertex_0 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->refVertex()[0];
      real_t rayVertex_0 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->rayVertex()[0];
      real_t thrVertex_0 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->thrVertex()[0];
      real_t refVertex_1 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->refVertex()[1];
      real_t rayVertex_1 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->rayVertex()[1];
      real_t thrVertex_1 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->thrVertex()[1];

      this->timingTree_->start("kernel");

      apply_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(

          _data_dst, _data_dstEdge, _data_dstVertex, _data_src, _data_srcEdge,
          _data_srcVertex, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_1comp0,
          macro_vertex_coord_id_1comp1, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float, radRayVertex, radRefVertex,
          rayVertex_0, rayVertex_1, refVertex_0, refVertex_1, thrVertex_0,
          thrVertex_1);

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
void P2PlusBubbleElementwiseMass_AnnulusMap_float64::toMatrix(
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
      WALBERLA_CHECK_NOT_NULLPTR(
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap()),
          "This operator requires the AnnulusMap to be registered as "
          "GeometryMap on every macro-cell.")
      real_t radRefVertex =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->radRefVertex();
      real_t radRayVertex =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->radRayVertex();
      real_t refVertex_0 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->refVertex()[0];
      real_t rayVertex_0 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->rayVertex()[0];
      real_t thrVertex_0 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->thrVertex()[0];
      real_t refVertex_1 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->refVertex()[1];
      real_t rayVertex_1 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->rayVertex()[1];
      real_t thrVertex_1 =
          std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
              ->thrVertex()[1];

      this->timingTree_->start("kernel");

      toMatrix_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(

          _data_dst, _data_dstEdge, _data_dstVertex, _data_src, _data_srcEdge,
          _data_srcVertex, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_1comp0,
          macro_vertex_coord_id_1comp1, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, mat, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float, radRayVertex, radRefVertex,
          rayVertex_0, rayVertex_1, refVertex_0, refVertex_1, thrVertex_0,
          thrVertex_1);

      this->timingTree_->stop("kernel");
    }
  }
  this->stopTiming("toMatrix");
}
void P2PlusBubbleElementwiseMass_AnnulusMap_float64::
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
        WALBERLA_CHECK_NOT_NULLPTR(
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap()),
            "This operator requires the AnnulusMap to be registered as "
            "GeometryMap on every macro-cell.")
        real_t radRefVertex =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->radRefVertex();
        real_t radRayVertex =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->radRayVertex();
        real_t refVertex_0 =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->refVertex()[0];
        real_t rayVertex_0 =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->rayVertex()[0];
        real_t thrVertex_0 =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->thrVertex()[0];
        real_t refVertex_1 =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->refVertex()[1];
        real_t rayVertex_1 =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->rayVertex()[1];
        real_t thrVertex_1 =
            std::dynamic_pointer_cast<AnnulusMap>(face.getGeometryMap())
                ->thrVertex()[1];

        this->timingTree_->start("kernel");

        computeInverseDiagonalOperatorValues_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(

            _data_invDiag_, _data_invDiag_Edge, _data_invDiag_Vertex,
            macro_vertex_coord_id_0comp0, macro_vertex_coord_id_0comp1,
            macro_vertex_coord_id_1comp0, macro_vertex_coord_id_1comp1,
            macro_vertex_coord_id_2comp0, macro_vertex_coord_id_2comp1,
            micro_edges_per_macro_edge, micro_edges_per_macro_edge_float,
            radRayVertex, radRefVertex, rayVertex_0, rayVertex_1, refVertex_0,
            refVertex_1, thrVertex_0, thrVertex_1);

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
P2PlusBubbleElementwiseMass_AnnulusMap_float64::getInverseDiagonalValues()
    const {
  return invDiag_;
}
void P2PlusBubbleElementwiseMass_AnnulusMap_float64::
    apply_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(
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
        walberla::float64 micro_edges_per_macro_edge_float,
        walberla::float64 radRayVertex, walberla::float64 radRefVertex,
        walberla::float64 rayVertex_0, walberla::float64 rayVertex_1,
        walberla::float64 refVertex_0, walberla::float64 refVertex_1,
        walberla::float64 thrVertex_0, walberla::float64 thrVertex_1) const {
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
          const walberla::float64 tmp_blending_op_0 =
              -rayVertex_0 + thrVertex_0;
          const walberla::float64 tmp_blending_op_1 =
              -p_affine_0_1 + p_affine_1_1;
          const walberla::float64 tmp_blending_op_2 =
              -p_affine_0_1 + p_affine_2_1;
          const walberla::float64 tmp_blending_op_3 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_4 =
              -p_affine_0_0 + p_affine_1_0;
          const walberla::float64 tmp_blending_op_5 =
              -p_affine_0_0 + p_affine_2_0;
          const walberla::float64 tmp_blending_op_6 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_7 =
              (tmp_blending_op_6 * tmp_blending_op_6);
          const walberla::float64 tmp_blending_op_8 =
              (tmp_blending_op_3 * tmp_blending_op_3);
          const walberla::float64 tmp_blending_op_9 =
              tmp_blending_op_7 + tmp_blending_op_8;
          const walberla::float64 tmp_blending_op_10 =
              -rayVertex_1 + thrVertex_1;
          const walberla::float64 tmp_blending_op_11 =
              (-radRayVertex + radRefVertex) * 1.0 /
              (-tmp_blending_op_0 * (-rayVertex_1 + refVertex_1) +
               tmp_blending_op_10 * (-rayVertex_0 + refVertex_0));
          const walberla::float64 tmp_blending_op_12 = tmp_blending_op_11 * 1.0;
          const walberla::float64 tmp_blending_op_13 =
              tmp_blending_op_12 * pow(tmp_blending_op_9, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_14 =
              tmp_blending_op_13 * tmp_blending_op_3;
          const walberla::float64 tmp_blending_op_15 =
              pow(tmp_blending_op_9, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_16 = -rayVertex_1;
          const walberla::float64 tmp_blending_op_17 = -rayVertex_0;
          const walberla::float64 tmp_blending_op_18 =
              radRayVertex + tmp_blending_op_11 *
                                 (-tmp_blending_op_0 *
                                      (tmp_blending_op_16 + tmp_blending_op_3) +
                                  tmp_blending_op_10 *
                                      (tmp_blending_op_17 + tmp_blending_op_6));
          const walberla::float64 tmp_blending_op_19 =
              tmp_blending_op_15 * tmp_blending_op_18 * 1.0;
          const walberla::float64 tmp_blending_op_20 =
              tmp_blending_op_13 * tmp_blending_op_6;
          const walberla::float64 tmp_blending_op_21 =
              p_affine_0_1 + tmp_blending_op_1 * 0.59999999999999998 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_22 =
              p_affine_0_0 + tmp_blending_op_4 * 0.59999999999999998 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_23 =
              (tmp_blending_op_22 * tmp_blending_op_22);
          const walberla::float64 tmp_blending_op_24 =
              (tmp_blending_op_21 * tmp_blending_op_21);
          const walberla::float64 tmp_blending_op_25 =
              tmp_blending_op_23 + tmp_blending_op_24;
          const walberla::float64 tmp_blending_op_26 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_25, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_27 =
              tmp_blending_op_21 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_28 =
              pow(tmp_blending_op_25, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_29 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_21) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_22));
          const walberla::float64 tmp_blending_op_30 =
              tmp_blending_op_28 * tmp_blending_op_29 * 1.0;
          const walberla::float64 tmp_blending_op_31 =
              tmp_blending_op_22 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_32 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_33 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_34 =
              (tmp_blending_op_33 * tmp_blending_op_33);
          const walberla::float64 tmp_blending_op_35 =
              (tmp_blending_op_32 * tmp_blending_op_32);
          const walberla::float64 tmp_blending_op_36 =
              tmp_blending_op_34 + tmp_blending_op_35;
          const walberla::float64 tmp_blending_op_37 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_36, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_38 =
              tmp_blending_op_32 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_39 =
              pow(tmp_blending_op_36, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_40 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_32) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_33));
          const walberla::float64 tmp_blending_op_41 =
              tmp_blending_op_39 * tmp_blending_op_40 * 1.0;
          const walberla::float64 tmp_blending_op_42 =
              tmp_blending_op_33 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_43 =
              p_affine_0_1 + tmp_blending_op_1 * 0.33333333333333331 +
              tmp_blending_op_2 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_44 =
              p_affine_0_0 + tmp_blending_op_4 * 0.33333333333333331 +
              tmp_blending_op_5 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_45 =
              (tmp_blending_op_44 * tmp_blending_op_44);
          const walberla::float64 tmp_blending_op_46 =
              (tmp_blending_op_43 * tmp_blending_op_43);
          const walberla::float64 tmp_blending_op_47 =
              tmp_blending_op_45 + tmp_blending_op_46;
          const walberla::float64 tmp_blending_op_48 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_47, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_49 =
              tmp_blending_op_43 * tmp_blending_op_48;
          const walberla::float64 tmp_blending_op_50 =
              pow(tmp_blending_op_47, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_51 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_43) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_44));
          const walberla::float64 tmp_blending_op_52 =
              tmp_blending_op_50 * tmp_blending_op_51 * 1.0;
          const walberla::float64 tmp_blending_op_53 =
              tmp_blending_op_44 * tmp_blending_op_48;
          const walberla::float64 jac_blending_1_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_14 +
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_7 * 1.0;
          const walberla::float64 jac_blending_1_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_14 -
              tmp_blending_op_19 * tmp_blending_op_3 * tmp_blending_op_6;
          const walberla::float64 jac_blending_0_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_20 -
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_3 *
                  tmp_blending_op_6;
          const walberla::float64 jac_blending_0_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_20 +
              tmp_blending_op_19 * tmp_blending_op_8;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_27 +
              tmp_blending_op_23 * tmp_blending_op_28 * tmp_blending_op_29 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_27 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_30;
          const walberla::float64 jac_blending_0_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_31 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_28 *
                  tmp_blending_op_29;
          const walberla::float64 jac_blending_0_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_31 +
              tmp_blending_op_24 * tmp_blending_op_30;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_38 +
              tmp_blending_op_34 * tmp_blending_op_39 * tmp_blending_op_40 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_38 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_41;
          const walberla::float64 jac_blending_0_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_42 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_39 *
                  tmp_blending_op_40;
          const walberla::float64 jac_blending_0_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_42 +
              tmp_blending_op_35 * tmp_blending_op_41;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_49 +
              tmp_blending_op_45 * tmp_blending_op_50 * tmp_blending_op_51 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_49 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_52;
          const walberla::float64 jac_blending_0_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_53 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_50 *
                  tmp_blending_op_51;
          const walberla::float64 jac_blending_0_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_53 +
              tmp_blending_op_46 * tmp_blending_op_52;
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
          const walberla::float64 tmp_kernel_op_0 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_1 =
              abs_det_jac_affine_GRAY * abs_det_jac_blending_q_0 * -0.28125;
          const walberla::float64 tmp_kernel_op_2 = 0.33333333333333343;
          const walberla::float64 tmp_kernel_op_3 = 0.66666666666666663;
          const walberla::float64 tmp_kernel_op_4 = -0.33333333333333337;
          const walberla::float64 tmp_kernel_op_5 =
              -tmp_kernel_op_3 - tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_6 =
              tmp_kernel_op_2 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_7 =
              tmp_kernel_op_1 * tmp_kernel_op_6 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_8 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_1 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_10 = 0.20000000000000007;
          const walberla::float64 tmp_kernel_op_11 = 1.2;
          const walberla::float64 tmp_kernel_op_12 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_15 =
              tmp_kernel_op_14 * tmp_kernel_op_9 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_16 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_17 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_2 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_18 = 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_19 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_20 = 0.19999999999999996;
          const walberla::float64 tmp_kernel_op_21 =
              -tmp_kernel_op_19 - tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_18 * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_17 * tmp_kernel_op_22 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_24 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_25 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_3 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_26 = 0.60000000000000009;
          const walberla::float64 tmp_kernel_op_27 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_28 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_29 =
              -tmp_kernel_op_27 - tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_26 * tmp_kernel_op_29;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_25 * tmp_kernel_op_30 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_0 * tmp_kernel_op_7 +
              tmp_kernel_op_15 * tmp_kernel_op_8 +
              tmp_kernel_op_16 * tmp_kernel_op_23 +
              tmp_kernel_op_24 * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_33 =
              tmp_kernel_op_1 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_34 =
              tmp_kernel_op_33 * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_35 =
              tmp_kernel_op_9 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_36 =
              tmp_kernel_op_14 * tmp_kernel_op_35;
          const walberla::float64 tmp_kernel_op_37 =
              tmp_kernel_op_17 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_22 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_39 =
              tmp_kernel_op_25 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_40 =
              tmp_kernel_op_30 * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_41 =
              tmp_kernel_op_12 * tmp_kernel_op_36 +
              tmp_kernel_op_20 * tmp_kernel_op_38 +
              tmp_kernel_op_28 * tmp_kernel_op_40 +
              tmp_kernel_op_34 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_42 = tmp_kernel_op_3 - 1.0;
          const walberla::float64 tmp_kernel_op_43 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_44 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_45 = tmp_kernel_op_27 - 1.0;
          const walberla::float64 tmp_kernel_op_46 =
              tmp_kernel_op_15 * tmp_kernel_op_43 +
              tmp_kernel_op_23 * tmp_kernel_op_44 +
              tmp_kernel_op_31 * tmp_kernel_op_45 +
              tmp_kernel_op_42 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_47 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_48 =
              -tmp_kernel_op_0 - tmp_kernel_op_47 + 4.0;
          const walberla::float64 tmp_kernel_op_49 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_50 =
              -tmp_kernel_op_49 - tmp_kernel_op_8 + 4.0;
          const walberla::float64 tmp_kernel_op_51 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_52 =
              -tmp_kernel_op_16 - tmp_kernel_op_51 + 4.0;
          const walberla::float64 tmp_kernel_op_53 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_54 =
              -tmp_kernel_op_24 - tmp_kernel_op_53 + 4.0;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_34 * tmp_kernel_op_48 +
              tmp_kernel_op_36 * tmp_kernel_op_50 +
              tmp_kernel_op_38 * tmp_kernel_op_52 +
              tmp_kernel_op_40 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_56 =
              tmp_kernel_op_15 * tmp_kernel_op_50 +
              tmp_kernel_op_23 * tmp_kernel_op_52 +
              tmp_kernel_op_31 * tmp_kernel_op_54 +
              tmp_kernel_op_48 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_57 = 9.0;
          const walberla::float64 tmp_kernel_op_58 =
              tmp_kernel_op_33 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_59 = 5.4000000000000021;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_35 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_61 = 5.4000000000000004;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_37 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_63 = 16.200000000000003;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_39 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_14 * tmp_kernel_op_59 * tmp_kernel_op_60 +
              tmp_kernel_op_22 * tmp_kernel_op_61 * tmp_kernel_op_62 +
              tmp_kernel_op_30 * tmp_kernel_op_63 * tmp_kernel_op_64 +
              tmp_kernel_op_57 * tmp_kernel_op_58 * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_1 * 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_9 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_17 * 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_25 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_4 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_12 * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_72 =
              tmp_kernel_op_20 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_28 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_74 =
              tmp_kernel_op_47 * tmp_kernel_op_70 +
              tmp_kernel_op_49 * tmp_kernel_op_71 +
              tmp_kernel_op_51 * tmp_kernel_op_72 +
              tmp_kernel_op_53 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_75 =
              tmp_kernel_op_4 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_12 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_77 =
              tmp_kernel_op_20 * tmp_kernel_op_62;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_28 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_42 * tmp_kernel_op_75 +
              tmp_kernel_op_43 * tmp_kernel_op_76 +
              tmp_kernel_op_44 * tmp_kernel_op_77 +
              tmp_kernel_op_45 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_48 * tmp_kernel_op_75 +
              tmp_kernel_op_50 * tmp_kernel_op_76 +
              tmp_kernel_op_52 * tmp_kernel_op_77 +
              tmp_kernel_op_54 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_81 =
              tmp_kernel_op_48 * tmp_kernel_op_70 +
              tmp_kernel_op_50 * tmp_kernel_op_71 +
              tmp_kernel_op_52 * tmp_kernel_op_72 +
              tmp_kernel_op_54 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_83 =
              tmp_kernel_op_59 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_61 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_70 * tmp_kernel_op_82 +
              tmp_kernel_op_71 * tmp_kernel_op_83 +
              tmp_kernel_op_72 * tmp_kernel_op_84 +
              tmp_kernel_op_73 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_87 = 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_1 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 = 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_89 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_91 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_17 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_93 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_25 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              tmp_kernel_op_42 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_43 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_44 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_45 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_99 =
              tmp_kernel_op_0 * tmp_kernel_op_95 +
              tmp_kernel_op_16 * tmp_kernel_op_97 +
              tmp_kernel_op_24 * tmp_kernel_op_98 +
              tmp_kernel_op_8 * tmp_kernel_op_96;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_42 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_43 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_102 =
              tmp_kernel_op_44 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_103 =
              tmp_kernel_op_45 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_104 =
              tmp_kernel_op_100 * tmp_kernel_op_58 +
              tmp_kernel_op_101 * tmp_kernel_op_60 +
              tmp_kernel_op_102 * tmp_kernel_op_62 +
              tmp_kernel_op_103 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_105 =
              tmp_kernel_op_100 * tmp_kernel_op_88 +
              tmp_kernel_op_101 * tmp_kernel_op_90 +
              tmp_kernel_op_102 * tmp_kernel_op_92 +
              tmp_kernel_op_103 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_107 =
              tmp_kernel_op_59 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_108 =
              tmp_kernel_op_61 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_109 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_106 * tmp_kernel_op_95 +
              tmp_kernel_op_107 * tmp_kernel_op_96 +
              tmp_kernel_op_108 * tmp_kernel_op_97 +
              tmp_kernel_op_109 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_111 =
              (tmp_kernel_op_48 * tmp_kernel_op_48);
          const walberla::float64 tmp_kernel_op_112 =
              (tmp_kernel_op_50 * tmp_kernel_op_50);
          const walberla::float64 tmp_kernel_op_113 =
              (tmp_kernel_op_52 * tmp_kernel_op_52);
          const walberla::float64 tmp_kernel_op_114 =
              (tmp_kernel_op_54 * tmp_kernel_op_54);
          const walberla::float64 tmp_kernel_op_115 =
              tmp_kernel_op_111 * tmp_kernel_op_58 +
              tmp_kernel_op_112 * tmp_kernel_op_60 +
              tmp_kernel_op_113 * tmp_kernel_op_62 +
              tmp_kernel_op_114 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_116 =
              tmp_kernel_op_48 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_117 =
              tmp_kernel_op_50 * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_118 =
              tmp_kernel_op_52 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_119 =
              tmp_kernel_op_54 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_120 =
              tmp_kernel_op_116 * tmp_kernel_op_47 +
              tmp_kernel_op_117 * tmp_kernel_op_49 +
              tmp_kernel_op_118 * tmp_kernel_op_51 +
              tmp_kernel_op_119 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_121 =
              tmp_kernel_op_116 * tmp_kernel_op_82 +
              tmp_kernel_op_117 * tmp_kernel_op_83 +
              tmp_kernel_op_118 * tmp_kernel_op_84 +
              tmp_kernel_op_119 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_122 =
              tmp_kernel_op_66 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_123 =
              tmp_kernel_op_67 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_124 =
              tmp_kernel_op_68 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_125 =
              tmp_kernel_op_69 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_126 =
              tmp_kernel_op_48 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_127 =
              tmp_kernel_op_50 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_128 =
              tmp_kernel_op_52 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_54 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_130 =
              tmp_kernel_op_0 * tmp_kernel_op_126 +
              tmp_kernel_op_127 * tmp_kernel_op_8 +
              tmp_kernel_op_128 * tmp_kernel_op_16 +
              tmp_kernel_op_129 * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_131 =
              tmp_kernel_op_122 * tmp_kernel_op_57 * 4.0 +
              tmp_kernel_op_123 * tmp_kernel_op_59 * 4.0 +
              tmp_kernel_op_124 * tmp_kernel_op_61 * 4.0 +
              tmp_kernel_op_125 * tmp_kernel_op_63 * 4.0;
          const walberla::float64 tmp_kernel_op_132 =
              tmp_kernel_op_106 * tmp_kernel_op_126 +
              tmp_kernel_op_107 * tmp_kernel_op_127 +
              tmp_kernel_op_108 * tmp_kernel_op_128 +
              tmp_kernel_op_109 * tmp_kernel_op_129;
          const walberla::float64 elMatVec_0 =
              src_dof_0 *
                  (tmp_kernel_op_1 * (tmp_kernel_op_2 * tmp_kernel_op_2) *
                       (tmp_kernel_op_5 * tmp_kernel_op_5) +
                   (tmp_kernel_op_10 * tmp_kernel_op_10) *
                       (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_9 +
                   tmp_kernel_op_17 * (tmp_kernel_op_18 * tmp_kernel_op_18) *
                       (tmp_kernel_op_21 * tmp_kernel_op_21) +
                   tmp_kernel_op_25 * (tmp_kernel_op_26 * tmp_kernel_op_26) *
                       (tmp_kernel_op_29 * tmp_kernel_op_29)) +
              src_dof_1 * tmp_kernel_op_41 + src_dof_2 * tmp_kernel_op_46 +
              src_dof_3 * tmp_kernel_op_55 + src_dof_4 * tmp_kernel_op_32 +
              src_dof_5 * tmp_kernel_op_56 + src_dof_6 * tmp_kernel_op_65;
          const walberla::float64 elMatVec_1 =
              src_dof_0 * tmp_kernel_op_41 +
              src_dof_1 *
                  ((tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_67 +
                   (tmp_kernel_op_20 * tmp_kernel_op_20) * tmp_kernel_op_68 +
                   (tmp_kernel_op_28 * tmp_kernel_op_28) * tmp_kernel_op_69 +
                   (tmp_kernel_op_4 * tmp_kernel_op_4) * tmp_kernel_op_66) +
              src_dof_2 * tmp_kernel_op_79 + src_dof_3 * tmp_kernel_op_81 +
              src_dof_4 * tmp_kernel_op_74 + src_dof_5 * tmp_kernel_op_80 +
              src_dof_6 * tmp_kernel_op_86;
          const walberla::float64 elMatVec_2 =
              src_dof_0 * tmp_kernel_op_46 + src_dof_1 * tmp_kernel_op_79 +
              src_dof_2 *
                  ((tmp_kernel_op_42 * tmp_kernel_op_42) * tmp_kernel_op_88 +
                   (tmp_kernel_op_43 * tmp_kernel_op_43) * tmp_kernel_op_90 +
                   (tmp_kernel_op_44 * tmp_kernel_op_44) * tmp_kernel_op_92 +
                   (tmp_kernel_op_45 * tmp_kernel_op_45) * tmp_kernel_op_94) +
              src_dof_3 * tmp_kernel_op_104 + src_dof_4 * tmp_kernel_op_99 +
              src_dof_5 * tmp_kernel_op_105 + src_dof_6 * tmp_kernel_op_110;
          const walberla::float64 elMatVec_3 =
              src_dof_0 * tmp_kernel_op_55 + src_dof_1 * tmp_kernel_op_81 +
              src_dof_2 * tmp_kernel_op_104 +
              src_dof_3 * (tmp_kernel_op_111 * tmp_kernel_op_66 +
                           tmp_kernel_op_112 * tmp_kernel_op_67 +
                           tmp_kernel_op_113 * tmp_kernel_op_68 +
                           tmp_kernel_op_114 * tmp_kernel_op_69) +
              src_dof_4 * tmp_kernel_op_120 + src_dof_5 * tmp_kernel_op_115 +
              src_dof_6 * tmp_kernel_op_121;
          const walberla::float64 elMatVec_4 =
              src_dof_0 * tmp_kernel_op_32 + src_dof_1 * tmp_kernel_op_74 +
              src_dof_2 * tmp_kernel_op_99 + src_dof_3 * tmp_kernel_op_120 +
              src_dof_4 *
                  (tmp_kernel_op_122 * 16.0 + tmp_kernel_op_123 * 16.0 +
                   tmp_kernel_op_124 * 16.0 + tmp_kernel_op_125 * 16.0) +
              src_dof_5 * tmp_kernel_op_130 + src_dof_6 * tmp_kernel_op_131;
          const walberla::float64 elMatVec_5 =
              src_dof_0 * tmp_kernel_op_56 + src_dof_1 * tmp_kernel_op_80 +
              src_dof_2 * tmp_kernel_op_105 + src_dof_3 * tmp_kernel_op_115 +
              src_dof_4 * tmp_kernel_op_130 +
              src_dof_5 * (tmp_kernel_op_111 * tmp_kernel_op_88 +
                           tmp_kernel_op_112 * tmp_kernel_op_90 +
                           tmp_kernel_op_113 * tmp_kernel_op_92 +
                           tmp_kernel_op_114 * tmp_kernel_op_94) +
              src_dof_6 * tmp_kernel_op_132;
          const walberla::float64 elMatVec_6 =
              src_dof_0 * tmp_kernel_op_65 + src_dof_1 * tmp_kernel_op_86 +
              src_dof_2 * tmp_kernel_op_110 + src_dof_3 * tmp_kernel_op_121 +
              src_dof_4 * tmp_kernel_op_131 + src_dof_5 * tmp_kernel_op_132 +
              src_dof_6 *
                  (tmp_kernel_op_122 * (tmp_kernel_op_57 * tmp_kernel_op_57) +
                   tmp_kernel_op_123 * (tmp_kernel_op_59 * tmp_kernel_op_59) +
                   tmp_kernel_op_124 * (tmp_kernel_op_61 * tmp_kernel_op_61) +
                   tmp_kernel_op_125 * (tmp_kernel_op_63 * tmp_kernel_op_63));
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
          const walberla::float64 tmp_blending_op_0 =
              -rayVertex_0 + thrVertex_0;
          const walberla::float64 tmp_blending_op_1 =
              -p_affine_0_1 + p_affine_1_1;
          const walberla::float64 tmp_blending_op_2 =
              -p_affine_0_1 + p_affine_2_1;
          const walberla::float64 tmp_blending_op_3 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_4 =
              -p_affine_0_0 + p_affine_1_0;
          const walberla::float64 tmp_blending_op_5 =
              -p_affine_0_0 + p_affine_2_0;
          const walberla::float64 tmp_blending_op_6 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_7 =
              (tmp_blending_op_6 * tmp_blending_op_6);
          const walberla::float64 tmp_blending_op_8 =
              (tmp_blending_op_3 * tmp_blending_op_3);
          const walberla::float64 tmp_blending_op_9 =
              tmp_blending_op_7 + tmp_blending_op_8;
          const walberla::float64 tmp_blending_op_10 =
              -rayVertex_1 + thrVertex_1;
          const walberla::float64 tmp_blending_op_11 =
              (-radRayVertex + radRefVertex) * 1.0 /
              (-tmp_blending_op_0 * (-rayVertex_1 + refVertex_1) +
               tmp_blending_op_10 * (-rayVertex_0 + refVertex_0));
          const walberla::float64 tmp_blending_op_12 = tmp_blending_op_11 * 1.0;
          const walberla::float64 tmp_blending_op_13 =
              tmp_blending_op_12 * pow(tmp_blending_op_9, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_14 =
              tmp_blending_op_13 * tmp_blending_op_3;
          const walberla::float64 tmp_blending_op_15 =
              pow(tmp_blending_op_9, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_16 = -rayVertex_1;
          const walberla::float64 tmp_blending_op_17 = -rayVertex_0;
          const walberla::float64 tmp_blending_op_18 =
              radRayVertex + tmp_blending_op_11 *
                                 (-tmp_blending_op_0 *
                                      (tmp_blending_op_16 + tmp_blending_op_3) +
                                  tmp_blending_op_10 *
                                      (tmp_blending_op_17 + tmp_blending_op_6));
          const walberla::float64 tmp_blending_op_19 =
              tmp_blending_op_15 * tmp_blending_op_18 * 1.0;
          const walberla::float64 tmp_blending_op_20 =
              tmp_blending_op_13 * tmp_blending_op_6;
          const walberla::float64 tmp_blending_op_21 =
              p_affine_0_1 + tmp_blending_op_1 * 0.59999999999999998 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_22 =
              p_affine_0_0 + tmp_blending_op_4 * 0.59999999999999998 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_23 =
              (tmp_blending_op_22 * tmp_blending_op_22);
          const walberla::float64 tmp_blending_op_24 =
              (tmp_blending_op_21 * tmp_blending_op_21);
          const walberla::float64 tmp_blending_op_25 =
              tmp_blending_op_23 + tmp_blending_op_24;
          const walberla::float64 tmp_blending_op_26 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_25, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_27 =
              tmp_blending_op_21 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_28 =
              pow(tmp_blending_op_25, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_29 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_21) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_22));
          const walberla::float64 tmp_blending_op_30 =
              tmp_blending_op_28 * tmp_blending_op_29 * 1.0;
          const walberla::float64 tmp_blending_op_31 =
              tmp_blending_op_22 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_32 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_33 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_34 =
              (tmp_blending_op_33 * tmp_blending_op_33);
          const walberla::float64 tmp_blending_op_35 =
              (tmp_blending_op_32 * tmp_blending_op_32);
          const walberla::float64 tmp_blending_op_36 =
              tmp_blending_op_34 + tmp_blending_op_35;
          const walberla::float64 tmp_blending_op_37 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_36, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_38 =
              tmp_blending_op_32 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_39 =
              pow(tmp_blending_op_36, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_40 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_32) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_33));
          const walberla::float64 tmp_blending_op_41 =
              tmp_blending_op_39 * tmp_blending_op_40 * 1.0;
          const walberla::float64 tmp_blending_op_42 =
              tmp_blending_op_33 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_43 =
              p_affine_0_1 + tmp_blending_op_1 * 0.33333333333333331 +
              tmp_blending_op_2 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_44 =
              p_affine_0_0 + tmp_blending_op_4 * 0.33333333333333331 +
              tmp_blending_op_5 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_45 =
              (tmp_blending_op_44 * tmp_blending_op_44);
          const walberla::float64 tmp_blending_op_46 =
              (tmp_blending_op_43 * tmp_blending_op_43);
          const walberla::float64 tmp_blending_op_47 =
              tmp_blending_op_45 + tmp_blending_op_46;
          const walberla::float64 tmp_blending_op_48 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_47, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_49 =
              tmp_blending_op_43 * tmp_blending_op_48;
          const walberla::float64 tmp_blending_op_50 =
              pow(tmp_blending_op_47, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_51 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_43) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_44));
          const walberla::float64 tmp_blending_op_52 =
              tmp_blending_op_50 * tmp_blending_op_51 * 1.0;
          const walberla::float64 tmp_blending_op_53 =
              tmp_blending_op_44 * tmp_blending_op_48;
          const walberla::float64 jac_blending_1_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_14 +
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_7 * 1.0;
          const walberla::float64 jac_blending_1_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_14 -
              tmp_blending_op_19 * tmp_blending_op_3 * tmp_blending_op_6;
          const walberla::float64 jac_blending_0_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_20 -
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_3 *
                  tmp_blending_op_6;
          const walberla::float64 jac_blending_0_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_20 +
              tmp_blending_op_19 * tmp_blending_op_8;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_27 +
              tmp_blending_op_23 * tmp_blending_op_28 * tmp_blending_op_29 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_27 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_30;
          const walberla::float64 jac_blending_0_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_31 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_28 *
                  tmp_blending_op_29;
          const walberla::float64 jac_blending_0_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_31 +
              tmp_blending_op_24 * tmp_blending_op_30;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_38 +
              tmp_blending_op_34 * tmp_blending_op_39 * tmp_blending_op_40 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_38 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_41;
          const walberla::float64 jac_blending_0_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_42 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_39 *
                  tmp_blending_op_40;
          const walberla::float64 jac_blending_0_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_42 +
              tmp_blending_op_35 * tmp_blending_op_41;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_49 +
              tmp_blending_op_45 * tmp_blending_op_50 * tmp_blending_op_51 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_49 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_52;
          const walberla::float64 jac_blending_0_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_53 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_50 *
                  tmp_blending_op_51;
          const walberla::float64 jac_blending_0_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_53 +
              tmp_blending_op_46 * tmp_blending_op_52;
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
          const walberla::float64 tmp_kernel_op_0 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_1 =
              abs_det_jac_affine_BLUE * abs_det_jac_blending_q_0 * -0.28125;
          const walberla::float64 tmp_kernel_op_2 = 0.33333333333333343;
          const walberla::float64 tmp_kernel_op_3 = 0.66666666666666663;
          const walberla::float64 tmp_kernel_op_4 = -0.33333333333333337;
          const walberla::float64 tmp_kernel_op_5 =
              -tmp_kernel_op_3 - tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_6 =
              tmp_kernel_op_2 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_7 =
              tmp_kernel_op_1 * tmp_kernel_op_6 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_8 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_1 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_10 = 0.20000000000000007;
          const walberla::float64 tmp_kernel_op_11 = 1.2;
          const walberla::float64 tmp_kernel_op_12 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_15 =
              tmp_kernel_op_14 * tmp_kernel_op_9 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_16 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_17 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_2 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_18 = 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_19 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_20 = 0.19999999999999996;
          const walberla::float64 tmp_kernel_op_21 =
              -tmp_kernel_op_19 - tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_18 * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_17 * tmp_kernel_op_22 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_24 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_25 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_3 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_26 = 0.60000000000000009;
          const walberla::float64 tmp_kernel_op_27 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_28 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_29 =
              -tmp_kernel_op_27 - tmp_kernel_op_28;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_26 * tmp_kernel_op_29;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_25 * tmp_kernel_op_30 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_0 * tmp_kernel_op_7 +
              tmp_kernel_op_15 * tmp_kernel_op_8 +
              tmp_kernel_op_16 * tmp_kernel_op_23 +
              tmp_kernel_op_24 * tmp_kernel_op_31;
          const walberla::float64 tmp_kernel_op_33 =
              tmp_kernel_op_1 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_34 =
              tmp_kernel_op_33 * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_35 =
              tmp_kernel_op_9 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_36 =
              tmp_kernel_op_14 * tmp_kernel_op_35;
          const walberla::float64 tmp_kernel_op_37 =
              tmp_kernel_op_17 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_22 * tmp_kernel_op_37;
          const walberla::float64 tmp_kernel_op_39 =
              tmp_kernel_op_25 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_40 =
              tmp_kernel_op_30 * tmp_kernel_op_39;
          const walberla::float64 tmp_kernel_op_41 =
              tmp_kernel_op_12 * tmp_kernel_op_36 +
              tmp_kernel_op_20 * tmp_kernel_op_38 +
              tmp_kernel_op_28 * tmp_kernel_op_40 +
              tmp_kernel_op_34 * tmp_kernel_op_4;
          const walberla::float64 tmp_kernel_op_42 = tmp_kernel_op_3 - 1.0;
          const walberla::float64 tmp_kernel_op_43 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_44 = tmp_kernel_op_19 - 1.0;
          const walberla::float64 tmp_kernel_op_45 = tmp_kernel_op_27 - 1.0;
          const walberla::float64 tmp_kernel_op_46 =
              tmp_kernel_op_15 * tmp_kernel_op_43 +
              tmp_kernel_op_23 * tmp_kernel_op_44 +
              tmp_kernel_op_31 * tmp_kernel_op_45 +
              tmp_kernel_op_42 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_47 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_48 =
              -tmp_kernel_op_0 - tmp_kernel_op_47 + 4.0;
          const walberla::float64 tmp_kernel_op_49 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_50 =
              -tmp_kernel_op_49 - tmp_kernel_op_8 + 4.0;
          const walberla::float64 tmp_kernel_op_51 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_52 =
              -tmp_kernel_op_16 - tmp_kernel_op_51 + 4.0;
          const walberla::float64 tmp_kernel_op_53 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_54 =
              -tmp_kernel_op_24 - tmp_kernel_op_53 + 4.0;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_34 * tmp_kernel_op_48 +
              tmp_kernel_op_36 * tmp_kernel_op_50 +
              tmp_kernel_op_38 * tmp_kernel_op_52 +
              tmp_kernel_op_40 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_56 =
              tmp_kernel_op_15 * tmp_kernel_op_50 +
              tmp_kernel_op_23 * tmp_kernel_op_52 +
              tmp_kernel_op_31 * tmp_kernel_op_54 +
              tmp_kernel_op_48 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_57 = 9.0;
          const walberla::float64 tmp_kernel_op_58 =
              tmp_kernel_op_33 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_59 = 5.4000000000000021;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_35 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_61 = 5.4000000000000004;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_37 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_63 = 16.200000000000003;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_39 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_14 * tmp_kernel_op_59 * tmp_kernel_op_60 +
              tmp_kernel_op_22 * tmp_kernel_op_61 * tmp_kernel_op_62 +
              tmp_kernel_op_30 * tmp_kernel_op_63 * tmp_kernel_op_64 +
              tmp_kernel_op_57 * tmp_kernel_op_58 * tmp_kernel_op_6;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_1 * 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_9 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_17 * 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_25 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_4 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_12 * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_72 =
              tmp_kernel_op_20 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_28 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_74 =
              tmp_kernel_op_47 * tmp_kernel_op_70 +
              tmp_kernel_op_49 * tmp_kernel_op_71 +
              tmp_kernel_op_51 * tmp_kernel_op_72 +
              tmp_kernel_op_53 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_75 =
              tmp_kernel_op_4 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_12 * tmp_kernel_op_60;
          const walberla::float64 tmp_kernel_op_77 =
              tmp_kernel_op_20 * tmp_kernel_op_62;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_28 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_42 * tmp_kernel_op_75 +
              tmp_kernel_op_43 * tmp_kernel_op_76 +
              tmp_kernel_op_44 * tmp_kernel_op_77 +
              tmp_kernel_op_45 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_48 * tmp_kernel_op_75 +
              tmp_kernel_op_50 * tmp_kernel_op_76 +
              tmp_kernel_op_52 * tmp_kernel_op_77 +
              tmp_kernel_op_54 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_81 =
              tmp_kernel_op_48 * tmp_kernel_op_70 +
              tmp_kernel_op_50 * tmp_kernel_op_71 +
              tmp_kernel_op_52 * tmp_kernel_op_72 +
              tmp_kernel_op_54 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_83 =
              tmp_kernel_op_59 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_61 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_70 * tmp_kernel_op_82 +
              tmp_kernel_op_71 * tmp_kernel_op_83 +
              tmp_kernel_op_72 * tmp_kernel_op_84 +
              tmp_kernel_op_73 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_87 = 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_1 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 = 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_89 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_91 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_17 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_93 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_25 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              tmp_kernel_op_42 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_43 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_44 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_45 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_99 =
              tmp_kernel_op_0 * tmp_kernel_op_95 +
              tmp_kernel_op_16 * tmp_kernel_op_97 +
              tmp_kernel_op_24 * tmp_kernel_op_98 +
              tmp_kernel_op_8 * tmp_kernel_op_96;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_42 * tmp_kernel_op_48;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_43 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_102 =
              tmp_kernel_op_44 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_103 =
              tmp_kernel_op_45 * tmp_kernel_op_54;
          const walberla::float64 tmp_kernel_op_104 =
              tmp_kernel_op_100 * tmp_kernel_op_58 +
              tmp_kernel_op_101 * tmp_kernel_op_60 +
              tmp_kernel_op_102 * tmp_kernel_op_62 +
              tmp_kernel_op_103 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_105 =
              tmp_kernel_op_100 * tmp_kernel_op_88 +
              tmp_kernel_op_101 * tmp_kernel_op_90 +
              tmp_kernel_op_102 * tmp_kernel_op_92 +
              tmp_kernel_op_103 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_107 =
              tmp_kernel_op_59 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_108 =
              tmp_kernel_op_61 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_109 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_106 * tmp_kernel_op_95 +
              tmp_kernel_op_107 * tmp_kernel_op_96 +
              tmp_kernel_op_108 * tmp_kernel_op_97 +
              tmp_kernel_op_109 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_111 =
              (tmp_kernel_op_48 * tmp_kernel_op_48);
          const walberla::float64 tmp_kernel_op_112 =
              (tmp_kernel_op_50 * tmp_kernel_op_50);
          const walberla::float64 tmp_kernel_op_113 =
              (tmp_kernel_op_52 * tmp_kernel_op_52);
          const walberla::float64 tmp_kernel_op_114 =
              (tmp_kernel_op_54 * tmp_kernel_op_54);
          const walberla::float64 tmp_kernel_op_115 =
              tmp_kernel_op_111 * tmp_kernel_op_58 +
              tmp_kernel_op_112 * tmp_kernel_op_60 +
              tmp_kernel_op_113 * tmp_kernel_op_62 +
              tmp_kernel_op_114 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_116 =
              tmp_kernel_op_48 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_117 =
              tmp_kernel_op_50 * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_118 =
              tmp_kernel_op_52 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_119 =
              tmp_kernel_op_54 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_120 =
              tmp_kernel_op_116 * tmp_kernel_op_47 +
              tmp_kernel_op_117 * tmp_kernel_op_49 +
              tmp_kernel_op_118 * tmp_kernel_op_51 +
              tmp_kernel_op_119 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_121 =
              tmp_kernel_op_116 * tmp_kernel_op_82 +
              tmp_kernel_op_117 * tmp_kernel_op_83 +
              tmp_kernel_op_118 * tmp_kernel_op_84 +
              tmp_kernel_op_119 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_122 =
              tmp_kernel_op_66 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_123 =
              tmp_kernel_op_67 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_124 =
              tmp_kernel_op_68 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_125 =
              tmp_kernel_op_69 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_126 =
              tmp_kernel_op_48 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_127 =
              tmp_kernel_op_50 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_128 =
              tmp_kernel_op_52 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_54 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_130 =
              tmp_kernel_op_0 * tmp_kernel_op_126 +
              tmp_kernel_op_127 * tmp_kernel_op_8 +
              tmp_kernel_op_128 * tmp_kernel_op_16 +
              tmp_kernel_op_129 * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_131 =
              tmp_kernel_op_122 * tmp_kernel_op_57 * 4.0 +
              tmp_kernel_op_123 * tmp_kernel_op_59 * 4.0 +
              tmp_kernel_op_124 * tmp_kernel_op_61 * 4.0 +
              tmp_kernel_op_125 * tmp_kernel_op_63 * 4.0;
          const walberla::float64 tmp_kernel_op_132 =
              tmp_kernel_op_106 * tmp_kernel_op_126 +
              tmp_kernel_op_107 * tmp_kernel_op_127 +
              tmp_kernel_op_108 * tmp_kernel_op_128 +
              tmp_kernel_op_109 * tmp_kernel_op_129;
          const walberla::float64 elMatVec_0 =
              src_dof_0 *
                  (tmp_kernel_op_1 * (tmp_kernel_op_2 * tmp_kernel_op_2) *
                       (tmp_kernel_op_5 * tmp_kernel_op_5) +
                   (tmp_kernel_op_10 * tmp_kernel_op_10) *
                       (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_9 +
                   tmp_kernel_op_17 * (tmp_kernel_op_18 * tmp_kernel_op_18) *
                       (tmp_kernel_op_21 * tmp_kernel_op_21) +
                   tmp_kernel_op_25 * (tmp_kernel_op_26 * tmp_kernel_op_26) *
                       (tmp_kernel_op_29 * tmp_kernel_op_29)) +
              src_dof_1 * tmp_kernel_op_41 + src_dof_2 * tmp_kernel_op_46 +
              src_dof_3 * tmp_kernel_op_55 + src_dof_4 * tmp_kernel_op_32 +
              src_dof_5 * tmp_kernel_op_56 + src_dof_6 * tmp_kernel_op_65;
          const walberla::float64 elMatVec_1 =
              src_dof_0 * tmp_kernel_op_41 +
              src_dof_1 *
                  ((tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_67 +
                   (tmp_kernel_op_20 * tmp_kernel_op_20) * tmp_kernel_op_68 +
                   (tmp_kernel_op_28 * tmp_kernel_op_28) * tmp_kernel_op_69 +
                   (tmp_kernel_op_4 * tmp_kernel_op_4) * tmp_kernel_op_66) +
              src_dof_2 * tmp_kernel_op_79 + src_dof_3 * tmp_kernel_op_81 +
              src_dof_4 * tmp_kernel_op_74 + src_dof_5 * tmp_kernel_op_80 +
              src_dof_6 * tmp_kernel_op_86;
          const walberla::float64 elMatVec_2 =
              src_dof_0 * tmp_kernel_op_46 + src_dof_1 * tmp_kernel_op_79 +
              src_dof_2 *
                  ((tmp_kernel_op_42 * tmp_kernel_op_42) * tmp_kernel_op_88 +
                   (tmp_kernel_op_43 * tmp_kernel_op_43) * tmp_kernel_op_90 +
                   (tmp_kernel_op_44 * tmp_kernel_op_44) * tmp_kernel_op_92 +
                   (tmp_kernel_op_45 * tmp_kernel_op_45) * tmp_kernel_op_94) +
              src_dof_3 * tmp_kernel_op_104 + src_dof_4 * tmp_kernel_op_99 +
              src_dof_5 * tmp_kernel_op_105 + src_dof_6 * tmp_kernel_op_110;
          const walberla::float64 elMatVec_3 =
              src_dof_0 * tmp_kernel_op_55 + src_dof_1 * tmp_kernel_op_81 +
              src_dof_2 * tmp_kernel_op_104 +
              src_dof_3 * (tmp_kernel_op_111 * tmp_kernel_op_66 +
                           tmp_kernel_op_112 * tmp_kernel_op_67 +
                           tmp_kernel_op_113 * tmp_kernel_op_68 +
                           tmp_kernel_op_114 * tmp_kernel_op_69) +
              src_dof_4 * tmp_kernel_op_120 + src_dof_5 * tmp_kernel_op_115 +
              src_dof_6 * tmp_kernel_op_121;
          const walberla::float64 elMatVec_4 =
              src_dof_0 * tmp_kernel_op_32 + src_dof_1 * tmp_kernel_op_74 +
              src_dof_2 * tmp_kernel_op_99 + src_dof_3 * tmp_kernel_op_120 +
              src_dof_4 *
                  (tmp_kernel_op_122 * 16.0 + tmp_kernel_op_123 * 16.0 +
                   tmp_kernel_op_124 * 16.0 + tmp_kernel_op_125 * 16.0) +
              src_dof_5 * tmp_kernel_op_130 + src_dof_6 * tmp_kernel_op_131;
          const walberla::float64 elMatVec_5 =
              src_dof_0 * tmp_kernel_op_56 + src_dof_1 * tmp_kernel_op_80 +
              src_dof_2 * tmp_kernel_op_105 + src_dof_3 * tmp_kernel_op_115 +
              src_dof_4 * tmp_kernel_op_130 +
              src_dof_5 * (tmp_kernel_op_111 * tmp_kernel_op_88 +
                           tmp_kernel_op_112 * tmp_kernel_op_90 +
                           tmp_kernel_op_113 * tmp_kernel_op_92 +
                           tmp_kernel_op_114 * tmp_kernel_op_94) +
              src_dof_6 * tmp_kernel_op_132;
          const walberla::float64 elMatVec_6 =
              src_dof_0 * tmp_kernel_op_65 + src_dof_1 * tmp_kernel_op_86 +
              src_dof_2 * tmp_kernel_op_110 + src_dof_3 * tmp_kernel_op_121 +
              src_dof_4 * tmp_kernel_op_131 + src_dof_5 * tmp_kernel_op_132 +
              src_dof_6 *
                  (tmp_kernel_op_122 * (tmp_kernel_op_57 * tmp_kernel_op_57) +
                   tmp_kernel_op_123 * (tmp_kernel_op_59 * tmp_kernel_op_59) +
                   tmp_kernel_op_124 * (tmp_kernel_op_61 * tmp_kernel_op_61) +
                   tmp_kernel_op_125 * (tmp_kernel_op_63 * tmp_kernel_op_63));
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
void P2PlusBubbleElementwiseMass_AnnulusMap_float64::
    toMatrix_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(
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
        walberla::float64 micro_edges_per_macro_edge_float,
        walberla::float64 radRayVertex, walberla::float64 radRefVertex,
        walberla::float64 rayVertex_0, walberla::float64 rayVertex_1,
        walberla::float64 refVertex_0, walberla::float64 refVertex_1,
        walberla::float64 thrVertex_0, walberla::float64 thrVertex_1) const {
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
          const walberla::float64 tmp_blending_op_0 =
              -rayVertex_0 + thrVertex_0;
          const walberla::float64 tmp_blending_op_1 =
              -p_affine_0_1 + p_affine_1_1;
          const walberla::float64 tmp_blending_op_2 =
              -p_affine_0_1 + p_affine_2_1;
          const walberla::float64 tmp_blending_op_3 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_4 =
              -p_affine_0_0 + p_affine_1_0;
          const walberla::float64 tmp_blending_op_5 =
              -p_affine_0_0 + p_affine_2_0;
          const walberla::float64 tmp_blending_op_6 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_7 =
              (tmp_blending_op_6 * tmp_blending_op_6);
          const walberla::float64 tmp_blending_op_8 =
              (tmp_blending_op_3 * tmp_blending_op_3);
          const walberla::float64 tmp_blending_op_9 =
              tmp_blending_op_7 + tmp_blending_op_8;
          const walberla::float64 tmp_blending_op_10 =
              -rayVertex_1 + thrVertex_1;
          const walberla::float64 tmp_blending_op_11 =
              (-radRayVertex + radRefVertex) * 1.0 /
              (-tmp_blending_op_0 * (-rayVertex_1 + refVertex_1) +
               tmp_blending_op_10 * (-rayVertex_0 + refVertex_0));
          const walberla::float64 tmp_blending_op_12 = tmp_blending_op_11 * 1.0;
          const walberla::float64 tmp_blending_op_13 =
              tmp_blending_op_12 * pow(tmp_blending_op_9, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_14 =
              tmp_blending_op_13 * tmp_blending_op_3;
          const walberla::float64 tmp_blending_op_15 =
              pow(tmp_blending_op_9, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_16 = -rayVertex_1;
          const walberla::float64 tmp_blending_op_17 = -rayVertex_0;
          const walberla::float64 tmp_blending_op_18 =
              radRayVertex + tmp_blending_op_11 *
                                 (-tmp_blending_op_0 *
                                      (tmp_blending_op_16 + tmp_blending_op_3) +
                                  tmp_blending_op_10 *
                                      (tmp_blending_op_17 + tmp_blending_op_6));
          const walberla::float64 tmp_blending_op_19 =
              tmp_blending_op_15 * tmp_blending_op_18 * 1.0;
          const walberla::float64 tmp_blending_op_20 =
              tmp_blending_op_13 * tmp_blending_op_6;
          const walberla::float64 tmp_blending_op_21 =
              p_affine_0_1 + tmp_blending_op_1 * 0.59999999999999998 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_22 =
              p_affine_0_0 + tmp_blending_op_4 * 0.59999999999999998 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_23 =
              (tmp_blending_op_22 * tmp_blending_op_22);
          const walberla::float64 tmp_blending_op_24 =
              (tmp_blending_op_21 * tmp_blending_op_21);
          const walberla::float64 tmp_blending_op_25 =
              tmp_blending_op_23 + tmp_blending_op_24;
          const walberla::float64 tmp_blending_op_26 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_25, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_27 =
              tmp_blending_op_21 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_28 =
              pow(tmp_blending_op_25, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_29 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_21) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_22));
          const walberla::float64 tmp_blending_op_30 =
              tmp_blending_op_28 * tmp_blending_op_29 * 1.0;
          const walberla::float64 tmp_blending_op_31 =
              tmp_blending_op_22 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_32 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_33 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_34 =
              (tmp_blending_op_33 * tmp_blending_op_33);
          const walberla::float64 tmp_blending_op_35 =
              (tmp_blending_op_32 * tmp_blending_op_32);
          const walberla::float64 tmp_blending_op_36 =
              tmp_blending_op_34 + tmp_blending_op_35;
          const walberla::float64 tmp_blending_op_37 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_36, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_38 =
              tmp_blending_op_32 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_39 =
              pow(tmp_blending_op_36, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_40 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_32) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_33));
          const walberla::float64 tmp_blending_op_41 =
              tmp_blending_op_39 * tmp_blending_op_40 * 1.0;
          const walberla::float64 tmp_blending_op_42 =
              tmp_blending_op_33 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_43 =
              p_affine_0_1 + tmp_blending_op_1 * 0.33333333333333331 +
              tmp_blending_op_2 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_44 =
              p_affine_0_0 + tmp_blending_op_4 * 0.33333333333333331 +
              tmp_blending_op_5 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_45 =
              (tmp_blending_op_44 * tmp_blending_op_44);
          const walberla::float64 tmp_blending_op_46 =
              (tmp_blending_op_43 * tmp_blending_op_43);
          const walberla::float64 tmp_blending_op_47 =
              tmp_blending_op_45 + tmp_blending_op_46;
          const walberla::float64 tmp_blending_op_48 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_47, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_49 =
              tmp_blending_op_43 * tmp_blending_op_48;
          const walberla::float64 tmp_blending_op_50 =
              pow(tmp_blending_op_47, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_51 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_43) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_44));
          const walberla::float64 tmp_blending_op_52 =
              tmp_blending_op_50 * tmp_blending_op_51 * 1.0;
          const walberla::float64 tmp_blending_op_53 =
              tmp_blending_op_44 * tmp_blending_op_48;
          const walberla::float64 jac_blending_1_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_14 +
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_7 * 1.0;
          const walberla::float64 jac_blending_1_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_14 -
              tmp_blending_op_19 * tmp_blending_op_3 * tmp_blending_op_6;
          const walberla::float64 jac_blending_0_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_20 -
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_3 *
                  tmp_blending_op_6;
          const walberla::float64 jac_blending_0_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_20 +
              tmp_blending_op_19 * tmp_blending_op_8;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_27 +
              tmp_blending_op_23 * tmp_blending_op_28 * tmp_blending_op_29 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_27 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_30;
          const walberla::float64 jac_blending_0_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_31 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_28 *
                  tmp_blending_op_29;
          const walberla::float64 jac_blending_0_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_31 +
              tmp_blending_op_24 * tmp_blending_op_30;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_38 +
              tmp_blending_op_34 * tmp_blending_op_39 * tmp_blending_op_40 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_38 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_41;
          const walberla::float64 jac_blending_0_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_42 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_39 *
                  tmp_blending_op_40;
          const walberla::float64 jac_blending_0_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_42 +
              tmp_blending_op_35 * tmp_blending_op_41;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_49 +
              tmp_blending_op_45 * tmp_blending_op_50 * tmp_blending_op_51 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_49 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_52;
          const walberla::float64 jac_blending_0_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_53 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_50 *
                  tmp_blending_op_51;
          const walberla::float64 jac_blending_0_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_53 +
              tmp_blending_op_46 * tmp_blending_op_52;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = 0.33333333333333343;
          const walberla::float64 tmp_kernel_op_1 = 0.66666666666666663;
          const walberla::float64 tmp_kernel_op_2 = -0.33333333333333337;
          const walberla::float64 tmp_kernel_op_3 =
              -tmp_kernel_op_1 - tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 =
              abs_det_jac_affine_GRAY * abs_det_jac_blending_q_0 * -0.28125;
          const walberla::float64 tmp_kernel_op_5 = 0.20000000000000007;
          const walberla::float64 tmp_kernel_op_6 = 1.2;
          const walberla::float64 tmp_kernel_op_7 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_8 =
              -tmp_kernel_op_6 - tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_1 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_10 = 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_11 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_12 = 0.19999999999999996;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_2 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_15 = 0.60000000000000009;
          const walberla::float64 tmp_kernel_op_16 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_17 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_18 =
              -tmp_kernel_op_16 - tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_19 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_3 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_20 =
              tmp_kernel_op_4 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_21 =
              tmp_kernel_op_0 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_20 * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_9 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_24 =
              tmp_kernel_op_5 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_25 =
              tmp_kernel_op_23 * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_26 =
              tmp_kernel_op_14 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_27 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_28 =
              tmp_kernel_op_26 * tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_29 =
              tmp_kernel_op_19 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_15 * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_29 * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_12 * tmp_kernel_op_28 +
              tmp_kernel_op_17 * tmp_kernel_op_31 +
              tmp_kernel_op_2 * tmp_kernel_op_22 +
              tmp_kernel_op_25 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_33 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_34 =
              tmp_kernel_op_21 * tmp_kernel_op_4 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_35 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_36 =
              tmp_kernel_op_24 * tmp_kernel_op_9 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_37 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_14 * tmp_kernel_op_27 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_39 = tmp_kernel_op_16 - 1.0;
          const walberla::float64 tmp_kernel_op_40 =
              tmp_kernel_op_19 * tmp_kernel_op_30 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_41 =
              tmp_kernel_op_33 * tmp_kernel_op_34 +
              tmp_kernel_op_35 * tmp_kernel_op_36 +
              tmp_kernel_op_37 * tmp_kernel_op_38 +
              tmp_kernel_op_39 * tmp_kernel_op_40;
          const walberla::float64 tmp_kernel_op_42 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_43 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_44 =
              -tmp_kernel_op_42 - tmp_kernel_op_43 + 4.0;
          const walberla::float64 tmp_kernel_op_45 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_46 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_47 =
              -tmp_kernel_op_45 - tmp_kernel_op_46 + 4.0;
          const walberla::float64 tmp_kernel_op_48 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_49 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_50 =
              -tmp_kernel_op_48 - tmp_kernel_op_49 + 4.0;
          const walberla::float64 tmp_kernel_op_51 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_52 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_53 =
              -tmp_kernel_op_51 - tmp_kernel_op_52 + 4.0;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_22 * tmp_kernel_op_44 +
              tmp_kernel_op_25 * tmp_kernel_op_47 +
              tmp_kernel_op_28 * tmp_kernel_op_50 +
              tmp_kernel_op_31 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_34 * tmp_kernel_op_42 +
              tmp_kernel_op_36 * tmp_kernel_op_45 +
              tmp_kernel_op_38 * tmp_kernel_op_48 +
              tmp_kernel_op_40 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_56 =
              tmp_kernel_op_34 * tmp_kernel_op_44 +
              tmp_kernel_op_36 * tmp_kernel_op_47 +
              tmp_kernel_op_38 * tmp_kernel_op_50 +
              tmp_kernel_op_40 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_57 = 9.0;
          const walberla::float64 tmp_kernel_op_58 =
              tmp_kernel_op_20 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_59 = 5.4000000000000021;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_23 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_61 = 5.4000000000000004;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_26 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_63 = 16.200000000000003;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_29 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_21 * tmp_kernel_op_57 * tmp_kernel_op_58 +
              tmp_kernel_op_24 * tmp_kernel_op_59 * tmp_kernel_op_60 +
              tmp_kernel_op_27 * tmp_kernel_op_61 * tmp_kernel_op_62 +
              tmp_kernel_op_30 * tmp_kernel_op_63 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_4 * 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_9 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_14 * 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_19 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_2 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_60 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_72 =
              tmp_kernel_op_12 * tmp_kernel_op_62;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_17 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_74 =
              tmp_kernel_op_33 * tmp_kernel_op_70 +
              tmp_kernel_op_35 * tmp_kernel_op_71 +
              tmp_kernel_op_37 * tmp_kernel_op_72 +
              tmp_kernel_op_39 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_75 =
              tmp_kernel_op_2 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_67 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_77 =
              tmp_kernel_op_12 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_17 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_44 * tmp_kernel_op_75 +
              tmp_kernel_op_47 * tmp_kernel_op_76 +
              tmp_kernel_op_50 * tmp_kernel_op_77 +
              tmp_kernel_op_53 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_43 * tmp_kernel_op_75 +
              tmp_kernel_op_46 * tmp_kernel_op_76 +
              tmp_kernel_op_49 * tmp_kernel_op_77 +
              tmp_kernel_op_52 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_81 =
              tmp_kernel_op_44 * tmp_kernel_op_70 +
              tmp_kernel_op_47 * tmp_kernel_op_71 +
              tmp_kernel_op_50 * tmp_kernel_op_72 +
              tmp_kernel_op_53 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_83 =
              tmp_kernel_op_59 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_61 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_75 * tmp_kernel_op_82 +
              tmp_kernel_op_76 * tmp_kernel_op_83 +
              tmp_kernel_op_77 * tmp_kernel_op_84 +
              tmp_kernel_op_78 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_87 = 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_4 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 = 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_89 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_91 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_14 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_93 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_19 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              tmp_kernel_op_33 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_35 * tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_37 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_39 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_99 =
              tmp_kernel_op_58 * tmp_kernel_op_95 +
              tmp_kernel_op_60 * tmp_kernel_op_96 +
              tmp_kernel_op_62 * tmp_kernel_op_97 +
              tmp_kernel_op_64 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_33 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_35 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_102 =
              tmp_kernel_op_37 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_103 =
              tmp_kernel_op_39 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_104 =
              tmp_kernel_op_100 * tmp_kernel_op_42 +
              tmp_kernel_op_101 * tmp_kernel_op_45 +
              tmp_kernel_op_102 * tmp_kernel_op_48 +
              tmp_kernel_op_103 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_105 =
              tmp_kernel_op_88 * tmp_kernel_op_95 +
              tmp_kernel_op_90 * tmp_kernel_op_96 +
              tmp_kernel_op_92 * tmp_kernel_op_97 +
              tmp_kernel_op_94 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_107 =
              tmp_kernel_op_59 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_108 =
              tmp_kernel_op_61 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_109 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_100 * tmp_kernel_op_106 +
              tmp_kernel_op_101 * tmp_kernel_op_107 +
              tmp_kernel_op_102 * tmp_kernel_op_108 +
              tmp_kernel_op_103 * tmp_kernel_op_109;
          const walberla::float64 tmp_kernel_op_111 =
              (tmp_kernel_op_44 * tmp_kernel_op_44);
          const walberla::float64 tmp_kernel_op_112 =
              (tmp_kernel_op_47 * tmp_kernel_op_47);
          const walberla::float64 tmp_kernel_op_113 =
              (tmp_kernel_op_50 * tmp_kernel_op_50);
          const walberla::float64 tmp_kernel_op_114 =
              (tmp_kernel_op_53 * tmp_kernel_op_53);
          const walberla::float64 tmp_kernel_op_115 =
              tmp_kernel_op_44 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_116 =
              tmp_kernel_op_47 * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_117 =
              tmp_kernel_op_50 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_118 =
              tmp_kernel_op_53 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_119 =
              tmp_kernel_op_115 * tmp_kernel_op_43 +
              tmp_kernel_op_116 * tmp_kernel_op_46 +
              tmp_kernel_op_117 * tmp_kernel_op_49 +
              tmp_kernel_op_118 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_120 =
              tmp_kernel_op_111 * tmp_kernel_op_58 +
              tmp_kernel_op_112 * tmp_kernel_op_60 +
              tmp_kernel_op_113 * tmp_kernel_op_62 +
              tmp_kernel_op_114 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_121 =
              tmp_kernel_op_115 * tmp_kernel_op_82 +
              tmp_kernel_op_116 * tmp_kernel_op_83 +
              tmp_kernel_op_117 * tmp_kernel_op_84 +
              tmp_kernel_op_118 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_122 =
              tmp_kernel_op_66 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_123 =
              tmp_kernel_op_67 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_124 =
              tmp_kernel_op_68 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_125 =
              tmp_kernel_op_69 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_126 =
              tmp_kernel_op_44 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_127 =
              tmp_kernel_op_47 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_128 =
              tmp_kernel_op_50 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_53 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_130 =
              tmp_kernel_op_126 * tmp_kernel_op_42 +
              tmp_kernel_op_127 * tmp_kernel_op_45 +
              tmp_kernel_op_128 * tmp_kernel_op_48 +
              tmp_kernel_op_129 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_131 =
              tmp_kernel_op_122 * tmp_kernel_op_57 * 4.0 +
              tmp_kernel_op_123 * tmp_kernel_op_59 * 4.0 +
              tmp_kernel_op_124 * tmp_kernel_op_61 * 4.0 +
              tmp_kernel_op_125 * tmp_kernel_op_63 * 4.0;
          const walberla::float64 tmp_kernel_op_132 =
              tmp_kernel_op_106 * tmp_kernel_op_126 +
              tmp_kernel_op_107 * tmp_kernel_op_127 +
              tmp_kernel_op_108 * tmp_kernel_op_128 +
              tmp_kernel_op_109 * tmp_kernel_op_129;
          const walberla::float64 elMat_0_0 =
              (tmp_kernel_op_0 * tmp_kernel_op_0) *
                  (tmp_kernel_op_3 * tmp_kernel_op_3) * tmp_kernel_op_4 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) *
                  (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_14 +
              (tmp_kernel_op_15 * tmp_kernel_op_15) *
                  (tmp_kernel_op_18 * tmp_kernel_op_18) * tmp_kernel_op_19 +
              (tmp_kernel_op_5 * tmp_kernel_op_5) *
                  (tmp_kernel_op_8 * tmp_kernel_op_8) * tmp_kernel_op_9;
          const walberla::float64 elMat_0_1 = tmp_kernel_op_32;
          const walberla::float64 elMat_0_2 = tmp_kernel_op_41;
          const walberla::float64 elMat_0_3 = tmp_kernel_op_54;
          const walberla::float64 elMat_0_4 = tmp_kernel_op_55;
          const walberla::float64 elMat_0_5 = tmp_kernel_op_56;
          const walberla::float64 elMat_0_6 = tmp_kernel_op_65;
          const walberla::float64 elMat_1_0 = tmp_kernel_op_32;
          const walberla::float64 elMat_1_1 =
              (tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_68 +
              (tmp_kernel_op_17 * tmp_kernel_op_17) * tmp_kernel_op_69 +
              (tmp_kernel_op_2 * tmp_kernel_op_2) * tmp_kernel_op_66 +
              tmp_kernel_op_67 * (tmp_kernel_op_7 * tmp_kernel_op_7);
          const walberla::float64 elMat_1_2 = tmp_kernel_op_74;
          const walberla::float64 elMat_1_3 = tmp_kernel_op_79;
          const walberla::float64 elMat_1_4 = tmp_kernel_op_80;
          const walberla::float64 elMat_1_5 = tmp_kernel_op_81;
          const walberla::float64 elMat_1_6 = tmp_kernel_op_86;
          const walberla::float64 elMat_2_0 = tmp_kernel_op_41;
          const walberla::float64 elMat_2_1 = tmp_kernel_op_74;
          const walberla::float64 elMat_2_2 =
              (tmp_kernel_op_33 * tmp_kernel_op_33) * tmp_kernel_op_88 +
              (tmp_kernel_op_35 * tmp_kernel_op_35) * tmp_kernel_op_90 +
              (tmp_kernel_op_37 * tmp_kernel_op_37) * tmp_kernel_op_92 +
              (tmp_kernel_op_39 * tmp_kernel_op_39) * tmp_kernel_op_94;
          const walberla::float64 elMat_2_3 = tmp_kernel_op_99;
          const walberla::float64 elMat_2_4 = tmp_kernel_op_104;
          const walberla::float64 elMat_2_5 = tmp_kernel_op_105;
          const walberla::float64 elMat_2_6 = tmp_kernel_op_110;
          const walberla::float64 elMat_3_0 = tmp_kernel_op_54;
          const walberla::float64 elMat_3_1 = tmp_kernel_op_79;
          const walberla::float64 elMat_3_2 = tmp_kernel_op_99;
          const walberla::float64 elMat_3_3 =
              tmp_kernel_op_111 * tmp_kernel_op_66 +
              tmp_kernel_op_112 * tmp_kernel_op_67 +
              tmp_kernel_op_113 * tmp_kernel_op_68 +
              tmp_kernel_op_114 * tmp_kernel_op_69;
          const walberla::float64 elMat_3_4 = tmp_kernel_op_119;
          const walberla::float64 elMat_3_5 = tmp_kernel_op_120;
          const walberla::float64 elMat_3_6 = tmp_kernel_op_121;
          const walberla::float64 elMat_4_0 = tmp_kernel_op_55;
          const walberla::float64 elMat_4_1 = tmp_kernel_op_80;
          const walberla::float64 elMat_4_2 = tmp_kernel_op_104;
          const walberla::float64 elMat_4_3 = tmp_kernel_op_119;
          const walberla::float64 elMat_4_4 =
              tmp_kernel_op_122 * 16.0 + tmp_kernel_op_123 * 16.0 +
              tmp_kernel_op_124 * 16.0 + tmp_kernel_op_125 * 16.0;
          const walberla::float64 elMat_4_5 = tmp_kernel_op_130;
          const walberla::float64 elMat_4_6 = tmp_kernel_op_131;
          const walberla::float64 elMat_5_0 = tmp_kernel_op_56;
          const walberla::float64 elMat_5_1 = tmp_kernel_op_81;
          const walberla::float64 elMat_5_2 = tmp_kernel_op_105;
          const walberla::float64 elMat_5_3 = tmp_kernel_op_120;
          const walberla::float64 elMat_5_4 = tmp_kernel_op_130;
          const walberla::float64 elMat_5_5 =
              tmp_kernel_op_111 * tmp_kernel_op_88 +
              tmp_kernel_op_112 * tmp_kernel_op_90 +
              tmp_kernel_op_113 * tmp_kernel_op_92 +
              tmp_kernel_op_114 * tmp_kernel_op_94;
          const walberla::float64 elMat_5_6 = tmp_kernel_op_132;
          const walberla::float64 elMat_6_0 = tmp_kernel_op_65;
          const walberla::float64 elMat_6_1 = tmp_kernel_op_86;
          const walberla::float64 elMat_6_2 = tmp_kernel_op_110;
          const walberla::float64 elMat_6_3 = tmp_kernel_op_121;
          const walberla::float64 elMat_6_4 = tmp_kernel_op_131;
          const walberla::float64 elMat_6_5 = tmp_kernel_op_132;
          const walberla::float64 elMat_6_6 =
              tmp_kernel_op_122 * (tmp_kernel_op_57 * tmp_kernel_op_57) +
              tmp_kernel_op_123 * (tmp_kernel_op_59 * tmp_kernel_op_59) +
              tmp_kernel_op_124 * (tmp_kernel_op_61 * tmp_kernel_op_61) +
              tmp_kernel_op_125 * (tmp_kernel_op_63 * tmp_kernel_op_63);

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
          const walberla::float64 tmp_blending_op_0 =
              -rayVertex_0 + thrVertex_0;
          const walberla::float64 tmp_blending_op_1 =
              -p_affine_0_1 + p_affine_1_1;
          const walberla::float64 tmp_blending_op_2 =
              -p_affine_0_1 + p_affine_2_1;
          const walberla::float64 tmp_blending_op_3 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_4 =
              -p_affine_0_0 + p_affine_1_0;
          const walberla::float64 tmp_blending_op_5 =
              -p_affine_0_0 + p_affine_2_0;
          const walberla::float64 tmp_blending_op_6 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_7 =
              (tmp_blending_op_6 * tmp_blending_op_6);
          const walberla::float64 tmp_blending_op_8 =
              (tmp_blending_op_3 * tmp_blending_op_3);
          const walberla::float64 tmp_blending_op_9 =
              tmp_blending_op_7 + tmp_blending_op_8;
          const walberla::float64 tmp_blending_op_10 =
              -rayVertex_1 + thrVertex_1;
          const walberla::float64 tmp_blending_op_11 =
              (-radRayVertex + radRefVertex) * 1.0 /
              (-tmp_blending_op_0 * (-rayVertex_1 + refVertex_1) +
               tmp_blending_op_10 * (-rayVertex_0 + refVertex_0));
          const walberla::float64 tmp_blending_op_12 = tmp_blending_op_11 * 1.0;
          const walberla::float64 tmp_blending_op_13 =
              tmp_blending_op_12 * pow(tmp_blending_op_9, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_14 =
              tmp_blending_op_13 * tmp_blending_op_3;
          const walberla::float64 tmp_blending_op_15 =
              pow(tmp_blending_op_9, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_16 = -rayVertex_1;
          const walberla::float64 tmp_blending_op_17 = -rayVertex_0;
          const walberla::float64 tmp_blending_op_18 =
              radRayVertex + tmp_blending_op_11 *
                                 (-tmp_blending_op_0 *
                                      (tmp_blending_op_16 + tmp_blending_op_3) +
                                  tmp_blending_op_10 *
                                      (tmp_blending_op_17 + tmp_blending_op_6));
          const walberla::float64 tmp_blending_op_19 =
              tmp_blending_op_15 * tmp_blending_op_18 * 1.0;
          const walberla::float64 tmp_blending_op_20 =
              tmp_blending_op_13 * tmp_blending_op_6;
          const walberla::float64 tmp_blending_op_21 =
              p_affine_0_1 + tmp_blending_op_1 * 0.59999999999999998 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_22 =
              p_affine_0_0 + tmp_blending_op_4 * 0.59999999999999998 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_23 =
              (tmp_blending_op_22 * tmp_blending_op_22);
          const walberla::float64 tmp_blending_op_24 =
              (tmp_blending_op_21 * tmp_blending_op_21);
          const walberla::float64 tmp_blending_op_25 =
              tmp_blending_op_23 + tmp_blending_op_24;
          const walberla::float64 tmp_blending_op_26 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_25, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_27 =
              tmp_blending_op_21 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_28 =
              pow(tmp_blending_op_25, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_29 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_21) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_22));
          const walberla::float64 tmp_blending_op_30 =
              tmp_blending_op_28 * tmp_blending_op_29 * 1.0;
          const walberla::float64 tmp_blending_op_31 =
              tmp_blending_op_22 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_32 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_33 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_34 =
              (tmp_blending_op_33 * tmp_blending_op_33);
          const walberla::float64 tmp_blending_op_35 =
              (tmp_blending_op_32 * tmp_blending_op_32);
          const walberla::float64 tmp_blending_op_36 =
              tmp_blending_op_34 + tmp_blending_op_35;
          const walberla::float64 tmp_blending_op_37 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_36, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_38 =
              tmp_blending_op_32 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_39 =
              pow(tmp_blending_op_36, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_40 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_32) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_33));
          const walberla::float64 tmp_blending_op_41 =
              tmp_blending_op_39 * tmp_blending_op_40 * 1.0;
          const walberla::float64 tmp_blending_op_42 =
              tmp_blending_op_33 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_43 =
              p_affine_0_1 + tmp_blending_op_1 * 0.33333333333333331 +
              tmp_blending_op_2 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_44 =
              p_affine_0_0 + tmp_blending_op_4 * 0.33333333333333331 +
              tmp_blending_op_5 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_45 =
              (tmp_blending_op_44 * tmp_blending_op_44);
          const walberla::float64 tmp_blending_op_46 =
              (tmp_blending_op_43 * tmp_blending_op_43);
          const walberla::float64 tmp_blending_op_47 =
              tmp_blending_op_45 + tmp_blending_op_46;
          const walberla::float64 tmp_blending_op_48 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_47, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_49 =
              tmp_blending_op_43 * tmp_blending_op_48;
          const walberla::float64 tmp_blending_op_50 =
              pow(tmp_blending_op_47, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_51 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_43) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_44));
          const walberla::float64 tmp_blending_op_52 =
              tmp_blending_op_50 * tmp_blending_op_51 * 1.0;
          const walberla::float64 tmp_blending_op_53 =
              tmp_blending_op_44 * tmp_blending_op_48;
          const walberla::float64 jac_blending_1_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_14 +
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_7 * 1.0;
          const walberla::float64 jac_blending_1_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_14 -
              tmp_blending_op_19 * tmp_blending_op_3 * tmp_blending_op_6;
          const walberla::float64 jac_blending_0_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_20 -
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_3 *
                  tmp_blending_op_6;
          const walberla::float64 jac_blending_0_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_20 +
              tmp_blending_op_19 * tmp_blending_op_8;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_27 +
              tmp_blending_op_23 * tmp_blending_op_28 * tmp_blending_op_29 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_27 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_30;
          const walberla::float64 jac_blending_0_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_31 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_28 *
                  tmp_blending_op_29;
          const walberla::float64 jac_blending_0_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_31 +
              tmp_blending_op_24 * tmp_blending_op_30;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_38 +
              tmp_blending_op_34 * tmp_blending_op_39 * tmp_blending_op_40 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_38 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_41;
          const walberla::float64 jac_blending_0_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_42 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_39 *
                  tmp_blending_op_40;
          const walberla::float64 jac_blending_0_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_42 +
              tmp_blending_op_35 * tmp_blending_op_41;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_49 +
              tmp_blending_op_45 * tmp_blending_op_50 * tmp_blending_op_51 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_49 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_52;
          const walberla::float64 jac_blending_0_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_53 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_50 *
                  tmp_blending_op_51;
          const walberla::float64 jac_blending_0_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_53 +
              tmp_blending_op_46 * tmp_blending_op_52;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = 0.33333333333333343;
          const walberla::float64 tmp_kernel_op_1 = 0.66666666666666663;
          const walberla::float64 tmp_kernel_op_2 = -0.33333333333333337;
          const walberla::float64 tmp_kernel_op_3 =
              -tmp_kernel_op_1 - tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_4 =
              abs_det_jac_affine_BLUE * abs_det_jac_blending_q_0 * -0.28125;
          const walberla::float64 tmp_kernel_op_5 = 0.20000000000000007;
          const walberla::float64 tmp_kernel_op_6 = 1.2;
          const walberla::float64 tmp_kernel_op_7 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_8 =
              -tmp_kernel_op_6 - tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_9 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_1 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_10 = 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_11 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_12 = 0.19999999999999996;
          const walberla::float64 tmp_kernel_op_13 =
              -tmp_kernel_op_11 - tmp_kernel_op_12;
          const walberla::float64 tmp_kernel_op_14 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_2 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_15 = 0.60000000000000009;
          const walberla::float64 tmp_kernel_op_16 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_17 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_18 =
              -tmp_kernel_op_16 - tmp_kernel_op_17;
          const walberla::float64 tmp_kernel_op_19 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_3 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_20 =
              tmp_kernel_op_4 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_21 =
              tmp_kernel_op_0 * tmp_kernel_op_3;
          const walberla::float64 tmp_kernel_op_22 =
              tmp_kernel_op_20 * tmp_kernel_op_21;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_9 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_24 =
              tmp_kernel_op_5 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_25 =
              tmp_kernel_op_23 * tmp_kernel_op_24;
          const walberla::float64 tmp_kernel_op_26 =
              tmp_kernel_op_14 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_27 =
              tmp_kernel_op_10 * tmp_kernel_op_13;
          const walberla::float64 tmp_kernel_op_28 =
              tmp_kernel_op_26 * tmp_kernel_op_27;
          const walberla::float64 tmp_kernel_op_29 =
              tmp_kernel_op_19 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_15 * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_29 * tmp_kernel_op_30;
          const walberla::float64 tmp_kernel_op_32 =
              tmp_kernel_op_12 * tmp_kernel_op_28 +
              tmp_kernel_op_17 * tmp_kernel_op_31 +
              tmp_kernel_op_2 * tmp_kernel_op_22 +
              tmp_kernel_op_25 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_33 = tmp_kernel_op_1 - 1.0;
          const walberla::float64 tmp_kernel_op_34 =
              tmp_kernel_op_21 * tmp_kernel_op_4 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_35 = tmp_kernel_op_6 - 1.0;
          const walberla::float64 tmp_kernel_op_36 =
              tmp_kernel_op_24 * tmp_kernel_op_9 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_37 = tmp_kernel_op_11 - 1.0;
          const walberla::float64 tmp_kernel_op_38 =
              tmp_kernel_op_14 * tmp_kernel_op_27 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_39 = tmp_kernel_op_16 - 1.0;
          const walberla::float64 tmp_kernel_op_40 =
              tmp_kernel_op_19 * tmp_kernel_op_30 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_41 =
              tmp_kernel_op_33 * tmp_kernel_op_34 +
              tmp_kernel_op_35 * tmp_kernel_op_36 +
              tmp_kernel_op_37 * tmp_kernel_op_38 +
              tmp_kernel_op_39 * tmp_kernel_op_40;
          const walberla::float64 tmp_kernel_op_42 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_43 = 1.3333333333333333;
          const walberla::float64 tmp_kernel_op_44 =
              -tmp_kernel_op_42 - tmp_kernel_op_43 + 4.0;
          const walberla::float64 tmp_kernel_op_45 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_46 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_47 =
              -tmp_kernel_op_45 - tmp_kernel_op_46 + 4.0;
          const walberla::float64 tmp_kernel_op_48 = 2.3999999999999999;
          const walberla::float64 tmp_kernel_op_49 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_50 =
              -tmp_kernel_op_48 - tmp_kernel_op_49 + 4.0;
          const walberla::float64 tmp_kernel_op_51 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_52 = 0.80000000000000004;
          const walberla::float64 tmp_kernel_op_53 =
              -tmp_kernel_op_51 - tmp_kernel_op_52 + 4.0;
          const walberla::float64 tmp_kernel_op_54 =
              tmp_kernel_op_22 * tmp_kernel_op_44 +
              tmp_kernel_op_25 * tmp_kernel_op_47 +
              tmp_kernel_op_28 * tmp_kernel_op_50 +
              tmp_kernel_op_31 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_55 =
              tmp_kernel_op_34 * tmp_kernel_op_42 +
              tmp_kernel_op_36 * tmp_kernel_op_45 +
              tmp_kernel_op_38 * tmp_kernel_op_48 +
              tmp_kernel_op_40 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_56 =
              tmp_kernel_op_34 * tmp_kernel_op_44 +
              tmp_kernel_op_36 * tmp_kernel_op_47 +
              tmp_kernel_op_38 * tmp_kernel_op_50 +
              tmp_kernel_op_40 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_57 = 9.0;
          const walberla::float64 tmp_kernel_op_58 =
              tmp_kernel_op_20 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_59 = 5.4000000000000021;
          const walberla::float64 tmp_kernel_op_60 =
              tmp_kernel_op_23 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_61 = 5.4000000000000004;
          const walberla::float64 tmp_kernel_op_62 =
              tmp_kernel_op_26 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_63 = 16.200000000000003;
          const walberla::float64 tmp_kernel_op_64 =
              tmp_kernel_op_29 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_65 =
              tmp_kernel_op_21 * tmp_kernel_op_57 * tmp_kernel_op_58 +
              tmp_kernel_op_24 * tmp_kernel_op_59 * tmp_kernel_op_60 +
              tmp_kernel_op_27 * tmp_kernel_op_61 * tmp_kernel_op_62 +
              tmp_kernel_op_30 * tmp_kernel_op_63 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_66 =
              tmp_kernel_op_4 * 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_67 =
              tmp_kernel_op_9 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_68 =
              tmp_kernel_op_14 * 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_69 =
              tmp_kernel_op_19 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_70 =
              tmp_kernel_op_2 * tmp_kernel_op_58;
          const walberla::float64 tmp_kernel_op_71 =
              tmp_kernel_op_60 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_72 =
              tmp_kernel_op_12 * tmp_kernel_op_62;
          const walberla::float64 tmp_kernel_op_73 =
              tmp_kernel_op_17 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_74 =
              tmp_kernel_op_33 * tmp_kernel_op_70 +
              tmp_kernel_op_35 * tmp_kernel_op_71 +
              tmp_kernel_op_37 * tmp_kernel_op_72 +
              tmp_kernel_op_39 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_75 =
              tmp_kernel_op_2 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_76 =
              tmp_kernel_op_67 * tmp_kernel_op_7;
          const walberla::float64 tmp_kernel_op_77 =
              tmp_kernel_op_12 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_78 =
              tmp_kernel_op_17 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_79 =
              tmp_kernel_op_44 * tmp_kernel_op_75 +
              tmp_kernel_op_47 * tmp_kernel_op_76 +
              tmp_kernel_op_50 * tmp_kernel_op_77 +
              tmp_kernel_op_53 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_80 =
              tmp_kernel_op_43 * tmp_kernel_op_75 +
              tmp_kernel_op_46 * tmp_kernel_op_76 +
              tmp_kernel_op_49 * tmp_kernel_op_77 +
              tmp_kernel_op_52 * tmp_kernel_op_78;
          const walberla::float64 tmp_kernel_op_81 =
              tmp_kernel_op_44 * tmp_kernel_op_70 +
              tmp_kernel_op_47 * tmp_kernel_op_71 +
              tmp_kernel_op_50 * tmp_kernel_op_72 +
              tmp_kernel_op_53 * tmp_kernel_op_73;
          const walberla::float64 tmp_kernel_op_82 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_83 =
              tmp_kernel_op_59 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_84 =
              tmp_kernel_op_61 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_85 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_86 =
              tmp_kernel_op_75 * tmp_kernel_op_82 +
              tmp_kernel_op_76 * tmp_kernel_op_83 +
              tmp_kernel_op_77 * tmp_kernel_op_84 +
              tmp_kernel_op_78 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_87 = 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_88 =
              tmp_kernel_op_4 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_89 = 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_90 =
              tmp_kernel_op_89 * tmp_kernel_op_9;
          const walberla::float64 tmp_kernel_op_91 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_92 =
              tmp_kernel_op_14 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_93 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_94 =
              tmp_kernel_op_19 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_95 =
              tmp_kernel_op_33 * tmp_kernel_op_44;
          const walberla::float64 tmp_kernel_op_96 =
              tmp_kernel_op_35 * tmp_kernel_op_47;
          const walberla::float64 tmp_kernel_op_97 =
              tmp_kernel_op_37 * tmp_kernel_op_50;
          const walberla::float64 tmp_kernel_op_98 =
              tmp_kernel_op_39 * tmp_kernel_op_53;
          const walberla::float64 tmp_kernel_op_99 =
              tmp_kernel_op_58 * tmp_kernel_op_95 +
              tmp_kernel_op_60 * tmp_kernel_op_96 +
              tmp_kernel_op_62 * tmp_kernel_op_97 +
              tmp_kernel_op_64 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_100 =
              tmp_kernel_op_33 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_101 =
              tmp_kernel_op_35 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_102 =
              tmp_kernel_op_37 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_103 =
              tmp_kernel_op_39 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_104 =
              tmp_kernel_op_100 * tmp_kernel_op_42 +
              tmp_kernel_op_101 * tmp_kernel_op_45 +
              tmp_kernel_op_102 * tmp_kernel_op_48 +
              tmp_kernel_op_103 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_105 =
              tmp_kernel_op_88 * tmp_kernel_op_95 +
              tmp_kernel_op_90 * tmp_kernel_op_96 +
              tmp_kernel_op_92 * tmp_kernel_op_97 +
              tmp_kernel_op_94 * tmp_kernel_op_98;
          const walberla::float64 tmp_kernel_op_106 =
              tmp_kernel_op_57 * 0.33333333333333331;
          const walberla::float64 tmp_kernel_op_107 =
              tmp_kernel_op_59 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_108 =
              tmp_kernel_op_61 * 0.59999999999999998;
          const walberla::float64 tmp_kernel_op_109 =
              tmp_kernel_op_63 * 0.20000000000000001;
          const walberla::float64 tmp_kernel_op_110 =
              tmp_kernel_op_100 * tmp_kernel_op_106 +
              tmp_kernel_op_101 * tmp_kernel_op_107 +
              tmp_kernel_op_102 * tmp_kernel_op_108 +
              tmp_kernel_op_103 * tmp_kernel_op_109;
          const walberla::float64 tmp_kernel_op_111 =
              (tmp_kernel_op_44 * tmp_kernel_op_44);
          const walberla::float64 tmp_kernel_op_112 =
              (tmp_kernel_op_47 * tmp_kernel_op_47);
          const walberla::float64 tmp_kernel_op_113 =
              (tmp_kernel_op_50 * tmp_kernel_op_50);
          const walberla::float64 tmp_kernel_op_114 =
              (tmp_kernel_op_53 * tmp_kernel_op_53);
          const walberla::float64 tmp_kernel_op_115 =
              tmp_kernel_op_44 * tmp_kernel_op_66;
          const walberla::float64 tmp_kernel_op_116 =
              tmp_kernel_op_47 * tmp_kernel_op_67;
          const walberla::float64 tmp_kernel_op_117 =
              tmp_kernel_op_50 * tmp_kernel_op_68;
          const walberla::float64 tmp_kernel_op_118 =
              tmp_kernel_op_53 * tmp_kernel_op_69;
          const walberla::float64 tmp_kernel_op_119 =
              tmp_kernel_op_115 * tmp_kernel_op_43 +
              tmp_kernel_op_116 * tmp_kernel_op_46 +
              tmp_kernel_op_117 * tmp_kernel_op_49 +
              tmp_kernel_op_118 * tmp_kernel_op_52;
          const walberla::float64 tmp_kernel_op_120 =
              tmp_kernel_op_111 * tmp_kernel_op_58 +
              tmp_kernel_op_112 * tmp_kernel_op_60 +
              tmp_kernel_op_113 * tmp_kernel_op_62 +
              tmp_kernel_op_114 * tmp_kernel_op_64;
          const walberla::float64 tmp_kernel_op_121 =
              tmp_kernel_op_115 * tmp_kernel_op_82 +
              tmp_kernel_op_116 * tmp_kernel_op_83 +
              tmp_kernel_op_117 * tmp_kernel_op_84 +
              tmp_kernel_op_118 * tmp_kernel_op_85;
          const walberla::float64 tmp_kernel_op_122 =
              tmp_kernel_op_66 * tmp_kernel_op_87;
          const walberla::float64 tmp_kernel_op_123 =
              tmp_kernel_op_67 * tmp_kernel_op_89;
          const walberla::float64 tmp_kernel_op_124 =
              tmp_kernel_op_68 * tmp_kernel_op_91;
          const walberla::float64 tmp_kernel_op_125 =
              tmp_kernel_op_69 * tmp_kernel_op_93;
          const walberla::float64 tmp_kernel_op_126 =
              tmp_kernel_op_44 * tmp_kernel_op_88;
          const walberla::float64 tmp_kernel_op_127 =
              tmp_kernel_op_47 * tmp_kernel_op_90;
          const walberla::float64 tmp_kernel_op_128 =
              tmp_kernel_op_50 * tmp_kernel_op_92;
          const walberla::float64 tmp_kernel_op_129 =
              tmp_kernel_op_53 * tmp_kernel_op_94;
          const walberla::float64 tmp_kernel_op_130 =
              tmp_kernel_op_126 * tmp_kernel_op_42 +
              tmp_kernel_op_127 * tmp_kernel_op_45 +
              tmp_kernel_op_128 * tmp_kernel_op_48 +
              tmp_kernel_op_129 * tmp_kernel_op_51;
          const walberla::float64 tmp_kernel_op_131 =
              tmp_kernel_op_122 * tmp_kernel_op_57 * 4.0 +
              tmp_kernel_op_123 * tmp_kernel_op_59 * 4.0 +
              tmp_kernel_op_124 * tmp_kernel_op_61 * 4.0 +
              tmp_kernel_op_125 * tmp_kernel_op_63 * 4.0;
          const walberla::float64 tmp_kernel_op_132 =
              tmp_kernel_op_106 * tmp_kernel_op_126 +
              tmp_kernel_op_107 * tmp_kernel_op_127 +
              tmp_kernel_op_108 * tmp_kernel_op_128 +
              tmp_kernel_op_109 * tmp_kernel_op_129;
          const walberla::float64 elMat_0_0 =
              (tmp_kernel_op_0 * tmp_kernel_op_0) *
                  (tmp_kernel_op_3 * tmp_kernel_op_3) * tmp_kernel_op_4 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) *
                  (tmp_kernel_op_13 * tmp_kernel_op_13) * tmp_kernel_op_14 +
              (tmp_kernel_op_15 * tmp_kernel_op_15) *
                  (tmp_kernel_op_18 * tmp_kernel_op_18) * tmp_kernel_op_19 +
              (tmp_kernel_op_5 * tmp_kernel_op_5) *
                  (tmp_kernel_op_8 * tmp_kernel_op_8) * tmp_kernel_op_9;
          const walberla::float64 elMat_0_1 = tmp_kernel_op_32;
          const walberla::float64 elMat_0_2 = tmp_kernel_op_41;
          const walberla::float64 elMat_0_3 = tmp_kernel_op_54;
          const walberla::float64 elMat_0_4 = tmp_kernel_op_55;
          const walberla::float64 elMat_0_5 = tmp_kernel_op_56;
          const walberla::float64 elMat_0_6 = tmp_kernel_op_65;
          const walberla::float64 elMat_1_0 = tmp_kernel_op_32;
          const walberla::float64 elMat_1_1 =
              (tmp_kernel_op_12 * tmp_kernel_op_12) * tmp_kernel_op_68 +
              (tmp_kernel_op_17 * tmp_kernel_op_17) * tmp_kernel_op_69 +
              (tmp_kernel_op_2 * tmp_kernel_op_2) * tmp_kernel_op_66 +
              tmp_kernel_op_67 * (tmp_kernel_op_7 * tmp_kernel_op_7);
          const walberla::float64 elMat_1_2 = tmp_kernel_op_74;
          const walberla::float64 elMat_1_3 = tmp_kernel_op_79;
          const walberla::float64 elMat_1_4 = tmp_kernel_op_80;
          const walberla::float64 elMat_1_5 = tmp_kernel_op_81;
          const walberla::float64 elMat_1_6 = tmp_kernel_op_86;
          const walberla::float64 elMat_2_0 = tmp_kernel_op_41;
          const walberla::float64 elMat_2_1 = tmp_kernel_op_74;
          const walberla::float64 elMat_2_2 =
              (tmp_kernel_op_33 * tmp_kernel_op_33) * tmp_kernel_op_88 +
              (tmp_kernel_op_35 * tmp_kernel_op_35) * tmp_kernel_op_90 +
              (tmp_kernel_op_37 * tmp_kernel_op_37) * tmp_kernel_op_92 +
              (tmp_kernel_op_39 * tmp_kernel_op_39) * tmp_kernel_op_94;
          const walberla::float64 elMat_2_3 = tmp_kernel_op_99;
          const walberla::float64 elMat_2_4 = tmp_kernel_op_104;
          const walberla::float64 elMat_2_5 = tmp_kernel_op_105;
          const walberla::float64 elMat_2_6 = tmp_kernel_op_110;
          const walberla::float64 elMat_3_0 = tmp_kernel_op_54;
          const walberla::float64 elMat_3_1 = tmp_kernel_op_79;
          const walberla::float64 elMat_3_2 = tmp_kernel_op_99;
          const walberla::float64 elMat_3_3 =
              tmp_kernel_op_111 * tmp_kernel_op_66 +
              tmp_kernel_op_112 * tmp_kernel_op_67 +
              tmp_kernel_op_113 * tmp_kernel_op_68 +
              tmp_kernel_op_114 * tmp_kernel_op_69;
          const walberla::float64 elMat_3_4 = tmp_kernel_op_119;
          const walberla::float64 elMat_3_5 = tmp_kernel_op_120;
          const walberla::float64 elMat_3_6 = tmp_kernel_op_121;
          const walberla::float64 elMat_4_0 = tmp_kernel_op_55;
          const walberla::float64 elMat_4_1 = tmp_kernel_op_80;
          const walberla::float64 elMat_4_2 = tmp_kernel_op_104;
          const walberla::float64 elMat_4_3 = tmp_kernel_op_119;
          const walberla::float64 elMat_4_4 =
              tmp_kernel_op_122 * 16.0 + tmp_kernel_op_123 * 16.0 +
              tmp_kernel_op_124 * 16.0 + tmp_kernel_op_125 * 16.0;
          const walberla::float64 elMat_4_5 = tmp_kernel_op_130;
          const walberla::float64 elMat_4_6 = tmp_kernel_op_131;
          const walberla::float64 elMat_5_0 = tmp_kernel_op_56;
          const walberla::float64 elMat_5_1 = tmp_kernel_op_81;
          const walberla::float64 elMat_5_2 = tmp_kernel_op_105;
          const walberla::float64 elMat_5_3 = tmp_kernel_op_120;
          const walberla::float64 elMat_5_4 = tmp_kernel_op_130;
          const walberla::float64 elMat_5_5 =
              tmp_kernel_op_111 * tmp_kernel_op_88 +
              tmp_kernel_op_112 * tmp_kernel_op_90 +
              tmp_kernel_op_113 * tmp_kernel_op_92 +
              tmp_kernel_op_114 * tmp_kernel_op_94;
          const walberla::float64 elMat_5_6 = tmp_kernel_op_132;
          const walberla::float64 elMat_6_0 = tmp_kernel_op_65;
          const walberla::float64 elMat_6_1 = tmp_kernel_op_86;
          const walberla::float64 elMat_6_2 = tmp_kernel_op_110;
          const walberla::float64 elMat_6_3 = tmp_kernel_op_121;
          const walberla::float64 elMat_6_4 = tmp_kernel_op_131;
          const walberla::float64 elMat_6_5 = tmp_kernel_op_132;
          const walberla::float64 elMat_6_6 =
              tmp_kernel_op_122 * (tmp_kernel_op_57 * tmp_kernel_op_57) +
              tmp_kernel_op_123 * (tmp_kernel_op_59 * tmp_kernel_op_59) +
              tmp_kernel_op_124 * (tmp_kernel_op_61 * tmp_kernel_op_61) +
              tmp_kernel_op_125 * (tmp_kernel_op_63 * tmp_kernel_op_63);

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
void P2PlusBubbleElementwiseMass_AnnulusMap_float64::
    computeInverseDiagonalOperatorValues_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(
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
        walberla::float64 micro_edges_per_macro_edge_float,
        walberla::float64 radRayVertex, walberla::float64 radRefVertex,
        walberla::float64 rayVertex_0, walberla::float64 rayVertex_1,
        walberla::float64 refVertex_0, walberla::float64 refVertex_1,
        walberla::float64 thrVertex_0, walberla::float64 thrVertex_1) const {
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
          const walberla::float64 tmp_blending_op_0 =
              -rayVertex_0 + thrVertex_0;
          const walberla::float64 tmp_blending_op_1 =
              -p_affine_0_1 + p_affine_1_1;
          const walberla::float64 tmp_blending_op_2 =
              -p_affine_0_1 + p_affine_2_1;
          const walberla::float64 tmp_blending_op_3 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_4 =
              -p_affine_0_0 + p_affine_1_0;
          const walberla::float64 tmp_blending_op_5 =
              -p_affine_0_0 + p_affine_2_0;
          const walberla::float64 tmp_blending_op_6 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_7 =
              (tmp_blending_op_6 * tmp_blending_op_6);
          const walberla::float64 tmp_blending_op_8 =
              (tmp_blending_op_3 * tmp_blending_op_3);
          const walberla::float64 tmp_blending_op_9 =
              tmp_blending_op_7 + tmp_blending_op_8;
          const walberla::float64 tmp_blending_op_10 =
              -rayVertex_1 + thrVertex_1;
          const walberla::float64 tmp_blending_op_11 =
              (-radRayVertex + radRefVertex) * 1.0 /
              (-tmp_blending_op_0 * (-rayVertex_1 + refVertex_1) +
               tmp_blending_op_10 * (-rayVertex_0 + refVertex_0));
          const walberla::float64 tmp_blending_op_12 = tmp_blending_op_11 * 1.0;
          const walberla::float64 tmp_blending_op_13 =
              tmp_blending_op_12 * pow(tmp_blending_op_9, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_14 =
              tmp_blending_op_13 * tmp_blending_op_3;
          const walberla::float64 tmp_blending_op_15 =
              pow(tmp_blending_op_9, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_16 = -rayVertex_1;
          const walberla::float64 tmp_blending_op_17 = -rayVertex_0;
          const walberla::float64 tmp_blending_op_18 =
              radRayVertex + tmp_blending_op_11 *
                                 (-tmp_blending_op_0 *
                                      (tmp_blending_op_16 + tmp_blending_op_3) +
                                  tmp_blending_op_10 *
                                      (tmp_blending_op_17 + tmp_blending_op_6));
          const walberla::float64 tmp_blending_op_19 =
              tmp_blending_op_15 * tmp_blending_op_18 * 1.0;
          const walberla::float64 tmp_blending_op_20 =
              tmp_blending_op_13 * tmp_blending_op_6;
          const walberla::float64 tmp_blending_op_21 =
              p_affine_0_1 + tmp_blending_op_1 * 0.59999999999999998 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_22 =
              p_affine_0_0 + tmp_blending_op_4 * 0.59999999999999998 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_23 =
              (tmp_blending_op_22 * tmp_blending_op_22);
          const walberla::float64 tmp_blending_op_24 =
              (tmp_blending_op_21 * tmp_blending_op_21);
          const walberla::float64 tmp_blending_op_25 =
              tmp_blending_op_23 + tmp_blending_op_24;
          const walberla::float64 tmp_blending_op_26 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_25, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_27 =
              tmp_blending_op_21 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_28 =
              pow(tmp_blending_op_25, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_29 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_21) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_22));
          const walberla::float64 tmp_blending_op_30 =
              tmp_blending_op_28 * tmp_blending_op_29 * 1.0;
          const walberla::float64 tmp_blending_op_31 =
              tmp_blending_op_22 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_32 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_33 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_34 =
              (tmp_blending_op_33 * tmp_blending_op_33);
          const walberla::float64 tmp_blending_op_35 =
              (tmp_blending_op_32 * tmp_blending_op_32);
          const walberla::float64 tmp_blending_op_36 =
              tmp_blending_op_34 + tmp_blending_op_35;
          const walberla::float64 tmp_blending_op_37 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_36, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_38 =
              tmp_blending_op_32 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_39 =
              pow(tmp_blending_op_36, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_40 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_32) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_33));
          const walberla::float64 tmp_blending_op_41 =
              tmp_blending_op_39 * tmp_blending_op_40 * 1.0;
          const walberla::float64 tmp_blending_op_42 =
              tmp_blending_op_33 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_43 =
              p_affine_0_1 + tmp_blending_op_1 * 0.33333333333333331 +
              tmp_blending_op_2 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_44 =
              p_affine_0_0 + tmp_blending_op_4 * 0.33333333333333331 +
              tmp_blending_op_5 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_45 =
              (tmp_blending_op_44 * tmp_blending_op_44);
          const walberla::float64 tmp_blending_op_46 =
              (tmp_blending_op_43 * tmp_blending_op_43);
          const walberla::float64 tmp_blending_op_47 =
              tmp_blending_op_45 + tmp_blending_op_46;
          const walberla::float64 tmp_blending_op_48 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_47, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_49 =
              tmp_blending_op_43 * tmp_blending_op_48;
          const walberla::float64 tmp_blending_op_50 =
              pow(tmp_blending_op_47, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_51 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_43) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_44));
          const walberla::float64 tmp_blending_op_52 =
              tmp_blending_op_50 * tmp_blending_op_51 * 1.0;
          const walberla::float64 tmp_blending_op_53 =
              tmp_blending_op_44 * tmp_blending_op_48;
          const walberla::float64 jac_blending_1_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_14 +
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_7 * 1.0;
          const walberla::float64 jac_blending_1_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_14 -
              tmp_blending_op_19 * tmp_blending_op_3 * tmp_blending_op_6;
          const walberla::float64 jac_blending_0_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_20 -
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_3 *
                  tmp_blending_op_6;
          const walberla::float64 jac_blending_0_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_20 +
              tmp_blending_op_19 * tmp_blending_op_8;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_27 +
              tmp_blending_op_23 * tmp_blending_op_28 * tmp_blending_op_29 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_27 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_30;
          const walberla::float64 jac_blending_0_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_31 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_28 *
                  tmp_blending_op_29;
          const walberla::float64 jac_blending_0_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_31 +
              tmp_blending_op_24 * tmp_blending_op_30;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_38 +
              tmp_blending_op_34 * tmp_blending_op_39 * tmp_blending_op_40 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_38 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_41;
          const walberla::float64 jac_blending_0_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_42 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_39 *
                  tmp_blending_op_40;
          const walberla::float64 jac_blending_0_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_42 +
              tmp_blending_op_35 * tmp_blending_op_41;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_49 +
              tmp_blending_op_45 * tmp_blending_op_50 * tmp_blending_op_51 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_49 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_52;
          const walberla::float64 jac_blending_0_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_53 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_50 *
                  tmp_blending_op_51;
          const walberla::float64 jac_blending_0_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_53 +
              tmp_blending_op_46 * tmp_blending_op_52;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = 0.66666666666666663;
          const walberla::float64 tmp_kernel_op_1 = -0.33333333333333337;
          const walberla::float64 tmp_kernel_op_2 =
              abs_det_jac_affine_GRAY * abs_det_jac_blending_q_0 * -0.28125;
          const walberla::float64 tmp_kernel_op_3 = 1.2;
          const walberla::float64 tmp_kernel_op_4 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_5 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_1 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_6 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_7 = 0.19999999999999996;
          const walberla::float64 tmp_kernel_op_8 = abs_det_jac_affine_GRAY *
                                                    abs_det_jac_blending_q_2 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_9 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_10 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_11 = abs_det_jac_affine_GRAY *
                                                     abs_det_jac_blending_q_3 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_12 =
              tmp_kernel_op_2 * 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_13 =
              tmp_kernel_op_5 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_8 * 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_15 =
              tmp_kernel_op_11 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_16 = 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_17 =
              tmp_kernel_op_16 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_18 = 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_19 =
              tmp_kernel_op_18 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_20 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_21 =
              tmp_kernel_op_20 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_22 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_11 * tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_24 = 1.7777777777777788;
          const walberla::float64 tmp_kernel_op_25 = 0.64000000000000046;
          const walberla::float64 tmp_kernel_op_26 = 0.64000000000000012;
          const walberla::float64 tmp_kernel_op_27 = 5.7600000000000016;
          const walberla::float64 tmp_kernel_op_28 =
              tmp_kernel_op_12 * tmp_kernel_op_16;
          const walberla::float64 tmp_kernel_op_29 =
              tmp_kernel_op_13 * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_14 * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_15 * tmp_kernel_op_22;
          const walberla::float64 elMatDiag_0 =
              tmp_kernel_op_11 *
                  ((-tmp_kernel_op_10 - tmp_kernel_op_9) *
                   (-tmp_kernel_op_10 - tmp_kernel_op_9)) *
                  0.3600000000000001 +
              tmp_kernel_op_2 *
                  ((-tmp_kernel_op_0 - tmp_kernel_op_1) *
                   (-tmp_kernel_op_0 - tmp_kernel_op_1)) *
                  0.11111111111111117 +
              tmp_kernel_op_5 *
                  ((-tmp_kernel_op_3 - tmp_kernel_op_4) *
                   (-tmp_kernel_op_3 - tmp_kernel_op_4)) *
                  0.040000000000000029 +
              tmp_kernel_op_8 *
                  ((-tmp_kernel_op_6 - tmp_kernel_op_7) *
                   (-tmp_kernel_op_6 - tmp_kernel_op_7)) *
                  0.040000000000000008;
          const walberla::float64 elMatDiag_1 =
              (tmp_kernel_op_1 * tmp_kernel_op_1) * tmp_kernel_op_12 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) * tmp_kernel_op_15 +
              tmp_kernel_op_13 * (tmp_kernel_op_4 * tmp_kernel_op_4) +
              tmp_kernel_op_14 * (tmp_kernel_op_7 * tmp_kernel_op_7);
          const walberla::float64 elMatDiag_2 =
              tmp_kernel_op_17 *
                  ((tmp_kernel_op_0 - 1.0) * (tmp_kernel_op_0 - 1.0)) +
              tmp_kernel_op_19 *
                  ((tmp_kernel_op_3 - 1.0) * (tmp_kernel_op_3 - 1.0)) +
              tmp_kernel_op_21 *
                  ((tmp_kernel_op_6 - 1.0) * (tmp_kernel_op_6 - 1.0)) +
              tmp_kernel_op_23 *
                  ((tmp_kernel_op_9 - 1.0) * (tmp_kernel_op_9 - 1.0));
          const walberla::float64 elMatDiag_3 =
              tmp_kernel_op_12 * tmp_kernel_op_24 +
              tmp_kernel_op_13 * tmp_kernel_op_25 +
              tmp_kernel_op_14 * tmp_kernel_op_26 +
              tmp_kernel_op_15 * tmp_kernel_op_27;
          const walberla::float64 elMatDiag_4 =
              tmp_kernel_op_28 * 16.0 + tmp_kernel_op_29 * 16.0 +
              tmp_kernel_op_30 * 16.0 + tmp_kernel_op_31 * 16.0;
          const walberla::float64 elMatDiag_5 =
              tmp_kernel_op_17 * tmp_kernel_op_24 +
              tmp_kernel_op_19 * tmp_kernel_op_25 +
              tmp_kernel_op_21 * tmp_kernel_op_26 +
              tmp_kernel_op_23 * tmp_kernel_op_27;
          const walberla::float64 elMatDiag_6 =
              tmp_kernel_op_28 * 81.0 + tmp_kernel_op_29 * 29.160000000000021 +
              tmp_kernel_op_30 * 29.160000000000004 +
              tmp_kernel_op_31 * 262.44000000000011;
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
          const walberla::float64 tmp_blending_op_0 =
              -rayVertex_0 + thrVertex_0;
          const walberla::float64 tmp_blending_op_1 =
              -p_affine_0_1 + p_affine_1_1;
          const walberla::float64 tmp_blending_op_2 =
              -p_affine_0_1 + p_affine_2_1;
          const walberla::float64 tmp_blending_op_3 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_4 =
              -p_affine_0_0 + p_affine_1_0;
          const walberla::float64 tmp_blending_op_5 =
              -p_affine_0_0 + p_affine_2_0;
          const walberla::float64 tmp_blending_op_6 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_7 =
              (tmp_blending_op_6 * tmp_blending_op_6);
          const walberla::float64 tmp_blending_op_8 =
              (tmp_blending_op_3 * tmp_blending_op_3);
          const walberla::float64 tmp_blending_op_9 =
              tmp_blending_op_7 + tmp_blending_op_8;
          const walberla::float64 tmp_blending_op_10 =
              -rayVertex_1 + thrVertex_1;
          const walberla::float64 tmp_blending_op_11 =
              (-radRayVertex + radRefVertex) * 1.0 /
              (-tmp_blending_op_0 * (-rayVertex_1 + refVertex_1) +
               tmp_blending_op_10 * (-rayVertex_0 + refVertex_0));
          const walberla::float64 tmp_blending_op_12 = tmp_blending_op_11 * 1.0;
          const walberla::float64 tmp_blending_op_13 =
              tmp_blending_op_12 * pow(tmp_blending_op_9, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_14 =
              tmp_blending_op_13 * tmp_blending_op_3;
          const walberla::float64 tmp_blending_op_15 =
              pow(tmp_blending_op_9, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_16 = -rayVertex_1;
          const walberla::float64 tmp_blending_op_17 = -rayVertex_0;
          const walberla::float64 tmp_blending_op_18 =
              radRayVertex + tmp_blending_op_11 *
                                 (-tmp_blending_op_0 *
                                      (tmp_blending_op_16 + tmp_blending_op_3) +
                                  tmp_blending_op_10 *
                                      (tmp_blending_op_17 + tmp_blending_op_6));
          const walberla::float64 tmp_blending_op_19 =
              tmp_blending_op_15 * tmp_blending_op_18 * 1.0;
          const walberla::float64 tmp_blending_op_20 =
              tmp_blending_op_13 * tmp_blending_op_6;
          const walberla::float64 tmp_blending_op_21 =
              p_affine_0_1 + tmp_blending_op_1 * 0.59999999999999998 +
              tmp_blending_op_2 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_22 =
              p_affine_0_0 + tmp_blending_op_4 * 0.59999999999999998 +
              tmp_blending_op_5 * 0.20000000000000001;
          const walberla::float64 tmp_blending_op_23 =
              (tmp_blending_op_22 * tmp_blending_op_22);
          const walberla::float64 tmp_blending_op_24 =
              (tmp_blending_op_21 * tmp_blending_op_21);
          const walberla::float64 tmp_blending_op_25 =
              tmp_blending_op_23 + tmp_blending_op_24;
          const walberla::float64 tmp_blending_op_26 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_25, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_27 =
              tmp_blending_op_21 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_28 =
              pow(tmp_blending_op_25, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_29 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_21) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_22));
          const walberla::float64 tmp_blending_op_30 =
              tmp_blending_op_28 * tmp_blending_op_29 * 1.0;
          const walberla::float64 tmp_blending_op_31 =
              tmp_blending_op_22 * tmp_blending_op_26;
          const walberla::float64 tmp_blending_op_32 =
              p_affine_0_1 + tmp_blending_op_1 * 0.20000000000000001 +
              tmp_blending_op_2 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_33 =
              p_affine_0_0 + tmp_blending_op_4 * 0.20000000000000001 +
              tmp_blending_op_5 * 0.59999999999999998;
          const walberla::float64 tmp_blending_op_34 =
              (tmp_blending_op_33 * tmp_blending_op_33);
          const walberla::float64 tmp_blending_op_35 =
              (tmp_blending_op_32 * tmp_blending_op_32);
          const walberla::float64 tmp_blending_op_36 =
              tmp_blending_op_34 + tmp_blending_op_35;
          const walberla::float64 tmp_blending_op_37 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_36, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_38 =
              tmp_blending_op_32 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_39 =
              pow(tmp_blending_op_36, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_40 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_32) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_33));
          const walberla::float64 tmp_blending_op_41 =
              tmp_blending_op_39 * tmp_blending_op_40 * 1.0;
          const walberla::float64 tmp_blending_op_42 =
              tmp_blending_op_33 * tmp_blending_op_37;
          const walberla::float64 tmp_blending_op_43 =
              p_affine_0_1 + tmp_blending_op_1 * 0.33333333333333331 +
              tmp_blending_op_2 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_44 =
              p_affine_0_0 + tmp_blending_op_4 * 0.33333333333333331 +
              tmp_blending_op_5 * 0.33333333333333331;
          const walberla::float64 tmp_blending_op_45 =
              (tmp_blending_op_44 * tmp_blending_op_44);
          const walberla::float64 tmp_blending_op_46 =
              (tmp_blending_op_43 * tmp_blending_op_43);
          const walberla::float64 tmp_blending_op_47 =
              tmp_blending_op_45 + tmp_blending_op_46;
          const walberla::float64 tmp_blending_op_48 =
              tmp_blending_op_12 *
              pow(tmp_blending_op_47, -0.50000000000000000);
          const walberla::float64 tmp_blending_op_49 =
              tmp_blending_op_43 * tmp_blending_op_48;
          const walberla::float64 tmp_blending_op_50 =
              pow(tmp_blending_op_47, -1.5000000000000000);
          const walberla::float64 tmp_blending_op_51 =
              radRayVertex +
              tmp_blending_op_11 * (-tmp_blending_op_0 * (tmp_blending_op_16 +
                                                          tmp_blending_op_43) +
                                    tmp_blending_op_10 * (tmp_blending_op_17 +
                                                          tmp_blending_op_44));
          const walberla::float64 tmp_blending_op_52 =
              tmp_blending_op_50 * tmp_blending_op_51 * 1.0;
          const walberla::float64 tmp_blending_op_53 =
              tmp_blending_op_44 * tmp_blending_op_48;
          const walberla::float64 jac_blending_1_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_14 +
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_7 * 1.0;
          const walberla::float64 jac_blending_1_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_14 -
              tmp_blending_op_19 * tmp_blending_op_3 * tmp_blending_op_6;
          const walberla::float64 jac_blending_0_1_q_3 =
              -tmp_blending_op_0 * tmp_blending_op_20 -
              tmp_blending_op_15 * tmp_blending_op_18 * tmp_blending_op_3 *
                  tmp_blending_op_6;
          const walberla::float64 jac_blending_0_0_q_3 =
              tmp_blending_op_10 * tmp_blending_op_20 +
              tmp_blending_op_19 * tmp_blending_op_8;
          const walberla::float64 abs_det_jac_blending_q_3 =
              jac_blending_0_0_q_3 * jac_blending_1_1_q_3 -
              jac_blending_0_1_q_3 * jac_blending_1_0_q_3;
          const walberla::float64 jac_blending_1_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_27 +
              tmp_blending_op_23 * tmp_blending_op_28 * tmp_blending_op_29 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_27 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_30;
          const walberla::float64 jac_blending_0_1_q_2 =
              -tmp_blending_op_0 * tmp_blending_op_31 -
              tmp_blending_op_21 * tmp_blending_op_22 * tmp_blending_op_28 *
                  tmp_blending_op_29;
          const walberla::float64 jac_blending_0_0_q_2 =
              tmp_blending_op_10 * tmp_blending_op_31 +
              tmp_blending_op_24 * tmp_blending_op_30;
          const walberla::float64 abs_det_jac_blending_q_2 =
              jac_blending_0_0_q_2 * jac_blending_1_1_q_2 -
              jac_blending_0_1_q_2 * jac_blending_1_0_q_2;
          const walberla::float64 jac_blending_1_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_38 +
              tmp_blending_op_34 * tmp_blending_op_39 * tmp_blending_op_40 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_38 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_41;
          const walberla::float64 jac_blending_0_1_q_1 =
              -tmp_blending_op_0 * tmp_blending_op_42 -
              tmp_blending_op_32 * tmp_blending_op_33 * tmp_blending_op_39 *
                  tmp_blending_op_40;
          const walberla::float64 jac_blending_0_0_q_1 =
              tmp_blending_op_10 * tmp_blending_op_42 +
              tmp_blending_op_35 * tmp_blending_op_41;
          const walberla::float64 abs_det_jac_blending_q_1 =
              jac_blending_0_0_q_1 * jac_blending_1_1_q_1 -
              jac_blending_0_1_q_1 * jac_blending_1_0_q_1;
          const walberla::float64 jac_blending_1_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_49 +
              tmp_blending_op_45 * tmp_blending_op_50 * tmp_blending_op_51 *
                  1.0;
          const walberla::float64 jac_blending_1_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_49 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_52;
          const walberla::float64 jac_blending_0_1_q_0 =
              -tmp_blending_op_0 * tmp_blending_op_53 -
              tmp_blending_op_43 * tmp_blending_op_44 * tmp_blending_op_50 *
                  tmp_blending_op_51;
          const walberla::float64 jac_blending_0_0_q_0 =
              tmp_blending_op_10 * tmp_blending_op_53 +
              tmp_blending_op_46 * tmp_blending_op_52;
          const walberla::float64 abs_det_jac_blending_q_0 =
              jac_blending_0_0_q_0 * jac_blending_1_1_q_0 -
              jac_blending_0_1_q_0 * jac_blending_1_0_q_0;
          const walberla::float64 tmp_kernel_op_0 = 0.66666666666666663;
          const walberla::float64 tmp_kernel_op_1 = -0.33333333333333337;
          const walberla::float64 tmp_kernel_op_2 =
              abs_det_jac_affine_BLUE * abs_det_jac_blending_q_0 * -0.28125;
          const walberla::float64 tmp_kernel_op_3 = 1.2;
          const walberla::float64 tmp_kernel_op_4 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_5 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_1 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_6 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_7 = 0.19999999999999996;
          const walberla::float64 tmp_kernel_op_8 = abs_det_jac_affine_BLUE *
                                                    abs_det_jac_blending_q_2 *
                                                    0.26041666666666669;
          const walberla::float64 tmp_kernel_op_9 = 0.40000000000000002;
          const walberla::float64 tmp_kernel_op_10 = -0.59999999999999998;
          const walberla::float64 tmp_kernel_op_11 = abs_det_jac_affine_BLUE *
                                                     abs_det_jac_blending_q_3 *
                                                     0.26041666666666669;
          const walberla::float64 tmp_kernel_op_12 =
              tmp_kernel_op_2 * 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_13 =
              tmp_kernel_op_5 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_14 =
              tmp_kernel_op_8 * 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_15 =
              tmp_kernel_op_11 * 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_16 = 0.1111111111111111;
          const walberla::float64 tmp_kernel_op_17 =
              tmp_kernel_op_16 * tmp_kernel_op_2;
          const walberla::float64 tmp_kernel_op_18 = 0.35999999999999999;
          const walberla::float64 tmp_kernel_op_19 =
              tmp_kernel_op_18 * tmp_kernel_op_5;
          const walberla::float64 tmp_kernel_op_20 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_21 =
              tmp_kernel_op_20 * tmp_kernel_op_8;
          const walberla::float64 tmp_kernel_op_22 = 0.040000000000000008;
          const walberla::float64 tmp_kernel_op_23 =
              tmp_kernel_op_11 * tmp_kernel_op_22;
          const walberla::float64 tmp_kernel_op_24 = 1.7777777777777788;
          const walberla::float64 tmp_kernel_op_25 = 0.64000000000000046;
          const walberla::float64 tmp_kernel_op_26 = 0.64000000000000012;
          const walberla::float64 tmp_kernel_op_27 = 5.7600000000000016;
          const walberla::float64 tmp_kernel_op_28 =
              tmp_kernel_op_12 * tmp_kernel_op_16;
          const walberla::float64 tmp_kernel_op_29 =
              tmp_kernel_op_13 * tmp_kernel_op_18;
          const walberla::float64 tmp_kernel_op_30 =
              tmp_kernel_op_14 * tmp_kernel_op_20;
          const walberla::float64 tmp_kernel_op_31 =
              tmp_kernel_op_15 * tmp_kernel_op_22;
          const walberla::float64 elMatDiag_0 =
              tmp_kernel_op_11 *
                  ((-tmp_kernel_op_10 - tmp_kernel_op_9) *
                   (-tmp_kernel_op_10 - tmp_kernel_op_9)) *
                  0.3600000000000001 +
              tmp_kernel_op_2 *
                  ((-tmp_kernel_op_0 - tmp_kernel_op_1) *
                   (-tmp_kernel_op_0 - tmp_kernel_op_1)) *
                  0.11111111111111117 +
              tmp_kernel_op_5 *
                  ((-tmp_kernel_op_3 - tmp_kernel_op_4) *
                   (-tmp_kernel_op_3 - tmp_kernel_op_4)) *
                  0.040000000000000029 +
              tmp_kernel_op_8 *
                  ((-tmp_kernel_op_6 - tmp_kernel_op_7) *
                   (-tmp_kernel_op_6 - tmp_kernel_op_7)) *
                  0.040000000000000008;
          const walberla::float64 elMatDiag_1 =
              (tmp_kernel_op_1 * tmp_kernel_op_1) * tmp_kernel_op_12 +
              (tmp_kernel_op_10 * tmp_kernel_op_10) * tmp_kernel_op_15 +
              tmp_kernel_op_13 * (tmp_kernel_op_4 * tmp_kernel_op_4) +
              tmp_kernel_op_14 * (tmp_kernel_op_7 * tmp_kernel_op_7);
          const walberla::float64 elMatDiag_2 =
              tmp_kernel_op_17 *
                  ((tmp_kernel_op_0 - 1.0) * (tmp_kernel_op_0 - 1.0)) +
              tmp_kernel_op_19 *
                  ((tmp_kernel_op_3 - 1.0) * (tmp_kernel_op_3 - 1.0)) +
              tmp_kernel_op_21 *
                  ((tmp_kernel_op_6 - 1.0) * (tmp_kernel_op_6 - 1.0)) +
              tmp_kernel_op_23 *
                  ((tmp_kernel_op_9 - 1.0) * (tmp_kernel_op_9 - 1.0));
          const walberla::float64 elMatDiag_3 =
              tmp_kernel_op_12 * tmp_kernel_op_24 +
              tmp_kernel_op_13 * tmp_kernel_op_25 +
              tmp_kernel_op_14 * tmp_kernel_op_26 +
              tmp_kernel_op_15 * tmp_kernel_op_27;
          const walberla::float64 elMatDiag_4 =
              tmp_kernel_op_28 * 16.0 + tmp_kernel_op_29 * 16.0 +
              tmp_kernel_op_30 * 16.0 + tmp_kernel_op_31 * 16.0;
          const walberla::float64 elMatDiag_5 =
              tmp_kernel_op_17 * tmp_kernel_op_24 +
              tmp_kernel_op_19 * tmp_kernel_op_25 +
              tmp_kernel_op_21 * tmp_kernel_op_26 +
              tmp_kernel_op_23 * tmp_kernel_op_27;
          const walberla::float64 elMatDiag_6 =
              tmp_kernel_op_28 * 81.0 + tmp_kernel_op_29 * 29.160000000000021 +
              tmp_kernel_op_30 * 29.160000000000004 +
              tmp_kernel_op_31 * 262.44000000000011;
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
