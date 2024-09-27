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

#include "P1ElementwiseDivKGrad_AnnulusMap_float64.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

P1ElementwiseDivKGrad_AnnulusMap_float64::
    P1ElementwiseDivKGrad_AnnulusMap_float64(
        const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
        size_t maxLevel, const P0Function<walberla::float64> &_k)
    : Operator(storage, minLevel, maxLevel), k(_k) {}

void P1ElementwiseDivKGrad_AnnulusMap_float64::apply(
    const P1Function<walberla::float64> &src,
    const P1Function<walberla::float64> &dst, uint_t level, DoFType flag,
    UpdateType updateType) const {
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
      walberla::float64 *_data_src =
          face.getData(src.getFaceDataID())->getPointer(level);
      walberla::float64 *_data_dst =
          face.getData(dst.getFaceDataID())->getPointer(level);
      walberla::float64 *_data_k =
          k.getDGFunction()->volumeDoFFunction()->dofMemory(it.first, level);

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

      apply_P1ElementwiseDivKGrad_AnnulusMap_float64_macro_2D(

          _data_dst, _data_k, _data_src, macro_vertex_coord_id_0comp0,
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
void P1ElementwiseDivKGrad_AnnulusMap_float64::toMatrix(
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

    WALBERLA_ABORT("Not implemented.");
  } else {
    this->timingTree_->start("pre-communication");

    this->timingTree_->stop("pre-communication");

    for (auto &it : storage_->getFaces()) {
      Face &face = *it.second;

      // get hold of the actual numerical data
      idx_t *_data_src = face.getData(src.getFaceDataID())->getPointer(level);
      idx_t *_data_dst = face.getData(dst.getFaceDataID())->getPointer(level);
      walberla::float64 *_data_k =
          k.getDGFunction()->volumeDoFFunction()->dofMemory(it.first, level);

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

      toMatrix_P1ElementwiseDivKGrad_AnnulusMap_float64_macro_2D(

          _data_dst, _data_k, _data_src, macro_vertex_coord_id_0comp0,
          macro_vertex_coord_id_0comp1, macro_vertex_coord_id_1comp0,
          macro_vertex_coord_id_1comp1, macro_vertex_coord_id_2comp0,
          macro_vertex_coord_id_2comp1, mat, micro_edges_per_macro_edge,
          micro_edges_per_macro_edge_float);

      this->timingTree_->stop("kernel");
    }
  }
  this->stopTiming("toMatrix");
}
void P1ElementwiseDivKGrad_AnnulusMap_float64::
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

      WALBERLA_ABORT("Not implemented.");
      (*invDiag_).invertElementwise(level);
    } else {
      this->timingTree_->start("pre-communication");

      this->timingTree_->stop("pre-communication");

      for (auto &it : storage_->getFaces()) {
        Face &face = *it.second;

        // get hold of the actual numerical data
        walberla::float64 *_data_invDiag_ =
            face.getData((*invDiag_).getFaceDataID())->getPointer(level);
        walberla::float64 *_data_k =
            k.getDGFunction()->volumeDoFFunction()->dofMemory(it.first, level);

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

        computeInverseDiagonalOperatorValues_P1ElementwiseDivKGrad_AnnulusMap_float64_macro_2D(

            _data_invDiag_, _data_k, macro_vertex_coord_id_0comp0,
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
      (*invDiag_).invertElementwise(level);
    }
  }

  this->stopTiming("computeInverseDiagonalOperatorValues");
}
std::shared_ptr<P1Function<walberla::float64>>
P1ElementwiseDivKGrad_AnnulusMap_float64::getInverseDiagonalValues() const {
  return invDiag_;
}
void P1ElementwiseDivKGrad_AnnulusMap_float64::
    apply_P1ElementwiseDivKGrad_AnnulusMap_float64_macro_2D(
        walberla::float64 *RESTRICT _data_dst,
        walberla::float64 *RESTRICT _data_k,
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
    {
      /* FaceType.GRAY */
      const walberla::float64 tmp_coords_jac_0_GRAY =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
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
              _data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 src_dof_1 =
              _data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          const walberla::float64 src_dof_2 =
              _data_src[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          const walberla::float64 k_dof_0 =
              _data_k[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 tmp_kernel_op_0 = k_dof_0 * 0.5;
          const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
          const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
          const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
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
              _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
        }
    }
    {
      /* FaceType.BLUE */
      const walberla::float64 tmp_coords_jac_0_BLUE =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
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
              _data_src[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          const walberla::float64 src_dof_1 =
              _data_src[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          const walberla::float64 src_dof_2 =
              _data_src[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
          const walberla::float64 k_dof_0 =
              _data_k[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      ((micro_edges_per_macro_edge *
                        (micro_edges_per_macro_edge + 1)) /
                       (2))];
          const walberla::float64 tmp_kernel_op_0 = k_dof_0 * 0.5;
          const walberla::float64 elMatVec_0 = src_dof_0 * tmp_kernel_op_0;
          const walberla::float64 elMatVec_1 = src_dof_1 * tmp_kernel_op_0;
          const walberla::float64 elMatVec_2 = src_dof_2 * tmp_kernel_op_0;
          _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                    ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
              elMatVec_0 +
              _data_dst[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                        ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatVec_1 +
              _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1] =
              elMatVec_2 +
              _data_dst[ctr_0 + (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                        (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
        }
    }
  }
}
void P1ElementwiseDivKGrad_AnnulusMap_float64::
    toMatrix_P1ElementwiseDivKGrad_AnnulusMap_float64_macro_2D(
        idx_t *RESTRICT _data_dst, walberla::float64 *RESTRICT _data_k,
        idx_t *RESTRICT _data_src,
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
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
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
          const walberla::float64 k_dof_0 =
              _data_k[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 tmp_kernel_op_0 = k_dof_0 * 0.5;
          const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
          const int64_t elMat_0_1 = 0;
          const int64_t elMat_0_2 = 0;
          const int64_t elMat_1_0 = 0;
          const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
          const int64_t elMat_1_2 = 0;
          const int64_t elMat_2_0 = 0;
          const int64_t elMat_2_1 = 0;
          const walberla::float64 elMat_2_2 = tmp_kernel_op_0;

          std::vector<uint_t> _data_rowIdx(3);
          std::vector<uint_t> _data_colIdx(3);
          std::vector<real_t> _data_mat(9);

          _data_rowIdx[0] =
              ((uint64_t)(_data_dst[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 2) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2))]));
          _data_rowIdx[1] =
              ((uint64_t)(_data_dst[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 2) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_rowIdx[2] =
              ((uint64_t)(_data_dst[ctr_0 +
                                    (ctr_1 + 1) *
                                        (micro_edges_per_macro_edge + 2) -
                                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_colIdx[0] =
              ((uint64_t)(_data_src[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 2) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2))]));
          _data_colIdx[1] =
              ((uint64_t)(_data_src[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 2) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_colIdx[2] =
              ((uint64_t)(_data_src[ctr_0 +
                                    (ctr_1 + 1) *
                                        (micro_edges_per_macro_edge + 2) -
                                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));

          /* Apply basis transformation */

          _data_mat[0] = ((real_t)(elMat_0_0));
          _data_mat[1] = ((real_t)(elMat_0_1));
          _data_mat[2] = ((real_t)(elMat_0_2));
          _data_mat[3] = ((real_t)(elMat_1_0));
          _data_mat[4] = ((real_t)(elMat_1_1));
          _data_mat[5] = ((real_t)(elMat_1_2));
          _data_mat[6] = ((real_t)(elMat_2_0));
          _data_mat[7] = ((real_t)(elMat_2_1));
          _data_mat[8] = ((real_t)(elMat_2_2));

          mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
        }
    }
    {
      /* FaceType.BLUE */
      const walberla::float64 tmp_coords_jac_0_BLUE =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
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
          const walberla::float64 k_dof_0 =
              _data_k[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      ((micro_edges_per_macro_edge *
                        (micro_edges_per_macro_edge + 1)) /
                       (2))];
          const walberla::float64 tmp_kernel_op_0 = k_dof_0 * 0.5;
          const walberla::float64 elMat_0_0 = tmp_kernel_op_0;
          const int64_t elMat_0_1 = 0;
          const int64_t elMat_0_2 = 0;
          const int64_t elMat_1_0 = 0;
          const walberla::float64 elMat_1_1 = tmp_kernel_op_0;
          const int64_t elMat_1_2 = 0;
          const int64_t elMat_2_0 = 0;
          const int64_t elMat_2_1 = 0;
          const walberla::float64 elMat_2_2 = tmp_kernel_op_0;

          std::vector<uint_t> _data_rowIdx(3);
          std::vector<uint_t> _data_colIdx(3);
          std::vector<real_t> _data_mat(9);

          _data_rowIdx[0] =
              ((uint64_t)(_data_dst[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 2) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_rowIdx[1] =
              ((uint64_t)(_data_dst[ctr_0 +
                                    (ctr_1 + 1) *
                                        (micro_edges_per_macro_edge + 2) -
                                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_rowIdx[2] =
              ((uint64_t)(_data_dst[ctr_0 +
                                    (ctr_1 + 1) *
                                        (micro_edges_per_macro_edge + 2) -
                                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1]));
          _data_colIdx[0] =
              ((uint64_t)(_data_src[ctr_0 +
                                    ctr_1 * (micro_edges_per_macro_edge + 2) -
                                    ((ctr_1 * (ctr_1 + 1)) / (2)) + 1]));
          _data_colIdx[1] =
              ((uint64_t)(_data_src[ctr_0 +
                                    (ctr_1 + 1) *
                                        (micro_edges_per_macro_edge + 2) -
                                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2))]));
          _data_colIdx[2] =
              ((uint64_t)(_data_src[ctr_0 +
                                    (ctr_1 + 1) *
                                        (micro_edges_per_macro_edge + 2) -
                                    (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1]));

          /* Apply basis transformation */

          _data_mat[0] = ((real_t)(elMat_0_0));
          _data_mat[1] = ((real_t)(elMat_0_1));
          _data_mat[2] = ((real_t)(elMat_0_2));
          _data_mat[3] = ((real_t)(elMat_1_0));
          _data_mat[4] = ((real_t)(elMat_1_1));
          _data_mat[5] = ((real_t)(elMat_1_2));
          _data_mat[6] = ((real_t)(elMat_2_0));
          _data_mat[7] = ((real_t)(elMat_2_1));
          _data_mat[8] = ((real_t)(elMat_2_2));

          mat->addValues(_data_rowIdx, _data_colIdx, _data_mat);
        }
    }
  }
}
void P1ElementwiseDivKGrad_AnnulusMap_float64::
    computeInverseDiagonalOperatorValues_P1ElementwiseDivKGrad_AnnulusMap_float64_macro_2D(
        walberla::float64 *RESTRICT _data_invDiag_,
        walberla::float64 *RESTRICT _data_k,
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
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
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
          const walberla::float64 k_dof_0 =
              _data_k[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 1) -
                      ((ctr_1 * (ctr_1 + 1)) / (2))];
          const walberla::float64 tmp_kernel_op_0 = k_dof_0 * 0.5;
          const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
          const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
          const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
          _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2))] =
              elMatDiag_0 +
              _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                             ((ctr_1 * (ctr_1 + 1)) / (2))];
          _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
              elMatDiag_1 +
              _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          _data_invDiag_[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatDiag_2 +
              _data_invDiag_[ctr_0 +
                             (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                             (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
        }
    }
    {
      /* FaceType.BLUE */
      const walberla::float64 tmp_coords_jac_0_BLUE =
          1.0 / (micro_edges_per_macro_edge_float) * 1.0;
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
          const walberla::float64 k_dof_0 =
              _data_k[ctr_0 + ctr_1 * micro_edges_per_macro_edge -
                      ((ctr_1 * (ctr_1 + 1)) / (2)) +
                      ((micro_edges_per_macro_edge *
                        (micro_edges_per_macro_edge + 1)) /
                       (2))];
          const walberla::float64 tmp_kernel_op_0 = k_dof_0 * 0.5;
          const walberla::float64 elMatDiag_0 = tmp_kernel_op_0;
          const walberla::float64 elMatDiag_1 = tmp_kernel_op_0;
          const walberla::float64 elMatDiag_2 = tmp_kernel_op_0;
          _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                         ((ctr_1 * (ctr_1 + 1)) / (2)) + 1] =
              elMatDiag_0 +
              _data_invDiag_[ctr_0 + ctr_1 * (micro_edges_per_macro_edge + 2) -
                             ((ctr_1 * (ctr_1 + 1)) / (2)) + 1];
          _data_invDiag_[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2))] =
              elMatDiag_1 +
              _data_invDiag_[ctr_0 +
                             (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                             (((ctr_1 + 1) * (ctr_1 + 2)) / (2))];
          _data_invDiag_[ctr_0 +
                         (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                         (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1] =
              elMatDiag_2 +
              _data_invDiag_[ctr_0 +
                             (ctr_1 + 1) * (micro_edges_per_macro_edge + 2) -
                             (((ctr_1 + 1) * (ctr_1 + 2)) / (2)) + 1];
        }
    }
  }
}

} // namespace operatorgeneration

} // namespace hyteg
