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

#pragma once

#include "core/DataTypes.h"

#include "hyteg/LikwidWrapper.hpp"

#include "hyteg/boundary/BoundaryConditions.hpp"

#include "hyteg/communication/Syncing.hpp"

#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"

#include "hyteg/geometry/AnnulusMap.hpp"

#include "hyteg/operators/Operator.hpp"

#include "hyteg/p2functionspace/P2PlusBubbleFunction.hpp"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "hyteg/solvers/Smoothables.hpp"

#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Mass operator.
///
/// Geometry map: AnnulusMap
///
/// Weak formulation
///
///     u: trial function (space: P2PlusBubble)
///     v: test function  (space: P2PlusBubble)
///
///     ∫ uv

class P2PlusBubbleElementwiseMass_AnnulusMap_float64
    : public Operator<P2PlusBubbleFunction<walberla::float64>,
                      P2PlusBubbleFunction<walberla::float64>>,
      public OperatorWithInverseDiagonal<
          P2PlusBubbleFunction<walberla::float64>> {
public:
  P2PlusBubbleElementwiseMass_AnnulusMap_float64(
      const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
      size_t maxLevel);

  void apply(const P2PlusBubbleFunction<walberla::float64> &src,
             const P2PlusBubbleFunction<walberla::float64> &dst, uint_t level,
             DoFType flag, UpdateType updateType = Replace) const;

  void toMatrix(const std::shared_ptr<SparseMatrixProxy> &mat,
                const P2PlusBubbleFunction<idx_t> &src,
                const P2PlusBubbleFunction<idx_t> &dst, uint_t level,
                DoFType flag) const;

  void computeInverseDiagonalOperatorValues();

  std::shared_ptr<P2PlusBubbleFunction<walberla::float64>>
  getInverseDiagonalValues() const;

protected:
private:
  /// Integral: P2PlusBubbleElementwiseMass_AnnulusMap_float64
  /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
  /// - kernel type:     apply
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Griener-Schmid 1 | points: 10, degree: 6
  /// - blending map:    AnnulusMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///    523     906      13      20      0              0                 0 0
  void apply_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(
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
      walberla::float64 thrVertex_0, walberla::float64 thrVertex_1) const;

  /// Integral: P2PlusBubbleElementwiseMass_AnnulusMap_float64
  /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
  /// - kernel type:     toMatrix
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Griener-Schmid 1 | points: 10, degree: 6
  /// - blending map:    AnnulusMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///    474     857      13      20      0              0                 0 3
  void toMatrix_P2PlusBubbleElementwiseMass_AnnulusMap_float64_macro_2D(
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
      walberla::float64 thrVertex_0, walberla::float64 thrVertex_1) const;

  /// Integral: P2PlusBubbleElementwiseMass_AnnulusMap_float64
  /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
  /// - kernel type:     computeInverseDiagonalOperatorValues
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Griener-Schmid 1 | points: 10, degree: 6
  /// - blending map:    AnnulusMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///    272     457      13      20      0              0                 0 0
  void
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
      walberla::float64 thrVertex_0, walberla::float64 thrVertex_1) const;

  std::shared_ptr<P2PlusBubbleFunction<walberla::float64>> invDiag_;
};

} // namespace operatorgeneration

} // namespace hyteg
