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

#include "hyteg/experimental/P2PlusBubbleFunction.hpp"

#include "hyteg/operators/Operator.hpp"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "hyteg/solvers/Smoothables.hpp"

#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Diffusion operator without coefficients.
///
/// Geometry map: IdentityMap
///
/// Weak formulation
///
///     u: trial function (space: P2PlusBubble)
///     v: test function  (space: P2PlusBubble)
///
///     ∫ ∇u : ∇v
///
///     Note that the double contraction (:) reduces to the dot product for
///     scalar function spaces, i.e. the form becomes
///
///     ∫ ∇u · ∇v

class P2PlusBubbleElementwiseDiffusion_float64
    : public Operator<P2PlusBubbleFunction<walberla::float64>,
                      P2PlusBubbleFunction<walberla::float64>>,
      public OperatorWithInverseDiagonal<
          P2PlusBubbleFunction<walberla::float64>> {
public:
  P2PlusBubbleElementwiseDiffusion_float64(
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
  /// Integral: P2PlusBubbleElementwiseDiffusion_float64
  /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
  /// - kernel type:     apply
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
  /// - blending map:    IdentityMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///    485     723      12       0      0              0                 0 0
  void apply_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(
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
      walberla::float64 micro_edges_per_macro_edge_float) const;

  /// Integral: P2PlusBubbleElementwiseDiffusion_float64
  /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
  /// - kernel type:     toMatrix
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
  /// - blending map:    IdentityMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///    436     674      12       0      0              0                 0 3
  void toMatrix_P2PlusBubbleElementwiseDiffusion_float64_macro_2D(
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
      walberla::float64 micro_edges_per_macro_edge_float) const;

  /// Integral: P2PlusBubbleElementwiseDiffusion_float64
  /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
  /// - kernel type:     computeInverseDiagonalOperatorValues
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
  /// - blending map:    IdentityMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///    212     272      12       0      0              0                 0 0
  void
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
      walberla::float64 micro_edges_per_macro_edge_float) const;

  std::shared_ptr<P2PlusBubbleFunction<walberla::float64>> invDiag_;
};

} // namespace operatorgeneration

} // namespace hyteg
