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

#pragma once

#include "core/DataTypes.h"

#include "hyteg/LikwidWrapper.hpp"

#include "hyteg/communication/Syncing.hpp"

#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"

#include "hyteg/operators/Operator.hpp"

#include "hyteg/p1functionspace/P1Function.hpp"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "hyteg/solvers/Smoothables.hpp"

#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Diffusion operator without coefficients.
///
/// Geometry map: IdentityMap
///
/// Weak formulation
///
///     u: trial function (space: Lagrange, degree: 1)
///     v: test function  (space: Lagrange, degree: 1)
///
///     ∫ ∇u · ∇v

class P1ElementwiseDiffusion_cubes_const_float64
    : public Operator<P1Function<walberla::float64>,
                      P1Function<walberla::float64>>,
      public OperatorWithInverseDiagonal<P1Function<walberla::float64>> {
public:
  P1ElementwiseDiffusion_cubes_const_float64(
      const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
      size_t maxLevel);

  void apply(const P1Function<walberla::float64> &src,
             const P1Function<walberla::float64> &dst, uint_t level,
             DoFType flag, UpdateType updateType = Replace) const;

  void toMatrix(const std::shared_ptr<SparseMatrixProxy> &mat,
                const P1Function<idx_t> &src, const P1Function<idx_t> &dst,
                uint_t level, DoFType flag) const;

  void computeInverseDiagonalOperatorValues();

  std::shared_ptr<P1Function<walberla::float64>>
  getInverseDiagonalValues() const;

protected:
private:
  /// Kernel type: apply
  /// - quadrature rule: Centroid rule | points: 1, degree: 1
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///     17      25       0       0      0              0                 0 0
  void apply_macro_2D(walberla::float64 *RESTRICT _data_dst,
                      walberla::float64 *RESTRICT _data_src,
                      walberla::float64 macro_vertex_coord_id_0comp0,
                      walberla::float64 macro_vertex_coord_id_0comp1,
                      walberla::float64 macro_vertex_coord_id_1comp0,
                      walberla::float64 macro_vertex_coord_id_1comp1,
                      walberla::float64 macro_vertex_coord_id_2comp0,
                      walberla::float64 macro_vertex_coord_id_2comp1,
                      int64_t micro_edges_per_macro_edge,
                      walberla::float64 micro_edges_per_macro_edge_float) const;
  /// Kernel type: apply
  /// - quadrature rule: Keast 0 | points: 1, degree: 1
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///     42      51       0       0      0              0                 0 0
  void apply_macro_3D(walberla::float64 *RESTRICT _data_dst,
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
                      walberla::float64 micro_edges_per_macro_edge_float) const;
  /// Kernel type: toMatrix
  /// - quadrature rule: Centroid rule | points: 1, degree: 1
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///      8      19       0       0      0              0                 0 3
  void
  toMatrix_macro_2D(idx_t *RESTRICT _data_dst, idx_t *RESTRICT _data_src,
                    walberla::float64 macro_vertex_coord_id_0comp0,
                    walberla::float64 macro_vertex_coord_id_0comp1,
                    walberla::float64 macro_vertex_coord_id_1comp0,
                    walberla::float64 macro_vertex_coord_id_1comp1,
                    walberla::float64 macro_vertex_coord_id_2comp0,
                    walberla::float64 macro_vertex_coord_id_2comp1,
                    std::shared_ptr<SparseMatrixProxy> mat,
                    int64_t micro_edges_per_macro_edge,
                    walberla::float64 micro_edges_per_macro_edge_float) const;
  /// Kernel type: toMatrix
  /// - quadrature rule: Keast 0 | points: 1, degree: 1
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///     26      41       0       0      0              0                 0 3
  void
  toMatrix_macro_3D(idx_t *RESTRICT _data_dst, idx_t *RESTRICT _data_src,
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
                    walberla::float64 micro_edges_per_macro_edge_float) const;
  /// Kernel type: computeInverseDiagonalOperatorValues
  /// - quadrature rule: Centroid rule | points: 1, degree: 1
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///      8      10       0       0      0              0                 0 0
  void computeInverseDiagonalOperatorValues_macro_2D(
      walberla::float64 *RESTRICT _data_invDiag_,
      walberla::float64 macro_vertex_coord_id_0comp0,
      walberla::float64 macro_vertex_coord_id_0comp1,
      walberla::float64 macro_vertex_coord_id_1comp0,
      walberla::float64 macro_vertex_coord_id_1comp1,
      walberla::float64 macro_vertex_coord_id_2comp0,
      walberla::float64 macro_vertex_coord_id_2comp1,
      int64_t micro_edges_per_macro_edge,
      walberla::float64 micro_edges_per_macro_edge_float) const;
  /// Kernel type: computeInverseDiagonalOperatorValues
  /// - quadrature rule: Keast 0 | points: 1, degree: 1
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///     18      17       0       0      0              0                 0 0
  void computeInverseDiagonalOperatorValues_macro_3D(
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
      walberla::float64 micro_edges_per_macro_edge_float) const;

  std::shared_ptr<P1Function<walberla::float64>> invDiag_;
};

} // namespace operatorgeneration

} // namespace hyteg
