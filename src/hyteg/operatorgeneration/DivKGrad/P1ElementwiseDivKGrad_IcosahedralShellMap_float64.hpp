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

#include "hyteg/geometry/IcosahedralShellMap.hpp"

#include "hyteg/operators/Operator.hpp"

#include "hyteg/p0functionspace/P0Function.hpp"

#include "hyteg/p1functionspace/P1Function.hpp"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "hyteg/solvers/Smoothables.hpp"

#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Diffusion operator with a scalar coefficient.
///
/// Geometry map: IcosahedralShellMap
///
/// Weak formulation
///
///     u: trial function (space: Lagrange, degree: 1)
///     v: test function  (space: Lagrange, degree: 1)
///     k: coefficient    (space: Lagrange, degree: 0)
///
///     ∫ k ∇u · ∇v

class P1ElementwiseDivKGrad_IcosahedralShellMap_float64
    : public Operator<P1Function<walberla::float64>,
                      P1Function<walberla::float64>>,
      public OperatorWithInverseDiagonal<P1Function<walberla::float64>> {
public:
  P1ElementwiseDivKGrad_IcosahedralShellMap_float64(
      const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
      size_t maxLevel, const P0Function<walberla::float64> &_k);

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
  /// Integral: P1ElementwiseDivKGrad_IcosahedralShellMap_float64
  /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
  /// - kernel type:     apply
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Hammer-Marlowe-Stroud 3 | points: 5, degree: 3
  /// - blending map:    IcosahedralShellMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///     97      77      36       0      0              0                 0 0
  void apply_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(
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
      walberla::float64 micro_edges_per_macro_edge_float) const;

  /// Integral: P1ElementwiseDivKGrad_IcosahedralShellMap_float64
  /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
  /// - kernel type:     toMatrix
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Hammer-Marlowe-Stroud 3 | points: 5, degree: 3
  /// - blending map:    IcosahedralShellMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///     93      73      36       0      0              0                 0 3
  void toMatrix_P1ElementwiseDivKGrad_IcosahedralShellMap_float64_macro_3D(
      idx_t *RESTRICT _data_dst, walberla::float64 *RESTRICT _data_k,
      idx_t *RESTRICT _data_src, walberla::float64 macro_vertex_coord_id_0comp0,
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

  /// Integral: P1ElementwiseDivKGrad_IcosahedralShellMap_float64
  /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
  /// - kernel type:     computeInverseDiagonalOperatorValues
  /// - loop strategy:   SAWTOOTH
  /// - quadrature rule: Hammer-Marlowe-Stroud 3 | points: 5, degree: 3
  /// - blending map:    IcosahedralShellMap
  /// - operations per element:
  ///   adds    muls    divs    pows    abs    assignments    function_calls
  ///   unknown_ops
  /// ------  ------  ------  ------  -----  -------------  ----------------
  /// -------------
  ///     97      73      36       0      0              0                 0 0
  void
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
      walberla::float64 micro_edges_per_macro_edge_float) const;

  std::shared_ptr<P1Function<walberla::float64>> invDiag_;
  P0Function<walberla::float64> k;
};

} // namespace operatorgeneration

} // namespace hyteg
