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
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Implements the fully coupled viscous operator of the Stokes problem.
/// The latter is the extension of the Epsilon operator to the case where
/// the velocity field need not be divergence-free. This is e.g. the case
/// in the (truncated) anelastic liquid approximation of mantle convection.
///
/// The strong representation of the operator is given by:
///
///    - div[ μ (grad(u)+grad(u)ᵀ) ] + 2/3 grad[ μ div(u) ]
///
/// Note that the factor 2/3 means that for 2D this is the pseudo-3D form
/// of the operator.
///
/// Component trial: 0
/// Component test:  0
/// Geometry map:    IdentityMap
///
/// Weak formulation
///
///     u: trial function (vectorial space: Lagrange, degree: 1)
///     v: test function  (vectorial space: Lagrange, degree: 1)
///     μ: coefficient    (scalar space:    Lagrange, degree: 1)
///
///     ∫ μ { ( 2 ε(u) : ε(v) ) - (2/3) [ ( ∇ · u ) · ( ∇ · v ) ] }
///
/// where
///
///     ε(w) := (1/2) (∇w + (∇w)ᵀ)

class P2EvaluateViscosityViscoplastic : public Operator< P1Function< walberla::float64 >, P1Function< walberla::float64 > >,
                                        public OperatorWithInverseDiagonal< P1Function< walberla::float64 > >
{
 public:
   P2EvaluateViscosityViscoplastic( const std::shared_ptr< PrimitiveStorage >& storage,
                                    size_t                                     minLevel,
                                    size_t                                     maxLevel,
                                    const P1Function< walberla::float64 >&     _mu_lin,
                                    const P1Function< walberla::float64 >&     _ux,
                                    const P1Function< walberla::float64 >&     _uy,
                                    walberla::float64                          mu_star_P2EvaluateViscosityViscoplastic,
                                    walberla::float64                          sigma_y_P2EvaluateViscosityViscoplastic );

   void apply( const P1Function< walberla::float64 >& src,
               const P1Function< walberla::float64 >& dst,
               uint_t                                 level,
               DoFType                                flag,
               UpdateType                             updateType = Replace ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< idx_t >&                  src,
                  const P1Function< idx_t >&                  dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const;

   void computeInverseDiagonalOperatorValues();

   std::shared_ptr< P1Function< walberla::float64 > > getInverseDiagonalValues() const;

 protected:
 private:
   /// Integral: P2EvaluateViscosityViscoplastic
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     apply
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Dunavant 2 | points: 3, degree: 2
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     53      46      19       1      0              0                 0              0
   void apply_P2EvaluateViscosityViscoplastic_macro_2D( walberla::float64* RESTRICT _data_dst,
                                                        walberla::float64* RESTRICT _data_mu_lin,
                                                        walberla::float64* RESTRICT _data_src,
                                                        walberla::float64* RESTRICT _data_ux,
                                                        walberla::float64* RESTRICT _data_uy,
                                                        walberla::float64           macro_vertex_coord_id_0comp0,
                                                        walberla::float64           macro_vertex_coord_id_0comp1,
                                                        walberla::float64           macro_vertex_coord_id_1comp0,
                                                        walberla::float64           macro_vertex_coord_id_1comp1,
                                                        walberla::float64           macro_vertex_coord_id_2comp0,
                                                        walberla::float64           macro_vertex_coord_id_2comp1,
                                                        int64_t                     micro_edges_per_macro_edge,
                                                        walberla::float64           micro_edges_per_macro_edge_float,
                                                        walberla::float64           mu_star,
                                                        walberla::float64           sigma_y ) const;

   /// Integral: P2EvaluateViscosityViscoplastic
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     toMatrix
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Dunavant 2 | points: 3, degree: 2
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     50      43      19       1      0              0                 0              3
   void toMatrix_P2EvaluateViscosityViscoplastic_macro_2D( idx_t* RESTRICT                      _data_dst,
                                                           walberla::float64* RESTRICT          _data_mu_lin,
                                                           idx_t* RESTRICT                      _data_src,
                                                           walberla::float64* RESTRICT          _data_ux,
                                                           walberla::float64* RESTRICT          _data_uy,
                                                           walberla::float64                    macro_vertex_coord_id_0comp0,
                                                           walberla::float64                    macro_vertex_coord_id_0comp1,
                                                           walberla::float64                    macro_vertex_coord_id_1comp0,
                                                           walberla::float64                    macro_vertex_coord_id_1comp1,
                                                           walberla::float64                    macro_vertex_coord_id_2comp0,
                                                           walberla::float64                    macro_vertex_coord_id_2comp1,
                                                           std::shared_ptr< SparseMatrixProxy > mat,
                                                           int64_t                              micro_edges_per_macro_edge,
                                                           walberla::float64                    micro_edges_per_macro_edge_float,
                                                           walberla::float64                    mu_star,
                                                           walberla::float64                    sigma_y ) const;

   /// Integral: P2EvaluateViscosityViscoplastic
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Dunavant 2 | points: 3, degree: 2
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     53      43      19       1      0              0                 0              0
   void computeInverseDiagonalOperatorValues_P2EvaluateViscosityViscoplastic_macro_2D(
       walberla::float64* RESTRICT _data_invDiag_,
       walberla::float64* RESTRICT _data_mu_lin,
       walberla::float64* RESTRICT _data_ux,
       walberla::float64* RESTRICT _data_uy,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float,
       walberla::float64           mu_star,
       walberla::float64           sigma_y ) const;

   std::shared_ptr< P1Function< walberla::float64 > > invDiag_;
   P1Function< walberla::float64 >                    mu_lin;
   P1Function< walberla::float64 >                    ux;
   P1Function< walberla::float64 >                    uy;
   walberla::float64                                  mu_star_P2EvaluateViscosityViscoplastic_;
   walberla::float64                                  sigma_y_P2EvaluateViscosityViscoplastic_;
};

} // namespace operatorgeneration

} // namespace hyteg
