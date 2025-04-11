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
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"
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
///     u: trial function (vectorial space: TensorialVectorSpace(Lagrange, degree: 2))
///     v: test function  (vectorial space: TensorialVectorSpace(Lagrange, degree: 2))
///     μ: coefficient    (scalar space:    Lagrange, degree: 0)
///
///     ∫ μ { ( 2 ε(u) : ε(v) ) - (2/3) [ ( ∇ · u ) · ( ∇ · v ) ] }
///
/// where
///
///     ε(w) := (1/2) (∇w + (∇w)ᵀ)
///
///
/// And the assembled FE matrix (K) is wrapped with a Rotation matrix (R) locally as below,
///
///     RKRᵀ uᵣ
///
/// where
///     R : Rotation matrix calculated with the normal vector (n̂) at the DoF
///     uᵣ: FE function but the components rotated at the boundaries according to the normal FE function passed
///
///     n̂ : normals (vectorial space: Lagrange, degree: 2)
///         * The passed normal vector must be normalized
///         * The radial component of the rotated vector will be pointing in the given normal direction
///         * If the normals are zero at a DoF, the rotation matrix is identity matrix
///
class P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator
: public Operator< P2VectorFunction< walberla::float64 >, P2VectorFunction< walberla::float64 > >,
  public OperatorWithInverseDiagonal< P2VectorFunction< walberla::float64 > >
{
 public:
   P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator(
       const std::shared_ptr< PrimitiveStorage >& storage,
       size_t                                     minLevel,
       size_t                                     maxLevel,
       const P0Function< walberla::float64 >&     _mu,
       const P2Function< walberla::float64 >&     _nx_rotation,
       const P2Function< walberla::float64 >&     _ny_rotation,
       const P2Function< walberla::float64 >&     _nz_rotation,
       walberla::float64 c_rot_penalty_P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator );

   void apply( const P2VectorFunction< walberla::float64 >& src,
               const P2VectorFunction< walberla::float64 >& dst,
               uint_t                                       level,
               DoFType                                      flag,
               UpdateType                                   updateType = Replace ) const;

   void computeInverseDiagonalOperatorValues();

   std::shared_ptr< P2VectorFunction< walberla::float64 > > getInverseDiagonalValues() const;

 protected:
 private:
   /// Integral: P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Jaśkowiec-Sukumar 04 | points: 11, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///  62439   94032     399     440    396              0                 0              1
   void apply_P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator_macro_3D(
       walberla::float64* RESTRICT _data_dst_edge_0,
       walberla::float64* RESTRICT _data_dst_edge_1,
       walberla::float64* RESTRICT _data_dst_edge_2,
       walberla::float64* RESTRICT _data_dst_vertex_0,
       walberla::float64* RESTRICT _data_dst_vertex_1,
       walberla::float64* RESTRICT _data_dst_vertex_2,
       walberla::float64* RESTRICT _data_mu,
       walberla::float64* RESTRICT _data_nx_rotationEdge,
       walberla::float64* RESTRICT _data_nx_rotationVertex,
       walberla::float64* RESTRICT _data_ny_rotationEdge,
       walberla::float64* RESTRICT _data_ny_rotationVertex,
       walberla::float64* RESTRICT _data_nz_rotationEdge,
       walberla::float64* RESTRICT _data_nz_rotationVertex,
       walberla::float64* RESTRICT _data_src_edge_0,
       walberla::float64* RESTRICT _data_src_edge_1,
       walberla::float64* RESTRICT _data_src_edge_2,
       walberla::float64* RESTRICT _data_src_vertex_0,
       walberla::float64* RESTRICT _data_src_vertex_1,
       walberla::float64* RESTRICT _data_src_vertex_2,
       walberla::float64           c_rot_penalty,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_0comp2,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_1comp2,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       walberla::float64           macro_vertex_coord_id_2comp2,
       walberla::float64           macro_vertex_coord_id_3comp0,
       walberla::float64           macro_vertex_coord_id_3comp1,
       walberla::float64           macro_vertex_coord_id_3comp2,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     apply
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Dunavant 3 | points: 4, degree: 3
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   2756    4856      12      24      0              0                 0              1
   void apply_P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator_macro_2D(
       walberla::float64* RESTRICT _data_dst_edge_0,
       walberla::float64* RESTRICT _data_dst_edge_1,
       walberla::float64* RESTRICT _data_dst_vertex_0,
       walberla::float64* RESTRICT _data_dst_vertex_1,
       walberla::float64* RESTRICT _data_mu,
       walberla::float64* RESTRICT _data_nx_rotationEdge,
       walberla::float64* RESTRICT _data_nx_rotationVertex,
       walberla::float64* RESTRICT _data_ny_rotationEdge,
       walberla::float64* RESTRICT _data_ny_rotationVertex,
       walberla::float64* RESTRICT _data_src_edge_0,
       walberla::float64* RESTRICT _data_src_edge_1,
       walberla::float64* RESTRICT _data_src_vertex_0,
       walberla::float64* RESTRICT _data_src_vertex_1,
       walberla::float64           c_rot_penalty,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Jaśkowiec-Sukumar 04 | points: 11, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   9044   16352     399     440    396              0                 0              1
   void computeInverseDiagonalOperatorValues_P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator_macro_3D(
       walberla::float64* RESTRICT _data_invDiag__edge_0,
       walberla::float64* RESTRICT _data_invDiag__edge_1,
       walberla::float64* RESTRICT _data_invDiag__edge_2,
       walberla::float64* RESTRICT _data_invDiag__vertex_0,
       walberla::float64* RESTRICT _data_invDiag__vertex_1,
       walberla::float64* RESTRICT _data_invDiag__vertex_2,
       walberla::float64* RESTRICT _data_mu,
       walberla::float64* RESTRICT _data_nx_rotationEdge,
       walberla::float64* RESTRICT _data_nx_rotationVertex,
       walberla::float64* RESTRICT _data_ny_rotationEdge,
       walberla::float64* RESTRICT _data_ny_rotationVertex,
       walberla::float64* RESTRICT _data_nz_rotationEdge,
       walberla::float64* RESTRICT _data_nz_rotationVertex,
       walberla::float64           c_rot_penalty,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_0comp2,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_1comp2,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       walberla::float64           macro_vertex_coord_id_2comp2,
       walberla::float64           macro_vertex_coord_id_3comp0,
       walberla::float64           macro_vertex_coord_id_3comp1,
       walberla::float64           macro_vertex_coord_id_3comp2,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Dunavant 3 | points: 4, degree: 3
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    688    1328      12      24      0              0                 0              1
   void computeInverseDiagonalOperatorValues_P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator_macro_2D(
       walberla::float64* RESTRICT _data_invDiag__edge_0,
       walberla::float64* RESTRICT _data_invDiag__edge_1,
       walberla::float64* RESTRICT _data_invDiag__vertex_0,
       walberla::float64* RESTRICT _data_invDiag__vertex_1,
       walberla::float64* RESTRICT _data_mu,
       walberla::float64* RESTRICT _data_nx_rotationEdge,
       walberla::float64* RESTRICT _data_nx_rotationVertex,
       walberla::float64* RESTRICT _data_ny_rotationEdge,
       walberla::float64* RESTRICT _data_ny_rotationVertex,
       walberla::float64           c_rot_penalty,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float ) const;

   std::shared_ptr< P2VectorFunction< walberla::float64 > > invDiag_;
   P0Function< walberla::float64 >                          mu;
   P2Function< walberla::float64 >                          nx_rotation;
   P2Function< walberla::float64 >                          ny_rotation;
   P2Function< walberla::float64 >                          nz_rotation;
   walberla::float64 c_rot_penalty_P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator_;
};

} // namespace operatorgeneration

} // namespace hyteg
