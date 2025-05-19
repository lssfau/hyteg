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
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// "Epsilon" operator.
///
/// Geometry map:    IcosahedralShellMap
///
/// Weak formulation
///
///     u: trial function (vectorial space: TensorialVectorSpace(Lagrange, degree: 2))
///     v: test function  (vectorial space: TensorialVectorSpace(Lagrange, degree: 2))
///     μ: coefficient    (scalar space:    Lagrange, degree: 2)
///
///     ∫ 2 μ ε(u) : ε(v)
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
class P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator
: public Operator< P2VectorFunction< walberla::float64 >, P2VectorFunction< walberla::float64 > >,
  public OperatorWithInverseDiagonal< P2VectorFunction< walberla::float64 > >
{
 public:
   P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator(
       const std::shared_ptr< PrimitiveStorage >& storage,
       size_t                                     minLevel,
       size_t                                     maxLevel,
       const P2Function< walberla::float64 >&     _mu,
       const P2Function< walberla::float64 >&     _nx_rotation,
       const P2Function< walberla::float64 >&     _ny_rotation,
       const P2Function< walberla::float64 >&     _nz_rotation,
       walberla::float64 c_rot_penalty_P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator );

   void apply( const P2VectorFunction< walberla::float64 >& src,
               const P2VectorFunction< walberla::float64 >& dst,
               uint_t                                       level,
               DoFType                                      flag,
               UpdateType                                   updateType = Replace ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P2VectorFunction< idx_t >&            dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const;

   void computeInverseDiagonalOperatorValues();

   std::shared_ptr< P2VectorFunction< walberla::float64 > > getInverseDiagonalValues() const;

 protected:
 private:
   /// Integral: P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Jaśkowiec-Sukumar 04 | points: 11, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///  49701   80887     454     473    385              0                 0              0
   void apply_P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator_macro_3D(
       walberla::float64* RESTRICT _data_dst_edge_0,
       walberla::float64* RESTRICT _data_dst_edge_1,
       walberla::float64* RESTRICT _data_dst_edge_2,
       walberla::float64* RESTRICT _data_dst_vertex_0,
       walberla::float64* RESTRICT _data_dst_vertex_1,
       walberla::float64* RESTRICT _data_dst_vertex_2,
       walberla::float64* RESTRICT _data_muEdge,
       walberla::float64* RESTRICT _data_muVertex,
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
       walberla::float64           forVertex_0,
       walberla::float64           forVertex_1,
       walberla::float64           forVertex_2,
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
       walberla::float64           micro_edges_per_macro_edge_float,
       walberla::float64           radRayVertex,
       walberla::float64           radRefVertex,
       walberla::float64           rayVertex_0,
       walberla::float64           rayVertex_1,
       walberla::float64           rayVertex_2,
       walberla::float64           refVertex_0,
       walberla::float64           refVertex_1,
       walberla::float64           refVertex_2,
       walberla::float64           thrVertex_0,
       walberla::float64           thrVertex_1,
       walberla::float64           thrVertex_2 ) const;

   /// Integral: P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Jaśkowiec-Sukumar 04 | points: 11, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///  48801   79987     454     473    385              0                 0              3
   void toMatrix_P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator_macro_3D(
       idx_t* RESTRICT                      _data_dst_edge_0,
       idx_t* RESTRICT                      _data_dst_edge_1,
       idx_t* RESTRICT                      _data_dst_edge_2,
       idx_t* RESTRICT                      _data_dst_vertex_0,
       idx_t* RESTRICT                      _data_dst_vertex_1,
       idx_t* RESTRICT                      _data_dst_vertex_2,
       walberla::float64* RESTRICT          _data_muEdge,
       walberla::float64* RESTRICT          _data_muVertex,
       walberla::float64* RESTRICT          _data_nx_rotationEdge,
       walberla::float64* RESTRICT          _data_nx_rotationVertex,
       walberla::float64* RESTRICT          _data_ny_rotationEdge,
       walberla::float64* RESTRICT          _data_ny_rotationVertex,
       walberla::float64* RESTRICT          _data_nz_rotationEdge,
       walberla::float64* RESTRICT          _data_nz_rotationVertex,
       idx_t* RESTRICT                      _data_src_edge_0,
       idx_t* RESTRICT                      _data_src_edge_1,
       idx_t* RESTRICT                      _data_src_edge_2,
       idx_t* RESTRICT                      _data_src_vertex_0,
       idx_t* RESTRICT                      _data_src_vertex_1,
       idx_t* RESTRICT                      _data_src_vertex_2,
       walberla::float64                    c_rot_penalty,
       walberla::float64                    forVertex_0,
       walberla::float64                    forVertex_1,
       walberla::float64                    forVertex_2,
       walberla::float64                    macro_vertex_coord_id_0comp0,
       walberla::float64                    macro_vertex_coord_id_0comp1,
       walberla::float64                    macro_vertex_coord_id_0comp2,
       walberla::float64                    macro_vertex_coord_id_1comp0,
       walberla::float64                    macro_vertex_coord_id_1comp1,
       walberla::float64                    macro_vertex_coord_id_1comp2,
       walberla::float64                    macro_vertex_coord_id_2comp0,
       walberla::float64                    macro_vertex_coord_id_2comp1,
       walberla::float64                    macro_vertex_coord_id_2comp2,
       walberla::float64                    macro_vertex_coord_id_3comp0,
       walberla::float64                    macro_vertex_coord_id_3comp1,
       walberla::float64                    macro_vertex_coord_id_3comp2,
       std::shared_ptr< SparseMatrixProxy > mat,
       int64_t                              micro_edges_per_macro_edge,
       walberla::float64                    micro_edges_per_macro_edge_float,
       walberla::float64                    radRayVertex,
       walberla::float64                    radRefVertex,
       walberla::float64                    rayVertex_0,
       walberla::float64                    rayVertex_1,
       walberla::float64                    rayVertex_2,
       walberla::float64                    refVertex_0,
       walberla::float64                    refVertex_1,
       walberla::float64                    refVertex_2,
       walberla::float64                    thrVertex_0,
       walberla::float64                    thrVertex_1,
       walberla::float64                    thrVertex_2 ) const;

   /// Integral: P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   SAWTOOTH
   /// - quadrature rule: Jaśkowiec-Sukumar 04 | points: 11, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///  10716   18596     454     473    385              0                 0              0
   void computeInverseDiagonalOperatorValues_P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator_macro_3D(
       walberla::float64* RESTRICT _data_invDiag__edge_0,
       walberla::float64* RESTRICT _data_invDiag__edge_1,
       walberla::float64* RESTRICT _data_invDiag__edge_2,
       walberla::float64* RESTRICT _data_invDiag__vertex_0,
       walberla::float64* RESTRICT _data_invDiag__vertex_1,
       walberla::float64* RESTRICT _data_invDiag__vertex_2,
       walberla::float64* RESTRICT _data_muEdge,
       walberla::float64* RESTRICT _data_muVertex,
       walberla::float64* RESTRICT _data_nx_rotationEdge,
       walberla::float64* RESTRICT _data_nx_rotationVertex,
       walberla::float64* RESTRICT _data_ny_rotationEdge,
       walberla::float64* RESTRICT _data_ny_rotationVertex,
       walberla::float64* RESTRICT _data_nz_rotationEdge,
       walberla::float64* RESTRICT _data_nz_rotationVertex,
       walberla::float64           c_rot_penalty,
       walberla::float64           forVertex_0,
       walberla::float64           forVertex_1,
       walberla::float64           forVertex_2,
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
       walberla::float64           micro_edges_per_macro_edge_float,
       walberla::float64           radRayVertex,
       walberla::float64           radRefVertex,
       walberla::float64           rayVertex_0,
       walberla::float64           rayVertex_1,
       walberla::float64           rayVertex_2,
       walberla::float64           refVertex_0,
       walberla::float64           refVertex_1,
       walberla::float64           refVertex_2,
       walberla::float64           thrVertex_0,
       walberla::float64           thrVertex_1,
       walberla::float64           thrVertex_2 ) const;

   std::shared_ptr< P2VectorFunction< walberla::float64 > > invDiag_;
   P2Function< walberla::float64 >                          mu;
   P2Function< walberla::float64 >                          nx_rotation;
   P2Function< walberla::float64 >                          ny_rotation;
   P2Function< walberla::float64 >                          nz_rotation;
   walberla::float64 c_rot_penalty_P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator_;
};

} // namespace operatorgeneration

} // namespace hyteg
