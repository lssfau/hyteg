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
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Mass operator.
///
/// Geometry map: IdentityMap
///
/// Weak formulation
///
///     u: trial function (space: Lagrange, degree: 2)
///     v: test function  (space: Lagrange, degree: 2)
///
///     ∫ uv ds

class P2ElementwiseBoundaryMass : public Operator< P2Function< walberla::float64 >, P2Function< walberla::float64 > >,
                                  public OperatorWithInverseDiagonal< P2Function< walberla::float64 > >
{
 public:
   P2ElementwiseBoundaryMass( const std::shared_ptr< PrimitiveStorage >& storage,
                              size_t                                     minLevel,
                              size_t                                     maxLevel,
                              BoundaryCondition                          boundaryCondition,
                              BoundaryUID                                P2ElementwiseBoundaryMass_boundary_uid );

   void apply( const P2Function< walberla::float64 >& src,
               const P2Function< walberla::float64 >& dst,
               uint_t                                 level,
               DoFType                                flag,
               UpdateType                             updateType = Replace ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  src,
                  const P2Function< idx_t >&                  dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const;

   void computeInverseDiagonalOperatorValues();

   std::shared_ptr< P2Function< walberla::float64 > > getInverseDiagonalValues() const;

 protected:
 private:
   /// Integral: P2ElementwiseBoundaryMass_facet_id_0
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    654     686      36       1      1              0                 0              0
   void apply_P2ElementwiseBoundaryMass_facet_id_0_macro_3D( walberla::float64* RESTRICT _data_dstEdge,
                                                             walberla::float64* RESTRICT _data_dstVertex,
                                                             walberla::float64* RESTRICT _data_srcEdge,
                                                             walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_1
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    651     686      36       1      1              0                 0              0
   void apply_P2ElementwiseBoundaryMass_facet_id_1_macro_3D( walberla::float64* RESTRICT _data_dstEdge,
                                                             walberla::float64* RESTRICT _data_dstVertex,
                                                             walberla::float64* RESTRICT _data_srcEdge,
                                                             walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_2
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    654     686      36       1      1              0                 0              0
   void apply_P2ElementwiseBoundaryMass_facet_id_2_macro_3D( walberla::float64* RESTRICT _data_dstEdge,
                                                             walberla::float64* RESTRICT _data_dstVertex,
                                                             walberla::float64* RESTRICT _data_srcEdge,
                                                             walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_3
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    657     686      36       1      1              0                 0              0
   void apply_P2ElementwiseBoundaryMass_facet_id_3_macro_3D( walberla::float64* RESTRICT _data_dstEdge,
                                                             walberla::float64* RESTRICT _data_dstVertex,
                                                             walberla::float64* RESTRICT _data_srcEdge,
                                                             walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_0
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    129     171      12       1      1              0                 0              0
   void apply_P2ElementwiseBoundaryMass_facet_id_0_macro_2D( walberla::float64* RESTRICT _data_dstEdge,
                                                             walberla::float64* RESTRICT _data_dstVertex,
                                                             walberla::float64* RESTRICT _data_srcEdge,
                                                             walberla::float64* RESTRICT _data_srcVertex,
                                                             walberla::float64           macro_vertex_coord_id_0comp0,
                                                             walberla::float64           macro_vertex_coord_id_0comp1,
                                                             walberla::float64           macro_vertex_coord_id_1comp0,
                                                             walberla::float64           macro_vertex_coord_id_1comp1,
                                                             walberla::float64           macro_vertex_coord_id_2comp0,
                                                             walberla::float64           macro_vertex_coord_id_2comp1,
                                                             int64_t                     micro_edges_per_macro_edge,
                                                             walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_1
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    129     171      12       1      1              0                 0              0
   void apply_P2ElementwiseBoundaryMass_facet_id_1_macro_2D( walberla::float64* RESTRICT _data_dstEdge,
                                                             walberla::float64* RESTRICT _data_dstVertex,
                                                             walberla::float64* RESTRICT _data_srcEdge,
                                                             walberla::float64* RESTRICT _data_srcVertex,
                                                             walberla::float64           macro_vertex_coord_id_0comp0,
                                                             walberla::float64           macro_vertex_coord_id_0comp1,
                                                             walberla::float64           macro_vertex_coord_id_1comp0,
                                                             walberla::float64           macro_vertex_coord_id_1comp1,
                                                             walberla::float64           macro_vertex_coord_id_2comp0,
                                                             walberla::float64           macro_vertex_coord_id_2comp1,
                                                             int64_t                     micro_edges_per_macro_edge,
                                                             walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_2
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    129     171      12       1      1              0                 0              0
   void apply_P2ElementwiseBoundaryMass_facet_id_2_macro_2D( walberla::float64* RESTRICT _data_dstEdge,
                                                             walberla::float64* RESTRICT _data_dstVertex,
                                                             walberla::float64* RESTRICT _data_srcEdge,
                                                             walberla::float64* RESTRICT _data_srcVertex,
                                                             walberla::float64           macro_vertex_coord_id_0comp0,
                                                             walberla::float64           macro_vertex_coord_id_0comp1,
                                                             walberla::float64           macro_vertex_coord_id_1comp0,
                                                             walberla::float64           macro_vertex_coord_id_1comp1,
                                                             walberla::float64           macro_vertex_coord_id_2comp0,
                                                             walberla::float64           macro_vertex_coord_id_2comp1,
                                                             int64_t                     micro_edges_per_macro_edge,
                                                             walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_0
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    554     586      36       1      1              0                 0              3
   void toMatrix_P2ElementwiseBoundaryMass_facet_id_0_macro_3D( idx_t* RESTRICT                      _data_dstEdge,
                                                                idx_t* RESTRICT                      _data_dstVertex,
                                                                idx_t* RESTRICT                      _data_srcEdge,
                                                                idx_t* RESTRICT                      _data_srcVertex,
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
                                                                walberla::float64 micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_1
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    551     586      36       1      1              0                 0              3
   void toMatrix_P2ElementwiseBoundaryMass_facet_id_1_macro_3D( idx_t* RESTRICT                      _data_dstEdge,
                                                                idx_t* RESTRICT                      _data_dstVertex,
                                                                idx_t* RESTRICT                      _data_srcEdge,
                                                                idx_t* RESTRICT                      _data_srcVertex,
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
                                                                walberla::float64 micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_2
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    554     586      36       1      1              0                 0              3
   void toMatrix_P2ElementwiseBoundaryMass_facet_id_2_macro_3D( idx_t* RESTRICT                      _data_dstEdge,
                                                                idx_t* RESTRICT                      _data_dstVertex,
                                                                idx_t* RESTRICT                      _data_srcEdge,
                                                                idx_t* RESTRICT                      _data_srcVertex,
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
                                                                walberla::float64 micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_3
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    557     586      36       1      1              0                 0              3
   void toMatrix_P2ElementwiseBoundaryMass_facet_id_3_macro_3D( idx_t* RESTRICT                      _data_dstEdge,
                                                                idx_t* RESTRICT                      _data_dstVertex,
                                                                idx_t* RESTRICT                      _data_srcEdge,
                                                                idx_t* RESTRICT                      _data_srcVertex,
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
                                                                walberla::float64 micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_0
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     93     135      12       1      1              0                 0              3
   void toMatrix_P2ElementwiseBoundaryMass_facet_id_0_macro_2D( idx_t* RESTRICT                      _data_dstEdge,
                                                                idx_t* RESTRICT                      _data_dstVertex,
                                                                idx_t* RESTRICT                      _data_srcEdge,
                                                                idx_t* RESTRICT                      _data_srcVertex,
                                                                walberla::float64                    macro_vertex_coord_id_0comp0,
                                                                walberla::float64                    macro_vertex_coord_id_0comp1,
                                                                walberla::float64                    macro_vertex_coord_id_1comp0,
                                                                walberla::float64                    macro_vertex_coord_id_1comp1,
                                                                walberla::float64                    macro_vertex_coord_id_2comp0,
                                                                walberla::float64                    macro_vertex_coord_id_2comp1,
                                                                std::shared_ptr< SparseMatrixProxy > mat,
                                                                int64_t                              micro_edges_per_macro_edge,
                                                                walberla::float64 micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_1
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     93     135      12       1      1              0                 0              3
   void toMatrix_P2ElementwiseBoundaryMass_facet_id_1_macro_2D( idx_t* RESTRICT                      _data_dstEdge,
                                                                idx_t* RESTRICT                      _data_dstVertex,
                                                                idx_t* RESTRICT                      _data_srcEdge,
                                                                idx_t* RESTRICT                      _data_srcVertex,
                                                                walberla::float64                    macro_vertex_coord_id_0comp0,
                                                                walberla::float64                    macro_vertex_coord_id_0comp1,
                                                                walberla::float64                    macro_vertex_coord_id_1comp0,
                                                                walberla::float64                    macro_vertex_coord_id_1comp1,
                                                                walberla::float64                    macro_vertex_coord_id_2comp0,
                                                                walberla::float64                    macro_vertex_coord_id_2comp1,
                                                                std::shared_ptr< SparseMatrixProxy > mat,
                                                                int64_t                              micro_edges_per_macro_edge,
                                                                walberla::float64 micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_2
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     93     135      12       1      1              0                 0              3
   void toMatrix_P2ElementwiseBoundaryMass_facet_id_2_macro_2D( idx_t* RESTRICT                      _data_dstEdge,
                                                                idx_t* RESTRICT                      _data_dstVertex,
                                                                idx_t* RESTRICT                      _data_srcEdge,
                                                                idx_t* RESTRICT                      _data_srcVertex,
                                                                walberla::float64                    macro_vertex_coord_id_0comp0,
                                                                walberla::float64                    macro_vertex_coord_id_0comp1,
                                                                walberla::float64                    macro_vertex_coord_id_1comp0,
                                                                walberla::float64                    macro_vertex_coord_id_1comp1,
                                                                walberla::float64                    macro_vertex_coord_id_2comp0,
                                                                walberla::float64                    macro_vertex_coord_id_2comp1,
                                                                std::shared_ptr< SparseMatrixProxy > mat,
                                                                int64_t                              micro_edges_per_macro_edge,
                                                                walberla::float64 micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_0
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    248     219      36       1      1              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseBoundaryMass_facet_id_0_macro_3D(
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_1
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    245     219      36       1      1              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseBoundaryMass_facet_id_1_macro_3D(
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_2
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    248     219      36       1      1              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseBoundaryMass_facet_id_2_macro_3D(
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_3
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///    251     219      36       1      1              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseBoundaryMass_facet_id_3_macro_3D(
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   /// Integral: P2ElementwiseBoundaryMass_facet_id_0
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     69      86      12       1      1              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseBoundaryMass_facet_id_0_macro_2D(
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_1
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     69      86      12       1      1              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseBoundaryMass_facet_id_1_macro_2D(
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float ) const;

   /// Integral: P2ElementwiseBoundaryMass_facet_id_2
   /// - volume element:  triangle, dim: 2, vertices: 3, spacedim: 2
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Gauss-Legendre | points: 3, degree: 5
   /// - blending map:    IdentityMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///     69      86      12       1      1              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseBoundaryMass_facet_id_2_macro_2D(
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
       walberla::float64           macro_vertex_coord_id_0comp0,
       walberla::float64           macro_vertex_coord_id_0comp1,
       walberla::float64           macro_vertex_coord_id_1comp0,
       walberla::float64           macro_vertex_coord_id_1comp1,
       walberla::float64           macro_vertex_coord_id_2comp0,
       walberla::float64           macro_vertex_coord_id_2comp1,
       int64_t                     micro_edges_per_macro_edge,
       walberla::float64           micro_edges_per_macro_edge_float ) const;

   std::shared_ptr< P2Function< walberla::float64 > > invDiag_;
   BoundaryCondition                                  boundaryCondition_;
   BoundaryUID                                        P2ElementwiseBoundaryMass_boundary_uid_;
};

} // namespace operatorgeneration

} // namespace hyteg
