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
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Radial Gradient Operator.
///
/// Geometry map: IcosahedralShellMap
///
/// Weak formulation
///
///     u: trial function (space: Lagrange, degree: 2)
///     v: test function  (space: Lagrange, degree: 2)
///
///     ∫ (∇𝑢 ⋅ n) v ds

class P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator
: public Operator< P2Function< walberla::float64 >, P2Function< walberla::float64 > >,
  public OperatorWithInverseDiagonal< P2Function< walberla::float64 > >
{
 public:
   P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator(
       const std::shared_ptr< PrimitiveStorage >& storage,
       size_t                                     minLevel,
       size_t                                     maxLevel,
       const P2Function< walberla::float64 >&     _T,
       BoundaryCondition                          boundaryCondition,
       BoundaryUID                                P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_boundary_uid );

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
   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_0
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1276    1743      50      30     36              0                 0              0
   void apply_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_0_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_dstEdge,
       walberla::float64* RESTRICT _data_dstVertex,
       walberla::float64* RESTRICT _data_srcEdge,
       walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_1
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1273    1743      50      30     36              0                 0              0
   void apply_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_1_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_dstEdge,
       walberla::float64* RESTRICT _data_dstVertex,
       walberla::float64* RESTRICT _data_srcEdge,
       walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_2
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1276    1743      50      30     36              0                 0              0
   void apply_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_2_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_dstEdge,
       walberla::float64* RESTRICT _data_dstVertex,
       walberla::float64* RESTRICT _data_srcEdge,
       walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_3
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1279    1743      50      30     36              0                 0              0
   void apply_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_3_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_dstEdge,
       walberla::float64* RESTRICT _data_dstVertex,
       walberla::float64* RESTRICT _data_srcEdge,
       walberla::float64* RESTRICT _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_0
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1236    1707      50      30     36              0                 0              3
   void toMatrix_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_0_macro_3D(
       walberla::float64* RESTRICT          _data_TEdge,
       walberla::float64* RESTRICT          _data_TVertex,
       idx_t* RESTRICT                      _data_dstEdge,
       idx_t* RESTRICT                      _data_dstVertex,
       idx_t* RESTRICT                      _data_srcEdge,
       idx_t* RESTRICT                      _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_1
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1233    1707      50      30     36              0                 0              3
   void toMatrix_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_1_macro_3D(
       walberla::float64* RESTRICT          _data_TEdge,
       walberla::float64* RESTRICT          _data_TVertex,
       idx_t* RESTRICT                      _data_dstEdge,
       idx_t* RESTRICT                      _data_dstVertex,
       idx_t* RESTRICT                      _data_srcEdge,
       idx_t* RESTRICT                      _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_2
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1236    1707      50      30     36              0                 0              3
   void toMatrix_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_2_macro_3D(
       walberla::float64* RESTRICT          _data_TEdge,
       walberla::float64* RESTRICT          _data_TVertex,
       idx_t* RESTRICT                      _data_dstEdge,
       idx_t* RESTRICT                      _data_dstVertex,
       idx_t* RESTRICT                      _data_srcEdge,
       idx_t* RESTRICT                      _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_3
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     toMatrix
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1239    1707      50      30     36              0                 0              3
   void toMatrix_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_3_macro_3D(
       walberla::float64* RESTRICT          _data_TEdge,
       walberla::float64* RESTRICT          _data_TVertex,
       idx_t* RESTRICT                      _data_dstEdge,
       idx_t* RESTRICT                      _data_dstVertex,
       idx_t* RESTRICT                      _data_srcEdge,
       idx_t* RESTRICT                      _data_srcVertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_0
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1117    1557      50      30     36              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_0_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_1
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1114    1557      50      30     36              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_1_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_2
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1117    1557      50      30     36              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_2_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   /// Integral: P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_3
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   BOUNDARY
   /// - quadrature rule: Dunavant 4 | points: 6, degree: 4
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   /// ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///   1120    1557      50      30     36              0                 0              0
   void computeInverseDiagonalOperatorValues_P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_facet_id_3_macro_3D(
       walberla::float64* RESTRICT _data_TEdge,
       walberla::float64* RESTRICT _data_TVertex,
       walberla::float64* RESTRICT _data_invDiag_Edge,
       walberla::float64* RESTRICT _data_invDiag_Vertex,
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

   std::shared_ptr< P2Function< walberla::float64 > > invDiag_;
   P2Function< walberla::float64 >                    T;
   BoundaryCondition                                  boundaryCondition_;
   BoundaryUID                                        P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator_boundary_uid_;
};

} // namespace operatorgeneration

} // namespace hyteg
