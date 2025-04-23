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
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/types/types.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Diffusion operator without coefficients.
///
/// Geometry map: IcosahedralShellMap
///
/// Weak formulation
///
///     u: trial function (space: Lagrange, degree: 1)
///     v: test function  (space: Lagrange, degree: 1)
///
///     ∫ ∇u : ∇v
///
///     Note that the double contraction (:) reduces to the dot product for scalar function spaces, i.e. the form becomes
///
///     ∫ ∇u · ∇v

class P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab
: public Operator< P1Function< double >, P1Function< double > >,
  public OperatorWithInverseDiagonal< P1Function< double > >
{
 public:
   P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab(
       const std::shared_ptr< PrimitiveStorage >& storage,
       size_t                                     minLevel,
       size_t                                     maxLevel );

   void apply( const P1Function< double >& src,
               const P1Function< double >& dst,
               uint_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const;

   void computeInverseDiagonalOperatorValues();

   std::shared_ptr< P1Function< double > > getInverseDiagonalValues() const;

 protected:
 private:
   /// Integral: P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     apply
   /// - loop strategy:   CUBES
   /// - quadrature rule: Keast 0 | points: 1, degree: 1
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   float_adds    float_muls    float_divs    int_adds    int_muls    int_divs    calls    branches    loops_with_dynamic_bounds
   /// ------------  ------------  ------------  ----------  ----------  ----------  -------  ----------  ---------------------------
   ///         1578          2910           234        1506         726         216       48           0                            0
   void apply_P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab_macro_3D(
       double* RESTRICT const _data_dst,
       double* RESTRICT const _data_src,
       const double           forVertex_0,
       const double           forVertex_1,
       const double           forVertex_2,
       const double           macro_vertex_coord_id_0comp0,
       const double           macro_vertex_coord_id_0comp1,
       const double           macro_vertex_coord_id_0comp2,
       const double           macro_vertex_coord_id_1comp0,
       const double           macro_vertex_coord_id_1comp1,
       const double           macro_vertex_coord_id_1comp2,
       const double           macro_vertex_coord_id_2comp0,
       const double           macro_vertex_coord_id_2comp1,
       const double           macro_vertex_coord_id_2comp2,
       const double           macro_vertex_coord_id_3comp0,
       const double           macro_vertex_coord_id_3comp1,
       const double           macro_vertex_coord_id_3comp2,
       const int64_t          micro_edges_per_macro_edge,
       const double           micro_edges_per_macro_edge_float,
       const double           radRayVertex,
       const double           radRefVertex,
       const double           rayVertex_0,
       const double           rayVertex_1,
       const double           rayVertex_2,
       const double           refVertex_0,
       const double           refVertex_1,
       const double           refVertex_2,
       const double           thrVertex_0,
       const double           thrVertex_1,
       const double           thrVertex_2 ) const;

   /// Integral: P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab
   /// - volume element:  tetrahedron, dim: 3, vertices: 4, spacedim: 3
   /// - kernel type:     computeInverseDiagonalOperatorValues
   /// - loop strategy:   CUBES
   /// - quadrature rule: Keast 0 | points: 1, degree: 1
   /// - blending map:    IcosahedralShellMap
   /// - operations per element:
   ///   float_adds    float_muls    float_divs    int_adds    int_muls    int_divs    calls    branches    loops_with_dynamic_bounds
   /// ------------  ------------  ------------  ----------  ----------  ----------  -------  ----------  ---------------------------
   ///         1542          2886           234        1260         774         144       48           0                            0
   void
       computeInverseDiagonalOperatorValues_P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab_macro_3D(
           double* RESTRICT const _data_invDiag_,
           const double           forVertex_0,
           const double           forVertex_1,
           const double           forVertex_2,
           const double           macro_vertex_coord_id_0comp0,
           const double           macro_vertex_coord_id_0comp1,
           const double           macro_vertex_coord_id_0comp2,
           const double           macro_vertex_coord_id_1comp0,
           const double           macro_vertex_coord_id_1comp1,
           const double           macro_vertex_coord_id_1comp2,
           const double           macro_vertex_coord_id_2comp0,
           const double           macro_vertex_coord_id_2comp1,
           const double           macro_vertex_coord_id_2comp2,
           const double           macro_vertex_coord_id_3comp0,
           const double           macro_vertex_coord_id_3comp1,
           const double           macro_vertex_coord_id_3comp2,
           const int64_t          micro_edges_per_macro_edge,
           const double           micro_edges_per_macro_edge_float,
           const double           radRayVertex,
           const double           radRefVertex,
           const double           rayVertex_0,
           const double           rayVertex_1,
           const double           rayVertex_2,
           const double           refVertex_0,
           const double           refVertex_1,
           const double           refVertex_2,
           const double           thrVertex_0,
           const double           thrVertex_1,
           const double           thrVertex_2 ) const;

   std::shared_ptr< P1Function< double > > invDiag_;
};

} // namespace operatorgeneration

} // namespace hyteg
