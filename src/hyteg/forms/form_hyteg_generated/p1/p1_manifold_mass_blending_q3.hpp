/*
 * Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_manifold_mass_blending_q3
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_manifold_mass_blending_q3 : public P1FormHyTeG
{



 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       embedded triangle, dim: 3, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Xiao-Gimbutas | points: 6, degree: 3
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                            10590   13004       0       6      6           2399                 6
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       embedded triangle, dim: 3, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Xiao-Gimbutas | points: 6, degree: 3
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                            10575   12986       0       6      6           2384                 6
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override
   {
      WALBERLA_UNUSED( coords );
      WALBERLA_UNUSED( elMat );
      WALBERLA_ABORT( "The manifold operator can only be used on 2D-manifolds in 3D-space." );
   };

 private:

   void Blending_DF_EmbeddedTriangle_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const;

};

} // namespace forms
} // namespace hyteg
