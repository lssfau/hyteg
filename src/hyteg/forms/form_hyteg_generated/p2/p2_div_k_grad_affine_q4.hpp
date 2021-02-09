/*
 * Copyright (c) 2017-2021 Nils Kohl.
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
 * Software:
 *
 * - quadpy version: 0.16.5
 *
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#pragma once

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"
#include "hyteg/forms/form_hyteg_base/P2FormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_div_k_grad_affine_q4
/// - description: 
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_div_k_grad_affine_q4 : public P2FormHyTeG
{

 public:

   p2_div_k_grad_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_div_k_grad_affine_q4( std::function< real_t ( const Point3D & ) > _callback2D, std::function< real_t ( const Point3D & ) > _callback3D )
   : callback2D(_callback2D)
   , callback3D(_callback3D)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback2D;
   std::function< real_t ( const Point3D & ) > callback3D;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    abs    assignments    function_calls
   ///                                           ------  ------  ------  -----  -------------  ----------------
   ///                                              589     999       2      1            473                 6
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    abs    assignments    function_calls
   ///                                           ------  ------  ------  -----  -------------  ----------------
   ///                                             4261    5928       2      1           2062                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

} // namespace forms
} // namespace hyteg
