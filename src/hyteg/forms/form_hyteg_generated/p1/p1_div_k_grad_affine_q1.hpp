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
/// - name:        p1_div_k_grad_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_div_k_grad_affine_q1 : public P1FormHyTeG
{

 public:

   p1_div_k_grad_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_div_k_grad_affine_q1( std::function< real_t ( const Point3D & ) > _callback3D, std::function< real_t ( const Point3D & ) > _callback2D )
   : callback3D(_callback3D)
   , callback2D(_callback2D)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback3D;
   std::function< real_t ( const Point3D & ) > callback2D;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Centroid rule | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               28      55       2       0      1             48                 1
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Centroid rule | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               23      32       1       0      1             29                 1
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 0 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               93     160       2       0      1            118                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 0 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               75     107       1       0      1             74                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const;

 private:

   void Scalar_Variable_Coefficient_2D( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

} // namespace forms
} // namespace hyteg
