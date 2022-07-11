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
 * Software:
 *
 * - quadpy version: 0.16.5
 *
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#pragma once

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P2FormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_invk_mass_affine_q6
/// - description: Implements a bilinear form of type (1/k(x) * phi, psi) where phi is a trial function, psi a test function and k = k(x) a scalar, spatially varying coefficient.
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_invk_mass_affine_q6 : public P2FormHyTeG
{

 public:

   p2_invk_mass_affine_q6() { WALBERLA_ABORT("Not implemented."); }

   p2_invk_mass_affine_q6( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_k, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_k )
   : callback_Scalar_Variable_Coefficient_3D_k(_callback_Scalar_Variable_Coefficient_3D_k)
   , callback_Scalar_Variable_Coefficient_2D_k(_callback_Scalar_Variable_Coefficient_2D_k)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_k;
   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_k;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Witherden-Vincent 6 | points: 12, degree: 6, test tolerance: 2.871e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              392     572      12       0      1            325                12
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Witherden-Vincent 6 | points: 12, degree: 6, test tolerance: 2.871e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              227     283      12       0      1            165                12
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             4568    4960      46       0      1           2343                46
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Witherden-Vincent 8 | points: 46, degree: 8, test tolerance: 1.332e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1991    1877      46       0      1            967                46
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

   bool assemble2D() const override { return true; }

   bool assembly2DDefined() const override { return true; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

} // namespace forms
} // namespace hyteg