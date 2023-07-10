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
#include "hyteg/forms/N1E1Form.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        n1e1_linear_form_blending_q6
/// - description: Implements a linear form of type: (k(x), psi) where psi a test function and k = k(x) a vectorial, external function.
/// - trial space: N1E1
/// - test space:  N1E1
///
class n1e1_linear_form_blending_q6 : public n1e1::N1E1Form
{

 public:

   n1e1_linear_form_blending_q6() { WALBERLA_ABORT("Not implemented."); }

   n1e1_linear_form_blending_q6( std::function< Point3D ( const Point3D & ) > _callback_Vector_Variable_Coefficient_3D_k )
   : callback_Vector_Variable_Coefficient_3D_k(_callback_Vector_Variable_Coefficient_3D_k)
   {}

 private:

   std::function< Point3D ( const Point3D & ) > callback_Vector_Variable_Coefficient_3D_k;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Xiao-Gimbutas 6 | points: 23, degree: 6, test tolerance: 7.137e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             2317    3234      24       0     24           1643                69
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

 private:

   void Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const;

   void Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const;

   void Vector_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const;

};

} // namespace forms
} // namespace hyteg
