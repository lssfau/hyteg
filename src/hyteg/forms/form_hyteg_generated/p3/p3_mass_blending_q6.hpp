/*
* Copyright (c) 2017-2026 Nils Kohl, Daniel Bauer, Fabian Böhm, Marcus Mohr.
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

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P3FormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p3_mass_blending_q6
/// - description: 
/// - trial space: Lagrange, degree: 3
/// - test space:  Lagrange, degree: 3
///
class p3_mass_blending_q6 : public P3FormHyTeG
{



 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3, spacedim: 2
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Griener-Schmid 1 | points: 10, degree: 6
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///                                              650    1522       0       0     11           1066                10              0
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

 private:

   void Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const;

};

} // namespace forms
} // namespace hyteg
