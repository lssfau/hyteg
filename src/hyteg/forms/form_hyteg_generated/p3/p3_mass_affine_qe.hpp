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
/// - name:        p3_mass_affine_qe
/// - description: 
/// - trial space: Lagrange, degree: 3
/// - test space:  Lagrange, degree: 3
///
class p3_mass_affine_qe : public P3FormHyTeG
{



 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3, spacedim: 2
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        exact
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///                                                5      44       0       0      1            247                 0              0
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

};

} // namespace forms
} // namespace hyteg
