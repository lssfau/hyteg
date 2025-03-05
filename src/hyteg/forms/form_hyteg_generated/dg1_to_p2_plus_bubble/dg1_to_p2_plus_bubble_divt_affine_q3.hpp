/*
* Copyright (c) 2025 Marcus Mohr.
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
* Interfaces needed slight modifications, as HOG support for DG1/P2+
* pair is not complete, yet.
*/

#pragma once

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/DG1ToP2PlusBubbleFormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        dg1_to_p2_plus_bubble_divt_0_affine_q3
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  P2PlusBubble
///
class dg1_to_p2_plus_bubble_divt_0_affine_q3 : public DG1ToP2PlusBubbleFormHyTeG
{



 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3, spacedim: 2
   /// - element matrix dimensions (rows, cols): (7, 3)
   /// - quadrature rule:                        Dunavant 3 | points: 4, degree: 3
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///                                              119     169       1       0      1            131                 0              0
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 7, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3, spacedim: 2
   /// - element matrix dimensions (rows, cols): (7, 3)
   /// - quadrature rule:                        Dunavant 3 | points: 4, degree: 3
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///                                               18      32       1       0      1             29                 0              0
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        dg1_to_p2_plus_bubble_divt_1_affine_q3
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  P2PlusBubble
///
class dg1_to_p2_plus_bubble_divt_1_affine_q3 : public DG1ToP2PlusBubbleFormHyTeG
{



 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3, spacedim: 2
   /// - element matrix dimensions (rows, cols): (7, 3)
   /// - quadrature rule:                        Dunavant 3 | points: 4, degree: 3
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///                                              119     169       1       0      1            131                 0              0
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 7, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3, spacedim: 2
   /// - element matrix dimensions (rows, cols): (7, 3)
   /// - quadrature rule:                        Dunavant 3 | points: 4, degree: 3
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls    unknown_ops
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------  -------------
   ///                                               18      32       1       0      1             29                 0              0
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        dg1_to_p2_plus_bubble_divt_2_affine_q3
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  P2PlusBubble
///
class dg1_to_p2_plus_bubble_divt_2_affine_q3 : public DG1ToP2PlusBubbleFormHyTeG
{



 public:

};

} // namespace forms
} // namespace hyteg
