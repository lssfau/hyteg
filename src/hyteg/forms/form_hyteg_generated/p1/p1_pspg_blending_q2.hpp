/*
* Copyright (c) 2023-2025 Andreas Burkhart
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

// This file has been generated with the AHFC. If buggy try fixing the generator itself.

#pragma once
#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"

namespace hyteg {
namespace forms {

using walberla::real_t;

/// PSPG stabilisation.
/// 
/// Weak formulation:
/// 
///     u: trial function (space: P1)
///     v: test function  (space: P1)
/// 
///     ∫ τ ∇u · ∇v
/// 
///     where tau is an element-wise factor given by
/// 
///     2D: 
///         τ = -CellVolume(triangle) / 5               
///     3D:
///         τ = -(CellVolume(tetrahedron))**(2/3) / 12
/// 
/// See e.g. 
/// 
///     Brezzi, F., & Douglas, J. (1988). 
///     Stabilized mixed methods for the Stokes problem. 
///     Numerische Mathematik, 53, 225-235.
/// 
/// or 
/// 
///     Hughes, T. J., Franca, L. P., & Balestra, M. (1986). 
///     A new finite element formulation for computational fluid dynamics: V. 
///     Circumventing the Babuška-Brezzi condition: A stable Petrov-Galerkin formulation of the Stokes problem accommodating
///     equal-order interpolations. 
///     Computer Methods in Applied Mechanics and Engineering, 59(1), 85-99.
/// 
/// for details.
/// 
/// Blending: True
/// Quadrature degree: 2
class p1_pspg_blending_q2 : public P1FormHyTeG {
public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override; 

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override; 

};

} // namespace forms
} // namespace hyteg