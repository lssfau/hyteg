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
#include "hyteg/forms/form_hyteg_base/P1ToP2FormHyTeG.hpp"

namespace hyteg {
namespace forms {

using walberla::real_t;

/// Mass operator with a scalar coefficient.
/// 
/// Weak formulation:
/// 
///     u: trial function (space: P1)
///     v: test function  (space: P2)
///     k: arbitrary scalar std::function coefficient
/// 
///     ∫ k uv
/// 
/// Blending: True
/// Quadrature degree: 5
class p1_to_p2_k_mass_blending_q5 : public P1ToP2FormHyTeG {
public:
   p1_to_p2_k_mass_blending_q5() { WALBERLA_ABORT("Not implemented."); }
   p1_to_p2_k_mass_blending_q5(std::function< real_t ( const Point3D& ) > k)
   : k_(k)
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 3 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override; 

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 4 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override; 

private:
   std::function< real_t ( const Point3D& ) > k_;
};

} // namespace forms
} // namespace hyteg