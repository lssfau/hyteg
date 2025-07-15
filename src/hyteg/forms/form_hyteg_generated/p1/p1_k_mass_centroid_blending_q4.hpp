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

/// Mass operator with a scalar coefficient.
/// 
/// Weak formulation:
/// 
///     u: trial function (space: P1)
///     v: test function  (space: P1)
///     k: arbitrary scalar std::function coefficient, provides the following additional inputs as parameters to the std::function:
///       element centroid position X, element centroid position Y, element centroid position Z
/// 
///     ∫ k uv
/// 
/// Blending: True
/// Quadrature degree: 4
class p1_k_mass_centroid_blending_q4 : public P1FormHyTeG {
public:
   p1_k_mass_centroid_blending_q4() { WALBERLA_ABORT("Not implemented."); }
   p1_k_mass_centroid_blending_q4(std::function< real_t ( const Point3D&, real_t, real_t, real_t ) > k)
   : k_(k)
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override; 

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override; 

private:
   std::function< real_t ( const Point3D&, real_t, real_t, real_t ) > k_;
};

} // namespace forms
} // namespace hyteg