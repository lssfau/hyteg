/*
* Copyright (c) 2023-2024 Andreas Burkhart
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
#include "hyteg/forms/form_hyteg_base/P2ToP1FormHyTeG.hpp"

namespace hyteg {
namespace forms {

using walberla::real_t;

/// Divergence.
/// 
/// Weak formulation:
/// 
///     u: trial function (vectorial space: P2)
///     v: test function  (scalar space:    P1)
/// 
///     ∫ - ( ∇ · u ) v
/// 
/// Blending: True
/// Quadrature degree: 6
class p2_to_p1_div_0_blending_q6 : public P2ToP1FormHyTeG {
public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override; 

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override; 

};

/// Divergence.
/// 
/// Weak formulation:
/// 
///     u: trial function (vectorial space: P2)
///     v: test function  (scalar space:    P1)
/// 
///     ∫ - ( ∇ · u ) v
/// 
/// Blending: True
/// Quadrature degree: 6
class p2_to_p1_div_1_blending_q6 : public P2ToP1FormHyTeG {
public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override; 

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override; 

};

/// Divergence.
/// 
/// Weak formulation:
/// 
///     u: trial function (vectorial space: P2)
///     v: test function  (scalar space:    P1)
/// 
///     ∫ - ( ∇ · u ) v
/// 
/// Blending: True
/// Quadrature degree: 6
class p2_to_p1_div_2_blending_q6 : public P2ToP1FormHyTeG {
public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override; 

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const override; 

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override; 

};

} // namespace forms
} // namespace hyteg