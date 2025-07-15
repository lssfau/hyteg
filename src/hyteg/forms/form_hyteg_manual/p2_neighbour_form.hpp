/*
 * Copyright (c) 2023-2024 Andreas Burkhart.
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

#pragma once

#include "hyteg/forms/form_hyteg_base/P2FormHyTeG.hpp"
#include "hyteg/geometry/GeometryMap.hpp"

namespace hyteg {
namespace forms {

/// Can be used to count the number of neighbours of a DoF
class p2_neighbour_form : public P2FormHyTeG
{
 public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
      elMat( 0, 3 ) = real_c( 0 );
      elMat( 0, 4 ) = real_c( 0 );
      elMat( 0, 5 ) = real_c( 0 );
      elMat( 1, 0 ) = real_c( 0 );
      elMat( 1, 1 ) = real_c( 1 );
      elMat( 1, 2 ) = real_c( 0 );
      elMat( 1, 3 ) = real_c( 0 );
      elMat( 1, 4 ) = real_c( 0 );
      elMat( 1, 5 ) = real_c( 0 );
      elMat( 2, 0 ) = real_c( 0 );
      elMat( 2, 1 ) = real_c( 0 );
      elMat( 2, 2 ) = real_c( 1 );
      elMat( 2, 3 ) = real_c( 0 );
      elMat( 2, 4 ) = real_c( 0 );
      elMat( 2, 5 ) = real_c( 0 );
      elMat( 3, 0 ) = real_c( 0 );
      elMat( 3, 1 ) = real_c( 0 );
      elMat( 3, 2 ) = real_c( 0 );
      elMat( 3, 3 ) = real_c( 1 );
      elMat( 3, 4 ) = real_c( 0 );
      elMat( 3, 5 ) = real_c( 0 );
      elMat( 4, 0 ) = real_c( 0 );
      elMat( 4, 1 ) = real_c( 0 );
      elMat( 4, 2 ) = real_c( 0 );
      elMat( 4, 3 ) = real_c( 0 );
      elMat( 4, 4 ) = real_c( 1 );
      elMat( 4, 5 ) = real_c( 0 );
      elMat( 5, 0 ) = real_c( 0 );
      elMat( 5, 1 ) = real_c( 0 );
      elMat( 5, 2 ) = real_c( 0 );
      elMat( 5, 3 ) = real_c( 0 );
      elMat( 5, 4 ) = real_c( 0 );
      elMat( 5, 5 ) = real_c( 1 );
   }

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
      elMat( 0, 3 ) = real_c( 0 );
      elMat( 0, 4 ) = real_c( 0 );
      elMat( 0, 5 ) = real_c( 0 );
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
      elMat( 0, 3 ) = real_c( 0 );
      elMat( 0, 4 ) = real_c( 0 );
      elMat( 0, 5 ) = real_c( 0 );
      elMat( 0, 6 ) = real_c( 0 );
      elMat( 0, 7 ) = real_c( 0 );
      elMat( 0, 8 ) = real_c( 0 );
      elMat( 0, 9 ) = real_c( 0 );
      elMat( 1, 0 ) = real_c( 0 );
      elMat( 1, 1 ) = real_c( 1 );
      elMat( 1, 2 ) = real_c( 0 );
      elMat( 1, 3 ) = real_c( 0 );
      elMat( 1, 4 ) = real_c( 0 );
      elMat( 1, 5 ) = real_c( 0 );
      elMat( 1, 6 ) = real_c( 0 );
      elMat( 1, 7 ) = real_c( 0 );
      elMat( 1, 8 ) = real_c( 0 );
      elMat( 1, 9 ) = real_c( 0 );
      elMat( 2, 0 ) = real_c( 0 );
      elMat( 2, 1 ) = real_c( 0 );
      elMat( 2, 2 ) = real_c( 1 );
      elMat( 2, 3 ) = real_c( 0 );
      elMat( 2, 4 ) = real_c( 0 );
      elMat( 2, 5 ) = real_c( 0 );
      elMat( 2, 6 ) = real_c( 0 );
      elMat( 2, 7 ) = real_c( 0 );
      elMat( 2, 8 ) = real_c( 0 );
      elMat( 2, 9 ) = real_c( 0 );
      elMat( 3, 0 ) = real_c( 0 );
      elMat( 3, 1 ) = real_c( 0 );
      elMat( 3, 2 ) = real_c( 0 );
      elMat( 3, 3 ) = real_c( 1 );
      elMat( 3, 4 ) = real_c( 0 );
      elMat( 3, 5 ) = real_c( 0 );
      elMat( 3, 6 ) = real_c( 0 );
      elMat( 3, 7 ) = real_c( 0 );
      elMat( 3, 8 ) = real_c( 0 );
      elMat( 3, 9 ) = real_c( 0 );
      elMat( 4, 0 ) = real_c( 0 );
      elMat( 4, 1 ) = real_c( 0 );
      elMat( 4, 2 ) = real_c( 0 );
      elMat( 4, 3 ) = real_c( 0 );
      elMat( 4, 4 ) = real_c( 1 );
      elMat( 4, 5 ) = real_c( 0 );
      elMat( 4, 6 ) = real_c( 0 );
      elMat( 4, 7 ) = real_c( 0 );
      elMat( 4, 8 ) = real_c( 0 );
      elMat( 4, 9 ) = real_c( 0 );
      elMat( 5, 0 ) = real_c( 0 );
      elMat( 5, 1 ) = real_c( 0 );
      elMat( 5, 2 ) = real_c( 0 );
      elMat( 5, 3 ) = real_c( 0 );
      elMat( 5, 4 ) = real_c( 0 );
      elMat( 5, 5 ) = real_c( 1 );
      elMat( 5, 6 ) = real_c( 0 );
      elMat( 5, 7 ) = real_c( 0 );
      elMat( 5, 8 ) = real_c( 0 );
      elMat( 5, 9 ) = real_c( 0 );
      elMat( 6, 0 ) = real_c( 0 );
      elMat( 6, 1 ) = real_c( 0 );
      elMat( 6, 2 ) = real_c( 0 );
      elMat( 6, 3 ) = real_c( 0 );
      elMat( 6, 4 ) = real_c( 0 );
      elMat( 6, 5 ) = real_c( 0 );
      elMat( 6, 6 ) = real_c( 1 );
      elMat( 6, 7 ) = real_c( 0 );
      elMat( 6, 8 ) = real_c( 0 );
      elMat( 6, 9 ) = real_c( 0 );
      elMat( 7, 0 ) = real_c( 0 );
      elMat( 7, 1 ) = real_c( 0 );
      elMat( 7, 2 ) = real_c( 0 );
      elMat( 7, 3 ) = real_c( 0 );
      elMat( 7, 4 ) = real_c( 0 );
      elMat( 7, 5 ) = real_c( 0 );
      elMat( 7, 6 ) = real_c( 0 );
      elMat( 7, 7 ) = real_c( 1 );
      elMat( 7, 8 ) = real_c( 0 );
      elMat( 7, 9 ) = real_c( 0 );
      elMat( 8, 0 ) = real_c( 0 );
      elMat( 8, 1 ) = real_c( 0 );
      elMat( 8, 2 ) = real_c( 0 );
      elMat( 8, 3 ) = real_c( 0 );
      elMat( 8, 4 ) = real_c( 0 );
      elMat( 8, 5 ) = real_c( 0 );
      elMat( 8, 6 ) = real_c( 0 );
      elMat( 8, 7 ) = real_c( 0 );
      elMat( 8, 8 ) = real_c( 1 );
      elMat( 8, 9 ) = real_c( 0 );
      elMat( 9, 0 ) = real_c( 0 );
      elMat( 9, 1 ) = real_c( 0 );
      elMat( 9, 2 ) = real_c( 0 );
      elMat( 9, 3 ) = real_c( 0 );
      elMat( 9, 4 ) = real_c( 0 );
      elMat( 9, 5 ) = real_c( 0 );
      elMat( 9, 6 ) = real_c( 0 );
      elMat( 9, 7 ) = real_c( 0 );
      elMat( 9, 8 ) = real_c( 0 );
      elMat( 9, 9 ) = real_c( 1 );
   }

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
      elMat( 0, 3 ) = real_c( 0 );
      elMat( 0, 4 ) = real_c( 0 );
      elMat( 0, 5 ) = real_c( 0 );
      elMat( 0, 6 ) = real_c( 0 );
      elMat( 0, 7 ) = real_c( 0 );
      elMat( 0, 8 ) = real_c( 0 );
      elMat( 0, 9 ) = real_c( 0 );
   }
};

} // namespace forms
} // namespace hyteg
