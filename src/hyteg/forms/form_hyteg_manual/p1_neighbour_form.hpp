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

#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"
#include "hyteg/geometry/GeometryMap.hpp"

namespace hyteg {
namespace forms {

/// Can be used to count the number of neighbours of a DoF   
class p1_neighbour_form : public P1FormHyTeG
{
 public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
      elMat( 1, 0 ) = real_c( 0 );
      elMat( 1, 1 ) = real_c( 1 );
      elMat( 1, 2 ) = real_c( 0 );
      elMat( 2, 0 ) = real_c( 0 );
      elMat( 2, 1 ) = real_c( 0 );
      elMat( 2, 2 ) = real_c( 1 );
   }

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
      elMat( 0, 3 ) = real_c( 0 );
      elMat( 1, 0 ) = real_c( 0 );
      elMat( 1, 1 ) = real_c( 1 );
      elMat( 1, 2 ) = real_c( 0 );
      elMat( 1, 3 ) = real_c( 0 );
      elMat( 2, 0 ) = real_c( 0 );
      elMat( 2, 1 ) = real_c( 0 );
      elMat( 2, 2 ) = real_c( 1 );
      elMat( 2, 3 ) = real_c( 0 );
      elMat( 3, 0 ) = real_c( 0 );
      elMat( 3, 1 ) = real_c( 0 );
      elMat( 3, 2 ) = real_c( 0 );
      elMat( 3, 3 ) = real_c( 1 );
   }

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override
   {
      elMat( 0, 0 ) = real_c( 1 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );
      elMat( 0, 3 ) = real_c( 0 );
   }
};

} // namespace forms
} // namespace hyteg
