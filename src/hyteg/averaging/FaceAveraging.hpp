/*
* Copyright (c) 2025 Ponsuganth Ilangovan.
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

#include "core/Abort.h"

#include "hyteg/types/Averaging.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

inline real_t evaluateSampledAverage( std::array< Point3D, 3 > microTriangles,
                                      std::array< real_t, 3 >  valueTriangles,
                                      hyteg::AveragingType     averagingMethod )
{
   Point3D microTet0 = microTriangles[0];
   Point3D microTet1 = microTriangles[1];
   Point3D microTet2 = microTriangles[2];

   real_t valueTet0 = valueTriangles[0];
   real_t valueTet1 = valueTriangles[1];
   real_t valueTet2 = valueTriangles[2];

   Point3D coordinates = ( microTet0 + microTet1 + microTet2 ) / 3.0;

   if ( averagingMethod == hyteg::AveragingType::ARITHMETIC )
   {
      return ( valueTet0 + valueTet1 + valueTet2 ) / 3.0;
   }
   else
   {
      WALBERLA_ABORT( "Not implemented" );
   }
}

} // namespace hyteg