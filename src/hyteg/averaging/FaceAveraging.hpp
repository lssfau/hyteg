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

#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/types/Averaging.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {
///
/// \brief Evaluate the average value depending on the \param averagingMethod
///
/// Computes the average using the function values at vertices and the centroid
///
/// \param microTriangles  Vertices of the triangle face
/// \param valueTets       FE function values at those vertices
/// \param averagingMethod Averaging method that should be used
///                        to compute a single average for the face
///
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

   Point3D centroid = ( microTet0 + microTet1 + microTet2 ) / 3.0;

   auto xLocal          = vertexdof::macroface::transformToLocalTri( microTet0, microTet1, microTet2, centroid );
   auto valueAtCentroid = valueTet0 * ( real_c( 1.0 ) - xLocal[0] - xLocal[1] ) + valueTet1 * xLocal[0] + valueTet2 * xLocal[1];

   if ( averagingMethod == hyteg::AveragingType::ARITHMETIC )
   {
      return ( valueAtCentroid + valueTet0 + valueTet1 + valueTet2 ) / 4.0;
   }
   else
   {
      WALBERLA_ABORT( "Not implemented" );
   }
}

} // namespace hyteg