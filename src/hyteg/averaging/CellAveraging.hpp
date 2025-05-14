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

#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/types/Averaging.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {
///
/// \brief Evaluate the average value depending on the \param averagingMethod
///
/// Possible to compute the averaging of the function values at the vertex and centroid.
/// Also a 5 point quadrature is hard coded 
/// where the function is evaluated and can be used for averaging.
///
/// \param microTets       Vertices of the tetrahedron
/// \param valueTets       FE function values at those vertices
/// \param averagingMethod Averaging method that should be used 
///                        to compute a single average for the cell
///
///
inline real_t evaluateSampledAverage( std::array< Point3D, 4 > microTets,
                                      std::array< real_t, 4 >  valueTets,
                                      hyteg::AveragingType     averagingMethod )
{
   Point3D microTet0 = microTets[0];
   Point3D microTet1 = microTets[1];
   Point3D microTet2 = microTets[2];
   Point3D microTet3 = microTets[3];

   real_t valueTet0 = valueTets[0];
   real_t valueTet1 = valueTets[1];
   real_t valueTet2 = valueTets[2];
   real_t valueTet3 = valueTets[3];

   Point3D centroid = ( microTet0 + microTet1 + microTet2 + microTet3 ) / 4.0;

   std::function< real_t( const Point3D& ) > locallyEvaluate = [&]( const Point3D& x ) {
      return valueTet0 * ( real_c( 1.0 ) - x[0] - x[1] - x[2] ) + valueTet1 * x[0] + valueTet2 * x[1] + valueTet3 * x[2];
   };

   // This 5 point quadrature is hard coded for now
   Point3D qp1 = Point3D( 0.25, 0.25, 0.25 );
   Point3D qp2 = Point3D( 0.16666667, 0.16666667, 0.5 );
   Point3D qp3 = Point3D( 0.16666667, 0.5, 0.16666667 );
   Point3D qp4 = Point3D( 0.5, 0.16666667, 0.16666667 );
   Point3D qp5 = Point3D( 0.16666667, 0.16666667, 0.16666667 );

   real_t val1 = locallyEvaluate( qp1 );
   real_t val2 = locallyEvaluate( qp2 );
   real_t val3 = locallyEvaluate( qp3 );
   real_t val4 = locallyEvaluate( qp4 );
   real_t val5 = locallyEvaluate( qp5 );

   auto xLocal = vertexdof::macrocell::detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, centroid );
   auto value  = valueTet0 * ( real_c( 1.0 ) - xLocal[0] - xLocal[1] - xLocal[2] ) + valueTet1 * xLocal[0] +
                valueTet2 * xLocal[1] + valueTet3 * xLocal[2];

   if ( averagingMethod == hyteg::AveragingType::ARITHMETIC )
   {
      return ( value + valueTet0 + valueTet1 + valueTet2 + valueTet3 ) / 5.0;
   }
   else if ( averagingMethod == hyteg::AveragingType::ARITHMETIC_QP )
   {
      return ( val1 + val2 + val3 + val4 + val5 ) / 5.0;
   }
   else if ( averagingMethod == hyteg::AveragingType::HARMONIC )
   {
      return 5.0 / ( ( 1.0 / value ) + ( 1.0 / valueTet0 ) + ( 1.0 / valueTet1 ) + ( 1.0 / valueTet2 ) + ( 1.0 / valueTet3 ) );
   }
   else if ( averagingMethod == hyteg::AveragingType::HARMONIC_QP )
   {
      return 5.0 / ( ( 1.0 / val1 ) + ( 1.0 / val2 ) + ( 1.0 / val3 ) + ( 1.0 / val4 ) + ( 1.0 / val5 ) );
   }
   else if ( averagingMethod == hyteg::AveragingType::GEOMETRIC )
   {
      return std::pow( value, 1.0 / 5.0 ) * std::pow( valueTet0, 1.0 / 5.0 ) * std::pow( valueTet1, 1.0 / 5.0 ) *
             std::pow( valueTet2, 1.0 / 5.0 ) * std::pow( valueTet3, 1.0 / 5.0 );
   }
   else if ( averagingMethod == hyteg::AveragingType::GEOMETRIC_QP )
   {
      return std::pow( val1, 1.0 / 5.0 ) * std::pow( val2, 1.0 / 5.0 ) * std::pow( val3, 1.0 / 5.0 ) *
             std::pow( val4, 1.0 / 5.0 ) * std::pow( val5, 1.0 / 5.0 );
   }
   else
   {
      WALBERLA_ABORT( "Not implemented" );
   }
}
} // namespace hyteg