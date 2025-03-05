/*
 * Copyright (c) 2017-2025 Marcus Mohr.
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

#include <sstream>
#include <string>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"

#include "hyteg/forms/form_hyteg_generated/dg1_to_p2_plus_bubble/dg1_to_p2_plus_bubble_divt_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p2_plus_bubble_to_dg1/p2_plus_bubble_to_dg1_div_affine_q3.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

void logSectionHeader( const std::string& header, const char* marker = "-" )
{
   size_t      len = header.length();
   std::string separator( len + 2, marker[0] );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << header << "\n" << separator );
}

std::array< Point3D, 3 > setupTriangle()
{
   std::array< Point3D, 3 > coords;

   coords[0][0] = real_c( 0 );
   coords[0][1] = real_c( 0 );
   coords[0][2] = real_c( 0 );

   coords[1][0] = real_c( 1 );
   coords[1][1] = real_c( 0 );
   coords[1][2] = real_c( 0 );

   coords[2][0] = real_c( 0 );
   coords[2][1] = real_c( 1 );
   coords[2][2] = real_c( 0 );

   // affinely shift unit simplex to another position
   Point3D shift{ real_c( 5 ), real_c( -1.2 ), real_c( 0 ) };

   coords[0] += shift;
   coords[1] += shift;
   coords[2] += shift;

   return coords;
}

Matrix< real_t, 7, 3 > computeGradientMatWithForm( int componentIdx )
{
   Matrix< real_t, 7, 3 >                        elMat;
   forms::dg1_to_p2_plus_bubble_divt_0_affine_q3 form0;
   forms::dg1_to_p2_plus_bubble_divt_1_affine_q3 form1;

   std::array< Point3D, 3 > coords = setupTriangle();
   switch ( componentIdx )
   {
   case 0:
      form0.integrateAll( coords, elMat );
      break;
   case 1:
      form1.integrateAll( coords, elMat );
      break;
   default:
      WALBERLA_ABORT( "componentIdx = " << componentIdx << ", but must be 0 or 1!" );
   }

   return elMat;
}

Matrix< real_t, 7, 3 > computeGradientMatWithMaple( int componentIdx )
{
   Matrix< real_t, 7, 3 > elMat;

   switch ( componentIdx )
   {
   case 0:
      elMat( 0, 0 ) = real_c( 20 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );

      elMat( 1, 0 ) = real_c( 0 );
      elMat( 1, 1 ) = real_c( -20 );
      elMat( 1, 2 ) = real_c( 0 );

      elMat( 2, 0 ) = real_c( 0 );
      elMat( 2, 1 ) = real_c( 0 );
      elMat( 2, 2 ) = real_c( 0 );

      elMat( 3, 0 ) = real_c( -20 );
      elMat( 3, 1 ) = real_c( -20 );
      elMat( 3, 2 ) = real_c( -40 );

      elMat( 4, 0 ) = real_c( 20 );
      elMat( 4, 1 ) = real_c( 20 );
      elMat( 4, 2 ) = real_c( 40 );

      elMat( 5, 0 ) = real_c( -20 );
      elMat( 5, 1 ) = real_c( 20 );
      elMat( 5, 2 ) = real_c( 0 );

      elMat( 6, 0 ) = real_c( -27 );
      elMat( 6, 1 ) = real_c( 27 );
      elMat( 6, 2 ) = real_c( 0 );
      break;

   case 1:
      elMat( 0, 0 ) = real_c( 20 );
      elMat( 0, 1 ) = real_c( 0 );
      elMat( 0, 2 ) = real_c( 0 );

      elMat( 1, 0 ) = real_c( 0 );
      elMat( 1, 1 ) = real_c( 0 );
      elMat( 1, 2 ) = real_c( 0 );

      elMat( 2, 0 ) = real_c( 0 );
      elMat( 2, 1 ) = real_c( 0 );
      elMat( 2, 2 ) = real_c( -20 );

      elMat( 3, 0 ) = real_c( -20 );
      elMat( 3, 1 ) = real_c( -40 );
      elMat( 3, 2 ) = real_c( -20 );

      elMat( 4, 0 ) = real_c( -20 );
      elMat( 4, 1 ) = real_c( 0 );
      elMat( 4, 2 ) = real_c( 20 );

      elMat( 5, 0 ) = real_c( 20 );
      elMat( 5, 1 ) = real_c( 40 );
      elMat( 5, 2 ) = real_c( 20 );

      elMat( 6, 0 ) = real_c( -27 );
      elMat( 6, 1 ) = real_c( 0 );
      elMat( 6, 2 ) = real_c( 27 );
      break;

   default:
      WALBERLA_ABORT( "componentIdx = " << componentIdx << ", but must be 0 or 1!" );
   }

   return elMat / real_c( 120 );
}

Matrix< double, 3, 7 > computeDivergenceMatWithForm( int componentIdx )
{
   Matrix< double, 3, 7 >                       elMat;
   forms::p2_plus_bubble_to_dg1_div_0_affine_q3 form0;
   forms::p2_plus_bubble_to_dg1_div_1_affine_q3 form1;

   std::array< Point3D, 3 > coords = setupTriangle();
   switch ( componentIdx )
   {
   case 0:
      form0.integrateAll( coords, elMat );
      break;
   case 1:
      form1.integrateAll( coords, elMat );
      break;
   default:
      WALBERLA_ABORT( "componentIdx = " << componentIdx << ", but must be 0 or 1!" );
   }

   return elMat;
}

template < typename T >
void showMat( const T& mat, double scalFac = 1.0, double threshold = 1.0e-15 )
{
   WALBERLA_LOG_INFO_ON_ROOT( " " << std::scientific << std::setprecision( 2 ) << std::showpos );

   std::stringstream sStr;

   sStr << " " << std::scientific << std::setprecision( 2 ) << std::showpos;
   for ( int i = 0; i < static_cast< int >( mat.rows() ); ++i )
   {
      for ( int j = 0; j < static_cast< int >( mat.cols() ); ++j )
      {
         if ( std::abs( mat( i, j ) ) >= threshold )
         {
            sStr << mat( i, j ) * scalFac << "  ";
         }
         else
         {
            sStr << "    *    "
                 << "  ";
         }
      }
      sStr << std::endl << " ";
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" << sStr.str() );
}

int main( int argc, char** argv )
{
   // environment stuff
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   // Gradient
   for ( int componentIdx = 0; componentIdx < 2; ++componentIdx )
   {
      std::stringstream sStr;
      sStr << "Gradient (componentIdx = " << componentIdx << ")";
      logSectionHeader( sStr.str() );

      Matrix< double, 7, 3 > elMat = computeGradientMatWithForm( componentIdx );

      WALBERLA_LOG_INFO_ON_ROOT( "Element matrix computed with HOG form:" );
      showMat( elMat );

      int scalingFactor = 120;
      WALBERLA_LOG_INFO_ON_ROOT( "Scaled by " << scalingFactor << ":" );
      showMat( elMat, static_cast< double >( scalingFactor ) );

      Matrix< real_t, 7, 3 > cntrl = computeGradientMatWithMaple( componentIdx );
      WALBERLA_LOG_INFO_ON_ROOT( "Element matrix computed with Maple:" );
      WALBERLA_LOG_INFO_ON_ROOT( "Scaled by " << scalingFactor << ":" );
      showMat( cntrl, static_cast< double >( scalingFactor ) );

      Matrix< real_t, 7, 3 > diff = cntrl - elMat;
      WALBERLA_LOG_INFO_ON_ROOT( "Frobenius norm of difference = " << diff.norm() );
      WALBERLA_CHECK_LESS( diff.norm(), real_c( 5e-16 ) );
   }

   // Divergence
   for ( int componentIdx = 0; componentIdx < 2; ++componentIdx )
   {
      std::stringstream sStr;
      sStr << "Gradient (componentIdx = " << componentIdx << ")";
      logSectionHeader( sStr.str() );

      Matrix< double, 3, 7 > elMat = computeDivergenceMatWithForm( componentIdx );

      WALBERLA_LOG_INFO_ON_ROOT( "Element matrix computed with HOG form:" );
      showMat( elMat );

      int scalingFactor = 120;
      WALBERLA_LOG_INFO_ON_ROOT( "Scaled by " << scalingFactor << ":" );
      showMat( elMat, static_cast< double >( scalingFactor ) );

      Matrix< real_t, 3, 7 > cntrl = computeGradientMatWithMaple( componentIdx ).transpose();
      Matrix< real_t, 3, 7 > diff  = cntrl - elMat;
      WALBERLA_LOG_INFO_ON_ROOT( "Frobenius norm of difference to control = " << diff.norm() );
      WALBERLA_CHECK_LESS( diff.norm(), real_c( 5e-16 ) );
   }

   return EXIT_SUCCESS;
}
