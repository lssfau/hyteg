/*
 * Copyright (c) 2017-2021 Marcus Mohr, Nils Kohl.
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

// compare element matrices from FEniCS forms and HyTeG forms
#include <cfenv>
#include <core/Environment.h>
#include <core/math/Constants.h>

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_affine_q2.hpp"
#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/IdentityMap.hpp"

using namespace hyteg;

void logSectionHeader( const char* header, const char* marker = "-" )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separator( len + 2, marker[0] );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << "\n" << separator );
}

template < uint_t nRows, uint_t nCols >
real_t
    normOfDifference( const Matrixr< nRows, nCols >& mat1, const Matrixr< nRows, nCols >& mat2, Matrixr< nRows, nCols >& diffMat )
{
   real_t        norm  = 0.0;
   const real_t* data1 = mat1.data();
   const real_t* data2 = mat2.data();
   real_t*       diff  = diffMat.data();
   for ( uint_t k = 0; k < nRows * nCols; k++ )
   {
      diff[k] = ( data1[k] - data2[k] );
      norm += diff[k] * diff[k];
   }
   return std::sqrt( norm );
}

template < typename Form, uint_t dim, uint_t rows, uint_t cols >
void compareRows( const Form& form, const std::array< Point3D, dim + 1 >& element, uint_t row, real_t tol )
{
   typedef Matrixr< rows, cols > MatType;
   typedef Matrixr< 1, cols >    RowType;

   MatType elMat;
   RowType elRow;
   RowType elRowOfFullMat;
   RowType elRowDifference;

   form.integrateAll( element, elMat );
   form.integrateRow( row, element, elRow );

   for ( uint_t col = 0; col < cols; col++ )
   {
      elRowOfFullMat( row, col ) = elRow( row, col );
   }

   real_t error = normOfDifference( elRow, elRowOfFullMat, elRowDifference );
   WALBERLA_LOG_INFO_ON_ROOT( " Difference: " << elRowDifference << "\n Frobenius norm: " << error << "\n Tolerance: " << tol );
   WALBERLA_CHECK_LESS_EQUAL( error, tol );
}

int main( int argc, char** argv )
{
#ifndef __APPLE__
   // abort in case of common floating-point exceptions
   feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#endif
   // environment stuff
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   std::array< Point3D, 3 > element2D = {
       Point3D( { 0.1, 0.345, 0 } ), Point3D( { 0.2, 0.083745, 0 } ), Point3D( { 0.985, 0.3, 0 } ) };

   forms::p1_diffusion_affine_q2 form_p1_diffusion_affine_q2;

   compareRows< forms::p1_diffusion_affine_q2, 2, 3, 3 >( form_p1_diffusion_affine_q2, element2D, 0, 1e-16 );

   return EXIT_SUCCESS;
}
