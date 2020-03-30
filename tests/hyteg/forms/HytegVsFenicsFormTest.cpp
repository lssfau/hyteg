/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/forms/form_hyteg_generated/P1FormLaplace.hpp"
#include "hyteg/forms/form_hyteg_generated/P1FormMass.hpp"
#include "hyteg/forms/form_hyteg_manual/P1FormMass3D.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormMass.hpp"
// #include "hyteg/forms/form_hyteg_manual/P2ToP1FormDiv.hpp"
#include "hyteg/geometry/IdentityMap.hpp"

using namespace hyteg;

void logSectionHeader( const char* header )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separator( len + 2, '-' );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << "\n" << separator );
}

template < uint_t nRows, uint_t nCols >
real_t normOfDifference( const Matrixr< nRows, nCols >& mat1,
                         const Matrixr< nRows, nCols >& mat2,
                         Matrixr< nRows, nCols >&       diffMat )
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

template < class FormFEniCS, class FormHyTeG, typename matType, uint_t dim >
void compareForms( const std::array< Point3D, dim+1 >& element, real_t tol )
{
   // setup our two forms
   FormFEniCS                     fenicsForm;
   FormHyTeG                      hytegForm;
   std::shared_ptr< GeometryMap > identMap( new IdentityMap );
   hytegForm.setGeometryMap( identMap );

   // assemble element matrices
   matType matFenics, matHyTeG;
   fenicsForm.integrateAll( element, matFenics );
   hytegForm.integrateAll( element, matHyTeG );

   WALBERLA_LOG_INFO_ON_ROOT( " FEniCS: " << matFenics );
   WALBERLA_LOG_INFO_ON_ROOT( " HyTeG:  " << matHyTeG );

   // compare results
   matType matDiff;
   real_t  error = normOfDifference( matFenics, matHyTeG, matDiff );
   WALBERLA_LOG_INFO_ON_ROOT( " Difference: " << matDiff << "\n Frobenius norm: " << error );
   WALBERLA_CHECK_LESS_EQUAL( error, tol );
}

int main( int argc, char** argv )
{
   // abort in case of common floating-point exceptions
   feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

   // environment stuff
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   // ------------
   //  2D Testing
   // ------------

   // define our test triangle
   std::array< Point3D, 3 > triangle{Point3D( {-0.7, -2.0, 0.0} ), Point3D( {1.0, 1.0, 0.0} ), Point3D( {-1.0, 0.5, 0.0} )};

   logSectionHeader( "P1 Mass Forms" );
   compareForms< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                 P1Form_mass,
                 Matrix3r,
                 2 >( triangle, 1e-15 );

   logSectionHeader( "P1 Diffusion Forms" );
   compareForms< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >,
                 P1Form_laplace,
                 Matrix3r,
                 2 >( triangle, 1e-15 );

   logSectionHeader( "P2 Mass Forms" );
   compareForms< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise >,
                 P2Form_mass,
                 Matrix6r,
                 2 >( triangle, 5e-14 );

   logSectionHeader( "P2 Laplace Form" );
   compareForms< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise >,
                 P2Form_laplace,
                 Matrix6r,
                 2 >( triangle, 5e-14 );

   logSectionHeader( "P2 DivKGrad Form" );
   compareForms< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise >,
                 P2Form_divKgrad,
                 Matrix6r,
                 2 >( triangle, 5e-14 );

   // logSectionHeader( "P2ToP1 DivX Forms" );
   // compareForms< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise >,
   //               P2ToP1Form_div<0>,
   //               Matrixr<3,6>,
   //               2 >( triangle, 1e-15 );

   // ------------
   //  3D Testing
   // ------------

   // define our test triangle
   std::array< Point3D, 4 > theTet{
       Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 1.0, 0.0} ), Point3D( {-1.0, 0.5, 0.0} ), Point3D( {0.3, 0.21, -1.2} )};
   // std::array<Point3D,4> theTet{ Point3D({0.0, 0.0, 0.0}), Point3D({1.0, 0.0, 0.0}), Point3D({0.0, 1.0, 0.0}), Point3D({0.0, 0.0, 1.0}) };

   logSectionHeader( "P1 Mass Forms (3D)" );
   compareForms< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                 P1Form_mass3D,
                 Matrix4r,
                 3 >( theTet, 1e-8 );
                 // 4 >( theTet, 1e-15 ); only works for lower-order quadrature rule in HyTeG form

   logSectionHeader( "P2 Mass Forms (3D)" );
   compareForms< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise >,
                 P2Form_mass,
                 Matrix10r,
                 3 >( theTet, 1e-8 ); // why the large difference? is our FEniCS form under-integrating?

   return EXIT_SUCCESS;
}
