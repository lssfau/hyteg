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
#include "hyteg/forms/form_fenics_base/P1ToP2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/forms/form_hyteg_generated/P1FormDiv.hpp"
#include "hyteg/forms/form_hyteg_generated/P1FormDivT.hpp"
#include "hyteg/forms/form_hyteg_generated/P1FormLaplace.hpp"
#include "hyteg/forms/form_hyteg_generated/P1FormMass.hpp"
#include "hyteg/forms/form_hyteg_manual/P1FormMass3D.hpp"
#include "hyteg/forms/form_hyteg_manual/P1ToP2FormDivT.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDiv.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivT.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormMass.hpp"
#include "hyteg/forms/form_hyteg_manual/P2ToP1FormDiv.hpp"
#include "hyteg/geometry/AffineMap2D.hpp"
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

template < class FormFEniCS, class FormHyTeG, typename matType, uint_t dim >
void compareForms( const std::array< Point3D, dim + 1 >& element,
                   real_t                                tol,
                   std::shared_ptr< GeometryMap >        map         = nullptr,
                   real_t                                scaleFactor = real_c( 1 ) )
{
   bool wBlending = true;

   // setup our two forms
   FormFEniCS fenicsForm;
   FormHyTeG  hytegForm;
   if ( map == nullptr )
   {
      map       = std::make_shared< IdentityMap >();
      wBlending = false;
      WALBERLA_LOG_INFO_ON_ROOT( " -> running with IDENTITY_MAP" );
   }
   hytegForm.setGeometryMap( map );

   // assemble element matrices
   matType matFenics, matHyTeG;
   fenicsForm.integrateAll( element, matFenics );
   hytegForm.integrateAll( element, matHyTeG );

   WALBERLA_LOG_INFO_ON_ROOT( " FEniCS: " << matFenics );
   if ( wBlending )
   {
      matFenics *= scaleFactor;
      WALBERLA_LOG_INFO_ON_ROOT( " FEniCS (scaled): " << matFenics );
   }
   WALBERLA_LOG_INFO_ON_ROOT( " HyTeG:  " << matHyTeG );

   // compare results
   matType matDiff;
   real_t  error = normOfDifference( matFenics, matHyTeG, matDiff );
   WALBERLA_LOG_INFO_ON_ROOT( " Difference: " << matDiff << "\n Frobenius norm: " << error );
   WALBERLA_CHECK_LESS_EQUAL( error, tol );
}

template < class FormFEniCS, class FormHyTeG, typename matType, uint_t dim >
void compareScaled( const std::array< Point3D, dim + 1 >& element,
                    real_t                                tol,
                    std::shared_ptr< GeometryMap >        map,
                    real_t                                scaleFactor )
{
   compareForms< FormFEniCS, FormHyTeG, matType, dim >( element, tol, map, scaleFactor );
}

void run2DTestsWithoutBlending()
{
   // define our test triangle
   std::array< Point3D, 3 > triangle{Point3D( {-0.7, -2.0, 0.0} ), Point3D( {1.0, 1.0, 0.0} ), Point3D( {-1.0, 0.5, 0.0} )};
   // std::array< Point3D, 3 > triangle{Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 0.0, 0.0} ), Point3D( {0.0, 1.0, 0.0} )};

   logSectionHeader( "P1 Mass Forms" );
   compareForms< P1FenicsForm< p1_mass_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_mass, Matrix3r, 2 >( triangle,
                                                                                                                    1e-15 );

   logSectionHeader( "P1 Diffusion Forms" );
   compareForms< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_laplace, Matrix3r, 2 >(
       triangle, 1e-15 );

   logSectionHeader( "P1 DivX Forms" );
   compareForms< P1FenicsForm< p1_div_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_div_1, Matrix3r, 2 >( triangle,
                                                                                                                    2e-15 );

   logSectionHeader( "P1 DivY Forms" );
   compareForms< P1FenicsForm< p1_div_cell_integral_1_otherwise, fenics::NoAssemble >, P1Form_div_2, Matrix3r, 2 >( triangle,
                                                                                                                    2e-15 );

   logSectionHeader( "P1 DivX^T Forms" );
   compareForms< P1FenicsForm< p1_divt_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_divT_1, Matrix3r, 2 >( triangle,
                                                                                                                      1.2e-14 );

   logSectionHeader( "P1 DivY^T Forms" );
   compareForms< P1FenicsForm< p1_divt_cell_integral_1_otherwise, fenics::NoAssemble >, P1Form_divT_2, Matrix3r, 2 >( triangle,
                                                                                                                      1.2e-14 );

   logSectionHeader( "P2 Mass Forms" );
   compareForms< P2FenicsForm< p2_mass_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_mass, Matrix6r, 2 >( triangle,
                                                                                                                    5e-14 );

   logSectionHeader( "P2 Laplace Form" );
   compareForms< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_laplace, Matrix6r, 2 >(
       triangle, 5e-14 );

   logSectionHeader( "P2 DivKGrad Form" );
   compareForms< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_divKgrad, Matrix6r, 2 >(
       triangle, 5e-14 );

   logSectionHeader( "P2 DivX Forms" );
   compareForms< P2FenicsForm< p2_div_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_div< 0 >, Matrix6r, 2 >( triangle,
                                                                                                                       1.5e-14 );

   logSectionHeader( "P2 DivY Forms" );
   compareForms< P2FenicsForm< p2_div_cell_integral_1_otherwise, fenics::NoAssemble >, P2Form_div< 1 >, Matrix6r, 2 >( triangle,
                                                                                                                       2e-14 );

   logSectionHeader( "P2 DivX^T Forms" );
   compareForms< P2FenicsForm< p2_divt_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_divt< 0 >, Matrix6r, 2 >(
       triangle, 1.5e-14 );

   logSectionHeader( "P2 DivY^T Forms" );
   compareForms< P2FenicsForm< p2_divt_cell_integral_1_otherwise, fenics::NoAssemble >, P2Form_divt< 1 >, Matrix6r, 2 >( triangle,
                                                                                                                         2e-14 );

   logSectionHeader( "P2ToP1 DivX Forms" );
   compareForms< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, fenics::NoAssemble >,
                 P2ToP1Form_div< 0 >,
                 Matrixr< 3, 6 >,
                 2 >( triangle, 1.2e-14 );

   logSectionHeader( "P2ToP1 DivY Forms" );
   compareForms< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, fenics::NoAssemble >,
                 P2ToP1Form_div< 1 >,
                 Matrixr< 3, 6 >,
                 2 >( triangle, 1.2e-14 );

   logSectionHeader( "P1ToP2 DivX^T Forms" );
   compareForms< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, fenics::NoAssemble >,
                 P1ToP2Form_divt< 0 >,
                 Matrixr< 6, 3 >,
                 2 >( triangle, 1.2e-14 );

   logSectionHeader( "P1ToP2 DivY^T Forms" );
   compareForms< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, fenics::NoAssemble >,
                 P1ToP2Form_divt< 1 >,
                 Matrixr< 6, 3 >,
                 2 >( triangle, 1.2e-14 );
}

void run2DTestsWithAffineMap()
{
   // define our affine map
   Matrix2r mat;
   mat( 0, 0 ) = real_c( std::exp( 1.0 ) );
   // mat(0,0) = real_c(4.0);
   mat( 1, 1 ) = real_c( 0.5 );
   Point2D vec( {-7.0, 3.0} );
   auto    map = std::make_shared< AffineMap2D >( mat, vec );

   // define our test triangle
   std::array< Point3D, 3 > triangle{Point3D( {-0.7, -2.0, 0.0} ), Point3D( {1.0, 1.0, 0.0} ), Point3D( {-1.0, 0.5, 0.0} )};
   // std::array< Point3D, 3 > triangle{Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 0.0, 0.0} ), Point3D( {0.0, 1.0, 0.0} )};

   logSectionHeader( "P1 Mass Forms" );
   compareScaled< P1FenicsForm< p1_mass_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_mass, Matrix3r, 2 >(
       triangle, 1e-15, map, mat( 0, 0 ) * mat( 1, 1 ) );

   logSectionHeader( "P2 Mass Forms" );
   compareScaled< P2FenicsForm< p2_mass_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_mass, Matrix6r, 2 >(
       triangle, 8e-14, map, mat( 0, 0 ) * mat( 1, 1 ) );

   logSectionHeader( "P2ToP1 DivX Forms" );
   compareScaled< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, fenics::NoAssemble >,
                  P2ToP1Form_div< 0 >,
                  Matrixr< 3, 6 >,
                  2 >( triangle, 2.5e-14, map, mat( 1, 1 ) );

   logSectionHeader( "P2ToP1 DivY Forms" );
   compareScaled< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, fenics::NoAssemble >,
                  P2ToP1Form_div< 1 >,
                  Matrixr< 3, 6 >,
                  2 >( triangle, 3e-14, map, mat( 0, 0 ) );

   logSectionHeader( "P1 DivX Forms" );
   compareScaled< P1FenicsForm< p1_div_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_div_1, Matrix3r, 2 >(
       triangle, 2e-15, map, mat( 1, 1 ) );

   logSectionHeader( "P1 DivY Forms" );
   compareScaled< P1FenicsForm< p1_div_cell_integral_1_otherwise, fenics::NoAssemble >, P1Form_div_2, Matrix3r, 2 >(
       triangle, 5e-15, map, mat( 0, 0 ) );

   logSectionHeader( "P2 DivX Forms" );
   compareScaled< P2FenicsForm< p2_div_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_div< 0 >, Matrix6r, 2 >(
       triangle, 2.5e-14, map, mat( 1, 1 ) );

   logSectionHeader( "P2 DivY Forms" );
   compareScaled< P2FenicsForm< p2_div_cell_integral_1_otherwise, fenics::NoAssemble >, P2Form_div< 1 >, Matrix6r, 2 >(
       triangle, 2.5e-14, map, mat( 0, 0 ) );

   logSectionHeader( "P2 DivX^T Forms" );
   compareScaled< P2FenicsForm< p2_divt_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_divt< 0 >, Matrix6r, 2 >(
       triangle, 2.5e-14, map, mat( 1, 1 ) );

   logSectionHeader( "P2 DivY^T Forms" );
   compareScaled< P2FenicsForm< p2_divt_cell_integral_1_otherwise, fenics::NoAssemble >, P2Form_divt< 1 >, Matrix6r, 2 >(
       triangle, 2.5e-14, map, mat( 0, 0 ) );

   logSectionHeader( "P1 DivX^T Forms" );
   compareScaled< P1FenicsForm< p1_divt_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_divT_1, Matrix3r, 2 >(
       triangle, 1.2e-14, map, mat( 1, 1 ) );

   logSectionHeader( "P1 DivY^T Forms" );
   compareScaled< P1FenicsForm< p1_divt_cell_integral_1_otherwise, fenics::NoAssemble >, P1Form_divT_2, Matrix3r, 2 >(
       triangle, 1.2e-14, map, mat( 0, 0 ) );

   logSectionHeader( "P1ToP2 DivX^T Forms" );
   compareScaled< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, fenics::NoAssemble >,
                  P1ToP2Form_divt< 0 >,
                  Matrixr< 6, 3 >,
                  2 >( triangle, 1.2e-14, map, mat( 1, 1 ) );

   logSectionHeader( "P1ToP2 DivY^T Forms" );
   compareScaled< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, fenics::NoAssemble >,
                  P1ToP2Form_divt< 1 >,
                  Matrixr< 6, 3 >,
                  2 >( triangle, 3e-14, map, mat( 0, 0 ) );

   // diffusion forms are not as straightforward; we just distort equally in both directions,
   // then we should get the same element matrices in 2D
   mat( 0, 0 ) = real_c( std::exp( 1.0 ) );
   mat( 1, 1 ) = mat( 0, 0 );
   map         = std::make_shared< AffineMap2D >( mat, vec );

   logSectionHeader( "P1 Diffusion Forms" );
   compareForms< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >, P1Form_laplace, Matrix3r, 2 >(
       triangle, 1e-15, map, real_c( 1 ) );

   logSectionHeader( "P2 Laplace Form" );
   compareForms< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_laplace, Matrix6r, 2 >(
       triangle, 5e-14, map, real_c( 1 ) );
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
   logSectionHeader( "2D TESTS W/O BLENDING", "=" );
   run2DTestsWithoutBlending();

   logSectionHeader( "2D TESTS W/ BLENDING", "=" );
   run2DTestsWithAffineMap();

   // ------------
   //  3D Testing
   // ------------
   logSectionHeader( "3D TESTS W/O BLENDING", "=" );

   // define our test tetrahedron
   std::array< Point3D, 4 > theTet{
       Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 1.0, 0.0} ), Point3D( {-1.0, 0.5, 0.0} ), Point3D( {0.3, 0.21, -1.2} )};
   // std::array<Point3D,4> theTet{ Point3D({0.0, 0.0, 0.0}), Point3D({1.0, 0.0, 0.0}), Point3D({0.0, 1.0, 0.0}), Point3D({0.0, 0.0, 1.0}) };

   logSectionHeader( "P1 Mass Forms (3D)" );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_mass_cell_integral_0_otherwise >, P1Form_mass3D, Matrix4r, 3 >( theTet,
                                                                                                                          1e-8 );
   // 4 >( theTet, 1e-15 ); only works for lower-order quadrature rule in HyTeG form

   logSectionHeader( "P2 Mass Forms (3D)" );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise >, P2Form_mass, Matrix10r, 3 >(
       theTet, 1e-8 ); // why the large difference? is our FEniCS form under-integrating?

   logSectionHeader( "P2 Laplace Forms (3D)" );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise >, P2Form_laplace, Matrix10r, 3 >(
       theTet, 2e-14 );

   return EXIT_SUCCESS;
}
