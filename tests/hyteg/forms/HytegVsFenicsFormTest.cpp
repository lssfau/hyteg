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
using walberla::math::pi;

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P1ToP2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_affine_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_divt_affine_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_divt_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_full_stokescc_affine_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_full_stokesvar_affine_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_affine_qe.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p1_to_p2/p1_to_p2_divt_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p1_to_p2/p1_to_p2_divt_blending_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_diffusion_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_full_stokescc_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_full_stokesvar_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_affine_qe.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_div_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_div_blending_q2.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/IdentityMap.hpp"

using namespace hyteg;

typedef std::function< real_t( const hyteg::Point3D& ) > func_t;

// parameter function for the callback with variable coeffient forms
std::shared_ptr< func_t > varCoeff = std::make_shared< func_t >( []( const hyteg::Point3D& ) { return real_c( 2 ); } );

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
                   std::shared_ptr< GeometryMap >        map               = nullptr,
                   real_t                                scaleFactor       = real_c( 1 ),
                   bool                                  applyMapToElement = false )
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

   // apply blending to element fed to FEniCS (useful for the tests with
   // an affine mapping)
   std::array< Point3D, dim + 1 > elementForFenics = element;
   if ( applyMapToElement )
   {
      for ( uint_t k = 0; k <= dim; ++k )
      {
         map->evalF( element[k], elementForFenics[k] );
      }
   }

   // assemble element matrices
   matType matFenics, matHyTeG;
   fenicsForm.integrateAll( elementForFenics, matFenics );
   hytegForm.integrateAll( element, matHyTeG );

   WALBERLA_LOG_INFO_ON_ROOT( " FEniCS: " << matFenics );
   if ( wBlending && !applyMapToElement )
   {
      matFenics *= scaleFactor;
      WALBERLA_LOG_INFO_ON_ROOT( " FEniCS (scaled): " << matFenics );
   }
   WALBERLA_LOG_INFO_ON_ROOT( " HyTeG:  " << matHyTeG );

   // compare results
   matType matDiff;
   real_t  error = normOfDifference( matFenics, matHyTeG, matDiff );
   WALBERLA_LOG_INFO_ON_ROOT( " Difference: " << matDiff << "\n Frobenius norm: " << error << "\n Tolerance: " << tol );
   WALBERLA_CHECK_LESS_EQUAL( error, tol );
}

// we need another version, since the "pure" 3D forms, i.e. the _2_2 ones, have a different constructor
// which accepts only a single callback
template < class FormFEniCS, class FormHyTeG, typename matType, uint_t dim >
void compareVarForms3D( const std::array< Point3D, dim + 1 >& element,
                        real_t                                tol,
                        std::shared_ptr< func_t >             callback,
                        std::shared_ptr< GeometryMap >        map,
                        real_t                                scaleFactor = real_c( 1 ) )
{
   // setup our two forms
   FormFEniCS fenicsForm;
   FormHyTeG  hytegForm( *callback );

   // if no map specified run with identity map
   if ( map == nullptr )
   {
      map = std::make_shared< IdentityMap >();
   }

   hytegForm.setGeometryMap( map );

   // apply blending to element fed to FEniCS (useful for the tests with
   // an affine mapping)
   std::array< Point3D, dim + 1 > elementForFenics = element;
   for ( uint_t k = 0; k <= dim; ++k )
   {
      map->evalF( element[k], elementForFenics[k] );
      WALBERLA_LOG_INFO_ON_ROOT( " -> running with IDENTITY_MAP" );
   }

   // assemble element matrices
   matType matFenics, matHyTeG;
   fenicsForm.integrateAll( elementForFenics, matFenics );
   hytegForm.integrateAll( element, matHyTeG );

   WALBERLA_LOG_INFO_ON_ROOT( " FEniCS: " << matFenics );
   matFenics *= scaleFactor;
   WALBERLA_LOG_INFO_ON_ROOT( " FEniCS (scaled): " << matFenics );
   WALBERLA_LOG_INFO_ON_ROOT( " HyTeG:  " << matHyTeG );

   // compare results
   matType matDiff;
   real_t  error = normOfDifference( matFenics, matHyTeG, matDiff );
   WALBERLA_LOG_INFO_ON_ROOT( " Difference: " << matDiff << "\n Frobenius norm: " << error << "\n Tolerance: " << tol );
   WALBERLA_CHECK_LESS_EQUAL( error, tol );
}

template < class FormFEniCS, class FormHyTeG, typename matType, uint_t dim >
void compareVarForms( const std::array< Point3D, dim + 1 >& element,
                      real_t                                tol,
                      std::shared_ptr< func_t >             callback,
                      std::shared_ptr< GeometryMap >        map,
                      real_t                                scaleFactor = real_c( 1 ) )
{
   // setup our two forms
   FormFEniCS fenicsForm;
   FormHyTeG  hytegForm( *callback, *callback );

   // if no map specified run with identity map
   if ( map == nullptr )
   {
      map = std::make_shared< IdentityMap >();
      WALBERLA_LOG_INFO_ON_ROOT( " -> running with IDENTITY_MAP" );
   }

   hytegForm.setGeometryMap( map );

   // apply blending to element fed to FEniCS (useful for the tests with
   // an affine mapping)
   std::array< Point3D, dim + 1 > elementForFenics = element;
   for ( uint_t k = 0; k <= dim; ++k )
   {
      map->evalF( element[k], elementForFenics[k] );
   }

   // assemble element matrices
   matType matFenics, matHyTeG;
   fenicsForm.integrateAll( elementForFenics, matFenics );
   hytegForm.integrateAll( element, matHyTeG );

   WALBERLA_LOG_INFO_ON_ROOT( " FEniCS: " << matFenics );
   matFenics *= scaleFactor;
   WALBERLA_LOG_INFO_ON_ROOT( " FEniCS (scaled): " << matFenics );
   WALBERLA_LOG_INFO_ON_ROOT( " HyTeG:  " << matHyTeG );

   // compare results
   matType matDiff;
   real_t  error = normOfDifference( matFenics, matHyTeG, matDiff );
   WALBERLA_LOG_INFO_ON_ROOT( " Difference: " << matDiff << "\n Frobenius norm: " << error << "\n Tolerance: " << tol );
   WALBERLA_CHECK_LESS_EQUAL( error, tol );
}

template < class FormFEniCS, class FormHyTeG, typename matType, uint_t dim >
void compareScaled( const std::array< Point3D, dim + 1 >& element,
                    real_t                                tol,
                    std::shared_ptr< GeometryMap >        map,
                    real_t                                scaleFactor )
{
   compareForms< FormFEniCS, FormHyTeG, matType, dim >( element, tol, map, scaleFactor, false );
}

template < class FormFEniCS, class FormHyTeG, typename matType, uint_t dim >
void compareUsingAffineMap( const std::array< Point3D, dim + 1 >& element, real_t tol, std::shared_ptr< GeometryMap > map )
{
   compareForms< FormFEniCS, FormHyTeG, matType, dim >( element, tol, map, 0.0, true );
}

void run2DTestsWithoutBlending()
{
   // define our test triangle
   std::array< Point3D, 3 > triangle{
       Point3D( { -0.7, -2.0, 0.0 } ), Point3D( { 1.0, 1.0, 0.0 } ), Point3D( { -1.0, 0.5, 0.0 } ) };
   // std::array< Point3D, 3 > triangle{Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 0.0, 0.0} ), Point3D( {0.0, 1.0, 0.0} )};

   logSectionHeader( "P1 DivX Forms" );
   compareForms< P1FenicsForm< p1_div_cell_integral_0_otherwise, fenics::NoAssemble >, forms::p1_div_0_affine_q1, Matrix3r, 2 >(
       triangle, 2e-15 );

   logSectionHeader( "P1 DivY Forms" );
   compareForms< P1FenicsForm< p1_div_cell_integral_1_otherwise, fenics::NoAssemble >, forms::p1_div_1_affine_q1, Matrix3r, 2 >(
       triangle, 2e-15 );

   logSectionHeader( "P1 DivX^T Forms" );
   compareForms< P1FenicsForm< p1_divt_cell_integral_0_otherwise, fenics::NoAssemble >, forms::p1_divt_0_affine_q1, Matrix3r, 2 >(
       triangle, 2e-15 );

   logSectionHeader( "P1 DivY^T Forms" );
   compareForms< P1FenicsForm< p1_divt_cell_integral_1_otherwise, fenics::NoAssemble >, forms::p1_divt_1_affine_q1, Matrix3r, 2 >(
       triangle, 2e-15 );

   logSectionHeader( "P2 Laplace Form" );
   compareForms< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_laplace, Matrix6r, 2 >(
       triangle, 5e-14 );

   logSectionHeader( "P2 DivKGrad Form" );
   compareForms< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >, P2Form_divKgrad, Matrix6r, 2 >(
       triangle, 5e-14 );

   logSectionHeader( "P2ToP1 DivX Forms" );
   compareForms< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, fenics::NoAssemble >,
                 forms::p2_to_p1_div_0_affine_q2,
                 Matrixr< 3, 6 >,
                 2 >( triangle, 1.2e-14 );

   logSectionHeader( "P2ToP1 DivY Forms" );
   compareForms< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, fenics::NoAssemble >,
                 forms::p2_to_p1_div_1_affine_q2,
                 Matrixr< 3, 6 >,
                 2 >( triangle, 1.2e-14 );

   logSectionHeader( "P1ToP2 DivX^T Forms" );
   compareForms< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, fenics::NoAssemble >,
                 forms::p1_to_p2_divt_0_affine_q2,
                 Matrixr< 6, 3 >,
                 2 >( triangle, 1.2e-14 );

   logSectionHeader( "P1ToP2 DivY^T Forms" );
   compareForms< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, fenics::NoAssemble >,
                 forms::p1_to_p2_divt_1_affine_q2,
                 Matrixr< 6, 3 >,
                 2 >( triangle, 1.2e-14 );

   // HyTeG form generator tests

   logSectionHeader( "P1 diffusion, 2D, no blending (HFG)" );
   compareForms< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >,
                 forms::p1_diffusion_affine_q2,
                 Matrix3r,
                 2 >( triangle, 1e-15 );

   logSectionHeader( "P1 mass, 2D, no blending (HFG)" );
   compareForms< P1FenicsForm< p1_mass_cell_integral_0_otherwise, fenics::NoAssemble >, forms::p1_mass_affine_qe, Matrix3r, 2 >(
       triangle, 1e-15 );

   logSectionHeader( "P1 Epsilon, 2D, no blending (HFG)" );
   compareForms< P1FenicsForm< p1_stokes_epsilon_cell_integral_0_otherwise, fenics::NoAssemble >,
                 forms::p1_epsiloncc_0_0_affine_q2,
                 Matrix3r,
                 2 >( triangle, 1e-15 );
   compareForms< P1FenicsForm< p1_stokes_epsilon_cell_integral_1_otherwise, fenics::NoAssemble >,
                 forms::p1_epsiloncc_0_1_affine_q2,
                 Matrix3r,
                 2 >( triangle, 1e-15 );
   compareForms< P1FenicsForm< p1_stokes_epsilon_cell_integral_2_otherwise, fenics::NoAssemble >,
                 forms::p1_epsiloncc_1_0_affine_q2,
                 Matrix3r,
                 2 >( triangle, 1e-15 );
   compareForms< P1FenicsForm< p1_stokes_epsilon_cell_integral_3_otherwise, fenics::NoAssemble >,
                 forms::p1_epsiloncc_1_1_affine_q2,
                 Matrix3r,
                 2 >( triangle, 1e-15 );

   logSectionHeader( "P1 Full Stokes, 2D, no blending (HFG)" );
   compareForms< P1FenicsForm< p1_stokes_full_cell_integral_0_otherwise, fenics::NoAssemble >,
                 forms::p1_full_stokescc_0_0_affine_q1,
                 Matrix3r,
                 2 >( triangle, 1e-15 );
   compareForms< P1FenicsForm< p1_stokes_full_cell_integral_1_otherwise, fenics::NoAssemble >,
                 forms::p1_full_stokescc_0_1_affine_q1,
                 Matrix3r,
                 2 >( triangle, 1e-15 );
   compareForms< P1FenicsForm< p1_stokes_full_cell_integral_2_otherwise, fenics::NoAssemble >,
                 forms::p1_full_stokescc_1_0_affine_q1,
                 Matrix3r,
                 2 >( triangle, 1e-15 );
   compareForms< P1FenicsForm< p1_stokes_full_cell_integral_3_otherwise, fenics::NoAssemble >,
                 forms::p1_full_stokescc_1_1_affine_q1,
                 Matrix3r,
                 2 >( triangle, 1e-15 );

   logSectionHeader( "P1 Full Stokes, 2D, no blending (HFG), variable coeff" );
   compareVarForms< P1FenicsForm< p1_stokes_full_cell_integral_0_otherwise, fenics::NoAssemble >,
                    forms::p1_full_stokesvar_0_0_affine_q1,
                    Matrix3r,
                    2 >( triangle, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms< P1FenicsForm< p1_stokes_full_cell_integral_1_otherwise, fenics::NoAssemble >,
                    forms::p1_full_stokesvar_0_1_affine_q1,
                    Matrix3r,
                    2 >( triangle, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms< P1FenicsForm< p1_stokes_full_cell_integral_2_otherwise, fenics::NoAssemble >,
                    forms::p1_full_stokesvar_1_0_affine_q1,
                    Matrix3r,
                    2 >( triangle, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms< P1FenicsForm< p1_stokes_full_cell_integral_3_otherwise, fenics::NoAssemble >,
                    forms::p1_full_stokesvar_1_1_affine_q1,
                    Matrix3r,
                    2 >( triangle, 5e-15, varCoeff, nullptr, real_c( 2 ) );

   logSectionHeader( "P2 mass, 2D, no blending (HFG)" );
   compareForms< P2FenicsForm< p2_mass_cell_integral_0_otherwise, fenics::NoAssemble >, forms::p2_mass_affine_qe, Matrix6r, 2 >(
       triangle, 1e-13 );

   logSectionHeader( "P2 Epsilon, 2D, no blending (HFG)" );
   compareForms< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, fenics::NoAssemble >,
                 forms::p2_epsiloncc_0_0_affine_q2,
                 Matrix6r,
                 2 >( triangle, 1e-13 );
   compareForms< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, fenics::NoAssemble >,
                 forms::p2_epsiloncc_0_1_affine_q2,
                 Matrix6r,
                 2 >( triangle, 1e-13 );
   compareForms< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, fenics::NoAssemble >,
                 forms::p2_epsiloncc_1_0_affine_q2,
                 Matrix6r,
                 2 >( triangle, 1e-13 );
   compareForms< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, fenics::NoAssemble >,
                 forms::p2_epsiloncc_1_1_affine_q2,
                 Matrix6r,
                 2 >( triangle, 1e-13 );

   logSectionHeader( "P2 Full Stokes, 2D, no blending (HFG)" );
   compareForms< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, fenics::NoAssemble >,
                 forms::p2_full_stokescc_0_0_affine_q3,
                 Matrix6r,
                 2 >( triangle, 1e-13 );
   compareForms< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, fenics::NoAssemble >,
                 forms::p2_full_stokescc_0_1_affine_q3,
                 Matrix6r,
                 2 >( triangle, 1e-13 );
   compareForms< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, fenics::NoAssemble >,
                 forms::p2_full_stokescc_1_0_affine_q3,
                 Matrix6r,
                 2 >( triangle, 1e-13 );
   compareForms< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, fenics::NoAssemble >,
                 forms::p2_full_stokescc_1_1_affine_q3,
                 Matrix6r,
                 2 >( triangle, 1e-13 );
}

void run2DTestsWithAffineMap()
{
   // define our affine map (rotation + scaling + shift)
   Matrix2r mat;
   real_t   phi = 2.0 / 9.0 * pi;
   mat( 0, 0 )  = +std::cos( phi );
   mat( 0, 1 )  = -std::sin( phi );
   mat( 1, 0 )  = +std::sin( phi ) * 2.25;
   mat( 1, 1 )  = +std::cos( phi ) * 2.25;
   Point2D vec( { -7.0, 3.0 } );
   auto    map = std::make_shared< AffineMap2D >( mat, vec );

   // define our test triangle
   std::array< Point3D, 3 > triangle{
       Point3D( { -0.7, -2.0, 0.0 } ), Point3D( { 1.0, 1.0, 0.0 } ), Point3D( { -1.0, 0.5, 0.0 } ) };
   // std::array< Point3D, 3 > triangle{Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 0.0, 0.0} ), Point3D( {0.0, 1.0, 0.0} )};

   logSectionHeader( "P2ToP1 DivX Forms" );
   compareUsingAffineMap< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p2_to_p1_div_0_blending_q2,
                          Matrixr< 3, 6 >,
                          2 >( triangle, 3.5e-14, map );

   logSectionHeader( "P2ToP1 DivY Forms" );
   compareUsingAffineMap< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, fenics::NoAssemble >,
                          forms::p2_to_p1_div_1_blending_q2,
                          Matrixr< 3, 6 >,
                          2 >( triangle, 3.5e-14, map );

   logSectionHeader( "P1 DivX Forms" );
   compareUsingAffineMap< P1FenicsForm< p1_div_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p1_div_0_blending_q1,
                          Matrix3r,
                          2 >( triangle, 5e-15, map );

   logSectionHeader( "P1 DivY Forms" );
   compareUsingAffineMap< P1FenicsForm< p1_div_cell_integral_1_otherwise, fenics::NoAssemble >,
                          forms::p1_div_1_blending_q1,
                          Matrix3r,
                          2 >( triangle, 5e-15, map );

   logSectionHeader( "P1 DivX^T Forms" );
   compareUsingAffineMap< P1FenicsForm< p1_divt_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p1_divt_0_blending_q1,
                          Matrix3r,
                          2 >( triangle, 5e-14, map );

   logSectionHeader( "P1 DivY^T Forms" );
   compareUsingAffineMap< P1FenicsForm< p1_divt_cell_integral_1_otherwise, fenics::NoAssemble >,
                          forms::p1_divt_1_blending_q1,
                          Matrix3r,
                          2 >( triangle, 5e-14, map );

   logSectionHeader( "P1ToP2 DivX^T Forms" );
   compareUsingAffineMap< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p1_to_p2_divt_0_blending_q2,
                          Matrixr< 6, 3 >,
                          2 >( triangle, 5e-14, map );

   logSectionHeader( "P1ToP2 DivY^T Forms" );
   compareUsingAffineMap< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, fenics::NoAssemble >,
                          forms::p1_to_p2_divt_1_blending_q2,
                          Matrixr< 6, 3 >,
                          2 >( triangle, 5e-14, map );

   logSectionHeader( "P2 Laplace Form" );
   compareUsingAffineMap< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >,
                          P2Form_laplace,
                          Matrix6r,
                          2 >( triangle, 1e-13, map );

   // HyTeG form generator test

   logSectionHeader( "P1 mass, 2D, with blending (HFG)" );
   compareUsingAffineMap< P1FenicsForm< p1_mass_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p1_mass_blending_q4,
                          Matrix3r,
                          2 >( triangle, 2.4e-15, map );

   logSectionHeader( "P2 mass, 2D, with blending (HFG)" );
   compareUsingAffineMap< P2FenicsForm< p2_mass_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p2_mass_blending_q4,
                          Matrix6r,
                          2 >( triangle, 8e-14, map );

   logSectionHeader( "P1 diffusion, 2D, with blending (HFG)" );
   compareUsingAffineMap< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p1_diffusion_blending_q3,
                          Matrix3r,
                          2 >( triangle, 3e-15, map );

   logSectionHeader( "P2 diffusion, 2D, with blending (HFG)" );
   compareUsingAffineMap< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p2_diffusion_blending_q3,
                          Matrix6r,
                          2 >( triangle, 1e-13, map );

   logSectionHeader( "P1 epsilon, 2D, with blending (HFG)" );
   compareUsingAffineMap< P1FenicsForm< p1_stokes_epsilon_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p1_epsiloncc_0_0_blending_q2,
                          Matrix3r,
                          2 >( triangle, 5e-15, map );
   compareUsingAffineMap< P1FenicsForm< p1_stokes_epsilon_cell_integral_1_otherwise, fenics::NoAssemble >,
                          forms::p1_epsiloncc_0_1_blending_q2,
                          Matrix3r,
                          2 >( triangle, 5e-15, map );
   compareUsingAffineMap< P1FenicsForm< p1_stokes_epsilon_cell_integral_2_otherwise, fenics::NoAssemble >,
                          forms::p1_epsiloncc_1_0_blending_q2,
                          Matrix3r,
                          2 >( triangle, 5e-15, map );
   compareUsingAffineMap< P1FenicsForm< p1_stokes_epsilon_cell_integral_3_otherwise, fenics::NoAssemble >,
                          forms::p1_epsiloncc_1_1_blending_q2,
                          Matrix3r,
                          2 >( triangle, 5e-15, map );

   logSectionHeader( "P2 epsilon, 2D, with blending (HFG), quad degree == 2" );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_0_0_blending_q2,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_0_1_blending_q2,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_1_0_blending_q2,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_1_1_blending_q2,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );

   logSectionHeader( "P2 epsilon, 2D, with blending (HFG), quad degree == 3" );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_0_0_blending_q3,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_0_1_blending_q3,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_1_0_blending_q3,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, fenics::NoAssemble >,
                          forms::p2_epsiloncc_1_1_blending_q3,
                          Matrix6r,
                          2 >( triangle, 5e-13, map );

   logSectionHeader( "P2 Full Stokes, 2D, with blending & var (HFG)" );
   compareVarForms< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, fenics::NoAssemble >,
                    forms::p2_full_stokesvar_0_0_blending_q3,
                    Matrix6r,
                    2 >( triangle, 5e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, fenics::NoAssemble >,
                    forms::p2_full_stokesvar_0_1_blending_q3,
                    Matrix6r,
                    2 >( triangle, 5e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, fenics::NoAssemble >,
                    forms::p2_full_stokesvar_1_0_blending_q3,
                    Matrix6r,
                    2 >( triangle, 5e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, fenics::NoAssemble >,
                    forms::p2_full_stokesvar_1_1_blending_q3,
                    Matrix6r,
                    2 >( triangle, 5e-13, varCoeff, map, real_c( 2 ) );
}

void run3DTestsWithoutBlending()
{
   // define our test tetrahedron
   // std::array< Point3D, 4 > theTet{
   //    Point3D( {0.0, 0.0, 0.0} ), Point3D( {1.0, 1.0, 0.0} ), Point3D( {-1.0, 0.5, 0.0} ), Point3D( {0.3, 0.21, -1.2} )};

   // std::array<Point3D,4> theTet{ Point3D({0.0, 0.0, 0.0}), Point3D({1.0, 0.0, 0.0}), Point3D({0.0, 1.0, 0.0}), Point3D({0.0, 0.0, 1.0}) };

   std::array< Point3D, 4 > theTet{ Point3D( { 1.80901699437495e-01, 1.31432778029783e-01, 8.61803398874989e-01 } ),
                                    Point3D( { 1.80901699437495e-01, -1.31432778029783e-01, 8.61803398874989e-01 } ),
                                    Point3D( { 1.80901699437495e-01, 1.31432778029783e-01, 1.11180339887499e+00 } ),
                                    Point3D( { 0.00000000000000e+00, 0.00000000000000e+00, 1.25000000000000e+00 } ) };

   logSectionHeader( "P2ToP1 DivX Forms (3D)" );
   compareForms< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise >,
                 forms::p2_to_p1_div_0_affine_q2,
                 Matrixr< 4, 10 >,
                 3 >( theTet, 1e-14 );

   logSectionHeader( "P2ToP1 DivY Forms (3D)" );
   compareForms< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise >,
                 forms::p2_to_p1_div_1_affine_q2,
                 Matrixr< 4, 10 >,
                 3 >( theTet, 1e-14 );

   logSectionHeader( "P2ToP1 DivZ Forms (3D)" );
   compareForms< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise >,
                 forms::p2_to_p1_div_2_affine_q2,
                 Matrixr< 4, 10 >,
                 3 >( theTet, 1e-14 );

   logSectionHeader( "P1ToP2 DivX^T Forms (3D)" );
   compareForms< P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise >,
                 forms::p1_to_p2_divt_0_affine_q2,
                 Matrixr< 10, 4 >,
                 3 >( theTet, 1e-14 );

   logSectionHeader( "P1ToP2 DivY^T Forms (3D)" );
   compareForms< P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise >,
                 forms::p1_to_p2_divt_1_affine_q2,
                 Matrixr< 10, 4 >,
                 3 >( theTet, 1e-14 );

   logSectionHeader( "P1ToP2 DivZ^T Forms (3D)" );
   compareForms< P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise >,
                 forms::p1_to_p2_divt_2_affine_q2,
                 Matrixr< 10, 4 >,
                 3 >( theTet, 1e-14 );

   logSectionHeader( "P2 Laplace Forms (3D)" );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise >, P2Form_laplace, Matrix10r, 3 >(
       theTet, 3e-14 );

   // HyTeG form generator test

   logSectionHeader( "P1 mass, 3D, no blending (HFG)" );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_mass_cell_integral_0_otherwise >,
                 forms::p1_mass_affine_qe,
                 Matrix4r,
                 3 >( theTet, 1e-15 );

   logSectionHeader( "P2 mass, 3D, no blending (HFG)" );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise >,
                 forms::p2_mass_affine_qe,
                 Matrix10r,
                 3 >( theTet, 1e-15 );

   logSectionHeader( "P1 diffusion, 3D, no blending (HFG)" );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_diffusion_cell_integral_0_otherwise >,
                 forms::p1_diffusion_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );

   logSectionHeader( "P2 diffusion, 3D, no blending (HFG)" );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise >,
                 forms::p2_diffusion_affine_q2,
                 Matrix10r,
                 3 >( theTet, 1e-13 );

   logSectionHeader( "P1 Epsilon, 3D, no blending (HFG)" );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_0_otherwise >,
                 forms::p1_epsiloncc_0_0_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_1_otherwise >,
                 forms::p1_epsiloncc_0_1_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_2_otherwise >,
                 forms::p1_epsiloncc_0_2_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_3_otherwise >,
                 forms::p1_epsiloncc_1_0_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_4_otherwise >,
                 forms::p1_epsiloncc_1_1_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_5_otherwise >,
                 forms::p1_epsiloncc_1_2_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_6_otherwise >,
                 forms::p1_epsiloncc_2_0_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_7_otherwise >,
                 forms::p1_epsiloncc_2_1_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_8_otherwise >,
                 forms::p1_epsiloncc_2_2_affine_q2,
                 Matrix4r,
                 3 >( theTet, 3e-14 );

   logSectionHeader( "P1 Full Stokes, 3D, no blending (HFG)" );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_0_otherwise >,
                 forms::p1_full_stokescc_0_0_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_1_otherwise >,
                 forms::p1_full_stokescc_0_1_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_2_otherwise >,
                 forms::p1_full_stokescc_0_2_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_3_otherwise >,
                 forms::p1_full_stokescc_1_0_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_4_otherwise >,
                 forms::p1_full_stokescc_1_1_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_5_otherwise >,
                 forms::p1_full_stokescc_1_2_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_6_otherwise >,
                 forms::p1_full_stokescc_2_0_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_7_otherwise >,
                 forms::p1_full_stokescc_2_1_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );
   compareForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_8_otherwise >,
                 forms::p1_full_stokescc_2_2_affine_q1,
                 Matrix4r,
                 3 >( theTet, 3e-14 );

   logSectionHeader( "P1 Full Stokes, 3D, no blending (HFG), var coeff" );
   compareVarForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_0_otherwise >,
                    forms::p1_full_stokesvar_0_0_affine_q1,
                    Matrix4r,
                    3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_1_otherwise >,
                    forms::p1_full_stokesvar_0_1_affine_q1,
                    Matrix4r,
                    3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms3D< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_2_otherwise >,
                      forms::p1_full_stokesvar_0_2_affine_q1,
                      Matrix4r,
                      3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_3_otherwise >,
                    forms::p1_full_stokesvar_1_0_affine_q1,
                    Matrix4r,
                    3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_4_otherwise >,
                    forms::p1_full_stokesvar_1_1_affine_q1,
                    Matrix4r,
                    3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms3D< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_5_otherwise >,
                      forms::p1_full_stokesvar_1_2_affine_q1,
                      Matrix4r,
                      3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms3D< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_6_otherwise >,
                      forms::p1_full_stokesvar_2_0_affine_q1,
                      Matrix4r,
                      3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms3D< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_7_otherwise >,
                      forms::p1_full_stokesvar_2_1_affine_q1,
                      Matrix4r,
                      3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );
   compareVarForms3D< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_full_tet_cell_integral_8_otherwise >,
                      forms::p1_full_stokesvar_2_2_affine_q1,
                      Matrix4r,
                      3 >( theTet, 5e-15, varCoeff, nullptr, real_c( 2 ) );

   logSectionHeader( "P2 Epsilon, 3D, no blending (HFG)" );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise >,
                 forms::p2_epsiloncc_0_0_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise >,
                 forms::p2_epsiloncc_0_1_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise >,
                 forms::p2_epsiloncc_0_2_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise >,
                 forms::p2_epsiloncc_1_0_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise >,
                 forms::p2_epsiloncc_1_1_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise >,
                 forms::p2_epsiloncc_1_2_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise >,
                 forms::p2_epsiloncc_2_0_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise >,
                 forms::p2_epsiloncc_2_1_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise >,
                 forms::p2_epsiloncc_2_2_affine_q2,
                 Matrix10r,
                 3 >( theTet, 3e-14 );

   logSectionHeader( "P2 Full Stokes, 3D, no blending (HFG)" );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_0_otherwise >,
                 forms::p2_full_stokescc_0_0_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_1_otherwise >,
                 forms::p2_full_stokescc_0_1_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_2_otherwise >,
                 forms::p2_full_stokescc_0_2_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_3_otherwise >,
                 forms::p2_full_stokescc_1_0_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_4_otherwise >,
                 forms::p2_full_stokescc_1_1_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_5_otherwise >,
                 forms::p2_full_stokescc_1_2_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_6_otherwise >,
                 forms::p2_full_stokescc_2_0_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_7_otherwise >,
                 forms::p2_full_stokescc_2_1_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
   compareForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_8_otherwise >,
                 forms::p2_full_stokescc_2_2_affine_q3,
                 Matrix10r,
                 3 >( theTet, 1e-13 );
}

void run3DTestsWithAffineMap()
{
   // define our affine map
   Matrix3r mat;

   // simple scaling
   mat( 0, 0 ) = real_c( 2.0 );
   mat( 1, 1 ) = real_c( 1.0 );
   mat( 2, 2 ) = real_c( std::exp( 1.0 ) );

   // more interesting variant (rotation around x-axis, y-axis and scaling of z-component)
#define CHALLENGING
#ifdef CHALLENGING
   mat( 0, 0 ) = +8.660254037844387e-01;
   mat( 0, 1 ) = -1.545084971874737e-01;
   mat( 0, 2 ) = -9.510565162951534e-01;
   mat( 1, 0 ) = +0.000000000000000e+00;
   mat( 1, 1 ) = +9.510565162951535e-01;
   mat( 1, 2 ) = -6.180339887498948e-01;
   mat( 2, 0 ) = +4.999999999999999e-01;
   mat( 2, 1 ) = +2.676165673298174e-01;
   mat( 2, 2 ) = +1.647278207092664e+00;
#endif
#undef CHALLENGING

   Point3D vec( { -7.0, 3.0, 2.0 } );
   auto    map = std::make_shared< AffineMap3D >( mat, vec );

   // define our test tetrahedrons
   Point3D                  v1( { 0.0, 0.00, 0.0 } );
   Point3D                  v2( { 1.0, 1.00, 0.0 } );
   Point3D                  v3( { -1.0, 0.50, 0.0 } );
   Point3D                  v4( { 0.3, 0.21, -1.2 } );
   std::array< Point3D, 4 > theTet{ v1, v2, v3, v4 };

   logSectionHeader( "P2 Laplace Forms (3D)" );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise >,
                          P2Form_laplace,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );

   logSectionHeader( "P2ToP1 DivX Forms (3D)" );
   compareUsingAffineMap< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise >,
                          forms::p2_to_p1_div_0_blending_q2,
                          Matrixr< 4, 10 >,
                          3 >( theTet, 1e-14, map );

   logSectionHeader( "P2ToP1 DivY Forms (3D)" );
   compareUsingAffineMap< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise >,
                          forms::p2_to_p1_div_1_blending_q2,
                          Matrixr< 4, 10 >,
                          3 >( theTet, 1e-14, map );

   logSectionHeader( "P2ToP1 DivZ Forms (3D)" );
   compareUsingAffineMap< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise >,
                          forms::p2_to_p1_div_2_blending_q2,
                          Matrixr< 4, 10 >,
                          3 >( theTet, 1e-14, map );

   // HyTeG form generator test

   logSectionHeader( "P1 mass, 3D, with blending (HFG)" );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_mass_cell_integral_0_otherwise >,
                          forms::p1_mass_blending_q4,
                          Matrix4r,
                          3 >( theTet, 5e-15, map );

   logSectionHeader( "P2 mass, 3D, with blending (HFG)" );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise >,
                          forms::p2_mass_blending_q4,
                          Matrix10r,
                          3 >( theTet, 1e-15, map );

   logSectionHeader( "P1 diffusion, 3D, with blending (HFG)" );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_diffusion_cell_integral_0_otherwise >,
                          forms::p1_diffusion_blending_q3,
                          Matrix4r,
                          3 >( theTet, 1e-15, map );

   logSectionHeader( "P2 diffusion, 3D, with blending (HFG)" );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise >,
                          forms::p2_diffusion_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );

   logSectionHeader( "P1 Epsilon, 3D, with blending (HFG)" );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_0_otherwise >,
                          forms::p1_epsiloncc_0_0_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_1_otherwise >,
                          forms::p1_epsiloncc_0_1_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_2_otherwise >,
                          forms::p1_epsiloncc_0_2_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_3_otherwise >,
                          forms::p1_epsiloncc_1_0_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_4_otherwise >,
                          forms::p1_epsiloncc_1_1_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_5_otherwise >,
                          forms::p1_epsiloncc_1_2_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_6_otherwise >,
                          forms::p1_epsiloncc_2_0_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_7_otherwise >,
                          forms::p1_epsiloncc_2_1_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );
   compareUsingAffineMap< P1FenicsForm< fenics::NoAssemble, p1_tet_stokes_epsilon_tet_cell_integral_8_otherwise >,
                          forms::p1_epsiloncc_2_2_blending_q2,
                          Matrix4r,
                          3 >( theTet, 3e-14, map );

   logSectionHeader( "P2 Epsilon, 3D, with blending (HFG), q2" );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise >,
                          forms::p2_epsiloncc_0_0_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise >,
                          forms::p2_epsiloncc_0_1_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise >,
                          forms::p2_epsiloncc_0_2_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise >,
                          forms::p2_epsiloncc_1_0_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise >,
                          forms::p2_epsiloncc_1_1_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise >,
                          forms::p2_epsiloncc_1_2_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise >,
                          forms::p2_epsiloncc_2_0_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise >,
                          forms::p2_epsiloncc_2_1_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise >,
                          forms::p2_epsiloncc_2_2_blending_q2,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );

   logSectionHeader( "P2 Epsilon, 3D, with blending (HFG), q3" );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise >,
                          forms::p2_epsiloncc_0_0_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise >,
                          forms::p2_epsiloncc_0_1_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise >,
                          forms::p2_epsiloncc_0_2_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise >,
                          forms::p2_epsiloncc_1_0_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-14, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise >,
                          forms::p2_epsiloncc_1_1_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise >,
                          forms::p2_epsiloncc_1_2_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise >,
                          forms::p2_epsiloncc_2_0_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise >,
                          forms::p2_epsiloncc_2_1_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );
   compareUsingAffineMap< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise >,
                          forms::p2_epsiloncc_2_2_blending_q3,
                          Matrix10r,
                          3 >( theTet, 5e-13, map );

   logSectionHeader( "P2 Epsilon, 3D, with blending & var (HFG)" );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise >,
                    forms::p2_epsilonvar_0_0_blending_q2,
                    Matrix10r,
                    3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise >,
                    forms::p2_epsilonvar_0_1_blending_q2,
                    Matrix10r,
                    3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise >,
                    forms::p2_epsilonvar_1_0_blending_q2,
                    Matrix10r,
                    3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise >,
                    forms::p2_epsilonvar_1_1_blending_q2,
                    Matrix10r,
                    3 >( theTet, 2e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise >,
                      forms::p2_epsilonvar_0_2_blending_q2,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise >,
                      forms::p2_epsilonvar_1_2_blending_q2,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise >,
                      forms::p2_epsilonvar_2_0_blending_q2,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise >,
                      forms::p2_epsilonvar_2_1_blending_q2,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise >,
                      forms::p2_epsilonvar_2_2_blending_q2,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );

   logSectionHeader( "P2 Full Stokes, 3D, with blending & var (HFG)" );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_0_otherwise >,
                    forms::p2_full_stokesvar_0_0_blending_q3,
                    Matrix10r,
                    3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_1_otherwise >,
                    forms::p2_full_stokesvar_0_1_blending_q3,
                    Matrix10r,
                    3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_3_otherwise >,
                    forms::p2_full_stokesvar_1_0_blending_q3,
                    Matrix10r,
                    3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_4_otherwise >,
                    forms::p2_full_stokesvar_1_1_blending_q3,
                    Matrix10r,
                    3 >( theTet, 2e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_2_otherwise >,
                      forms::p2_full_stokesvar_0_2_blending_q3,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_5_otherwise >,
                      forms::p2_full_stokesvar_1_2_blending_q3,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_6_otherwise >,
                      forms::p2_full_stokesvar_2_0_blending_q3,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_7_otherwise >,
                      forms::p2_full_stokesvar_2_1_blending_q3,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
   compareVarForms3D< P2FenicsForm< fenics::NoAssemble, p2_tet_stokes_full_tet_cell_integral_8_otherwise >,
                      forms::p2_full_stokesvar_2_2_blending_q3,
                      Matrix10r,
                      3 >( theTet, 1e-13, varCoeff, map, real_c( 2 ) );
}

int main( int argc, char** argv )
{
#ifndef __APPLE__
   // clang 9 seams to produce a problem related to vectorized division
   // https://stackoverflow.com/questions/63125919/how-to-avoid-floating-point-exceptions-in-unused-simd-lanes
#if ( defined( NDEBUG ) && defined( __clang__ ) ) || ( !defined( NDEBUG ) && defined( __INTEL_LLVM_COMPILER ) )
   // clang 10 has problems with some of the forms in Release mode (see issue #147)
   // intel llvm has problems in debug mode
#else
// abort in case of common floating-point exceptions
#ifndef _MSC_VER
   feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#endif
#endif
#endif
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
   run3DTestsWithoutBlending();

   logSectionHeader( "3D TESTS W/ BLENDING", "=" );
   run3DTestsWithAffineMap();

   return EXIT_SUCCESS;
}
