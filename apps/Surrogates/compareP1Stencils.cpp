/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1SurrogateOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator_new.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

typedef std::function< real_t( const hyteg::Point3D& ) > expression;

typedef struct
{
   uint_t degree;
   uint_t evalLevel;
   uint_t lsqLevel;
   uint_t weightIdx;
} param_t;

void showUsage()
{
   std::cout << "\n --------\n  USAGE:\n --------\n\n"
             << " show_mesh demonstrates generation of a MeshInfo object by one of the\n following methods:\n\n"
             << " This is steered by choosing one of the options below:\n\n"
             << "  --help ............ show this usage info\n"
             << "  --degree [n] ...... degree of interpolating polynomial\n"
             << "  --level [n] ....... refinement level for stencil evaluation\n"
             << "  --lsqLevel [n] .... sampling level for fitting surrogate polynomial\n"
             << "  --weight [n] ...... index of stencil weight to examine\n"
             << "                      (see enum stencilDirection for numeric values)\n\n"
             << " For [n] please supply an integer value\n"
             << std::endl;
}

bool getParamValue( int argc, char* argv[], const char* param, uint_t& value )
{
   bool found = false;
   for ( uint_t k = 1; k < argc; k += 2 )
   {
      if ( strcmp( argv[k], param ) == 0 )
      {
         value = (uint_t) atoi( argv[k + 1] );
         found = true;
      }
   }
   return found;
}

void echoParams( param_t& params )
{
   WALBERLA_LOG_INFO_ON_ROOT( " -> degree ...... " << params.degree );
   WALBERLA_LOG_INFO_ON_ROOT( " -> level ....... " << params.evalLevel );
   WALBERLA_LOG_INFO_ON_ROOT( " -> lsqLevel .... " << params.lsqLevel );
   WALBERLA_LOG_INFO_ON_ROOT( " -> weight ...... " << params.weightIdx );
}

void processCLI( int argc, char* argv[], param_t& params )
{
   bool weAreFine = true;

   if ( argc == 2 )
   {
      if ( strcmp( argv[1], "--help" ) == 0 )
      {
         showUsage();
         exit( EXIT_SUCCESS );
      }
      weAreFine = false;
   }

   weAreFine = weAreFine && getParamValue( argc, argv, "--degree", params.degree );
   weAreFine = weAreFine && getParamValue( argc, argv, "--level", params.evalLevel );
   weAreFine = weAreFine && getParamValue( argc, argv, "--lsqLevel", params.lsqLevel );
   weAreFine = weAreFine && getParamValue( argc, argv, "--weight", params.weightIdx );

   if ( !weAreFine )
   {
      showUsage();
      exit( EXIT_FAILURE );
   }
}

void logMessage( const char* msg )
{
   WALBERLA_LOG_INFO_ON_ROOT( " " << msg );
}

void logSectionHeader( const char* header )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separatorBot( len + 2, '-' );
   std::string separatorTop( 50, '=' );
   WALBERLA_LOG_INFO_ON_ROOT( separatorTop << "\n " << hdr << "\n" << separatorBot );
}

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   if ( argc == 1 )
   {
      showUsage();
      exit( EXIT_SUCCESS );
   }

   logSectionHeader( "General Setup" );

   logMessage( "Processing command-line" );
   param_t params;
   processCLI( argc, argv, params );
   echoParams( params );

   logMessage( "Generating single triangle mesh" );
   MeshInfo meshInfo = MeshInfo::singleTriangle( Point2D( {0.0, 0.0} ), Point2D( {0.0, 1.0} ), Point2D( {1.0, 0.0} ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t maxLevel = params.evalLevel;
   uint_t minLevel = maxLevel;

   logMessage( "Selecting parameter function" );
   real_t omg1 = real_c( 0 );
   real_t omg2 = real_c( 0 );

   expression parameterFunction = [omg1, omg2]( const hyteg::Point3D& x ) {
      return real_c( 1.5 ) + std::sin( omg1 * pi * x[0] ) * std::sin( omg2 * pi * x[1] );
   };

   P1Function< real_t > kappa( "demo", storage, minLevel, maxLevel );
   kappa.interpolate( parameterFunction, maxLevel, All );

   logMessage( "Setting up P1Form" );
   forms::p1_div_k_grad_affine_q3 form( parameterFunction, parameterFunction );

   // ===================
   //  Standard FEM Part
   // ===================
   logSectionHeader( "Standard FEM Part" );

   logMessage( "Preparing P1VariableOperator_new" );
   P1AffineDivkGradOperator_new varOp( storage, minLevel, maxLevel, form );

   logMessage( "Preparing Stencil Generation" );
   std::vector< PrimitiveID > faceIDs = storage->getFaceIDs();
   WALBERLA_ASSERT_EQUAL( faceIDs.size(), 1 );

   Face&                                                   face        = *storage->getFace( faceIDs[0] );
   const PrimitiveDataID< StencilMemory< real_t >, Face >& stencilID   = varOp.getFaceStencilID();
   uint_t                                                  stencilSize = face.getData( stencilID )->getSize( maxLevel );

   std::vector< real_t > stencil( stencilSize );

   real_t h = real_c( 1 ) / ( walberla::real_c( levelinfo::num_microvertices_per_edge( maxLevel ) - 1 ) );

   stencil::Directions2D sDir = stencil::Directions2D( h, face );

   WALBERLA_LOG_INFO_ON_ROOT( "-> test level  = " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "-> mesh width  = " << h );
   WALBERLA_LOG_INFO_ON_ROOT( "-> stencilSize = " << stencilSize );

   // select stencil weight
   // auto direction = stencilDirection::VERTEX_C;
   auto direction = static_cast< stencilDirection >( params.weightIdx ); // makes me shiver

   // defines a function to evaluate the variable stencils and return one of the weights
   expression weight = [&stencil, &form, sDir, direction]( const hyteg::Point3D& x ) {
      for ( size_t k = 0; k < stencil.size(); k++ )
      {
         stencil[k] = real_c( 0 );
      }

      // add contributions from six elements surrounding node to stencil weights
      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form, {x, x + sDir.W, x + sDir.S}, P1Elements::P1Elements2D::elementSW, stencil.data() );

      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form, {x, x + sDir.S, x + sDir.SE}, P1Elements::P1Elements2D::elementS, stencil.data() );

      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form, {x, x + sDir.SE, x + sDir.E}, P1Elements::P1Elements2D::elementSE, stencil.data() );

      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form, {x, x + sDir.E, x + sDir.N}, P1Elements::P1Elements2D::elementNE, stencil.data() );

      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form, {x, x + sDir.N, x + sDir.NW}, P1Elements::P1Elements2D::elementN, stencil.data() );

      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form, {x, x + sDir.NW, x + sDir.W}, P1Elements::P1Elements2D::elementNW, stencil.data() );

      return stencil[vertexdof::stencilIndexFromVertex( direction )];
   };

   P1Function< real_t > stencilWeight( "stencil weights", storage, minLevel, maxLevel );
   stencilWeight.interpolate( weight, maxLevel, All );

   // =================
   //  Surrogates Part
   // =================
   logSectionHeader( "Surrogates Part" );

   logMessage( "Preparing P1SurrogateOperator" );
   P1SurrogateOperator sOp( storage, minLevel, maxLevel, form );

   logMessage( "Interpolating surrogate polynomials" );
   sOp.interpolateStencils( params.degree, params.lsqLevel );

   logMessage( "Extracting polynomial" );
   typedef Polynomial< 2, Point2D, MonomialBasis2D > Poly2D;

   auto                   polyID      = sOp.getFacePolyID();
   std::vector< Poly2D >& stencilPoly = face.getData( polyID )->getData( maxLevel );
   Poly2D&                poly        = stencilPoly[params.weightIdx];

   for ( uint_t k = 0; k <= params.degree; k++ )
   {
     WALBERLA_LOG_INFO_ON_ROOT( " -> coeff[" << k << "] = " << poly.getCoefficient( k ) );
   }

   logMessage( "Evaluating surrogate polynomial" );
   expression surrogate = [&poly]( const hyteg::Point3D& x ) {
      Point2D p( {x[0], x[1]} );
      return poly.eval( p );
   };

   P1Function< real_t > surrogateWeight( "surrogate weights", storage, minLevel, maxLevel );
   surrogateWeight.interpolate( surrogate, maxLevel, All );

   // =================
   //  Postprocessing
   // =================
   logSectionHeader( "Postprocessing" );

   logMessage( "Exporting weight functions" );
   std::string oFile   = "surrogates";
   std::string oFolder = "../../output";
   WALBERLA_LOG_INFO_ON_ROOT( " -> directory = " << oFolder );
   WALBERLA_LOG_INFO_ON_ROOT( " -> baseName  = " << oFile );
   VTKOutput vtkOutput( oFolder, oFile, storage );
   vtkOutput.add( kappa );
   vtkOutput.add( stencilWeight );
   vtkOutput.add( surrogateWeight );
   vtkOutput.write( maxLevel );

   return 0;
}
