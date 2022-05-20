/*
 * Copyright (c) 2021 Marcus Mohr.
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
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1SurrogateOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
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
             << " compareP1Stencils: For the div-k-grad problem we evaluate a specific stencil weight\n"
             << " for a standard FE approach and for a corresponding surrogate polynomial. The results\n"
             << " are obtained for a single macro-face, compared and stored for visualisation.\n\n"
             << " Required arguments are:\n\n"
             << "  --degree [n] ...... degree of interpolating polynomial\n"
             << "  --level [n] ....... refinement level for stencil evaluation\n"
             << "  --lsqLevel [n] .... sampling level for fitting surrogate polynomial\n"
             << "  --weight [n] ...... index of stencil weight to examine\n"
             << "                      (see enum stencilDirection for numeric values)\n\n"
             << " For [n] please supply an integer value\n\n"
             << " or:\n\n"
             << "  --help ............ show this usage info\n"
             << std::endl;
}

bool getParamValue( int argc, char* argv[], const char* param, uint_t& value )
{
   bool found = false;
   for ( int k = 1; k < argc; k += 2 )
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

uint_t getNumInnerMicroVertices( uint_t level )
{
   return levelinfo::num_microvertices_per_face( level ) - 3 * levelinfo::num_microvertices_per_edge( level ) + 3;
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

Matrix< real_t, 2, 2 > computeMappingFunction( std::vector< Point3D >& vertices )
{
   Matrix< real_t, 2, 2 > M;

   real_t x0 = vertices[0][0];
   real_t x1 = vertices[1][0];
   real_t x2 = vertices[2][0];

   real_t y0 = vertices[0][1];
   real_t y1 = vertices[1][1];
   real_t y2 = vertices[2][1];

   real_t det = ( x1 - x0 ) * ( y2 - y0 ) - ( x2 - x0 ) * ( y1 - y0 );

   M( 0, 0 ) = ( y2 - y0 ) / det;
   M( 0, 1 ) = ( x0 - x2 ) / det;
   M( 1, 0 ) = ( y0 - y1 ) / det;
   M( 1, 1 ) = ( x1 - x0 ) / det;

   return M;
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
   // MeshInfo meshInfo = MeshInfo::singleTriangle( Point2D( {0.0, 0.0} ), Point2D( {1.0, 0.0} ), Point2D( {0.0, 1.0} ) );
   MeshInfo meshInfo = MeshInfo::singleTriangle( Point2D( {0.0, 0.0} ), Point2D( {2.0, 0.0} ), Point2D( {2.0, 1.0} ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t maxLevel = params.evalLevel;
   uint_t minLevel = maxLevel;

   logMessage( "Selecting parameter function" );
   real_t omg1 = real_c( 0.25 );
   real_t omg2 = real_c( 0.50 );

   expression parameterFunction1 = [omg1, omg2]( const hyteg::Point3D& x ) {
      real_t damping = real_c( 0.1 );
      return real_c( 0.2 ) + damping * std::cos( omg1 * pi * x[0] ) * std::sin( omg2 * pi * ( x[1] - 0.5 ) );
   };

   expression parameterFunction = []( const hyteg::Point3D& x ) {
      real_t value = real_c( 1.2 ) + std::sin( pi * ( x[1] * 0.5 + 0.1 ) ) + std::cos( 0.5 * pi * x[0] );
      return 0.1 * value;
   };

   expression parameterFunction2 = []( const hyteg::Point3D& x ) { return 1.0 * x[0] + 0.0 * x[1] + 1.0; };

   P1Function< real_t > kappa( "kappa", storage, minLevel, maxLevel );
   kappa.interpolate( parameterFunction, maxLevel, All );
   real_t kappaMin = kappa.getMinValue( maxLevel, All );
   WALBERLA_LOG_INFO_ON_ROOT( " -> kappa (minVal) = " << kappaMin );
   WALBERLA_LOG_INFO_ON_ROOT( " -> kappa (maxVal) = " << kappa.getMaxValue( maxLevel, All ) );
   if ( kappaMin < real_c( 0 ) )
   {
      WALBERLA_ABORT( "Expecting positive values of kappa!" );
   }

   logMessage( "Setting up P1Form" );
   forms::p1_div_k_grad_affine_q3 form( parameterFunction, parameterFunction );

   logMessage( "Performing sanity check" );
   uint_t nSamples = getNumInnerMicroVertices( params.lsqLevel );
   uint_t nCoeffs  = ( params.degree + 2 ) * ( params.degree + 1 ) / 2;
   WALBERLA_ASSERT_EQUAL( nCoeffs, Polynomial2D< MonomialBasis2D >::getNumCoefficients( params.degree ) );
   WALBERLA_LOG_INFO_ON_ROOT( " -> Sampling on level " << params.lsqLevel << " gives " << nSamples << " samples" );
   WALBERLA_LOG_INFO_ON_ROOT( " -> Polynomial of degree " << params.degree << " has " << nCoeffs << " coefficients" );
   if ( nSamples < nCoeffs )
   {
      WALBERLA_ABORT( " -> Cowardly refusing to work with these parameters!" );
   }

   logMessage( "Extracting face and vertex info" );
   std::vector< PrimitiveID > faceIDs = storage->getFaceIDs();
   WALBERLA_ASSERT_EQUAL( faceIDs.size(), 1 );
   Face& face = *storage->getFace( faceIDs[0] );

   std::vector< Point3D > vertices;
   for ( auto vertexID : face.neighborVertices() )
   {
      Vertex& vert = *storage->getVertex( vertexID );
      vertices.push_back( vert.getCoordinates() );
   }

   logMessage( "Preparing affine mapping function" );
   Matrix< real_t, 2, 2 > Phi = computeMappingFunction( vertices );

   // ===================
   //  Standard FEM Part
   // ===================
   logSectionHeader( "Standard FEM Part" );

   logMessage( "Preparing P1VariableOperator" );
   P1AffineDivkGradOperator varOp( storage, minLevel, maxLevel, form );

   logMessage( "Preparing Stencil Generation" );
   const PrimitiveDataID< StencilMemory< real_t >, Face >& stencilID   = varOp.getFaceStencilID();
   uint_t                                                  stencilSize = face.getData( stencilID )->getSize( maxLevel );
   std::vector< real_t >                                   stencil( stencilSize );

   real_t delta = real_c( 1 ) / ( walberla::real_c( levelinfo::num_microvertices_per_edge( maxLevel ) - 1 ) );

   stencil::Directions2D sDir = stencil::Directions2D( delta, face );

   WALBERLA_LOG_INFO_ON_ROOT( "-> test level  = " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "-> mesh width  = " << delta );
   WALBERLA_LOG_INFO_ON_ROOT( "-> stencilSize = " << stencilSize );

   // select stencil weight
   auto direction = static_cast< stencilDirection >( params.weightIdx ); // makes me shiver a little ;-)

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
   auto                                            polyID      = sOp.getFacePolyID();
   std::vector< Polynomial2D< MonomialBasis2D > >& stencilPoly = face.getData( polyID )->getData( maxLevel );
   Polynomial2D< MonomialBasis2D >&                poly        = stencilPoly[vertexdof::stencilIndexFromVertex( direction )];

   for ( uint_t k = 0; k < nCoeffs; k++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " -> coeff[" << k << "] = " << poly.getCoefficient( k ) );
   }

   logMessage( "Evaluating surrogate polynomial" );
   expression surrogate = [&poly, &Phi, &vertices]( const hyteg::Point3D& x ) {
      // polynomial is given in barycentric coordinates, so we must map x
      // to the reference element here
      Point2D aux( {x[0] - vertices[0][0], x[1] - vertices[0][1]} );
      Point2D pnt = Phi.mul( aux );
      return poly.eval( pnt );
   };

   P1Function< real_t > surrogateWeight( "surrogate weights", storage, minLevel, maxLevel );
   surrogateWeight.interpolate( surrogate, maxLevel, All );

   // =================
   //  Postprocessing
   // =================
   logSectionHeader( "Postprocessing" );

   logMessage( "Checking differences:" );
   P1Function< real_t > difference( "approximation error", storage, minLevel, maxLevel );
   difference.assign( {1.0, -1.0}, {stencilWeight, surrogateWeight}, maxLevel );

   // compute error norms
   real_t errorNormMax = difference.getMaxMagnitude( maxLevel, Inner );
   real_t errorNormL2discr =
       std::sqrt( difference.dotGlobal( difference, maxLevel, Inner ) ) / real_c( getNumInnerMicroVertices( params.lsqLevel ) );

   P1ConstantMassOperator mass( storage, minLevel, maxLevel );
   P1Function< real_t >   relErr( "relative error", storage, minLevel, maxLevel );
   P1Function< real_t >   aux( "auxilliary", storage, minLevel, maxLevel );
   mass.apply( difference, aux, maxLevel, Inner );
   real_t errorNormL2 = std::sqrt( aux.dotGlobal( difference, maxLevel, Inner ) );

   relErr.assign( {1.0}, {stencilWeight}, maxLevel, All );
   relErr.invertElementwise( maxLevel, All );
   relErr.multElementwise( {difference, relErr}, maxLevel, All );
   real_t errorMaxRel = relErr.getMaxMagnitude( maxLevel, Inner );

   WALBERLA_LOG_INFO_ON_ROOT( " -> maximum norm ........ " << std::scientific << errorNormMax );
   WALBERLA_LOG_INFO_ON_ROOT( " -> discrete L2 norm .... " << std::scientific << errorNormL2discr );
   WALBERLA_LOG_INFO_ON_ROOT( " -> H0 norm ............. " << std::scientific << errorNormL2 );
   WALBERLA_LOG_INFO_ON_ROOT( " -> maximum relative .... " << std::scientific << errorMaxRel );

   logMessage( "Exporting weight functions" );
   std::stringstream stream;
   stream << "surrogates_level-" << params.evalLevel << "_lsq-" << params.lsqLevel << "_degree-" << params.degree;
   std::string oFile   = stream.str();
   std::string oFolder = "../../output";
   WALBERLA_LOG_INFO_ON_ROOT( " -> directory = " << oFolder );
   WALBERLA_LOG_INFO_ON_ROOT( " -> baseName  = " << oFile );
   VTKOutput vtkOutput( oFolder, oFile, storage );
   vtkOutput.add( kappa );
   vtkOutput.add( stencilWeight );
   vtkOutput.add( surrogateWeight );
   vtkOutput.add( difference );
   vtkOutput.add( relErr );

   // set stencil functions to zero on boundaries
   // stencilWeight.interpolate( real_c(0), maxLevel, Boundary );
   // surrogateWeight.interpolate( real_c(0), maxLevel, Boundary );
   // difference.interpolate( real_c(0), maxLevel, Boundary );

   vtkOutput.write( maxLevel );

   std::string separatorTop( 50, '=' );
   WALBERLA_LOG_INFO_ON_ROOT( "" << separatorTop );

   return 0;
}
