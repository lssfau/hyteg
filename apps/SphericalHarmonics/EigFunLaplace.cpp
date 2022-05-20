/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

// Setup an analytic eigenfunction for the Laplacian on a ball of fixed
// radius with Dirichlet boundary conditions. Apply Laplace and mass
// operator to it and compare (scaled) results.

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>

#ifndef __cpp_lib_math_special_functions
#if !( __STDCPP_MATH_SPEC_FUNCS__ >= 201003L )
#error "Sorry, app needs some special math functions!"
#endif
#endif

#include <array>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

using terraneo::SphericalHarmonicsTool;
using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // ============
   //  Parameters
   // ============

   // check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./EigFunLaplace.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle params = cfg->getBlock( "Parameters" );
   if ( walberla::MPIManager::instance()->worldRank() == 0 )
   {
      params.listParameters();
   }

   // =========
   //  Meshing
   // =========

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Preparing mesh ..." )

   const uint_t level = params.getParameter< uint_t >( "level" );
   const uint_t nRad  = params.getParameter< uint_t >( "nRad" );
   const uint_t nTan  = params.getParameter< uint_t >( "nTan" );

   real_t outerRad = 2.0;
   real_t innerRad = 1.0;

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, innerRad, outerRad );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   IcosahedralShellMap::setMap( setupStorage );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   // ===============
   //  Eigenfunction
   // ===============

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Setting up eigenfunction ..." )

   uint_t degree = params.getParameter< uint_t >( "degree" );
   int    order  = params.getParameter< int >( "order" );
   uint_t nroot  = params.getParameter< uint_t >( "nroot" );

   // check validity of inputs
   if ( degree >= 5 )
   {
      WALBERLA_ABORT( "Please choose a degree from {0,1,2,3,4}." );
   }
   if ( nroot < 1 || nroot > 4 )
   {
      WALBERLA_ABORT( "Please choose nroot from {1,2,3,4}." );
   }
   if ( order < 0 || walberla::uint_c( std::abs( order ) ) > degree )
   {
      WALBERLA_ABORT( "Order cannot exceed degree." );
   }

   // set first roots of first spherical Bessel functions
   std::array< std::array< double, 4 >, 5 > roots;
   roots[0] = {3.141592653592e+00, 6.283185307185e+00, 9.424777960777e+00, 1.256637061437e+01};
   roots[1] = {4.493409457908e+00, 7.725251836935e+00, 1.090412165944e+01, 1.406619391284e+01};
   roots[2] = {5.763459196896e+00, 9.095011330470e+00, 1.232294097056e+01, 1.551460301090e+01};
   roots[3] = {6.987932000507e+00, 1.041711854738e+01, 1.369802315324e+01, 1.692362128522e+01};
   roots[4] = {8.182561452563e+00, 1.170490715457e+01, 1.503966470762e+01, 1.830125595953e+01};

   // prepare computation of spherical harmonics
   uint_t                                    lmax    = degree;
   std::shared_ptr< SphericalHarmonicsTool > sphTool = std::make_shared< SphericalHarmonicsTool >( lmax );

   // describe eigenfunction
   real_t                                    besselRoot    = roots[degree][nroot - 1];
   std::function< real_t( const Point3D& ) > eigenFunction = [sphTool, degree, order, besselRoot, outerRad]( const Point3D& x ) {
      real_t sph = sphTool->shconvert_eval( degree, order, x[0], x[1], x[2] );
      real_t rad = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      real_t arg = besselRoot * rad / outerRad;
      real_t scl = std::sqrt( real_c( degree ) * outerRad * 0.5 / besselRoot );
      return scl * std::sph_bessel( (unsigned int) degree, arg ) * sph;
   };

   P2Function< real_t > eigfunc( "eigenfunction", storage, level, level );
   eigfunc.interpolate( eigenFunction, level, All );

   // compute analytic eigenvalue
   real_t eigenval = ( besselRoot / outerRad ) * ( besselRoot / outerRad );
   WALBERLA_LOG_INFO_ON_ROOT( "Eigenvalue is " << std::scientific << eigenval )

   // ============
   //  Comparison
   // ============

   // apply Laplace operator
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Applying Laplace operator ..." )
   P2ElementwiseBlendingLaplaceOperator lapOp( storage, level, level );
   P2Function< real_t >                 lapApplied( "Laplace operator applied", storage, level, level );
   lapOp.apply( eigfunc, lapApplied, level, Inner );

   // apply Mass operator and scale result be eigenvalue
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Applying Mass operator ..." )
   P2ElementwiseBlendingMassOperator massOp( storage, level, level );
   P2Function< real_t >              massApplied( "Mass operator applied", storage, level, level );
   massOp.apply( eigfunc, massApplied, level, Inner );
   massApplied.assign( {eigenval}, {massApplied}, level );

   // compute difference
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Checking difference ..." )
   P2Function< real_t > diff( "Difference", storage, level, level );
   diff.assign( {1.0, -1.0}, {massApplied, lapApplied}, level );
   uint_t ndofs = numberOfGlobalInnerDoFs< hyteg::P2FunctionTag >( *storage, level );
   real_t norm  = std::sqrt( diff.dotGlobal( diff, level ) ) / real_c( ndofs );
   WALBERLA_LOG_INFO_ON_ROOT( "Discrete L_2 norm of difference = " << norm );

   // store results
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Exporting results to VTK ..." )
   hyteg::VTKOutput vtkOutput( "./output", "EigFunLaplace", storage );
   vtkOutput.add( eigfunc );
   vtkOutput.add( lapApplied );
   vtkOutput.add( massApplied );
   vtkOutput.write( level, 0 );

   return EXIT_SUCCESS;
}
