/*
* Copyright (c) 2022 Daniel Bauer.
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

#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormCurlCurl.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormMass.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::math::pi;

/// Returns the approximate L2 error.
real_t test( uint_t                                             level,
             MeshInfo                                           meshInfo,
             std::function< Eigen::Vector3r( const Point3D& ) > solFunc,
             std::function< Eigen::Vector3r( const Point3D& ) > rhsFunc,
             bool                                               writeVTK = false )
{
   using namespace n1e1;
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   n1e1::N1E1Form_curl_curl curlCurlForm;
   n1e1::N1E1Form_mass      massForm;

   N1E1ElementwiseMassOperator              M( storage, level, level );
   N1E1ElementwiseLinearCombinationOperator A( storage, level, level, { { 1.0, 1.0 }, { &curlCurlForm, &massForm } } );

   N1E1VectorFunction< real_t > u( "u", storage, level, level );
   N1E1VectorFunction< real_t > f( "f", storage, level, level );
   N1E1VectorFunction< real_t > sol( "sol", storage, level, level );
   N1E1VectorFunction< real_t > err( "err", storage, level, level );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, level, level );

   const uint_t nDoFs = numberOfGlobalDoFs( u, level );
   WALBERLA_LOG_INFO_ON_ROOT( "dofs on level " << level << ": " << nDoFs );

   // Assemble RHS.
   tmp.interpolate( rhsFunc, level );
   M.apply( tmp, f, level, DoFType::All );

   // Boundary conditions: homogeneous tangential trace
   u.interpolate( Eigen::Vector3r{ 0.0, 0.0, 0.0 }, level, DoFType::Boundary );

   // Interpolate solution
   sol.interpolate( solFunc, level );

   // Solve system.
#ifdef HYTEG_BUILD_WITH_PETSC
   WALBERLA_LOG_INFO_ON_ROOT( "Using PETSc solver" )
   auto solverA = PETScCGSolver< n1e1::N1E1ElementwiseLinearCombinationOperator >( storage, level );
#else
   WALBERLA_LOG_INFO_ON_ROOT( "Using HyTeG solver" )
   auto solverA = CGSolver< n1e1::N1E1ElementwiseLinearCombinationOperator >( storage, level, level, 10000, 1e-12 );
#endif
   solverA.solve( A, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   M.apply( err, tmp, level, DoFType::All );
   const real_t discrL2 = std::sqrt( err.dotGlobal( tmp, level ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "N1E1CurlCurlConvergenceTest", storage );
      vtk.add( u );
      vtk.add( f );
      vtk.add( sol );
      vtk.add( err );
      vtk.write( level );
   }

   return discrL2;
}

void runTest( uint_t                                                    minLevel,
              uint_t                                                    maxLevel,
              const MeshInfo&                                           meshInfo,
              const std::function< Eigen::Vector3r( const Point3D& ) >& solFunc,
              const std::function< Eigen::Vector3r( const Point3D& ) >& rhsFunc )
{
   const real_t l2ConvRate  = 1.0 / 4.0;
   const real_t convRateEps = l2ConvRate * 0.1;
   real_t       err         = test( minLevel, meshInfo, solFunc, rhsFunc );

   WALBERLA_LOG_INFO_ON_ROOT( "expected L2 rate: " << l2ConvRate << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );

   for ( uint_t level = minLevel + 1; level <= maxLevel; level++ )
   {
      const real_t errFiner     = test( level, meshInfo, solFunc, rhsFunc );
      const real_t computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << level << ": " << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << level << " / " << level - 1 << ": " << computedRate );

      WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << level << " vs level " << level - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );
      err = errFiner;
   }
}

int main( int argc, char** argv )
{
   using std::sin;

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, polynomial ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );

      std::function< Eigen::Vector3r( const Point3D& ) > solFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ y * z * ( x + y + z - 1 ), x * z * ( x + y + z - 1 ), x * y * ( x + y + z - 1 ) };
      };

      std::function< Eigen::Vector3r( const Point3D& ) > rhsFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{
             y * z * ( x + y + z - 1 ) - y - z, x * z * ( x + y + z - 1 ) - x - z, x * y * ( x + y + z - 1 ) - x - y };
      };

      runTest( 3, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, sinusoidal ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );

      std::function< Eigen::Vector3r( const Point3D& ) > solFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ sin( y ) * sin( z ) * sin( x + y + z - 1 ),
                                 sin( x ) * sin( z ) * sin( x + y + z - 1 ),
                                 sin( x ) * sin( y ) * sin( x + y + z - 1 ) };
      };

      std::function< Eigen::Vector3r( const Point3D& ) > rhsFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ -sin( x ) * sin( y ) * sin( x + y + z - 1 ) - sin( x ) * sin( z ) * sin( x + y + z - 1 ) +
                                     5 * sin( y ) * sin( z ) * sin( x + y + z - 1 ) + sin( y ) * cos( x ) * cos( x + y + z - 1 ) -
                                     2 * sin( y ) * cos( z ) * cos( x + y + z - 1 ) + sin( z ) * cos( x ) * cos( x + y + z - 1 ) -
                                     2 * sin( z ) * cos( y ) * cos( x + y + z - 1 ),
                                 -sin( x ) * sin( y ) * sin( x + y + z - 1 ) + 5 * sin( x ) * sin( z ) * sin( x + y + z - 1 ) +
                                     sin( x ) * cos( y ) * cos( x + y + z - 1 ) - 2 * sin( x ) * cos( z ) * cos( x + y + z - 1 ) -
                                     sin( y ) * sin( z ) * sin( x + y + z - 1 ) - 2 * sin( z ) * cos( x ) * cos( x + y + z - 1 ) +
                                     sin( z ) * cos( y ) * cos( x + y + z - 1 ),
                                 5 * sin( x ) * sin( y ) * sin( x + y + z - 1 ) - sin( x ) * sin( z ) * sin( x + y + z - 1 ) -
                                     2 * sin( x ) * cos( y ) * cos( x + y + z - 1 ) + sin( x ) * cos( z ) * cos( x + y + z - 1 ) -
                                     sin( y ) * sin( z ) * sin( x + y + z - 1 ) - 2 * sin( y ) * cos( x ) * cos( x + y + z - 1 ) +
                                     sin( y ) * cos( z ) * cos( x + y + z - 1 ) };
      };

      runTest( 4, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, polynomial ###" );

      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

      std::function< Eigen::Vector3r( const Point3D& ) > solFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ y * ( 1 - y ) * z * ( 1 - z ), x * ( 1 - x ) * z * ( 1 - z ), x * ( 1 - x ) * y * ( 1 - y ) };
      };

      std::function< Eigen::Vector3r( const Point3D& ) > rhsFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ 2 * ( y * ( 1 - y ) + z * ( 1 - z ) ) + y * ( 1 - y ) * z * ( 1 - z ),
                                 2 * ( x * ( 1 - x ) + z * ( 1 - z ) ) + x * ( 1 - x ) * z * ( 1 - z ),
                                 2 * ( x * ( 1 - x ) + y * ( 1 - y ) ) + x * ( 1 - x ) * y * ( 1 - y ) };
      };

      runTest( 4, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, sinusoidal ###" );

      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

      std::function< Eigen::Vector3r( const Point3D& ) > solFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ sin( 2 * pi * y ) * sin( 2 * pi * z ),
                                 sin( 2 * pi * x ) * sin( 2 * pi * z ),
                                 sin( 2 * pi * x ) * sin( 2 * pi * y ) };
      };

      std::function< Eigen::Vector3r( const Point3D& ) > rhsFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ sin( 2 * pi * y ) * sin( 2 * pi * z ) + 8 * pi * pi * sin( 2 * pi * y ) * sin( 2 * pi * z ),
                                 sin( 2 * pi * x ) * sin( 2 * pi * z ) + 8 * pi * pi * sin( 2 * pi * x ) * sin( 2 * pi * z ),
                                 sin( 2 * pi * x ) * sin( 2 * pi * y ) + 8 * pi * pi * sin( 2 * pi * x ) * sin( 2 * pi * y ) };
      };

      runTest( 4, 5, meshInfo, solFunc, rhsFunc );
   }

   // TODO inhomogeneous BCs?

   return EXIT_SUCCESS;
}
