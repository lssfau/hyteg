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
#include "hyteg/gridtransferoperators/N1E1toN1E1Prolongation.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Restriction.hpp"
#include "hyteg/n1e1functionspace/HybridSmoother.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::math::pi;

/// Returns the approximate L2 error.
real_t test( const uint_t                                             maxLevel,
             const MeshInfo                                           meshInfo,
             const std::function< Eigen::Vector3r( const Point3D& ) > solFunc,
             const std::function< Eigen::Vector3r( const Point3D& ) > rhsFunc,
             const bool                                               writeVTK = false )
{
   using namespace n1e1;

   using P1LaplaceOperator = P1ConstantLaplaceOperator;
   using N1E1Smoother      = ChebyshevSmoother< N1E1ElementwiseLinearCombinationOperator >;
   using P1Smoother        = GaussSeidelSmoother< P1LaplaceOperator >;

   const uint_t minLevel                = 0;
   const uint_t spectralRadiusEstLevel  = std::min( uint_c( 3 ), maxLevel );
   const int    numSpectralRadiusEstIts = 40;
   const int    numVCycles              = 8;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   N1E1Form_curl_curl curlCurlForm;
   N1E1Form_mass      massForm;

   N1E1ElementwiseMassOperator              M( storage, minLevel, maxLevel );
   N1E1ElementwiseLinearCombinationOperator A( storage, minLevel, maxLevel, { { 1.0, 1.0 }, { &curlCurlForm, &massForm } } );

   N1E1VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > sol( "sol", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > err( "err", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   const uint_t nDoFs = numberOfGlobalDoFs( u, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "dofs on level " << maxLevel << ": " << nDoFs );

   // Assemble RHS.
   tmp.interpolate( rhsFunc, maxLevel );
   M.apply( tmp, f, maxLevel, DoFType::All );

   // Boundary conditions: homogeneous tangential trace
   u.interpolate( Eigen::Vector3r{ 0.0, 0.0, 0.0 }, maxLevel, DoFType::Boundary );

   // Hybrid smoother
   auto p1LaplaceOperator = std::make_shared< P1LaplaceOperator >( storage, minLevel, maxLevel );

   auto chebyshevSmoother = std::make_shared< N1E1Smoother >( storage, minLevel, maxLevel );
   auto p1Smoother        = std::make_shared< P1Smoother >();

   sol.interpolate( solFunc, spectralRadiusEstLevel );
   const real_t spectralRadius =
       chebyshev::estimateRadius( A, spectralRadiusEstLevel, numSpectralRadiusEstIts, storage, sol, tmp );
   chebyshevSmoother->setupCoefficients( 4, spectralRadius );
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( spectralRadius );

   auto hybridSmoother = std::make_shared< HybridSmoother< N1E1ElementwiseLinearCombinationOperator, P1LaplaceOperator > >(
       storage, p1LaplaceOperator, chebyshevSmoother, p1Smoother, minLevel, maxLevel );

   // GMG solver
#ifdef HYTEG_BUILD_WITH_PETSC
   WALBERLA_LOG_INFO_ON_ROOT( "Using PETSc solver" )
   auto coarseGridSolver = std::make_shared< PETScCGSolver< N1E1ElementwiseLinearCombinationOperator > >( storage, minLevel );
#else
   WALBERLA_LOG_INFO_ON_ROOT( "Using HyTeG solver" )
   auto coarseGridSolver =
       std::make_shared< CGSolver< N1E1ElementwiseLinearCombinationOperator > >( storage, minLevel, minLevel, 10000, 1e-12 );
#endif
   auto restrictionOperator  = std::make_shared< N1E1toN1E1Restriction >();
   auto prolongationOperator = std::make_shared< N1E1toN1E1Prolongation >();

   auto gmgSolver = GeometricMultigridSolver< N1E1ElementwiseLinearCombinationOperator >(
       storage, hybridSmoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

   // Interpolate solution
   sol.interpolate( solFunc, maxLevel );

   // Solve system.
   real_t discrL2 = 0.0;

   for ( int i = 0; i < numVCycles; ++i )
   {
      gmgSolver.solve( A, u, f, maxLevel );

      err.assign( { 1.0, -1.0 }, { u, sol }, maxLevel );
      M.apply( err, tmp, maxLevel, DoFType::All );
      discrL2 = std::sqrt( err.dotGlobal( tmp, maxLevel ) );
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( discrL2 )
   }

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "N1E1GMGConvergenceTest", storage );
      vtk.add( u );
      vtk.add( f );
      vtk.add( sol );
      vtk.add( err );
      vtk.write( maxLevel );
   }

   return discrL2;
}

// TODO extract to common header file
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
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << std::scientific << err );

   for ( uint_t level = minLevel + 1; level <= maxLevel; level++ )
   {
      const real_t errFiner     = test( level, meshInfo, solFunc, rhsFunc );
      const real_t computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << level << ": " << std::scientific << errFiner );
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

      runTest( 4, 6, meshInfo, solFunc, rhsFunc );
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

      runTest( 4, 6, meshInfo, solFunc, rhsFunc );
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

      runTest( 4, 6, meshInfo, solFunc, rhsFunc );
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

      runTest( 4, 6, meshInfo, solFunc, rhsFunc );
   }

   return EXIT_SUCCESS;
}
