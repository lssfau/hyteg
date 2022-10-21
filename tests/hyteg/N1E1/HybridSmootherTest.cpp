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

#include "hyteg/n1e1functionspace/HybridSmoother.hpp"

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormCurlCurl.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormMass.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"

using namespace hyteg;

void test( const uint_t                                       level,
           const MeshInfo                                     meshInfo,
           std::function< Eigen::Vector3r( const Point3D& ) > analyticSolution,
           std::function< Eigen::Vector3r( const Point3D& ) > rhs,
           const uint_t                                       nIterations,
           const real_t                                       expectedError,
           const bool                                         writeVTK = false )
{
   using namespace n1e1;

   using P1LaplaceOperator = P1ConstantLaplaceOperator;
   using N1E1Smoother      = ChebyshevSmoother< N1E1ElementwiseLinearCombinationOperator >;
   using P1Smoother        = GaussSeidelSmoother< P1LaplaceOperator >;

   const uint_t minLevel = level;
   const uint_t maxLevel = level;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "HybridSmootherTest_partitioning" );

   N1E1VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > sol( "sol", storage, level, level );
   N1E1VectorFunction< real_t > err( "err", storage, level, level );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, level, level );

   // operators
   N1E1Form_curl_curl                       curlCurlForm;
   N1E1Form_mass                            massForm;
   N1E1ElementwiseMassOperator              M( storage, level, level );
   N1E1ElementwiseLinearCombinationOperator A( storage, minLevel, maxLevel, { { 1.0, 1.0 }, { &curlCurlForm, &massForm } } );

   tmp.interpolate( rhs, level );
   M.apply( tmp, f, level, DoFType::All );
   sol.interpolate( analyticSolution, level );

   auto p1LaplaceOperator = std::make_shared< P1LaplaceOperator >( storage, minLevel, maxLevel );

   // smoothers
   auto chebyshevSmoother = std::make_shared< N1E1Smoother >( storage, minLevel, maxLevel );

   const real_t spectralRadius = chebyshev::estimateRadius( A, level, 100, storage, sol, tmp );
   sol.interpolate( analyticSolution, level );
   chebyshevSmoother->setupCoefficients( 4, spectralRadius );
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( spectralRadius );

   std::shared_ptr< Solver< N1E1ElementwiseLinearCombinationOperator > > n1e1Smoother = chebyshevSmoother;
   std::shared_ptr< Solver< P1LaplaceOperator > >                        p1Smoother   = std::make_shared< P1Smoother >();

   // solve Au = f with hybrid smoother
   HybridSmoother hybridSmoother( storage, p1LaplaceOperator, n1e1Smoother, p1Smoother, minLevel, maxLevel );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   M.apply( err, tmp, level, DoFType::All );
   double prevDiscL2 = std::sqrt( err.dotGlobal( tmp, level ) );

   for ( uint_t i = 0; i < nIterations; ++i )
   {
      hybridSmoother.solve( A, u, f, level );

      // determine error
      err.assign( { 1.0, -1.0 }, { u, sol }, level );
      M.apply( err, tmp, level, DoFType::All );
      const real_t discrL2 = std::sqrt( err.dotGlobal( tmp, level ) );

      WALBERLA_CHECK_LESS( discrL2, prevDiscL2, "hybrid smoother, does not decrease error" );
      prevDiscL2 = discrL2;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Approximate L2 norm of error: " << prevDiscL2 );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "HybridSmootherTest", storage );
      vtk.add( u );
      vtk.add( f );
      vtk.add( sol );
      vtk.add( err );
      vtk.write( level );
   }

   // the calculated error should up to 10% match the expected error
   WALBERLA_CHECK_LESS( 0.9 * expectedError, prevDiscL2, "error for hybrid smoother has changed " );
   WALBERLA_CHECK_GREATER( 1.1 * expectedError, prevDiscL2, "error for hybrid smoother has changed " );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Single Macro ###" );

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

      test( 4, meshInfo, solFunc, rhsFunc, 100, 0.000111476 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Pyramid ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" );

      std::function< Eigen::Vector3r( const Point3D& ) > solFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{ z * ( z - 2 * x ) * ( z + 2 * x - 2 ) * ( z - 2 * y ) * ( z + 2 * y - 2 ),
                                 z * ( z - 2 * x ) * ( z + 2 * x - 2 ) * ( z - 2 * y ) * ( z + 2 * y - 2 ),
                                 1 * ( z - 2 * x ) * ( z + 2 * x - 2 ) * ( z - 2 * y ) * ( z + 2 * y - 2 ) };
      };

      // yes, I know it's ugly
      std::function< Eigen::Vector3r( const Point3D& ) > rhsFunc = []( const Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return Eigen::Vector3r{
             16 * x * x * y * y * z - 16 * x * x * y * z - 4 * x * x * std::pow( z, 3 ) + 8 * x * x * z * z - 8 * x * x * z -
                 16 * x * x - 16 * x * y * y * z + 80 * x * y * z + 4 * x * std::pow( z, 3 ) - 8 * x * z * z - 40 * x * z +
                 32 * x - 4 * y * y * std::pow( z, 3 ) + 8 * y * y * z * z + 24 * y * y * z - 16 * y * y +
                 4 * y * std::pow( z, 3 ) - 8 * y * z * z - 56 * y * z + 16 * y + std::pow( z, 5 ) - 4 * std::pow( z, 4 ) -
                 8 * std::pow( z, 3 ) + 32 * z * z - 8,
             16 * x * x * y * y * z - 16 * x * x * y * z - 4 * x * x * std::pow( z, 3 ) + 8 * x * x * z * z + 24 * x * x * z -
                 16 * x * x - 16 * x * y * y * z + 80 * x * y * z + 4 * x * std::pow( z, 3 ) - 8 * x * z * z - 56 * x * z +
                 16 * x - 4 * y * y * std::pow( z, 3 ) + 8 * y * y * z * z - 8 * y * y * z - 16 * y * y +
                 4 * y * std::pow( z, 3 ) - 8 * y * z * z - 40 * y * z + 32 * y + std::pow( z, 5 ) - 4 * std::pow( z, 4 ) -
                 8 * std::pow( z, 3 ) + 32 * z * z - 8,
             16 * x * x * y * y + 16 * x * x * y - 4 * x * x * z * z + 8 * x * x * z - 48 * x * x + 16 * x * y * y - 48 * x * y -
                 20 * x * z * z + 24 * x * z + 48 * x - 4 * y * y * z * z + 8 * y * y * z - 48 * y * y - 20 * y * z * z +
                 24 * y * z + 48 * y + std::pow( z, 4 ) - 4 * std::pow( z, 3 ) + 44 * z * z - 64 * z };
      };

      test( 4, meshInfo, solFunc, rhsFunc, 100, 0.0269902 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Cube ###" );

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

      test( 4, meshInfo, solFunc, rhsFunc, 100, 0.0349963 );
   }

   return EXIT_SUCCESS;
}
