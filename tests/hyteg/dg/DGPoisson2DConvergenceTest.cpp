/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

/// Returns the scaled L2 error.
real_t test( uint_t                                    level,
             uint_t                                    degree,
             MeshInfo                                  meshInfo,
             std::function< real_t( const Point3D& ) > solFunc,
             std::function< real_t( const Point3D& ) > rhsFunc,
             bool                                      writeVTK = false )
{
   using namespace dg;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   real_t beta_0 = storage->hasGlobalCells() ? 0.5 : 1.0;

   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm = std::make_shared< DGDiffusionForm_Example >( beta_0, solFunc, solFunc );
   auto massForm    = std::make_shared< DGMassForm_Example >();

   DGFunction< real_t > u( "u", storage, level, level, basis, degree );
   DGFunction< real_t > f( "f", storage, level, level, basis, degree );
   DGFunction< real_t > sol( "sol", storage, level, level, basis, degree );
   DGFunction< real_t > tmp( "tmp", storage, level, level, basis, degree );
   DGFunction< real_t > err( "err", storage, level, level, basis, degree );

   WALBERLA_LOG_INFO_ON_ROOT( "dofs: " << u.getNumberOfGlobalDoFs( level ) );

   DGFunction< idx_t > numerator( "numerator", storage, level, level, basis, degree );
   numerator.enumerate( level );

   DGOperator A( storage, level, level, laplaceForm );
   DGOperator M( storage, level, level, massForm );

   // Assemble RHS.
   f.evaluateLinearFunctional( rhsFunc, level );
   f.applyDirichletBoundaryConditions( laplaceForm, level );

   // Interpolate solution
   tmp.evaluateLinearFunctional( solFunc, level );
   PETScCGSolver< DGOperator > solverM( storage, level, numerator );
   solverM.solve( M, sol, tmp, level );

   // Solve system.
   PETScCGSolver< DGOperator > solverA( storage, level, numerator, 1e-12, 1e-12, 10000 );
   solverA.solve( A, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   auto discrL2 = sqrt( err.dotGlobal( err, level ) / real_c( numberOfGlobalDoFs( u, level ) ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "DGPoisson3DConvergenceTest", storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( f );
      vtk.write( level );
   }

   return discrL2;
}

void runTest( uint_t                                           minLevel,
              uint_t                                           maxLevel,
              uint_t                                           minDegree,
              uint_t                                           maxDegree,
              const MeshInfo&                                  meshInfo,
              const std::function< real_t( const Point3D& ) >& solFunc,
              const std::function< real_t( const Point3D& ) >& rhsFunc )
{
   for ( uint_t degree = minDegree; degree <= maxDegree; degree++ )
   {
      auto l2ConvRate  = std::pow( 2, -( int( degree ) + 1 ) );
      auto convRateEps = l2ConvRate * 0.1;
      auto err         = hyteg::test( minLevel, degree, meshInfo, solFunc, rhsFunc );
      WALBERLA_LOG_INFO_ON_ROOT( "degree " << degree << ", expected L2 rate: " << l2ConvRate
                                           << ", threshold: " << l2ConvRate + convRateEps );
      WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );
      for ( uint_t l = minLevel + 1; l <= maxLevel; l++ )
      {
         auto errFiner     = hyteg::test( l, degree, meshInfo, solFunc, rhsFunc );
         auto computedRate = errFiner / err;

         WALBERLA_LOG_INFO_ON_ROOT( "error level " << l << ": " << errFiner );
         WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << l << " / " << l - 1 << ": " << computedRate );

         WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                    l2ConvRate + convRateEps,
                                    "Convergence L2 rate level " << l << " vs level " << l - 1
                                                                 << " not sufficiently small (computed: " << computedRate
                                                                 << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );
         err = errFiner;
      }
   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   using hyteg::MeshInfo;
   using hyteg::Point2D;
   using hyteg::Point3D;
   using walberla::real_t;
   using walberla::math::pi;

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * ( x[0] + x[1] - 1 ) );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 4 * pi * pi * ( -2 * sin( 4 * pi * ( x[0] + x[1] ) ) + sin( 4 * pi * x[0] ) + sin( 4 * pi * x[1] ) );
      };

      hyteg::runTest( 3, 6, 1, 1, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 8 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) );
      };

      hyteg::runTest( 3, 6, 1, 1, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, inhom. BC, rhs = 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };

      hyteg::runTest( 2, 6, 1, 1, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, inhom. BC, rhs != 0 ###" );

      MeshInfo meshInfo =
          hyteg::MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), hyteg::MeshInfo::CRISS, 2, 2 );

      std::function< real_t( const hyteg::Point3D& ) > solFunc = []( const hyteg::Point3D& x ) {
         return std::exp( -x[0] - ( x[1] * x[1] ) );
      };
      std::function< real_t( const hyteg::Point3D& ) > rhsFunc = []( const hyteg::Point3D& x ) {
         return -( 4 * x[1] * x[1] - 1 ) * std::exp( -x[0] - ( x[1] * x[1] ) );
      };

      hyteg::runTest( 3, 6, 1, 1, meshInfo, solFunc, rhsFunc );
   }

   return EXIT_SUCCESS;
}
