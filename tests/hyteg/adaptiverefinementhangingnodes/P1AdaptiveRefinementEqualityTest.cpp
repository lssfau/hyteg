/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

/// Checks if the computed solution differs when comparing globally uniform macro-refinement with micro-refinement.
real_t cmpMacroMicroRefinementTest( uint_t microRefinementLevel, uint_t macroRefinementLevel, const MeshInfo& meshInfo )
{
   const bool writeVTK = false;

   std::function< real_t( const hyteg::Point3D& ) > solFunc = []( const hyteg::Point3D& x ) {
      return std::exp( -x[0] - ( x[1] * x[1] ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhsFunc = []( const hyteg::Point3D& x ) {
      return -( 4 * x[1] * x[1] - 1 ) * std::exp( -x[0] - ( x[1] * x[1] ) );
   };

   const auto level      = microRefinementLevel;
   const auto refinement = macroRefinementLevel;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // Refine storage
   for ( uint_t r = 0; r < refinement; r++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Refining (step " << r + 1 << "/" << refinement << ") ..." );
      auto faceIDs = storage->getFaceIDs();
      storage->refinementAndCoarseningHanging( faceIDs, {} );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Done refining." );

   WALBERLA_LOG_INFO_ON_ROOT( "Allocation ..." );

   P1Function< real_t > u( "u", storage, level, level );
   P1Function< real_t > f( "f", storage, level, level );
   P1Function< real_t > sol( "sol", storage, level, level );
   P1Function< real_t > tmp( "tmp", storage, level, level );
   P1Function< real_t > Merr( "Merr", storage, level, level );
   P1Function< real_t > err( "err", storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "DoFs: " << numberOfGlobalDoFs( u, level ) );

   P1Function< idx_t > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   P1ConstantLaplaceOperator A( storage, level, level );
   P1ConstantMassOperator    M( storage, level, level );

   tmp.interpolate( rhsFunc, level );
   M.apply( tmp, f, level, All, Replace );

   // Interpolate solution
   u.interpolate( solFunc, level, DirichletBoundary );
   sol.interpolate( solFunc, level );

   // Solve system.
   PETScCGSolver< P1ConstantLaplaceOperator > solverA( storage, level, numerator, 1e-12, 1e-12, 10000 );
   solverA.solve( A, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );

   M.apply( err, Merr, level, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, level, Inner ) );

   WALBERLA_LOG_INFO_ON_ROOT( "L2 error: " << discrL2 );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/",
                     "P1AdaptiveRefinementEqualityTest_MacroRef_" + std::to_string( refinement ) + "_MicroRef_" +
                         std::to_string( level ),
                     storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( f );
      vtk.write( level );
   }

   return discrL2;
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
   using hyteg::prependHyTeGMeshDir;
   using walberla::real_t;

   const uint_t sum = 5;

   MeshInfo meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/quad_16el.msh" ) );

   for ( uint_t a = 1; a <= sum / 2; a++ )
   {
      const uint_t b = sum - a;

      WALBERLA_LOG_INFO_ON_ROOT( "### Comparing " << a << " macro refinements w/ " << b << " levels and vice versa." );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      const auto e1 = hyteg::cmpMacroMicroRefinementTest( a, b, meshInfo );
      const auto e2 = hyteg::cmpMacroMicroRefinementTest( b, a, meshInfo );

      const auto diff = std::abs( e1 - e2 );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "### Diff: " << diff );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      WALBERLA_CHECK_LESS( diff, 1e-15 );
   }

   return EXIT_SUCCESS;
}
