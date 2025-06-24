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
#include "hyteg/forms/form_hyteg_dg/DGDiffusionForm_Example.hpp"
#include "hyteg/forms/form_hyteg_dg/DGMassForm_Example.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

/// Successively refines an entire macro-face to see if global refinement works.
/// Re-allocates all data in each run.
///
/// Should produce correct error convergence rates across either refinement route
/// (i.e. regardless of whether micro- or macro-refinement is chosen)
///
/// \param minLevel min micro-refinement level
/// \param maxLevel max micro-refinement level
/// \param degree   element poly degree
void singleMacroFaceGlobalRefinement( uint_t minLevel, uint_t maxLevel, uint_t degree )
{
   using namespace dg;

   const uint_t numGlobalRefinements = 4;
   const bool   writeVTK             = true;

   // macro-level -> micro-level -> error
   std::map< uint_t, std::map< uint_t, real_t > > l2Errors;

   std::function< real_t( const hyteg::Point3D& ) > solFunc = []( const hyteg::Point3D& x ) {
      return std::exp( -x[0] - ( x[1] * x[1] ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhsFunc = []( const hyteg::Point3D& x ) {
      return -( 4 * x[1] * x[1] - 1 ) * std::exp( -x[0] - ( x[1] * x[1] ) );
   };

   for ( uint_t refinement = 0; refinement < numGlobalRefinements; refinement++ )
   {
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         MeshInfo              meshInfo = MeshInfo::meshFaceChain( 1 );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
         setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
         std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

         // Refine storage
         for ( uint_t r = 0; r < refinement; r++ )
         {
            auto faceIDs = storage->getFaceIDs();
            storage->refinementAndCoarseningHanging( faceIDs, {} );
         }

         real_t beta_0 = storage->hasGlobalCells() ? 0.5 : 1.0;

         auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
         auto laplaceForm = std::make_shared< DGDiffusionForm_Example >( beta_0, solFunc, solFunc );
         auto massForm    = std::make_shared< DGMassForm_Example >();

         DGFunction< real_t > u( "u", storage, level, level, basis, degree );
         DGFunction< real_t > f( "f", storage, level, level, basis, degree );
         DGFunction< real_t > sol( "sol", storage, level, level, basis, degree );
         DGFunction< real_t > tmp( "tmp", storage, level, level, basis, degree );
         DGFunction< real_t > Merr( "Merr", storage, level, level, basis, degree );
         DGFunction< real_t > err( "err", storage, level, level, basis, degree );

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

         M.apply( err, Merr, level, All, Replace );
         auto discrL2 = sqrt( err.dotGlobal( Merr, level, Inner ) );

         l2Errors[refinement][level] = discrL2;

         WALBERLA_LOG_INFO_ON_ROOT( "error (macro-ref: " << refinement << " | level: " << level << "): " << discrL2 );

         if ( writeVTK )
         {
            VTKOutput vtk( "../../output/", "DGGlobalRefinementTest_SingleFace_" + std::to_string( refinement ), storage );
            vtk.add( u );
            vtk.add( sol );
            vtk.add( err );
            vtk.add( f );
            vtk.write( level );
         }
      }

      WALBERLA_LOG_INFO_ON_ROOT( "" );
   }

   // Comparing convergence rates across micro-refinement levels
   for ( uint_t i = 0; i < numGlobalRefinements; i++ )
   {
      auto coarseError = l2Errors[i][minLevel];
      for ( uint_t ii = minLevel + 1; ii <= maxLevel; ii++ )
      {
         auto fineError = l2Errors[i][ii];
         auto rateMicro = fineError / coarseError;
         WALBERLA_LOG_INFO_ON_ROOT( "micro conv rate: " << rateMicro );
         WALBERLA_CHECK_LESS( rateMicro, 0.275 );
         coarseError = fineError;
      }
   }

   // Comparing convergence rates across macro-refinement levels
   for ( uint_t i = minLevel; i <= maxLevel; i++ )
   {
      auto coarseError = l2Errors[0][i];
      for ( uint_t ii = 1; ii < numGlobalRefinements; ii++ )
      {
         auto fineError = l2Errors[ii][i];
         auto rateMacro = fineError / coarseError;
         WALBERLA_LOG_INFO_ON_ROOT( "macro conv rate: " << rateMacro );
         WALBERLA_CHECK_LESS( rateMacro, 0.275 );
         coarseError = fineError;
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

   hyteg::singleMacroFaceGlobalRefinement( 3, 5, 1 );

   return EXIT_SUCCESS;
}
