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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormCurlCurl.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormMass.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

#include "common.hpp"

using namespace hyteg;
using walberla::real_t;

/// Returns the approximate L2 error.
real_t test( const uint_t level, const n1e1::System& system, const bool writeVTK = false )
{
   using namespace n1e1;
   SetupPrimitiveStorage setupStorage( system.domain_, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
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
   tmp.interpolate( system.rhs_, level );
   M.apply( tmp, f, level, DoFType::All );

   // Boundary conditions: homogeneous tangential trace
   u.interpolate( Eigen::Vector3r{ 0.0, 0.0, 0.0 }, level, DoFType::Boundary );

   // Interpolate solution
   sol.interpolate( system.analyticalSol_, level );

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

int main( int argc, char** argv )
{
   using std::sin;

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, polynomial ###" );
   n1e1::L2ConvergenceTest( 3, 5, n1e1::System::polynomialOnSingleTet(), test );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, sinusoidal ###" );
   n1e1::L2ConvergenceTest( 4, 5, n1e1::System::sinusoidalOnSingleTet(), test );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, polynomial ###" );
   n1e1::L2ConvergenceTest( 4, 5, n1e1::System::polynomialOnCube(), test );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, sinusoidal ###" );
   n1e1::L2ConvergenceTest( 4, 5, n1e1::System::sinusoidalOnCube(), test );

   // TODO inhomogeneous BCs?

   return EXIT_SUCCESS;
}
