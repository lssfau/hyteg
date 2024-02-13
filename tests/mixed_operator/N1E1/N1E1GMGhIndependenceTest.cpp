/*
* Copyright (c) 2022-2023 Daniel Bauer.
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

#include "core/mpi/Environment.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_curl_curl_affine_q0.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_mass_affine_qe.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Prolongation.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Restriction.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

#include "../../hyteg/N1E1/common.hpp"
#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "mixed_operator/HybridSmoother.hpp"

using namespace hyteg;

/// returns the number of iterations to reduce residual by 1e-6.
uint_t test( const uint_t maxLevel, const uint_t numMaxVCycles, const n1e1::System& system, const bool writeVTK = false )
{
   using namespace n1e1;

   using P1LaplaceOperator = P1ConstantLaplaceOperator;
   using N1E1Smoother      = ChebyshevSmoother< N1E1ElementwiseLinearCombinationOperator >;
   using P1Smoother        = GaussSeidelSmoother< P1LaplaceOperator >;

   const uint_t minLevel                = 0;
   const uint_t spectralRadiusEstLevel  = std::min( uint_c( 3 ), maxLevel );
   const uint_t numSpectralRadiusEstIts = 40;

   SetupPrimitiveStorage setupStorage( system.domain_, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   forms::n1e1_curl_curl_affine_q0 curlCurlForm;
   forms::n1e1_mass_affine_qe      massForm;

   N1E1ElementwiseMassOperator              M( storage, minLevel, maxLevel );
   N1E1ElementwiseLinearCombinationOperator A( storage, minLevel, maxLevel, { { 1.0, 1.0 }, { &curlCurlForm, &massForm } } );

   N1E1VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > sol( "sol", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > res( "res", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   // Assemble RHS.
   assembleLinearForm< forms::n1e1_linear_form_affine_q6 >( maxLevel, maxLevel, { system.rhs_ }, f );

   // Boundary conditions: homogeneous tangential trace
   u.interpolate( Point3D{ 0.0, 0.0, 0.0 }, maxLevel, DoFType::Boundary );

   // Hybrid smoother
   auto p1LaplaceOperator = std::make_shared< P1LaplaceOperator >( storage, minLevel, maxLevel );

   auto chebyshevSmoother = std::make_shared< N1E1Smoother >( storage, minLevel, maxLevel );
   auto p1Smoother        = std::make_shared< P1Smoother >();

   sol.interpolate( system.analyticalSol_, spectralRadiusEstLevel );
   const real_t spectralRadius =
       chebyshev::estimateRadius( A, spectralRadiusEstLevel, numSpectralRadiusEstIts, storage, sol, tmp );
   chebyshevSmoother->setupCoefficients( 4, spectralRadius );
   // WALBERLA_LOG_DEVEL_VAR_ON_ROOT( spectralRadius );

   auto hybridSmoother = std::make_shared< HybridSmoother< N1E1ElementwiseLinearCombinationOperator, P1LaplaceOperator > >(
       storage, p1LaplaceOperator, chebyshevSmoother, p1Smoother, minLevel, maxLevel );

   // GMG solver
#ifdef HYTEG_BUILD_WITH_PETSC
   auto coarseGridSolver = std::make_shared< PETScCGSolver< N1E1ElementwiseLinearCombinationOperator > >( storage, minLevel );
#else
   auto coarseGridSolver =
       std::make_shared< CGSolver< N1E1ElementwiseLinearCombinationOperator > >( storage, minLevel, minLevel, 10000, 1e-12 );
#endif
   auto restrictionOperator  = std::make_shared< N1E1toN1E1Restriction >();
   auto prolongationOperator = std::make_shared< N1E1toN1E1Prolongation >();

   auto gmgSolver = GeometricMultigridSolver< N1E1ElementwiseLinearCombinationOperator >(
       storage, hybridSmoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

   // Interpolate solution
   sol.interpolate( system.analyticalSol_, maxLevel );

   // Determine initial residual norm
   A.apply( u, tmp, maxLevel, DoFType::All );
   res.assign( { 1.0, -1.0 }, { f, tmp }, maxLevel, DoFType::Inner );
   const real_t initResNorm = std::sqrt( res.dotGlobal( res, maxLevel, DoFType::Inner ) );

   // Solve system.
   real_t residualNorm = 1.0;
   uint_t its          = 0;
   for ( ; its < numMaxVCycles && residualNorm / initResNorm > 1e-6; ++its )
   {
      gmgSolver.solve( A, u, f, maxLevel );

      A.apply( u, tmp, maxLevel, DoFType::All );
      res.assign( { 1.0, -1.0 }, { f, tmp }, maxLevel, DoFType::Inner );
      residualNorm = std::sqrt( res.dotGlobal( res, maxLevel, DoFType::Inner ) );
      // WALBERLA_LOG_DEVEL_VAR_ON_ROOT( residualNorm / initResNorm )
   }

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "N1E1hIndependenceTest", storage );
      vtk.add( u );
      vtk.add( f );
      vtk.add( sol );
      vtk.add( res );
      vtk.write( maxLevel );
   }

   return its;
}

void hIndependenceTest( const uint_t minLevel, const uint_t maxLevel, const n1e1::System& system )
{
   const uint_t maxIts         = 10;
   const uint_t maxAbsIncrease = 1;
   const real_t maxRelIncrease = real_c( 0.2 );

   uint_t nItsCoarse = test( minLevel, maxIts, system );
   WALBERLA_LOG_INFO_ON_ROOT( "level " << minLevel << ": " << nItsCoarse )
   WALBERLA_CHECK_LESS( nItsCoarse, maxIts );

   for ( uint_t fineLevel = minLevel + 1; fineLevel <= maxLevel; ++fineLevel )
   {
      const uint_t nItsFine = test( fineLevel, maxIts, system );
      WALBERLA_LOG_INFO_ON_ROOT( "level " << fineLevel << ": " << nItsFine )
      WALBERLA_CHECK_LESS( nItsFine, maxIts );

      WALBERLA_CHECK_LESS_EQUAL( nItsFine, nItsCoarse + maxAbsIncrease )
      WALBERLA_CHECK_LESS_EQUAL( real_c( nItsFine ) / real_c( nItsCoarse ), 1.0 + maxRelIncrease )

      nItsCoarse = nItsFine;
   }
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   hIndependenceTest( 3, 6, n1e1::System::polynomialOnSingleTet() );
   hIndependenceTest( 2, 5, n1e1::System::sinusoidalOnCube() );
}
