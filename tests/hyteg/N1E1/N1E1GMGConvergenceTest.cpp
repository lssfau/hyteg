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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_curl_curl_affine_qe.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_curl_curl_blending_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_blending_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_mass_affine_qe.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_mass_blending_q2.hpp"
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
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

#include "common.hpp"

using namespace hyteg;
using walberla::real_t;

/// Returns the approximate L2 error.
template < class N1E1CurlCurlForm,
           class N1E1MassForm,
           class N1E1LinearForm,
           class N1E1MassOperator,
           class P1LaplaceOperator,
           class P1Smoother >
real_t test( const uint_t maxLevel, const n1e1::System& system, const bool writeVTK = false )
{
   using namespace n1e1;

   using N1E1Smoother = ChebyshevSmoother< N1E1ElementwiseLinearCombinationOperator >;

   const uint_t minLevel                = 0;
   const uint_t spectralRadiusEstLevel  = std::min( uint_c( 3 ), maxLevel );
   const int    numSpectralRadiusEstIts = 40;
   const int    numVCycles              = 8;

   SetupPrimitiveStorage setupStorage( system.domain_, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   system.setMap_( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   N1E1CurlCurlForm curlCurlForm;
   N1E1MassForm     massForm;

   N1E1MassOperator                         M( storage, minLevel, maxLevel );
   N1E1ElementwiseLinearCombinationOperator A( storage, minLevel, maxLevel, { { 1.0, 1.0 }, { &curlCurlForm, &massForm } } );

   N1E1VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > sol( "sol", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > err( "err", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   const uint_t nDoFs = numberOfGlobalDoFs( u, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "dofs on level " << maxLevel << ": " << nDoFs );

   // Assemble RHS.
   assembleLinearForm< N1E1LinearForm >( maxLevel, maxLevel, { system.rhs_ }, f );

   // Boundary conditions: homogeneous tangential trace
   u.interpolate( Point3D{ 0.0, 0.0, 0.0 }, maxLevel, DoFType::Boundary );

   // Hybrid smoother
   auto p1LaplaceOperator = std::make_shared< P1LaplaceOperator >( storage, minLevel, maxLevel );
   auto chebyshevSmoother = std::make_shared< N1E1Smoother >( storage, minLevel, maxLevel );

   std::shared_ptr< P1Smoother > p1Smoother;
   if constexpr ( std::is_same< P1Smoother, WeightedJacobiSmoother< P1LaplaceOperator > >::value )
   {
      p1Smoother = std::make_shared< P1Smoother >( storage, minLevel, maxLevel, 2.0 / 3.0 );
   }
   else
   {
      p1Smoother = std::make_shared< P1Smoother >();
   }

   sol.interpolate( system.analyticalSol_, spectralRadiusEstLevel );
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
   sol.interpolate( system.analyticalSol_, maxLevel );

   // Solve system.
   real_t discrL2 = 0.0;

   for ( int i = 0; i < numVCycles; ++i )
   {
      gmgSolver.solve( A, u, f, maxLevel );

      err.assign( { 1.0, -1.0 }, { u, sol }, maxLevel );
      M.apply( err, tmp, maxLevel, DoFType::All );
      discrL2 = std::sqrt( err.dotGlobal( tmp, maxLevel ) );
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( discrL2 )

      // determine residual
      A.apply( u, tmp, maxLevel, DoFType::Inner );
      tmp.assign( { 1.0, -1.0 }, { f, tmp }, maxLevel, DoFType::Inner );
      const real_t residual = std::sqrt( tmp.dotGlobal( tmp, maxLevel, DoFType::Inner ) );
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( residual )
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

real_t testNoBlending( const uint_t maxLevel, const n1e1::System& system, const bool writeVTK = false )
{
   return test< forms::n1e1_curl_curl_affine_qe,
                forms::n1e1_mass_affine_qe,
                forms::n1e1_linear_form_affine_q6,
                n1e1::N1E1ElementwiseMassOperator,
                P1ConstantLaplaceOperator,
                GaussSeidelSmoother< P1ConstantLaplaceOperator > >( maxLevel, system, writeVTK );
}

real_t testBlending( const uint_t maxLevel, const n1e1::System& system, const bool writeVTK = false )
{
   return test< forms::n1e1_curl_curl_blending_q2,
                forms::n1e1_mass_blending_q2,
                forms::n1e1_linear_form_blending_q6,
                n1e1::N1E1ElementwiseBlendingMassOperatorQ2,
                P1ElementwiseBlendingLaplaceOperator,
                WeightedJacobiSmoother< P1ElementwiseBlendingLaplaceOperator > >( maxLevel, system, writeVTK );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, polynomial ###" );
   n1e1::L2ConvergenceTest( 4, 6, n1e1::System::polynomialOnSingleTet(), testNoBlending );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, sinusoidal ###" );
   n1e1::L2ConvergenceTest( 4, 6, n1e1::System::sinusoidalOnSingleTet(), testNoBlending );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, polynomial ###" );
   n1e1::L2ConvergenceTest( 4, 5, n1e1::System::polynomialOnCube(), testNoBlending );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, sinusoidal ###" );
   n1e1::L2ConvergenceTest( 4, 5, n1e1::System::sinusoidalOnCube(), testNoBlending );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Test on slice of solid torus, blending ###" );
   n1e1::L2ConvergenceTest( 4, 5, n1e1::System::onToroidalSlice(), testBlending );

   return EXIT_SUCCESS;
}
